"""Independent checks of the recording ingestion/calibration contract."""

import importlib
import importlib.util
import hashlib
import json
import subprocess
import sys
import warnings

import numpy as np
import pytest
import soundfile as sf


def pipeline():
    assert importlib.util.find_spec("oceansignatures.pipeline") is not None, (
        "recording pipeline has not been implemented"
    )
    return importlib.import_module("oceansignatures.pipeline")


def test_scalar_adc_calibration_uses_voltage_gain_and_micropascal_reference():
    module = pipeline()
    samples = np.array([-0.5, 0.0, 0.5])
    calibration = module.Calibration(
        sensitivity_db_v_per_upa=-160.0, gain_db=20.0, full_scale_voltage=2.0
    )
    # V = 0.5 * 2 / 10; sensitivity = 1e-8 V/uPa; 1e6 uPa = 1 Pa.
    result = module.calibrate_samples(samples, calibration)
    np.testing.assert_allclose(result, [-10.0, 0.0, 10.0], rtol=1e-14, atol=0)


@pytest.mark.parametrize(
    "field,value",
    [
        ("full_scale_voltage", 0.0),
        ("full_scale_voltage", -1.0),
        ("gain_db", float("nan")),
        ("sensitivity_db_v_per_upa", float("inf")),
        ("sensitivity_db_v_per_upa", -10000.0),
    ],
)
def test_invalid_calibration_is_rejected(field, value):
    values = dict(sensitivity_db_v_per_upa=-160.0, gain_db=20.0, full_scale_voltage=2.0)
    values[field] = value
    with pytest.raises(ValueError):
        pipeline().Calibration(**values)


@pytest.mark.parametrize("coefficients", [(-600.0, -600.0, 1e300), (600.0, 600.0, 1e-300)])
@pytest.mark.parametrize("samples", [[0.0, 0.0], [-0.5, 0.5]])
def test_unrepresentable_combined_calibration_is_rejected_without_warnings(
    tmp_path, coefficients, samples
):
    module = pipeline()
    calibration = module.Calibration(*coefficients)
    path = tmp_path / "combined-scale.wav"
    sf.write(path, samples, 8000, subtype="DOUBLE")
    with warnings.catch_warnings(record=True) as seen:
        warnings.simplefilter("always")
        with pytest.raises(ValueError, match="calibration.*scale.*representable"):
            module.calibrate_samples(samples, calibration)
        with pytest.raises(ValueError, match="calibration.*scale.*representable"):
            module.analyse_recording(path, calibration)
    assert not seen, [str(w.message) for w in seen]


def test_nonzero_pressure_conversion_cannot_underflow_to_silence(tmp_path):
    module = pipeline()
    calibration = module.Calibration(600.0, 600.0, 1e-254)
    samples = np.array([1e-4, -1e-4])
    path = tmp_path / "pressure-underflow.wav"
    sf.write(path, samples, 8000, subtype="DOUBLE")
    with warnings.catch_warnings(record=True) as seen:
        warnings.simplefilter("always")
        with pytest.raises(ValueError, match="Calibrated pressure.*representable"):
            module.calibrate_samples(samples, calibration)
        with pytest.raises(ValueError, match="Calibrated pressure.*representable"):
            module.analyse_recording(path, calibration)
    assert not seen, [str(w.message) for w in seen]


@pytest.mark.parametrize("coefficients", [(-600.0, -600.0, 1e200), (600.0, 600.0, 1e-100)])
def test_unrepresentable_pressure_square_is_rejected_without_warnings(tmp_path, coefficients):
    module = pipeline()
    samples = 0.2 * np.sin(2 * np.pi * 3000 * np.arange(48000) / 48000)
    path = tmp_path / "pressure-square.wav"
    sf.write(path, samples, 48000, subtype="DOUBLE")
    calibration = module.Calibration(*coefficients)
    with warnings.catch_warnings(record=True) as seen:
        warnings.simplefilter("always")
        pressure = module.calibrate_samples(samples, calibration)
        assert np.isfinite(pressure).all() and np.any(pressure != 0)
        with pytest.raises(ValueError, match="pressure-square.*representable"):
            module.analyse_recording(path, calibration)
    assert not seen, [str(w.message) for w in seen]


def test_unrepresentable_psd_is_rejected_without_warnings(tmp_path):
    module = pipeline()
    samples = np.ones(2048)
    path = tmp_path / "psd-overflow.wav"
    sf.write(path, samples, 1, subtype="DOUBLE")
    calibration = module.Calibration(-120.0, 0.0, 1e153)
    # This calibration is 1e153 Pa/fraction, so each pressure-square is 1e306 Pa².
    # Periodic Hann DC density: (sum(w))² / (fs * sum(w²)) = 2*N/(3*fs).
    assert np.isfinite(1e306)
    assert 2 * len(samples) / 3 > np.finfo(float).max / 1e306
    with warnings.catch_warnings(record=True) as seen:
        warnings.simplefilter("always")
        with pytest.raises(ValueError, match="PSD.*representable"):
            module.analyse_recording(path, calibration)
    assert not seen, [str(w.message) for w in seen]


def test_nonzero_psd_cannot_underflow_to_zero_density(tmp_path):
    module = pipeline()
    path = tmp_path / "psd-underflow.wav"
    sf.write(path, np.ones(2), 48000, subtype="DOUBLE")
    calibration = module.Calibration(-120.0, 0.0, 1e-161)
    # Each pressure-square is nonzero float64, but the two-sample Hann density
    # is pressure-square / fs: smaller than the least positive float64.
    assert 1e-161 * 1e-161 > 0
    assert (1e-161 * 1e-161) / 48000 == 0
    with warnings.catch_warnings(record=True) as seen:
        warnings.simplefilter("always")
        with pytest.raises(ValueError, match="PSD.*representable"):
            module.analyse_recording(path, calibration)
    assert not seen, [str(w.message) for w in seen]


@pytest.mark.parametrize("pressure_scale,fs", [(1e153, 8000), (1e154, 8000), (1e153, 10)])
def test_representable_large_metrics_use_stable_reductions(tmp_path, pressure_scale, fs):
    module = pipeline()
    path = tmp_path / "large-representable.wav"
    sf.write(path, np.ones(4096), fs, subtype="DOUBLE")
    calibration = module.Calibration(-120.0, 0.0, pressure_scale)
    with warnings.catch_warnings(record=True) as seen:
        warnings.simplefilter("always")
        result = module.analyse_recording(path, calibration)
        module.write_analysis(result, tmp_path / "products")
    assert not seen, [str(w.message) for w in seen]
    assert result["rms"] == pytest.approx(pressure_scale, rel=1e-14)
    assert result["spl_db_re_1upa"] == pytest.approx(20 * (np.log10(pressure_scale) + 6))
    # Exact periodic-Hann DC density for the 2048-sample Welch segments.
    assert result["psd"][0] / pressure_scale / pressure_scale == pytest.approx(
        2 * 2048 / (3 * fs), rel=1e-14
    )
    for key in ("time_s", "signal", "frequency_hz", "psd"):
        assert np.isfinite(result[key]).all()


def test_recording_analysis_distinguishes_unknown_and_declared_calibration(tmp_path):
    module = pipeline()
    assert callable(getattr(module, "analyse_recording", None)), "missing recording analysis"
    fs = 8000
    waveform = 0.5 * np.cos(2 * np.pi * 1000 * np.arange(fs) / fs)
    path = tmp_path / "tone.wav"
    sf.write(path, waveform, fs, subtype="DOUBLE")
    unknown = module.analyse_recording(path)
    assert unknown["unit"] == "normalized digital amplitude"
    assert unknown["spl_db_re_1upa"] is None
    assert unknown["calibration_status"] == "unknown"
    assert unknown["sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert unknown["sample_count"] == fs
    assert unknown["rms"] == pytest.approx(np.sqrt(0.125))
    frequency = unknown["frequency_hz"]
    psd = unknown["psd"]
    assert frequency[np.argmax(psd)] == pytest.approx(1000.0)
    assert np.sum(psd) * (frequency[1] - frequency[0]) == pytest.approx(0.125, rel=1e-10)
    calibration = module.Calibration(-160.0, 20.0, 2.0)
    calibrated = module.analyse_recording(path, calibration=calibration)
    assert calibrated["unit"] == "Pa"
    assert calibrated["rms"] == pytest.approx(np.sqrt(50))
    assert calibrated["spl_db_re_1upa"] == pytest.approx(10 * np.log10(50 / 1e-12))
    assert "scalar" in calibrated["calibration_status"]


@pytest.mark.parametrize(
    "coefficients,expected_scale",
    [(None, None), ((-180.0, 20.0, 2.0), 200.0), ((-160.0, 20.0, 2.0), 20.0)],
)
def test_analysis_and_export_retain_declared_calibration_provenance(
    tmp_path, coefficients, expected_scale
):
    module = pipeline()
    samples = np.array([0.0, 0.25, -0.5, 1.0])
    path = tmp_path / "provenance.wav"
    sf.write(path, samples, 8000, subtype="DOUBLE")
    calibration = module.Calibration(*coefficients) if coefficients is not None else None
    analysis = module.analyse_recording(path, calibration)
    output = tmp_path / "products"
    module.write_analysis(analysis, output)
    summary = json.loads((output / "summary.json").read_text())
    expected = {
        "calibration": (
            dict(zip(("sensitivity_db_v_per_upa", "gain_db", "full_scale_voltage"), coefficients))
            if coefficients is not None
            else None
        ),
        "pressure_scale_pa_per_digital_fraction": expected_scale,
        "calibration_assumption": (
            "flat scalar response over the analysed band" if calibration is not None else None
        ),
        "adc_full_scale_convention": (
            "digital fractions +/-1 correspond to +/-full_scale_voltage V at the ADC input "
            "(peak, not RMS or peak-to-peak)"
            if calibration is not None
            else None
        ),
    }
    for result in (analysis, summary):
        for key, value in expected.items():
            assert key in result, f"Missing calibration provenance: {key}"
            if key == "pressure_scale_pa_per_digital_fraction" and value is not None:
                assert result[key] == pytest.approx(value, rel=1e-14)
            else:
                assert result[key] == value
        if calibration is None:
            assert result["unit"] == "normalized digital amplitude"
            assert result["spl_db_re_1upa"] is None
    np.testing.assert_allclose(
        analysis["signal"],
        samples if calibration is None else samples * expected_scale,
        rtol=1e-14,
        atol=0,
    )


def test_qc_reports_potential_full_scale_samples_without_claiming_hardware_clipping(tmp_path):
    path = tmp_path / "qc.wav"
    sf.write(path, [0.0, 1.0, -1.0, 0.5], 8000, subtype="DOUBLE")
    result = pipeline().analyse_recording(path)
    assert "qc" in result, "QC status has not been implemented"
    assert result["qc"]["near_full_scale_samples"] == 2
    assert result["qc"]["timing_continuity"] == "not established from a single audio file"


@pytest.mark.parametrize("samples", [[], [float("nan"), 0.1], [1.1, 0.0]])
def test_unsupported_audio_values_are_rejected(tmp_path, samples):
    path = tmp_path / "invalid.wav"
    sf.write(path, np.asarray(samples), 8000, subtype="DOUBLE")
    with pytest.raises(ValueError):
        pipeline().analyse_recording(path)


def test_channel_and_file_size_guard_are_explicit(tmp_path):
    path = tmp_path / "stereo.wav"
    sf.write(path, np.zeros((32, 2)), 8000, subtype="DOUBLE")
    with pytest.raises(ValueError):
        pipeline().analyse_recording(path, channel=-1)
    with pytest.raises(ValueError):
        pipeline().analyse_recording(path, max_bytes=1)


def test_export_is_reproducible_and_preserves_unrelated_existing_outputs(tmp_path):
    module = pipeline()
    assert callable(getattr(module, "write_analysis", None)), "missing export pipeline"
    audio = tmp_path / "source.wav"
    sf.write(audio, np.linspace(-0.5, 0.5, 128), 8000, subtype="DOUBLE")
    result = module.analyse_recording(audio)
    directory = tmp_path / "output"
    paths = module.write_analysis(result, directory)
    first = {path.name: path.read_bytes() for path in paths}
    assert set(first) == {"summary.json", "spectrum.csv", "signal.npz"}
    summary = json.loads(first["summary.json"])
    assert summary["spl_db_re_1upa"] is None
    assert summary["sha256"] == hashlib.sha256(audio.read_bytes()).hexdigest()
    module.write_analysis(result, directory)
    assert first == {path.name: path.read_bytes() for path in paths}
    (directory / "summary.json").write_text("unrelated")
    with pytest.raises(FileExistsError):
        module.write_analysis(result, directory)
    assert (directory / "summary.json").read_text() == "unrelated"


def test_cli_exercises_the_actual_file_to_products_path(tmp_path):
    assert callable(getattr(pipeline(), "main", None)), "missing executable pipeline"
    path = tmp_path / "input.wav"
    sf.write(path, np.zeros(800), 8000, subtype="PCM_16")
    destination = tmp_path / "products"
    result = subprocess.run(
        [sys.executable, "-m", "oceansignatures.pipeline", str(path), str(destination)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    summary = json.loads((destination / "summary.json").read_text())
    assert summary["sample_count"] == 800
    assert summary["spl_db_re_1upa"] is None
    assert summary["calibration_status"] == "unknown"


@pytest.mark.parametrize("coefficients", [None, (-180.0, 20.0, 2.0), (-600.0, -600.0, 1e200)])
def test_silence_has_finite_zero_products_without_a_pressure_level(tmp_path, coefficients):
    module = pipeline()
    path = tmp_path / "silence.wav"
    sf.write(path, np.zeros(32), 8000, subtype="DOUBLE")
    calibration = module.Calibration(*coefficients) if coefficients is not None else None
    result = module.analyse_recording(path, calibration)
    assert result["rms"] == 0
    assert result["spl_db_re_1upa"] is None
    assert np.array_equal(result["signal"], np.zeros(32))
    assert np.all(result["psd"] == 0)
    assert result["qc"]["finite_samples"] is True
    module.write_analysis(result, tmp_path / "products")


@pytest.mark.parametrize(
    "coefficients,expected_scale",
    [((600.0, -600.0, 1e308), 1e302), ((-600.0, 600.0, 1e-300), 1e-306)],
)
def test_combined_conversion_avoids_unrepresentable_intermediate_voltages(
    coefficients, expected_scale
):
    module = pipeline()
    samples = np.array([-0.5, 0.0, 0.5])
    pressure = module.calibrate_samples(samples, module.Calibration(*coefficients))
    np.testing.assert_allclose(pressure / expected_scale, samples, rtol=1e-14, atol=0)


@pytest.mark.parametrize("size", [3, 127, 128, 5001])
def test_calibrated_density_matches_independent_hann_dft(tmp_path, size):
    module = pipeline()
    fs = 8000
    samples = np.random.default_rng(793).uniform(-0.3, 0.7, size)
    path = tmp_path / "independent-psd.wav"
    sf.write(path, samples, fs, subtype="DOUBLE")
    result = module.analyse_recording(path, module.Calibration(-180.0, 20.0, 2.0))
    pressure = samples * 200
    n = min(2048, size)
    window = 0.5 - 0.5 * np.cos(2 * np.pi * np.arange(n) / n)
    starts = range(0, size - n + 1, n - n // 2)
    transforms = np.fft.rfft(np.array([pressure[i : i + n] * window for i in starts]))
    density = np.abs(transforms) ** 2 / (fs * np.sum(window**2))
    density[:, 1 : (-1 if n % 2 == 0 else None)] *= 2
    np.testing.assert_allclose(result["frequency_hz"], np.fft.rfftfreq(n, 1 / fs))
    np.testing.assert_allclose(result["psd"], density.mean(axis=0), rtol=1e-12, atol=1e-14)
    assert result["rms"] == pytest.approx(np.sqrt(np.mean(pressure**2)), rel=1e-14)


@pytest.mark.parametrize(
    "coefficients,valid", [((-180.0, 20.0, 2.0), True), ((-600.0, -600.0, 1e200), False)]
)
def test_calibrated_cli_preserves_provenance_or_rejects_before_writing(
    tmp_path, coefficients, valid
):
    path = tmp_path / "input.wav"
    sf.write(path, np.full(32, 0.2), 8000, subtype="DOUBLE")
    values = dict(zip(("sensitivity_db_v_per_upa", "gain_db", "full_scale_voltage"), coefficients))
    calibration_path = tmp_path / "calibration.json"
    calibration_path.write_text(json.dumps(values))
    destination = tmp_path / "products"
    result = subprocess.run(
        [
            sys.executable,
            "-W",
            "error",
            "-m",
            "oceansignatures.pipeline",
            str(path),
            str(destination),
            "--calibration-json",
            str(calibration_path),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert "RuntimeWarning" not in result.stderr
    if valid:
        assert result.returncode == 0, result.stderr
        summary = json.loads((destination / "summary.json").read_text())
        assert summary["calibration"] == values
        assert summary["pressure_scale_pa_per_digital_fraction"] == pytest.approx(200)
        assert summary["unit"] == "Pa"
        assert summary["rms"] == pytest.approx(40)
        assert (destination / "spectrum.csv").read_text().splitlines()[
            0
        ] == "frequency_hz,psd_Pa2_per_Hz"
    else:
        assert result.returncode != 0
        assert (
            "ValueError: Signal pressure-square/amplitude-square is not representable"
            in result.stderr
        )
        assert not destination.exists()
