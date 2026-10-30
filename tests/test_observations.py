"""Independent numerical and clock controls for offline observations."""

from datetime import datetime, timedelta, timezone
import hashlib
import importlib
import importlib.util
import json
import csv

import numpy as np
import pytest
import soundfile as sf

from oceansignatures.pipeline import Calibration

EPOCH = datetime(2026, 1, 1, tzinfo=timezone.utc)


def observations():
    assert importlib.util.find_spec("oceansignatures.observations") is not None, (
        "timestamped observation API is missing"
    )
    return importlib.import_module("oceansignatures.observations")


def chunk(tmp_path, samples, *, start=0.0, fs=8192, name="record.wav", **kwargs):
    module = observations()
    path = tmp_path / name
    sf.write(path, samples, fs, subtype="DOUBLE")
    return module.RecordingChunk(
        path=path,
        start_utc=EPOCH + timedelta(seconds=start),
        evidence_kind="synthetic",
        **kwargs,
    )


def test_timestamped_tone_reuses_calibration_and_keeps_source_identity(tmp_path):
    module = observations()
    fs = 8192
    x = 0.5 * np.cos(2 * np.pi * 512 * np.arange(2 * fs) / fs)
    recording = chunk(
        tmp_path,
        np.column_stack((np.zeros_like(x), x)),
        channel=1,
        calibration=Calibration(-160, 20, 2),
    )
    before = recording.path.read_bytes()
    result = module.analyse_chunks([recording], bands_hz={"tone": (400, 600)})
    assert result["unit"] == "Pa"
    assert result["calibration_status"] == "declared scalar calibration"
    assert result["psd"].shape == (1, 513)
    assert result["valid"].tolist() == [True]
    assert result["window_start_sample"].tolist() == [0]
    assert result["window_stop_sample"].tolist() == [2 * fs]
    assert result["window_start_utc"].tolist() == ["2026-01-01T00:00:00.000000Z"]
    assert result["window_stop_s"].tolist() == [2.0]
    assert result["band_power"][0, 0] == pytest.approx(50, rel=1e-12)
    assert result["band_spl_db_re_1upa"][0, 0] == pytest.approx(10 * np.log10(50) + 120)
    provenance = result["metadata"]["chunks"][0]
    assert provenance["sha256"] == hashlib.sha256(before).hexdigest()
    assert provenance["channel"] == 1
    assert provenance["evidence_kind"] == "synthetic"
    assert provenance["calibration"] == {
        "sensitivity_db_v_per_upa": -160,
        "gain_db": 20,
        "full_scale_voltage": 2,
    }
    assert provenance["pressure_scale_pa_per_digital_fraction"] == pytest.approx(20)
    assert "peak, not RMS" in provenance["adc_full_scale_convention"]
    assert recording.path.read_bytes() == before


def test_non_aligned_chunk_boundaries_preserve_anchored_windows(tmp_path):
    module = observations()
    fs = 8192
    x = np.random.default_rng(18).uniform(-0.4, 0.4, 4 * fs)
    whole = module.analyse_chunks([chunk(tmp_path, x, name="whole.wav")])
    cuts = [0, 3173, 19111, len(x)]
    pieces = [
        chunk(tmp_path, x[a:b], start=a / fs, name=f"part-{a}.wav")
        for a, b in zip(cuts[:-1], cuts[1:])
    ]
    split = module.analyse_chunks(pieces[::-1], anchor_utc=EPOCH)
    for key in ("window_start_sample", "window_stop_sample", "psd", "stft", "band_power"):
        np.testing.assert_array_equal(split[key], whole[key])
    assert split["band_spl_db_re_1upa"] is None
    assert split["metadata"]["anchor_utc"] == "2026-01-01T00:00:00.000000Z"
    assert len(split["metadata"]["chunks"]) == 3
    assert split["metadata"]["processing"]["detrend"] is False
    assert split["metadata"]["processing"]["window"] == "hann"
    assert split["metadata"]["processing"]["noverlap"] == 512


def test_gaps_partial_support_and_active_contact_are_not_silence(tmp_path):
    module = observations()
    a = chunk(tmp_path, np.ones(8192), start=0, name="a.wav")
    b = chunk(tmp_path, np.full(3 * 8192, 0.2), start=3, name="b.wav")
    c = chunk(tmp_path, np.ones(8192 // 2), start=6, name="c.wav")
    result = module.analyse_chunks(
        [a, b, c], active_intervals=[(EPOCH + timedelta(seconds=4), EPOCH + timedelta(seconds=5))]
    )
    assert result["window_start_s"].tolist() == [0, 2, 4, 6]
    assert result["valid"].tolist() == [False, False, True, False]
    assert result["coverage_fraction"].tolist() == [0.5, 0.5, 1.0, 0.25]
    assert np.isnan(result["psd"][[0, 1, 3]]).all()
    assert np.isnan(result["stft"][[0, 1, 3]]).all()
    assert np.isnan(result["band_power"][[0, 1, 3]]).all()
    assert result["active_transmission"].tolist() == [False, False, True, False]
    assert result["near_full_scale_samples"].tolist() == [8192, 0, 0, 4096]
    assert result["qc_ok"].tolist() == [False, False, False, False]
    assert result["partial_window"].tolist() == [False, False, False, True]
    assert result["metadata"]["gaps"] == [{"start_sample": 8192, "stop_sample": 24576}]


def test_duplicate_and_overlap_chunks_are_rejected_not_overwritten(tmp_path):
    module = observations()
    a = chunk(tmp_path, np.ones(16384), name="a.wav")
    with pytest.raises(ValueError, match="overlap|duplicate"):
        module.analyse_chunks([a, a])
    b = chunk(tmp_path, np.ones(8192), start=1, name="b.wav")
    with pytest.raises(ValueError, match="overlap|duplicate"):
        module.analyse_chunks([a, b])


@pytest.mark.parametrize("n", [15, 16])
@pytest.mark.parametrize("window", ["hann", "boxcar"])
@pytest.mark.parametrize("detrend", [False, "constant"])
def test_independent_dft_parseval_and_declared_detrend(n, window, detrend):
    # Verification-only controls for the initial shared kernel.
    module = observations()
    x = np.arange(n, dtype=float) / n + 0.2 * (-1.0) ** np.arange(n)
    y = x - x.mean() if detrend == "constant" else x
    w = 0.5 - 0.5 * np.cos(2 * np.pi * np.arange(n) / n) if window == "hann" else np.ones(n)
    k = np.arange(n // 2 + 1)
    direct = np.exp(-2j * np.pi * k[:, None] * np.arange(n) / n) @ (y * w)
    expected = np.abs(direct) ** 2 / (100 * np.sum(w * w))
    expected[1 : -1 if n % 2 == 0 else None] *= 2
    result = module.frame_spectra(x, 100, nperseg=n, noverlap=0, window=window, detrend=detrend)
    np.testing.assert_allclose(result["psd"], expected, atol=1e-16, rtol=1e-12)
    assert result["psd"].sum() * 100 / n == pytest.approx(np.sum((y * w) ** 2) / np.sum(w * w))
    np.testing.assert_allclose(np.abs(result["stft"][0]) ** 2, expected, atol=1e-16)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"detrend": "linear"},
        {"detrend": True},
        {"window": "bad"},
        {"nperseg": 2.5},
        {"nperseg": True},
        {"nperseg": 17},
        {"noverlap": -1},
        {"noverlap": 16},
        {"noverlap": True},
    ],
)
def test_spectral_contract_rejects_implicit_changes(kwargs):
    params = {"nperseg": 16, **kwargs}
    with pytest.raises(ValueError):
        observations().frame_spectra(np.ones(16), 100, **params)


@pytest.mark.parametrize("fs", [1, 44100, 384000])
def test_weak_bin_survives_frame_dynamic_range(fs):
    # Separate DC/Nyquist segments: no subtraction or spectral leakage in the oracle.
    result = observations().frame_spectra(
        [1e153, 1e153, 1e-10, -1e-10], fs, nperseg=2, noverlap=0, window="boxcar"
    )
    np.testing.assert_allclose(result["psd"], [1e306 / fs, 1e-20 / fs], rtol=1e-13, atol=0)
    np.testing.assert_allclose(
        result["psd"], np.mean(np.abs(result["stft"]) ** 2, axis=0), rtol=1e-13, atol=0
    )
    assert result["stft"][0, 1] == result["stft"][1, 0] == 0


@pytest.mark.parametrize("n", [2, 1024])
def test_calibrated_wav_weak_nyquist_band_is_not_silence(tmp_path, n):
    fs = 44100
    large = chunk(
        tmp_path, np.full(n, 0.5), fs=fs, name="large.wav", calibration=Calibration(-120, 0, 2e153)
    )
    small = chunk(
        tmp_path,
        1e-10 * (-1.0) ** np.arange(n),
        fs=fs,
        start=n / fs,
        name="small.wav",
        calibration=Calibration(-120, 0, 1),
    )
    result = observations().analyse_chunks(
        [large, small],
        frame_seconds=2 * n / fs,
        nperseg=n,
        noverlap=0,
        window="boxcar",
        bands_hz={"nyquist": (fs / 2 - fs / (2 * n), fs / 2)},
    )
    assert result["valid"].tolist() == result["qc_ok"].tolist() == [True]
    assert result["unit"] == "Pa"
    assert np.isfinite(result["stft"]).all()
    assert result["band_power"][0, 0] == pytest.approx(5e-21, rel=1e-13, abs=0)
    assert result["band_spl_db_re_1upa"][0, 0] == pytest.approx(10 * np.log10(5e-21) + 120)
    np.testing.assert_allclose(
        result["psd"], np.mean(np.abs(result["stft"]) ** 2, axis=1), rtol=1e-13, atol=0
    )


@pytest.mark.parametrize("weak", [1e-158, 1e-161])
def test_weak_bin_float64_boundary_is_not_false_silence(weak):
    params = {"nperseg": 2, "noverlap": 0, "window": "boxcar"}
    x = [1e153, 1e153, weak, -weak]
    if weak == 1e-161:
        with pytest.raises(ValueError, match="representable"):
            observations().frame_spectra(x, 384000, **params)
    else:
        # Representable subnormal PSD: permit two float64 subnormal rounding steps.
        expected = weak**2 / 384000
        result = observations().frame_spectra(x, 384000, **params)
        assert result["psd"][1] > 0
        assert result["psd"][1] == pytest.approx(expected, rel=0, abs=2 * np.nextafter(0.0, 1.0))


@pytest.mark.parametrize("samples,fs", [(np.full(16, 1e154), 2), ([1e154, 1e154, 0, 0], 1)])
def test_large_finite_mean_power_does_not_require_finite_unscaled_sum(samples, fs):
    result = observations().frame_spectra(samples, fs, nperseg=2, noverlap=0, window="boxcar")
    assert result["psd"][0] == pytest.approx(1e308, rel=1e-13)
    assert result["psd"][1] == 0


def test_weak_transform_underflow_before_physical_restoration_is_rejected():
    from decimal import Decimal, localcontext

    with localcontext() as context:
        context.prec = 100
        amplitude = Decimal("4e-162") / Decimal("2e20").sqrt()
        density = amplitude**2 / 2
    assert float(amplitude) > 0 and density > 0 and float(density) == 0
    # Both raw sample squares survive; only the final weak PSD is unrepresentable.
    assert np.float64(2e-162) ** 2 > 0
    with pytest.raises(ValueError, match="representable"):
        observations().frame_spectra(
            [1e153, 1e153, 2e-162, -2e-162], 1e20, nperseg=2, noverlap=0, window="boxcar"
        )


@pytest.mark.parametrize("amplitude,fs", [(1e154, 1e308), (1e-154, 1e-309), (0.0, 1e308)])
def test_standalone_rate_extremes_preserve_representable_spectrum(amplitude, fs):
    from decimal import Decimal, localcontext

    with localcontext() as context:
        context.prec = 100
        expected_psd = float(2 * Decimal.from_float(amplitude) ** 2 / Decimal.from_float(fs))
        expected_stft = float(
            2 * Decimal.from_float(amplitude) / (2 * Decimal.from_float(fs)).sqrt()
        )
    result = observations().frame_spectra(
        [amplitude, amplitude], fs, nperseg=2, noverlap=0, window="boxcar"
    )
    np.testing.assert_allclose(result["psd"], [expected_psd, 0], rtol=1e-13, atol=0)
    np.testing.assert_allclose(result["stft"][0], [expected_stft, 0], rtol=1e-13, atol=0)
    np.testing.assert_array_equal(result["frequency_hz"], [0, fs / 2])


@pytest.mark.parametrize("fs", [1, 1e-308])
def test_unrepresentable_derived_spectrum_rejected_without_warnings(fs):
    # Respect the raw-square guard while the final physical PSD overflows.
    with pytest.raises(ValueError, match="representable"):
        observations().frame_spectra(np.full(4, 1e154), fs, nperseg=2, window="boxcar")


@pytest.mark.parametrize("detrend", [False, "constant"])
def test_exact_spectral_zeros_remain_valid(detrend):
    result = observations().frame_spectra(
        np.full(16, 1e153), 44100, nperseg=2, noverlap=0, window="boxcar", detrend=detrend
    )
    assert result["psd"][1] == 0
    assert np.all(result["stft"][:, 1] == 0)
    if detrend == "constant":
        assert np.all(result["psd"] == 0)
        assert np.all(result["stft"] == 0)


@pytest.mark.parametrize("n", [3, 5, 7, 16, 31, 128, 1024])
@pytest.mark.parametrize("window", ["boxcar", "hann"])
@pytest.mark.parametrize("amplitude", [1e-154, 0.1, 1e153, -0.1, 0.0])
def test_exact_constant_detrend_is_true_silence(n, window, amplitude):
    # Exact repeated input has no residual, irrespective of sum/mean roundoff.
    result = observations().frame_spectra(
        np.full(2 * n, amplitude),
        44100,
        nperseg=n,
        noverlap=0,
        window=window,
        detrend="constant",
    )
    np.testing.assert_array_equal(result["segment_start_sample"], [0, n])
    assert np.all(result["stft"] == 0)
    assert np.all(result["psd"] == 0)


@pytest.mark.parametrize("n", [3, 7, 31, 128, 1024])
@pytest.mark.parametrize("window", ["boxcar", "hann"])
def test_nonconstant_detrend_matches_independent_complex_dft(n, window):
    # Verification-only: scalar compensated mean and direct DFT, not the kernel.
    import math

    fs = 44100
    x = np.random.default_rng(83).uniform(-0.2, 0.8, 3 * n + 1)
    hop = max(1, n // 3)
    starts = list(range(0, len(x) - n + 1, hop))
    w = np.ones(n) if window == "boxcar" else 0.5 - 0.5 * np.cos(2 * np.pi * np.arange(n) / n)
    d = np.full(n // 2 + 1, 2.0)
    d[0] = 1
    if n % 2 == 0:
        d[-1] = 1
    direct = np.exp(-2j * np.pi * np.arange(n // 2 + 1)[:, None] * np.arange(n) / n)
    expected = []
    for start in starts:
        row = x[start : start + n]
        residual = row - math.fsum(row) / n
        expected.append((direct @ (residual * w)) * np.sqrt(d / (fs * math.fsum(w * w))))
    expected = np.asarray(expected)
    expected_psd = np.mean(np.abs(expected) ** 2, axis=0)
    result = observations().frame_spectra(
        x, fs, nperseg=n, noverlap=n - hop, window=window, detrend="constant"
    )
    np.testing.assert_array_equal(result["segment_start_sample"], starts)
    np.testing.assert_allclose(result["stft"], expected, rtol=3e-12, atol=3e-15)
    np.testing.assert_allclose(result["psd"], expected_psd, rtol=3e-12, atol=1e-18)


@pytest.mark.parametrize("exponent", [-400, -3, 400])
@pytest.mark.parametrize("window", ["boxcar", "hann"])
def test_one_ulp_nonconstant_detrend_is_not_silenced(exponent, window):
    # Exact binary residuals. An approximate-equality constant test would erase them.
    base = float(2.0**exponent)
    step = np.spacing(base)
    residual = np.array([0.0, 1.0, -1.0, 0.0])
    x = base + step * residual
    w = np.ones(4) if window == "boxcar" else np.array([0.0, 0.5, 1.0, 0.5])
    direct = np.exp(-2j * np.pi * np.arange(3)[:, None] * np.arange(4) / 4)
    expected_unit = (direct @ (residual * w)) * np.sqrt(np.array([1, 2, 1]) / (44100 * sum(w * w)))
    actual = observations().frame_spectra(
        x, 44100, nperseg=4, noverlap=0, window=window, detrend="constant"
    )
    assert np.any(actual["psd"] > 0)
    np.testing.assert_allclose(actual["stft"][0] / step, expected_unit, rtol=1e-13, atol=1e-17)
    np.testing.assert_allclose(
        actual["psd"] / step / step, np.abs(expected_unit) ** 2, rtol=1e-13, atol=1e-30
    )


def test_constant_rows_do_not_erase_a_separate_weak_nonconstant_row():
    n = 31
    weak = 1e-10 * np.where(np.arange(n) % 2, -1.0, 1.0)
    x = np.concatenate([np.full(n, 1e153), np.full(n, 1e-154), weak])
    params = dict(nperseg=n, noverlap=0, window="boxcar", detrend="constant")
    separate = observations().frame_spectra(weak, 44100, **params)
    together = observations().frame_spectra(x, 44100, **params)
    assert np.all(together["stft"][:2] == 0)
    np.testing.assert_array_equal(together["stft"][2], separate["stft"][0])
    np.testing.assert_allclose(together["psd"], separate["psd"] / 3, rtol=1e-14, atol=0)
    assert together["psd"][1] > 0


@pytest.mark.parametrize("scale", [1e200, 1e-200])
def test_constant_detrend_does_not_bypass_raw_square_guard(scale):
    with pytest.raises(ValueError, match="representable"):
        observations().frame_spectra(np.full(3, scale), 44100, nperseg=3, detrend="constant")


@pytest.mark.parametrize("scale", [1e200, 1e-200])
def test_unrepresentable_squares_rejected_without_warnings(scale):
    with pytest.raises(ValueError, match="representable"):
        observations().frame_spectra(np.full(16, scale), 100, nperseg=16)


def test_large_representable_spectrum_and_silent_level(tmp_path):
    module = observations()
    large = module.frame_spectra(np.full(16, 1e153), 8000, nperseg=16)
    assert large["psd"][0] / 1e306 == pytest.approx(2 * 16 / (3 * 8000))
    silent = module.analyse_chunks(
        [chunk(tmp_path, np.zeros(16384), calibration=Calibration(-120, 0, 1))]
    )
    assert silent["valid"].tolist() == [True]
    assert np.all(silent["psd"] == 0)
    assert np.isneginf(silent["band_spl_db_re_1upa"]).all()


@pytest.mark.parametrize(
    "kwargs",
    [
        {"max_chunk_bytes": 1},
        {"max_chunk_bytes": True},
        {"max_decode_samples": 0},
        {"max_frame_samples": 100},
        {"max_output_cells": 10},
        {"max_windows": 1, "frame_seconds": 1},
        {"frame_seconds": 0},
        {"frame_seconds": float("nan")},
        {"frame_seconds": 2.00001},
        {"nperseg": 32768},
        {"bands_hz": {"bad": (100, 100)}},
        {"bands_hz": {"bad": (-1, 100)}},
        {"bands_hz": {"bad": (0, 5000)}},
        {"bands_hz": {}},
        {"bands_hz": {"empty": (0.1, 0.2)}},
        {"active_intervals": [(EPOCH, EPOCH)]},
    ],
)
def test_invalid_analysis_configuration_fails_closed(tmp_path, kwargs):
    recording = chunk(tmp_path, np.zeros(16384))
    with pytest.raises(ValueError):
        observations().analyse_chunks([recording], **kwargs)


@pytest.mark.parametrize("bad_samples", [[], [float("nan"), 0], [1.1, 0]])
def test_invalid_wav_samples_fail_even_in_partial_window(tmp_path, bad_samples):
    recording = chunk(tmp_path, np.asarray(bad_samples))
    with pytest.raises(ValueError):
        observations().analyse_chunks([recording])


@pytest.mark.parametrize("channel", [-1, 1, True, 0.1])
def test_invalid_channel(tmp_path, channel):
    with pytest.raises(ValueError):
        observations().analyse_chunks([chunk(tmp_path, np.zeros(16384), channel=channel)])


def test_utc_clock_and_mixed_units_are_not_guessed(tmp_path):
    module = observations()
    a = chunk(tmp_path, np.zeros(8192), name="a.wav")
    b = chunk(tmp_path, np.zeros(8192), start=1.00001, name="b.wav")
    with pytest.raises(ValueError, match="sample clock"):
        module.analyse_chunks([a, b])
    with pytest.raises(ValueError, match="UTC"):
        module.analyse_chunks(
            [module.RecordingChunk(a.path, datetime(2026, 1, 1), evidence_kind="synthetic")]
        )
    c = chunk(tmp_path, np.zeros(8192), start=1, name="c.wav", calibration=Calibration(-120, 0, 1))
    with pytest.raises(ValueError, match="calibration"):
        module.analyse_chunks([a, c])
    d = chunk(tmp_path, np.zeros(8000), fs=8000, start=1, name="d.wav")
    with pytest.raises(ValueError, match="sample rate"):
        module.analyse_chunks([a, d])
    with pytest.raises(ValueError):
        module.analyse_chunks([])


def test_hash_pin_and_bounded_decoder(tmp_path, monkeypatch):
    module = observations()
    a = chunk(tmp_path, np.zeros((16384, 2)))
    pinned = module.RecordingChunk(
        a.path, EPOCH, evidence_kind="synthetic", expected_sha256="0" * 64
    )
    with pytest.raises(ValueError, match="SHA256"):
        module.analyse_chunks([pinned])
    original = sf.SoundFile.read
    sizes = []

    def read(self, frames=-1, *args, **kwargs):
        sizes.append(frames * self.channels)
        assert 0 < frames * self.channels <= 100
        return original(self, frames, *args, **kwargs)

    monkeypatch.setattr(sf.SoundFile, "read", read)
    result = module.analyse_chunks([a], max_decode_samples=100)
    assert result["valid"].all()
    assert sizes


def test_export_roundtrip_determinism_nulls_and_no_overwrite(tmp_path):
    module = observations()
    assert callable(getattr(module, "write_observations", None)), (
        "protected observation exporter missing"
    )
    a = chunk(tmp_path, np.zeros(20000), calibration=Calibration(-120, 0, 1))
    result = module.analyse_chunks([a])
    output = tmp_path / "products"
    paths = module.write_observations(result, output)
    bodies = {p.name: p.read_bytes() for p in paths}
    assert set(bodies) == {"observations.npz", "windows.csv", "summary.json"}
    with np.load(output / "observations.npz", allow_pickle=False) as archive:
        for key, value in result.items():
            if isinstance(value, np.ndarray):
                np.testing.assert_array_equal(archive[key], value)
    summary = json.loads(
        (output / "summary.json").read_text(), parse_constant=lambda x: pytest.fail(x)
    )
    assert summary["band_spl_db_re_1upa"] == [[None], [None]]
    assert summary["band_power"] == [[0.0], [None]]
    assert summary["valid"] == [True, False]
    assert "negative infinity" in summary["serialization"]["nonfinite"]
    with (output / "windows.csv").open() as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 2
    assert rows[0]["band_power"] == "0"
    assert rows[0]["band_spl_db_re_1upa"] == "-inf"
    assert rows[1]["band_power"] == ""
    module.write_observations(result, output)
    other = tmp_path / "replay"
    module.write_observations(result, other)
    assert {p.name: p.read_bytes() for p in other.iterdir()} == bodies
    changed = dict(result, unit="changed")
    with pytest.raises(FileExistsError):
        module.write_observations(changed, output)
    assert {p.name: p.read_bytes() for p in output.iterdir()} == bodies


def test_export_rejects_symlink_destination_without_outside_writes(tmp_path):
    module = observations()
    result = module.analyse_chunks([chunk(tmp_path, np.zeros(16384))])
    outside = tmp_path / "outside"
    outside.mkdir()
    link = tmp_path / "link"
    link.symlink_to(outside, target_is_directory=True)
    with pytest.raises(ValueError, match="symlink"):
        module.write_observations(result, link / "products")
    assert list(outside.iterdir()) == []
    output = tmp_path / "existing"
    output.mkdir()
    (output / "summary.json").symlink_to(tmp_path / "record.wav")
    before = (tmp_path / "record.wav").read_bytes()
    with pytest.raises(ValueError, match="symlink"):
        module.write_observations(result, output)
    assert (tmp_path / "record.wav").read_bytes() == before


def test_anchor_phase_and_half_open_band_and_activity_boundaries(tmp_path):
    module = observations()
    fs = 8192
    x = 0.25 * np.cos(2 * np.pi * 512 * np.arange(4 * fs) / fs)
    record = chunk(tmp_path, x, start=0.5)
    result = module.analyse_chunks(
        [record],
        anchor_utc=EPOCH,
        nperseg=2048,
        window="boxcar",
        noverlap=0,
        bands_hz={"below": (0, 512), "above": (512, 4096)},
        active_intervals=[(EPOCH, EPOCH + timedelta(seconds=2))],
    )
    assert result["valid"].tolist() == [False, True, False]
    assert result["window_start_s"].tolist() == [0, 2, 4]
    assert result["active_transmission"].tolist() == [True, False, False]
    assert result["band_power"][1, 0] < 1e-25
    assert result["band_power"][1, 1] == pytest.approx(0.25**2 / 2)
    assert result["band_level_db"][1, 1] == pytest.approx(10 * np.log10(0.25**2 / 2))


def test_source_change_after_inspection_fails_hash_gate(tmp_path, monkeypatch):
    module = observations()
    record = chunk(tmp_path, np.zeros(16384))
    original = module._snapshot
    calls = 0

    def snapshot(path, limit, expected=None):
        nonlocal calls
        calls += 1
        if calls == 2:
            sf.write(path, np.ones(16384), 8192, subtype="DOUBLE")
        return original(path, limit, expected)

    monkeypatch.setattr(module, "_snapshot", snapshot)
    with pytest.raises(ValueError, match="SHA256"):
        module.analyse_chunks([record])


@pytest.mark.parametrize("fs", [True, 0, float("nan"), float("inf")])
def test_invalid_shared_kernel_sample_rate(fs):
    with pytest.raises(ValueError):
        observations().frame_spectra(np.ones(16), fs, nperseg=16)


def test_high_precision_utc_is_not_silently_truncated(tmp_path):
    module = observations()
    record = chunk(tmp_path, np.zeros(16384))
    imprecise = module.RecordingChunk(
        record.path, "2026-01-01T00:00:00.0000009Z", evidence_kind="synthetic"
    )
    with pytest.raises(ValueError, match="microsecond"):
        module.analyse_chunks([imprecise])


def test_shared_kernel_rejects_excess_spectral_work_before_allocation():
    with pytest.raises(ValueError, match="bound"):
        observations().frame_spectra(np.zeros(10000), 8192, nperseg=8192, noverlap=8191)


@pytest.mark.parametrize("fs", [1e300, 1e308])
@pytest.mark.parametrize("sign", [-1, 1])
def test_normalization_rejects_lost_nonconstant_pair_power(fs, sign):
    from decimal import Decimal, localcontext

    x = sign * np.array([1e-154, np.nextafter(1e-154, np.inf)])
    assert np.all(x * x > 0)
    with localcontext() as context:
        context.prec = 110
        # Nyquist ignores the mean: exact difference of the stored samples.
        delta = Decimal.from_float(float(x[0])) - Decimal.from_float(float(x[1]))
        power = delta**2 / (2 * Decimal.from_float(fs))
    assert power > 0 and float(power) == 0
    with pytest.raises(ValueError, match="representable"):
        observations().frame_spectra(
            x, fs, nperseg=2, noverlap=0, window="boxcar", detrend="constant"
        )


@pytest.mark.parametrize("fs", [1e300, 1e308])
@pytest.mark.parametrize("n", [8192, 16384])
def test_normalization_rejects_lost_window_supported_impulse_power(fs, n):
    from decimal import Decimal, localcontext
    import math

    x = np.zeros(n)
    x[1] = 2e-162
    assert x[1] * x[1] > 0
    with localcontext() as context:
        context.prec = 110
        weight = Decimal.from_float(math.sin(math.pi / n) ** 2)
        # Independent DC power; analytic energy of the periodic Hann.
        power = (Decimal.from_float(float(x[1])) * weight) ** 2 / (
            Decimal.from_float(fs) * 3 * n / 8
        )
    assert power > 0 and float(power) == 0
    with pytest.raises(ValueError, match="representable"):
        observations().frame_spectra(x, fs, nperseg=n, noverlap=0, window="hann")


@pytest.mark.parametrize("kind", ["pair", "impulse"])
def test_normalization_accepts_representable_mean_with_lost_weak_segment(kind):
    from decimal import Decimal, localcontext

    fs = 1e308
    if kind == "pair":
        strong = np.array([1e153, np.nextafter(1e153, np.inf)])
        weak = np.array([1e-154, np.nextafter(1e-154, np.inf)])
        params = dict(nperseg=2, noverlap=0, window="boxcar", detrend="constant")
        # Rounded subtraction leaves one nonzero sample in each pair; both bins
        # have the same magnitude. Mean power includes two segments, not one.
        assert np.count_nonzero(strong - strong.mean()) == 1
        with localcontext() as context:
            context.prec = 110
            delta = Decimal.from_float(float(strong[1])) - Decimal.from_float(float(strong[0]))
            expected = float(delta**2 / (4 * Decimal.from_float(fs)))
        expected_power = np.full(2, expected)
    else:
        n = 8192
        strong, weak = np.zeros(n), np.zeros(n)
        strong[n // 2], weak[1] = 1e154, 2e-162
        params = dict(nperseg=n, noverlap=0, window="hann", detrend=False)
        with localcontext() as context:
            context.prec = 110
            expected = float(
                Decimal.from_float(1e154) ** 2 / (2 * Decimal.from_float(fs) * 3 * n / 8)
            )
        expected_power = np.full(n // 2 + 1, 2 * expected)
        expected_power[[0, -1]] = expected
    result = observations().frame_spectra(np.r_[strong, weak], fs, **params)
    assert np.all(result["stft"][1] == 0)
    assert np.all(result["psd"] > 0)
    np.testing.assert_allclose(result["psd"], expected_power, rtol=2e-14, atol=0)
    np.testing.assert_allclose(
        result["psd"], np.abs(result["stft"][0]) ** 2 / 2, rtol=2e-14, atol=0
    )


@pytest.mark.parametrize("fs", [np.nextafter(0.0, 1.0), 44100, 1e308])
def test_normalization_keeps_window_annihilated_input_silent(fs):
    x = np.zeros(8192)
    x[0] = 2e-162  # Hann is exactly zero here: nonzero input is not spectral evidence.
    result = observations().frame_spectra(x, fs, nperseg=len(x), noverlap=0, window="hann")
    assert np.all(result["stft"] == 0) and np.all(result["psd"] == 0)


def test_normalization_does_not_invent_content_below_fft_accuracy():
    # In exact arithmetic, DC contains the tiny impulse; the float64 FFT
    # cannot resolve it amid the cancelling strong values. This is not evidence
    # of a *normalization* loss and must not be inferred from nonzero samples.
    x = np.array([1e153, -1e153, 2e-162, 0.0])
    result = observations().frame_spectra(x, 1e308, nperseg=4, noverlap=0, window="boxcar")
    assert result["psd"][0] == 0
    assert result["stft"][0, 0] == 0
    np.testing.assert_allclose(result["psd"][1:], [0.01, 0.01], rtol=2e-14, atol=0)
