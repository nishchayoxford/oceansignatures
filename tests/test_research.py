"""Independent full-chain checks for the public synthetic research export."""

import importlib
import json
from pathlib import Path

import numpy as np
import pytest
import soundfile as sf
from scipy.signal import welch


def research():
    return importlib.import_module("oceansignatures.research")


def test_full_chain_recovers_gain_and_retains_excluded_windows(tmp_path):
    module = research()
    output = tmp_path / "run"
    report = module.export_research(output, seed=1729, gain_db=6.0)
    assert report["kind"] == "predictive approximations + synthetic verification"
    arrays = np.load(output / "research.npz", allow_pickle=False)
    assert arrays["observed_psd_pa2_hz"].shape == (120, 2049)
    assert np.array_equal(arrays["window_start_s"], np.arange(120) * 2)
    assert np.flatnonzero(~arrays["valid"]).tolist() == [40]
    assert np.flatnonzero(arrays["active_transmission"]).tolist() == [75, 76]
    assert np.flatnonzero(~arrays["qc_ok"]).tolist() == [40, 75, 76]
    assert np.isnan(arrays["observed_psd_pa2_hz"][40]).all()
    assert np.isfinite(arrays["observed_psd_pa2_hz"][[75, 76]]).all()
    assert np.isnan(arrays["matched_residual_db"][[40, 75, 76]]).all()
    np.testing.assert_allclose(arrays["matched_residual_db"][arrays["qc_ok"]], 6, atol=2e-12)
    summary = report["numericalsummary"]
    assert summary["compared_windows"] == 117
    assert summary["excluded_windows"] == 3
    assert summary["gain_bias_db"] == pytest.approx(6, abs=1e-12)
    assert summary["observation_output_cells"] <= 2_000_000
    # Independently decode the first 2 seconds of channel one, apply the
    # declared ADC algebra and scipy's Welch (not the production frame kernel).
    metadata = json.loads((output / "metadata.json").read_text())
    first = metadata["observed"]["chunks"][0]
    digital, fs = sf.read(output / first["path"], always_2d=True)
    assert fs == 8192 and digital.shape[1] == 2 and first["channel"] == 1
    cal = first["calibration"]
    scale = (
        cal["full_scale_voltage"]
        * 10 ** (-(cal["gain_db"] + cal["sensitivity_db_v_per_upa"]) / 20)
        * 1e-6
    )
    f, expected = welch(
        digital[:16384, 1] * scale, fs=fs, nperseg=4096, noverlap=2048, detrend=False, window="hann"
    )
    np.testing.assert_array_equal(arrays["frequency_hz"], f)
    np.testing.assert_allclose(arrays["observed_psd_pa2_hz"][0], expected, rtol=2e-12, atol=1e-22)
    np.testing.assert_allclose(
        arrays["observed_psd_pa2_hz"], arrays["digital_psd"] / 1 * scale**2, rtol=2e-12, atol=1e-22
    )
    assert metadata["digital_control"]["psd_units"] == "digital^2/Hz"
    assert metadata["digital_control"]["chunks"][0]["calibration"] is None
    assert len(list(output.glob("wav/**/*.wav"))) > 2
    assert (output / "windows.csv").is_file()
    assert str(output) not in (output / "metadata.json").read_text()


def test_adapter_requires_calibration_and_preserves_exact_support():
    module = research()
    from oceansignatures.comparison import compare_prediction_observation

    analysis = dict(
        unit="Pa",
        calibration_status="declared scalar calibration",
        metadata=dict(psd_units="Pa^2/Hz", evidence_kind="synthetic", anchor_utc=module.ORIGIN),
        frequency_hz=np.array([0.0, 100.0, 102.0, 4000.0, 4096.0]),
        window_start_s=np.array([0.0, 2.0]),
        window_stop_s=np.array([2.0, 4.0]),
        psd=np.ones((2, 5)),
        valid=np.array([True, True]),
        qc_ok=np.array([True, False]),
    )
    batch = module.comparison_batch(analysis)
    assert batch["evidence_kind"] == "synthetic_observation"
    assert batch["psd_unit"] == "Pa2/Hz"
    np.testing.assert_array_equal(batch["valid_window"], [True, False])
    np.testing.assert_array_equal(batch["frequency_hz"], [100.0, 102.0, 4000.0])
    predicted = module.comparison_batch(analysis, as_prediction=True)
    result = compare_prediction_observation(predicted, batch)
    assert result["summary"]["excluded_windows"] == 1
    for key, changed in (
        ("unit", "normalized digital amplitude"),
        ("calibration_status", "unknown"),
    ):
        with pytest.raises(ValueError, match="physical|calibration"):
            module.comparison_batch({**analysis, key: changed})
    with pytest.raises(ValueError, match="physical|units"):
        module.comparison_batch(
            {**analysis, "metadata": {**analysis["metadata"], "psd_units": "digital^2/Hz"}}
        )
    for key in ("frequency_hz", "window_start_s", "window_stop_s"):
        with pytest.raises(ValueError, match="exactly aligned"):
            compare_prediction_observation(predicted, {**batch, key: batch[key] + 0.001})
    with pytest.raises(ValueError, match="UTC origin"):
        compare_prediction_observation(
            predicted, {**batch, "time_origin_utc": "2000-01-02T00:00:00Z"}
        )


def test_expected_welch_matches_independent_finite_fourier_oracle(tmp_path):
    module = research()
    # A finite circular Gaussian ensemble with power only in one Fourier bin.
    # Direct covariance/DFT calculation, independent of production convolution.
    fs, n, length = 16, 16, 8
    one_sided = np.zeros(9)
    one_sided[3] = 2.5
    w = 0.5 - 0.5 * np.cos(2 * np.pi * np.arange(length) / length)
    lag = np.arange(length)[:, None] - np.arange(length)[None, :]
    covariance = 2.5 * (fs / n) * np.cos(2 * np.pi * 3 * lag / n)
    expected = []
    for k in range(length // 2 + 1):
        v = w * np.exp(-2j * np.pi * k * np.arange(length) / length)
        expected.append(
            (v @ covariance @ v.conj()).real
            / (fs * np.sum(w * w))
            * (1 if k in (0, length // 2) else 2)
        )
    actual = module.expected_welch_psd(one_sided, fs=fs, nperseg=length)
    np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-15)
    for bad in (one_sided - 1, one_sided * np.nan):
        with pytest.raises(ValueError):
            module.expected_welch_psd(bad, fs=fs, nperseg=length)
    output = tmp_path / "expected"
    report = module.export_research(output)
    arrays = np.load(output / "research.npz")
    assert arrays["expected_welch_psd_pa2_hz"].shape == (120, 2049)
    good = arrays["qc_ok"]
    mask = (arrays["frequency_hz"] >= 100) & (arrays["frequency_hz"] < 4000)
    # Broad-band finite-realization error is measured, never forced away by RMS scaling.
    expected_power = arrays["expected_welch_psd_pa2_hz"][:, mask].sum(axis=1) * 2
    actual_power = arrays["baseline_psd_pa2_hz"][:, mask].sum(axis=1) * 2
    err = 10 * np.log10(actual_power[good] / expected_power[good])
    assert np.sqrt(np.mean(err**2)) < 0.5
    assert report["numericalsummary"]["sampling_broadband_rmse_db"] == pytest.approx(
        np.sqrt(np.mean(err**2))
    )
    assert not np.array_equal(
        arrays["raw_expected_psd_pa2_hz"], arrays["expected_welch_psd_pa2_hz"]
    )
    assert report["numericalsummary"]["welch_smoothing_max_abs_db"] > 0.1
    assert arrays["source_psd_pa2_m2_hz"].shape == arrays["comparison_frequency_hz"].shape
    assert np.nanmax(np.abs(arrays["matched_band_residual_db"] - 6)) < 1e-11


@pytest.mark.parametrize("scale", [1.0, 1e307, 1e308, np.nextafter(0.0, 1.0)])
def test_expected_welch_uniform_float64_boundary(scale):
    from decimal import Decimal, localcontext

    # Direct circular-covariance / eight-sample Hann DFT, not FFT convolution:
    # a uniform one-sided density has exact ratios 2/3, 13/12, 1, 13/12, 2/3.
    # Multiply in Decimal before rounding, so halving a subnormal cannot erase it.
    with localcontext() as context:
        context.prec = 100
        ratios = [Decimal(2) / 3, Decimal(13) / 12, Decimal(1), Decimal(13) / 12, Decimal(2) / 3]
        expected = np.array([float(Decimal.from_float(float(scale)) * r) for r in ratios])
    actual = research().expected_welch_psd(np.full(9, scale), fs=16, nperseg=8)
    assert np.isfinite(actual).all()
    assert np.all(actual > 0)
    np.testing.assert_allclose(actual / scale, expected / scale, rtol=1e-13, atol=0)


def test_expected_welch_boundary_rows_are_independent_and_silence_is_zero():
    scales = np.array([0.0, np.nextafter(0.0, 1.0), 1.0, 1e308])
    result = research().expected_welch_psd(np.repeat(scales[:, None], 9, axis=1), nperseg=8)
    np.testing.assert_array_equal(result[0], 0)
    for i in range(1, len(scales)):
        assert np.all(result[i] > 0)
        np.testing.assert_allclose(
            result[i] / scales[i],
            np.ones(5) if i == 1 else [2 / 3, 13 / 12, 1, 13 / 12, 2 / 3],
            rtol=1e-13,
            atol=0,
        )


@pytest.mark.parametrize("kind", ["overflow", "underflow"])
def test_expected_welch_unrepresentable_output_rejected(kind):
    if kind == "overflow":
        # 13/12 times float64 max is genuinely beyond the output range.
        density = np.full(9, np.finfo(float).max)
    else:
        # A DC-only ensemble gives DC expectation p[0]/3, below half a subnormal.
        density = np.zeros(9)
        density[0] = np.nextafter(0.0, 1.0)
    with pytest.raises(ValueError, match="representable"):
        research().expected_welch_psd(density, fs=16, nperseg=8)


def test_export_figures_and_one_at_a_time_sensitivity_are_computed(tmp_path):
    output = tmp_path / "figures"
    report = research().export_research(output)
    assert len(report["figures"]) == 6
    assert {f["path"] for f in report["figures"]} == {
        "figures/spectra.png",
        "figures/passage.png",
        "figures/propagation.png",
        "figures/comparison.png",
        "figures/residuals.png",
        "figures/sensitivity.png",
    }
    for item in report["figures"]:
        assert (output / item["path"]).read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
        assert item["caption"] and item["title"]
    cases = json.loads((output / "sensitivity.json").read_text())["cases"]
    sensitivity = np.load(output / "sensitivity.npz")
    assert len(cases) == 19
    assert all(len(c["changes"]) == 1 for c in cases[1:])
    assert {next(iter(c["changes"])) for c in cases[1:]} == {
        "source_offset_db",
        "wind_offset_m_s",
        "speed_kn",
        "closest_approach_m",
        "receiver_depth_m",
        "temperature_c",
        "salinity_ppt",
        "depth_m",
        "ph",
    }
    data = np.load(output / "research.npz")
    f = data["comparison_frequency_hz"]
    band = (f >= 100) & (f < 500)
    for offset in (-6, 6):
        index = next(i for i, c in enumerate(cases) if c["changes"] == {"source_offset_db": offset})
        expected = (
            data["vessel_psd_pa2_hz"][:, band] * 10 ** (offset / 10)
            + data["wind_psd_pa2_hz"][:, band]
        ).sum(axis=1) * 2
        np.testing.assert_allclose(sensitivity["band_power_pa2"][index, :, 0], expected, rtol=2e-14)
    # Travel-delay and retarded geometry oracle, separate from predictor internals.
    speed = 12 * 1852 / 3600
    distance = np.sqrt((speed * (data["emission_time_s"] - 120)) ** 2 + 250**2 + (50 - 8) ** 2)
    np.testing.assert_allclose(data["travel_time_s"], distance / 1500, rtol=1e-12)
    np.testing.assert_allclose(
        data["receiver_time_s"], data["emission_time_s"] + distance / 1500, atol=1e-12
    )
    assert report["numericalsummary"]["quasi_stationary_edge_max_abs_db"] < 0.3
    assert report["configuration"]["water"] == {
        "temperature_c": 10.0,
        "salinity_ppt": 35.0,
        "depth_m": 50.0,
        "ph": 8.0,
    }
    assert report["source_refs"]


def test_export_is_root_independent_and_refuses_differing_replay(tmp_path):
    import hashlib

    module = research()
    first, other = tmp_path / "first", tmp_path / "elsewhere" / "other"
    report = module.export_research(first)
    assert module.export_research(first) == report
    assert module.export_research(other) == report

    def inventory(root):
        return {
            p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in root.rglob("*")
            if p.is_file()
        }

    original = inventory(first)
    assert original == inventory(other)
    manifest = json.loads((first / "manifest.json").read_text())
    assert manifest["files"] == {k: v for k, v in original.items() if k != "manifest.json"}
    for item in report["downloads"] + report["figures"]:
        assert not Path(item["path"]).is_absolute() and ".." not in Path(item["path"]).parts
        assert (first / item["path"]).is_file()
    for path in first.rglob("*.json"):
        assert str(tmp_path) not in path.read_text()
        assert "/home/" not in path.read_text()
    with pytest.raises(FileExistsError):
        module.export_research(first, seed=1234)
    assert inventory(first) == original
    (first / "unrecognized.txt").write_text("do not replace")
    with pytest.raises(FileExistsError):
        module.export_research(first)
    assert (first / "unrecognized.txt").read_text() == "do not replace"
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(FileExistsError):
        module.export_research(empty)
    linked = tmp_path / "linked"
    linked.symlink_to(other, target_is_directory=True)
    with pytest.raises(ValueError, match="symlink"):
        module.export_research(linked)
    with pytest.raises(ValueError):
        module.export_research(tmp_path / "parent" / ".." / "bad")
    for seed in (-1, True, 1.5):
        with pytest.raises(ValueError):
            module.export_research(tmp_path / "invalid", seed=seed)
    for gain in (float("nan"), float("inf"), True, 1000):
        with pytest.raises(ValueError):
            module.export_research(tmp_path / "invalid", gain_db=gain)


def test_public_cli_executes_and_reports_real_counts(tmp_path):
    import hashlib
    import subprocess
    import sys

    root = Path(__file__).resolve().parents[1]
    output = tmp_path / "cli"
    command = [
        sys.executable,
        str(root / "scripts" / "build_research.py"),
        "--output",
        str(output),
        "--seed",
        "1729",
        "--gain-db",
        "6",
    ]
    result = subprocess.run(command, cwd=root, text=True, capture_output=True)
    assert result.returncode == 0, result.stderr
    report = json.loads(result.stdout)
    assert report["files_total"] == len([p for p in output.rglob("*") if p.is_file()])
    assert report["figures"] == 6
    assert (
        report["manifest_sha256"]
        == hashlib.sha256((output / "manifest.json").read_bytes()).hexdigest()
    )
    assert report["numericalsummary"]["gain_bias_db"] == pytest.approx(6)
    replay = subprocess.run(command, cwd=root, text=True, capture_output=True)
    assert replay.returncode == 0, replay.stderr
    assert replay.stdout == result.stdout
    refused = subprocess.run(command[:-1] + ["5"], cwd=root, text=True, capture_output=True)
    assert refused.returncode != 0 and "Refusing differing" in refused.stderr


def test_exports_keep_band_order_and_missing_plot_support(tmp_path, monkeypatch):
    import csv
    from matplotlib.axes import Axes

    captured = []
    original = Axes.plot

    def record_plot(self, *args, **kwargs):
        if kwargs.get("label") == "Matched band residual − injected gain":
            captured.append(np.asarray(args[1]))
        return original(self, *args, **kwargs)

    monkeypatch.setattr(Axes, "plot", record_plot)
    output = tmp_path / "quality"
    report = research().export_research(output)
    assert len(captured) == 1
    assert np.isnan(captured[0][[40, 75, 76]]).all()
    arrays = np.load(output / "research.npz")
    assert arrays["band_names"].tolist() == list(research().BANDS)
    assert report["configuration"]["band_names"] == list(research().BANDS)
    with (output / "windows.csv").open() as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 360
    gap = next(r for r in rows if r["start_s"] == "80.0")
    active = next(r for r in rows if r["start_s"] == "150.0")
    assert gap["observed_pa2"] == "" and float(gap["expected_pa2"]) > 0
    assert active["qc_ok"] == "False" and float(active["observed_pa2"]) > 0
    assert active["observed_minus_expected_db"] == ""
    for i, row in enumerate(rows[:3]):
        assert float(row["observed_pa2"]) == pytest.approx(arrays["observed_band_power_pa2"][0, i])
        assert float(row["expected_pa2"]) == pytest.approx(arrays["expected_band_power_pa2"][0, i])
    assert not arrays["near_full_scale_samples"].any()
