"""Offline tests: generated WAV fixtures are unit-test signals, not reference data."""

from __future__ import annotations

import hashlib
import importlib
import importlib.util
import warnings
from pathlib import Path

import numpy as np
import pytest
import soundfile as sf


def reference_module():
    try:
        return importlib.import_module("oceansignatures.reference")
    except ModuleNotFoundError as exc:
        pytest.fail(f"reference loader feature is missing: {exc}")


def wav_fixture(tmp_path, samples=None, fs=8000):
    if samples is None:
        samples = np.array([[0.0, 0.25], [-0.5, 0.5], [0.125, -0.125]])
    path = tmp_path / "test.wav"
    sf.write(path, samples, fs, subtype="PCM_16")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    return path, digest


def test_loader_preserves_channels_and_checks_pcm_metadata(tmp_path):
    path, digest = wav_fixture(tmp_path)
    result = reference_module().load_reference(path, expected_sha256=digest)
    np.testing.assert_array_equal(result["samples"], [[0, 0.25], [-0.5, 0.5], [0.125, -0.125]])
    assert result["fs"] == 8000
    metadata = result["metadata"]
    assert metadata["frames"] == 3
    assert metadata["channels"] == 2
    assert metadata["sample_count"] == 6
    assert metadata["duration_s"] == 3 / 8000
    assert metadata["sha256"] == digest
    assert metadata["magic"] == "RIFF/WAVE"
    assert metadata["subtype"] == "PCM_16"
    assert metadata["amplitude_units"] == "normalized digital units"
    assert metadata["calibration_status"] == "uncalibrated; conversion to Pa unavailable"


@pytest.mark.parametrize("failure", ["hash", "magic", "size", "truncated", "nonfinite", "empty"])
def test_loader_rejects_untrusted_or_invalid_wav(tmp_path, failure):
    path, digest = wav_fixture(tmp_path)
    limit = 20_000_000
    match = ""
    if failure == "hash":
        digest, match = "0" * 64, "SHA256"
    elif failure == "magic":
        sf.write(path, np.zeros(10), 8000, format="FLAC")
        match = "RIFF/WAVE"
    elif failure == "size":
        limit, match = 16, "size limit"
    elif failure == "truncated":
        path.write_bytes(path.read_bytes()[:-2])
        match = "truncated"
    elif failure == "nonfinite":
        sf.write(path, np.array([0, np.nan, np.inf]), 8000, subtype="FLOAT")
        match = "non-finite"
    else:
        sf.write(path, np.zeros(0), 8000)
        match = "empty"
    if failure != "hash":
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
    with pytest.raises(ValueError, match=match):
        reference_module().load_reference(path, expected_sha256=digest, max_bytes=limit)


def test_analyser_returns_relative_psd_spectrogram_and_qc():
    module = reference_module()
    assert hasattr(module, "analyse_reference"), "reference analyser feature is missing"
    fs = 1024
    t = np.arange(4096) / fs
    selected = 0.125 + 0.5 * np.sin(2 * np.pi * 128 * t)
    recording = {
        "samples": np.column_stack([np.zeros(len(t)), selected]),
        "fs": fs,
        "metadata": {"amplitude_units": "normalized digital units"},
    }
    result = module.analyse_reference(recording, channel=1, nperseg=256, spectrogram_nperseg=128)
    np.testing.assert_array_equal(result["time_s"], t)
    np.testing.assert_array_equal(result["waveform"], selected)
    f, psd = result["frequency_hz"], result["psd"]
    assert f[np.argmax(psd)] == 128
    assert np.sum(psd) * (f[1] - f[0]) == pytest.approx(0.5**2 / 2)
    spec = result["spectrogram_psd"]
    assert spec.shape == (
        len(result["spectrogram_frequency_hz"]),
        len(result["spectrogram_time_s"]),
    )
    np.testing.assert_allclose(np.sum(spec, axis=0) * fs / 128, 0.5**2 / 2)
    assert result["metadata"]["psd_units"] == "digital units squared per Hz"
    assert result["metadata"]["selected_channel"] == 1
    assert result["metadata"]["welch"]["detrend"] == "constant"
    assert result["qc"]["mean"] == pytest.approx(0.125)
    assert result["qc"]["rms"] == pytest.approx(np.sqrt(0.125**2 + 0.5**2 / 2))
    assert result["qc"]["peak_absolute"] == pytest.approx(0.625)
    assert result["qc"]["near_full_scale_fraction"] == 0


@pytest.mark.parametrize(
    "overrides",
    [
        {"samples": np.zeros(4)},
        {"samples": np.zeros((0, 1))},
        {"samples": np.zeros((4, 0))},
        {"samples": np.array([[0.0], [np.nan]])},
        {"fs": 0},
        {"fs": np.inf},
        {"channel": -1},
        {"channel": 1},
        {"channel": 0.5},
        {"nperseg": 1},
        {"nperseg": 2.5},
        {"spectrogram_nperseg": 0},
    ],
)
def test_analyser_rejects_invalid_inputs(overrides):
    recording = {"samples": np.zeros((4, 1)), "fs": 8000}
    kwargs = {}
    for key, value in overrides.items():
        (recording if key in recording else kwargs)[key] = value
    with pytest.raises(ValueError):
        reference_module().analyse_reference(recording, **kwargs)


def test_silent_short_record_has_finite_zero_density():
    result = reference_module().analyse_reference({"samples": np.zeros((3, 1)), "fs": 100})
    assert result["qc"]["silent"] is True
    assert result["qc"]["rms"] == 0
    assert np.all(result["psd"] == 0)
    assert np.all(result["spectrogram_psd"] == 0)
    assert result["metadata"]["welch"]["nperseg"] == 3


@pytest.mark.parametrize("amplitude", [1e-200, 1e200])
def test_double_wav_rejects_unrepresentable_square_without_warnings(tmp_path, amplitude):
    samples = np.tile([amplitude, -amplitude], 2048)
    path = tmp_path / "extreme.wav"
    sf.write(path, samples, 8000, subtype="DOUBLE")
    recording = reference_module().load_reference(
        path, expected_sha256=hashlib.sha256(path.read_bytes()).hexdigest()
    )
    np.testing.assert_array_equal(recording["samples"][:, 0], samples)
    with warnings.catch_warnings(record=True) as seen:
        warnings.simplefilter("always")
        with pytest.raises(ValueError, match="representable"):
            reference_module().analyse_reference(recording)
    assert not seen


@pytest.mark.parametrize("amplitude", [0.25, 1e153, 1e-150])
def test_double_wav_representable_scale_retains_waveform_and_density(tmp_path, amplitude):
    samples = np.tile([amplitude, -amplitude], 2048)
    path = tmp_path / "supported.wav"
    sf.write(path, samples, 8, subtype="DOUBLE")
    recording = reference_module().load_reference(
        path, expected_sha256=hashlib.sha256(path.read_bytes()).hexdigest()
    )
    with warnings.catch_warnings(record=True) as seen:
        warnings.simplefilter("always")
        result = reference_module().analyse_reference(recording, nperseg=4, spectrogram_nperseg=4)
    assert not seen
    np.testing.assert_array_equal(result["waveform"], samples)
    assert result["qc"]["finite"] is True
    assert result["qc"]["rms"] == pytest.approx(amplitude, rel=1e-14, abs=0)
    assert result["qc"]["mean"] == 0
    assert not result["qc"]["silent"]
    # Periodic Hann: Nyquist density is a**2 / 3 at fs=8, nperseg=4.
    assert result["psd"][-1] == pytest.approx(amplitude**2 / 3, rel=1e-14, abs=0)
    np.testing.assert_allclose(result["spectrogram_psd"][-1], amplitude**2 / 3, rtol=1e-14, atol=0)
    assert np.isfinite(result["psd"]).all()
    assert np.isfinite(result["spectrogram_psd"]).all()


@pytest.mark.parametrize(
    "amplitude,fs", [(1e153, 1e-10), (1e-155, 1e20), (0.0, 1e-320), (0.25, 1e308)]
)
def test_analyser_rejects_unrepresentable_density_or_axes_without_warnings(amplitude, fs):
    recording = {"samples": np.tile([amplitude, -amplitude], 32)[:, None], "fs": fs}
    with warnings.catch_warnings(record=True) as seen:
        warnings.simplefilter("always")
        with pytest.raises(ValueError, match="representable"):
            reference_module().analyse_reference(recording, nperseg=4, spectrogram_nperseg=4)
    assert not seen


def test_ordinary_scipy_density_convention_is_unchanged():
    from scipy import signal

    waveform = np.random.default_rng(67).normal(0.125, 0.25, 513)
    result = reference_module().analyse_reference(
        {"samples": waveform[:, None], "fs": 8000}, nperseg=128, spectrogram_nperseg=64
    )
    shared = dict(window="hann", detrend="constant", return_onesided=True, scaling="density")
    f, psd = signal.welch(waveform, fs=8000, nperseg=128, noverlap=64, average="mean", **shared)
    sfreq, stime, spec = signal.spectrogram(
        waveform, fs=8000, nperseg=64, noverlap=48, mode="psd", **shared
    )
    np.testing.assert_array_equal(result["waveform"], waveform)
    np.testing.assert_array_equal(result["frequency_hz"], f)
    np.testing.assert_array_equal(result["spectrogram_frequency_hz"], sfreq)
    np.testing.assert_array_equal(result["spectrogram_time_s"], stime)
    np.testing.assert_allclose(result["psd"], psd, rtol=1e-13, atol=0)
    np.testing.assert_allclose(result["spectrogram_psd"], spec, rtol=1e-12, atol=0)
    assert result["qc"]["mean"] == pytest.approx(np.mean(waveform))
    assert result["qc"]["rms"] == pytest.approx(np.sqrt(np.mean(waveform**2)))


def fetch_module():
    path = Path(__file__).resolve().parents[1] / "scripts" / "fetch_reference.py"
    assert path.exists(), "bounded reference fetch feature is missing"
    spec = importlib.util.spec_from_file_location("reference_fetch_test", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_fetch_pins_bytes_and_verifies_decoding_offline(tmp_path, monkeypatch):
    import io

    module = fetch_module()
    path, digest = wav_fixture(tmp_path)
    payload = path.read_bytes()
    monkeypatch.setattr(module, "EXPECTED_SHA256", digest)

    class Response(io.BytesIO):
        headers = {"Content-Length": str(len(payload)), "Content-Type": "audio/x-wav"}

        def geturl(self):
            return module.SOURCE_URL

    def urlopen(request, *, timeout):
        assert request.full_url == module.SOURCE_URL
        assert 0 < timeout <= 20
        return Response(payload)

    monkeypatch.setattr(module.urllib.request, "urlopen", urlopen)
    destination = tmp_path / "download" / "reference.wav"
    result = module.fetch_reference(destination)
    assert destination.read_bytes() == payload
    assert result["sha256"] == digest
    assert result["frames"] == 3
    assert result["channels"] == 2


@pytest.mark.parametrize(
    "failure", ["declared_size", "streamed_size", "hash", "redirect", "deadline"]
)
def test_fetch_rejects_bad_responses_without_replacing_destination(tmp_path, monkeypatch, failure):
    import io

    module = fetch_module()
    path, digest = wav_fixture(tmp_path)
    payload = path.read_bytes()
    monkeypatch.setattr(module, "EXPECTED_SHA256", digest if failure != "hash" else "0" * 64)
    max_bytes = len(payload) - 1 if "size" in failure else 20_000_000

    class Response(io.BytesIO):
        headers = {"Content-Length": str(len(payload))} if failure == "declared_size" else {}

        def geturl(self):
            return (
                "https://not-noaa.example/surprise.wav"
                if failure == "redirect"
                else module.SOURCE_URL
            )

        def read1(self, size=-1):
            assert size >= 0, "unbounded network read"
            assert failure != "declared_size", (
                "oversized Content-Length was not checked before read"
            )
            return super().read1(size)

    monkeypatch.setattr(module.urllib.request, "urlopen", lambda *args, **kwargs: Response(payload))
    if failure == "deadline":
        assert hasattr(module, "time"), "total transfer deadline is missing"
        ticks = iter([0, 61, 62])
        monkeypatch.setattr(module.time, "monotonic", lambda: next(ticks))
    destination = tmp_path / "existing.wav"
    destination.write_bytes(b"preserve existing file")
    with pytest.raises(
        ValueError,
        match={
            "declared_size": "size limit",
            "streamed_size": "size limit",
            "hash": "SHA256",
            "redirect": "redirect",
            "deadline": "deadline",
        }[failure],
    ):
        module.fetch_reference(destination, max_bytes=max_bytes)
    assert destination.read_bytes() == b"preserve existing file"
    assert not list(tmp_path.glob("*.part"))


def test_cli_verify_only_analyses_local_file_without_network(tmp_path, monkeypatch, capsys):
    import json

    module = fetch_module()
    assert hasattr(module, "main"), "fetch and verify CLI feature is missing"
    path, digest = wav_fixture(tmp_path)
    monkeypatch.setattr(module, "EXPECTED_SHA256", digest)

    def forbidden(*args, **kwargs):
        pytest.fail("verify-only must not use the network")

    monkeypatch.setattr(module.urllib.request, "urlopen", forbidden)
    assert module.main(["--verify-only", "--destination", str(path)]) == 0
    summary = json.loads(capsys.readouterr().out)
    assert summary["metadata"]["sha256"] == digest
    assert summary["metadata"]["channels"] == 2
    assert summary["qc"]["finite"] is True
    assert summary["array_shapes"]["waveform"] == [3]


def test_cli_reports_missing_file(tmp_path, capsys):
    module = fetch_module()
    assert hasattr(module, "main"), "fetch and verify CLI feature is missing"
    assert module.main(["--verify-only", "--destination", str(tmp_path / "absent.wav")]) == 1
    assert "reference error" in capsys.readouterr().err


@pytest.mark.parametrize("max_bytes", [-1, 0, 20_000_001, 2.5, True])
def test_byte_bound_configuration_must_remain_bounded(tmp_path, monkeypatch, max_bytes):
    path, digest = wav_fixture(tmp_path)
    with pytest.raises(ValueError, match="max_bytes"):
        reference_module().load_reference(path, expected_sha256=digest, max_bytes=max_bytes)
    module = fetch_module()
    monkeypatch.setattr(
        module.urllib.request,
        "urlopen",
        lambda *args, **kwargs: pytest.fail("invalid bounds must fail before network"),
    )
    with pytest.raises(ValueError, match="max_bytes"):
        module.fetch_reference(tmp_path / "out.wav", max_bytes=max_bytes)


@pytest.mark.parametrize("digest", ["not-a-hash", "A" * 64, "0" * 63])
def test_expected_digest_requires_literal_lowercase_sha256(tmp_path, digest):
    path, _ = wav_fixture(tmp_path)
    with pytest.raises(ValueError, match="64 lowercase"):
        reference_module().load_reference(path, expected_sha256=digest)


def test_qc_flags_rail_samples_and_dc_without_claiming_physical_clipping():
    recording = {"samples": np.array([[-1.0], [-1.0], [0.999969482421875], [0.5]]), "fs": 8000}
    result = reference_module().analyse_reference(recording)
    assert result["qc"]["near_full_scale_fraction"] == 0.75
    assert "possible digital clipping" in result["qc"].get("warnings", [])
    assert "substantial DC offset relative to digital RMS" in result["qc"]["warnings"]


@pytest.mark.parametrize(
    "kwargs",
    [
        {"timeout_s": None},
        {"timeout_s": 0},
        {"timeout_s": np.nan},
        {"timeout_s": 21},
        {"total_timeout_s": np.inf},
        {"total_timeout_s": 61},
    ],
)
def test_fetch_rejects_unbounded_timeouts_before_network(tmp_path, monkeypatch, kwargs):
    module = fetch_module()
    monkeypatch.setattr(
        module.urllib.request,
        "urlopen",
        lambda *args, **kw: pytest.fail("invalid timeout reached network"),
    )
    with pytest.raises(ValueError, match="timeout"):
        module.fetch_reference(tmp_path / "unused.wav", **kwargs)


def test_fetch_uses_single_buffer_reads_so_deadline_can_be_checked(tmp_path, monkeypatch):
    import io

    module = fetch_module()
    path, digest = wav_fixture(tmp_path)
    payload = path.read_bytes()
    monkeypatch.setattr(module, "EXPECTED_SHA256", digest)

    class Response(io.BytesIO):
        headers = {}

        def geturl(self):
            return module.SOURCE_URL

        def read(self, size=-1):
            pytest.fail("bulk read can hide a slow trickle from deadline checks")

        def read1(self, size=-1):
            return super().read(min(size, 10))

    monkeypatch.setattr(module.urllib.request, "urlopen", lambda *args, **kwargs: Response(payload))
    result = module.fetch_reference(tmp_path / "single-read.wav")
    assert result["sha256"] == digest
