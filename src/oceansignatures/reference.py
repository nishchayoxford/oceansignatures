"""Uncalibrated public reference audio; digital units are never acoustic pressure."""

from __future__ import annotations

import hashlib
import io
import re
from pathlib import Path

import numpy as np
import soundfile as sf
from scipy import signal

MAX_BYTES = 20_000_000


def load_reference(path, *, expected_sha256, max_bytes=MAX_BYTES):
    """Decode a SHA256-pinned RIFF/WAVE without downmixing or normalization.

    PCM integers are converted by libsndfile to digital full-scale fractions;
    float WAVs retain their stored scale. Neither representation is in Pa.
    The byte bound is enforced before decoding, including a bounded read.
    """
    if type(max_bytes) is not int or not 0 < max_bytes <= MAX_BYTES:
        raise ValueError(f"max_bytes must be an integer in [1, {MAX_BYTES}]")
    if not isinstance(expected_sha256, str) or not re.fullmatch(r"[0-9a-f]{64}", expected_sha256):
        raise ValueError("expected_sha256 must be exactly 64 lowercase hexadecimal characters")
    path = Path(path)
    if path.stat().st_size > max_bytes:
        raise ValueError("WAV exceeds size limit")
    with path.open("rb") as stream:
        payload = stream.read(max_bytes + 1)
    if len(payload) > max_bytes:
        raise ValueError("WAV exceeds size limit")
    digest = hashlib.sha256(payload).hexdigest()
    if digest != expected_sha256:
        raise ValueError("SHA256 mismatch; refusing changed reference")
    if len(payload) < 12 or payload[:4] != b"RIFF" or payload[8:12] != b"WAVE":
        raise ValueError("Expected RIFF/WAVE magic")
    if int.from_bytes(payload[4:8], "little") + 8 != len(payload):
        raise ValueError("RIFF size mismatch: truncated or trailing content")
    with sf.SoundFile(io.BytesIO(payload)) as stream:
        info = sf.info(io.BytesIO(payload))
        samples = stream.read(dtype="float64", always_2d=True)
        fs = stream.samplerate
    if not len(samples):
        raise ValueError("WAV is empty")
    if not np.isfinite(samples).all():
        raise ValueError("WAV contains non-finite samples")
    return {
        "samples": samples,
        "fs": fs,
        "metadata": {
            "filename": path.name,
            "sha256": hashlib.sha256(payload).hexdigest(),
            "bytes": len(payload),
            "magic": "RIFF/WAVE",
            "format": info.format,
            "subtype": info.subtype,
            "frames": len(samples),
            "channels": info.channels,
            "sample_count": samples.size,
            "sample_rate_hz": fs,
            "duration_s": len(samples) / fs,
            "amplitude_units": "normalized digital units",
            "calibration_status": "uncalibrated; conversion to Pa unavailable",
        },
    }


def analyse_reference(recording, *, channel=0, nperseg=2048, spectrogram_nperseg=1024):
    """Return unchanged waveform, relative density spectra, and digital QC.

    One selected channel, periodic Hann windows, constant detrending per segment,
    one-sided density scaling, Welch mean averaging with 50% overlap, spectrogram
    with 75% overlap and no padding. Time/frequency refer to the distributed WAV,
    not necessarily to the original acoustic event. No species classification.
    """
    samples = np.asarray(recording["samples"], dtype=np.float64)
    fs = float(recording["fs"])
    if samples.ndim != 2 or samples.shape[0] < 2 or samples.shape[1] < 1:
        raise ValueError("samples must be a frames-by-channels array with at least two frames")
    if not np.isfinite(samples).all() or not np.isfinite(fs) or fs <= 0:
        raise ValueError("samples and positive fs must be finite")
    if not isinstance(channel, (int, np.integer)) or not 0 <= channel < samples.shape[1]:
        raise ValueError("channel must be an in-range nonnegative integer")
    for value in (nperseg, spectrogram_nperseg):
        if not isinstance(value, (int, np.integer)) or value < 2:
            raise ValueError("segment lengths must be integers >= 2")
    waveform = samples[:, channel].copy()
    square_error = "Signal amplitude-square is not representable as float64"
    try:
        with np.errstate(over="raise", invalid="raise", under="ignore"):
            squared = waveform**2
    except FloatingPointError as exc:
        raise ValueError(square_error) from exc
    if np.any((waveform != 0) & (squared == 0)):
        raise ValueError(square_error)
    n = len(waveform)
    welch_n = min(nperseg, n)
    spec_n = min(spectrogram_nperseg, n)
    shared = {
        "window": "hann",
        "detrend": "constant",
        "return_onesided": True,
        "scaling": "density",
    }
    # Scaling is internal only: restore squared digital units after SciPy's
    # FFT and averaging, which can otherwise overflow before density scaling.
    normalizer = float(np.max(np.abs(waveform))) or 1.0
    numeric_error = "Reference density or time/frequency axes are not representable as float64"
    # SciPy may form fs * window energy using Python arithmetic (no warning).
    # Bound fs * segment length, which also bounds Hann window energy.
    if fs > np.finfo(np.float64).max / max(welch_n, spec_n):
        raise ValueError(numeric_error)
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise", under="ignore"):
            normalized = waveform / normalizer
            time_s = np.arange(n) / fs
            f, unit_psd = signal.welch(
                normalized, fs=fs, nperseg=welch_n, noverlap=welch_n // 2, average="mean", **shared
            )
            sfreq, stime, unit_spec = signal.spectrogram(
                normalized, fs=fs, nperseg=spec_n, noverlap=3 * spec_n // 4, mode="psd", **shared
            )
            psd = unit_psd * normalizer**2
            spec = unit_spec * normalizer**2
    except (FloatingPointError, OverflowError, ZeroDivisionError) as exc:
        raise ValueError(numeric_error) from exc
    if (
        not all(np.isfinite(a).all() for a in (time_s, f, stime, sfreq, psd, spec))
        or np.any((unit_psd != 0) & (psd == 0))
        or np.any((unit_spec != 0) & (spec == 0))
    ):
        raise ValueError(numeric_error)
    metadata = {
        **recording.get("metadata", {}),
        "selected_channel": channel,
        "amplitude_units": "normalized digital units",
        "calibration_status": "uncalibrated; conversion to Pa unavailable",
        "psd_units": "digital units squared per Hz",
        "time_frequency_basis": "distributed WAV sample rate; original event timing unverified",
        "waveform_processing": "none; no downmixing, gain normalization, or resampling",
        "welch": {
            **shared,
            "nperseg": welch_n,
            "noverlap": welch_n // 2,
            "average": "mean",
            "frequency_bin_hz": fs / welch_n,
        },
        "spectrogram": {
            **shared,
            "nperseg": spec_n,
            "noverlap": 3 * spec_n // 4,
            "mode": "psd",
            "frequency_bin_hz": fs / spec_n,
            "time_bin_s": (spec_n - 3 * spec_n // 4) / fs,
        },
    }
    qc = {
        "finite": True,
        "mean": float(normalizer * np.mean(normalized)),
        "rms": float(normalizer * np.sqrt(np.mean(normalized**2))),
        "peak_absolute": float(np.max(np.abs(waveform))),
        "zero_fraction": float(np.mean(waveform == 0)),
        "silent": bool(np.all(waveform == 0)),
        "near_full_scale_threshold": 0.999,
        "near_full_scale_fraction": float(np.mean(np.abs(waveform) >= 0.999)),
        "clipping_caveat": "near-full-scale samples flag possible clipping; not proof",
    }
    qc["warnings"] = []
    if qc["near_full_scale_fraction"] > 0:
        qc["warnings"].append("possible digital clipping")
    if abs(qc["mean"]) > 0.1 * qc["rms"]:
        qc["warnings"].append("substantial DC offset relative to digital RMS")
    return {
        "fs": fs,
        "time_s": time_s,
        "waveform": waveform,
        "frequency_hz": f,
        "psd": psd,
        "spectrogram_time_s": stime,
        "spectrogram_frequency_hz": sfreq,
        "spectrogram_psd": spec,
        "metadata": metadata,
        "qc": qc,
    }
