"""Real, uniformly sampled pressure metrics; no listening normalization."""

import numpy as np
from scipy import signal as sp_signal


def _signal(signal):
    raw = np.asarray(signal)
    if np.iscomplexobj(raw):
        raise ValueError("signal must be real")
    x = np.asarray(raw, dtype=float)
    if x.ndim != 1 or x.size == 0 or not np.all(np.isfinite(x)):
        raise ValueError("signal must be a nonempty finite 1-D real array")
    return x


def rms(signal):
    """Unweighted full-array RMS, including DC, in the input amplitude units."""
    x = _signal(signal)
    scale = float(np.max(np.abs(x)))
    return scale * float(np.sqrt(np.mean((x / scale) ** 2))) if scale else 0.0


def _spectral_inputs(signal, fs, nperseg):
    x = _signal(signal)
    if not np.isscalar(fs) or not np.isfinite(fs) or fs <= 0:
        raise ValueError("fs must be finite and positive")
    if isinstance(nperseg, (bool, np.bool_)) or not isinstance(nperseg, (int, np.integer)):
        raise ValueError("nperseg must be an integer")
    if not 2 <= nperseg <= x.size:
        raise ValueError("nperseg must be between 2 and signal length; no implicit shortening")
    return x


def welch_psd(signal, fs, nperseg=2048):
    """(Hz, Pa²/Hz) for Pa input: periodic Hann, segment-mean removal,
    50% overlap (floor), mean average, one-sided density. No padding.
    Incomplete final segments are discarded; nperseg > length is rejected.
    """
    x = _spectral_inputs(signal, fs, nperseg)
    return sp_signal.welch(
        x,
        fs=fs,
        window="hann",
        nperseg=nperseg,
        noverlap=nperseg // 2,
        detrend="constant",
        return_onesided=True,
        scaling="density",
        average="mean",
    )


def spectrogram_psd(signal, fs, nperseg=512, noverlap=384):
    """(frame-centre seconds, Hz, frequency-by-time Pa²/Hz).

    Periodic Hann, per-frame mean removal, one-sided density; no boundary
    extension or padding. Frames cover [start, start+nperseg)/fs and any
    incomplete final frame is discarded. Window-weighted Parseval applies,
    not equality to the raw full-record RMS (especially for transients).
    """
    x = _spectral_inputs(signal, fs, nperseg)
    if (
        isinstance(noverlap, (bool, np.bool_))
        or not isinstance(noverlap, (int, np.integer))
        or not 0 <= noverlap < nperseg
    ):
        raise ValueError("noverlap must be an integer in [0, nperseg)")
    f, t, p = sp_signal.spectrogram(
        x,
        fs=fs,
        window="hann",
        nperseg=nperseg,
        noverlap=noverlap,
        detrend="constant",
        return_onesided=True,
        scaling="density",
        mode="psd",
    )
    return t, f, p
