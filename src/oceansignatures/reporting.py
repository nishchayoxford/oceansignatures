"""Faithful display reductions and independent analytic channel references."""

import numpy as np


def envelope(time_s, signal, max_bins=1200):
    """Return block-midpoint time and actual extrema, retaining transient peaks."""
    time_s = np.asarray(time_s, dtype=float)
    signal = np.asarray(signal, dtype=float)
    if time_s.ndim != 1 or signal.shape != time_s.shape or not signal.size:
        raise ValueError("Time and signal must be equal nonempty one-dimensional arrays")
    if not np.all(np.isfinite(time_s)) or not np.all(np.isfinite(signal)) or max_bins < 1:
        raise ValueError("Display inputs must be finite with a positive bin count")
    edges = np.linspace(0, signal.size, min(max_bins, signal.size) + 1, dtype=int)
    low = np.array([np.min(signal[a:b]) for a, b in zip(edges[:-1], edges[1:])])
    high = np.array([np.max(signal[a:b]) for a, b in zip(edges[:-1], edges[1:])])
    centres = np.array([(time_s[a] + time_s[b - 1]) / 2 for a, b in zip(edges[:-1], edges[1:])])
    return centres, low, high


def channel_transfer(
    frequency_hz,
    range_m,
    source_depth_m=8.0,
    receiver_depth_m=12.0,
    sound_speed=1500.0,
    surface=True,
):
    """Continuous two-ray H(f), relative to an ideal 1 m reference pressure.

    This curve is analytic and is not a discretized waveform propagation result.
    """
    frequency = np.asarray(frequency_hz, dtype=float)
    direct = np.sqrt(np.asarray(range_m) ** 2 + (source_depth_m - receiver_depth_m) ** 2)
    reflected = np.sqrt(np.asarray(range_m) ** 2 + (source_depth_m + receiver_depth_m) ** 2)
    h = np.exp(-2j * np.pi * frequency * direct / sound_speed) / direct
    if surface:
        h -= np.exp(-2j * np.pi * frequency * reflected / sound_speed) / reflected
    return h
