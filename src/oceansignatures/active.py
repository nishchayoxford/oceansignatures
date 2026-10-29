"""Stationary scalar point-target benchmark in homogeneous free space."""

from numbers import Integral, Real

import numpy as np
from scipy.signal import correlate

_SOUND_SPEED = 1500.0
_PING_DURATION = 0.05
_F_START, _F_END = 3000.0, 9000.0
_REFERENCE_DISTANCE = 1.0
_TARGET_GAIN = 1.0
_SEARCH_DURATION = 1.0


def _finite_ping(time_s):
    """Unscaled finite 50 ms Hann-tapered real 3–9 kHz LFM."""
    duration_s = _PING_DURATION
    u = np.asarray(time_s)
    inside = (u >= 0.0) & (u < duration_s)
    values = np.zeros_like(u, dtype=float)
    phase = (
        2
        * np.pi
        * (_F_START * u[inside] + 0.5 * ((_F_END - _F_START) / duration_s) * u[inside] ** 2)
    )
    taper = 0.5 * (1 - np.cos(2 * np.pi * u[inside] / duration_s))
    values[inside] = taper * np.cos(phase)
    return values


def active_case(range_m=150.0, snr_db=-12.0, fs=48000, seed=2026):
    """Return a stationary synthetic two-way echo and its range estimate.

    ``range_m`` is horizontal separation at equal depth, in [1, 750] m;
    ``snr_db`` is clean echo-support mean square / expected noise variance,
    in [-120, 120] dB. ``fs`` is an integer in [24000, 384000] Hz and
    ``seed`` a nonnegative integer. The 1 Pa RMS reference ping is a finite
    50 ms real Hann-tapered 3–9 kHz LFM. The chosen target scalar is 1.

    Time starts at transmission. Echo delay is exactly 2*range_m/1500 s;
    the continuous waveform is evaluated at fractional sample offsets.
    Returned ``lag_s`` contains valid-correlation *template-start* lags
    from 0 to 1 s. ``matched_filter`` is signed correlation divided by
    noise_std*sqrt(template_energy), hence dimensionless unit noise SD at
    each lag (not independent across lags). Range is the integer-lag global
    absolute maximum times 1500/2: an estimator, not a detection decision.

    All pressure arrays share ``time_s``. ``received_clean_pa`` and
    ``noise_pa`` supplement the frozen API for reproducibility. Metadata
    records computation-derived quantities and assumptions. See ACTIVE.md.
    """
    if (
        isinstance(range_m, bool)
        or not isinstance(range_m, Real)
        or not np.isfinite(range_m)
        or not 1 <= range_m <= 750
    ):
        raise ValueError("range_m must be finite and between 1 and 750 m")
    if (
        isinstance(snr_db, bool)
        or not isinstance(snr_db, Real)
        or not np.isfinite(snr_db)
        or not -120 <= snr_db <= 120
    ):
        raise ValueError("snr_db must be finite and between -120 and 120 dB")
    if isinstance(fs, bool) or not isinstance(fs, Integral) or not 24000 <= fs <= 384000:
        raise ValueError("fs must be an integer between 24000 and 384000 Hz")
    if isinstance(seed, bool) or not isinstance(seed, Integral) or seed < 0:
        raise ValueError("seed must be a nonnegative integer")
    range_m, snr_db, fs, seed = float(range_m), float(snr_db), int(fs), int(seed)
    ping_time = np.arange(int(np.ceil(_PING_DURATION * fs))) / fs
    raw_ping = _finite_ping(ping_time)
    scale = 1 / np.sqrt(np.mean(raw_ping**2))
    ping = scale * raw_ping
    time_s = np.arange(int(_SEARCH_DURATION * fs) + len(ping)) / fs
    transmitted = np.zeros_like(time_s)
    transmitted[: len(ping)] = ping
    delay_s = 2 * range_m / _SOUND_SPEED
    echo_sample_offsets = np.arange(len(time_s)) - delay_s * fs
    echo_time = echo_sample_offsets / fs
    gain = _TARGET_GAIN * (_REFERENCE_DISTANCE / range_m) ** 2
    clean_echo = scale * _finite_ping(echo_time) * gain
    echo_support = (echo_sample_offsets >= 0) & (echo_sample_offsets < _PING_DURATION * fs)
    support_indices = np.flatnonzero(echo_support)
    echo_power = float(np.mean(clean_echo[echo_support] ** 2))
    noise_std = np.sqrt(echo_power / 10 ** (snr_db / 10))
    noise = np.random.default_rng(seed).normal(0, noise_std, len(time_s))
    received = clean_echo + noise
    normalization = noise_std * np.sqrt(np.sum(ping**2))
    matched_filter = correlate(received, ping, mode="valid", method="fft") / normalization
    lag_s = np.arange(len(matched_filter)) / fs
    peak_index = int(np.argmax(np.abs(matched_filter)))
    estimated_range = lag_s[peak_index] * _SOUND_SPEED / 2
    realized_noise_power = float(np.mean(noise[echo_support] ** 2))
    return {
        "time_s": time_s,
        "transmitted_pa": transmitted,
        "received_pa": received,
        "received_clean_pa": clean_echo,
        "noise_pa": noise,
        "lag_s": lag_s,
        "matched_filter": matched_filter,
        "estimated_range_m": float(estimated_range),
        "true_range_m": float(range_m),
        "fs": fs,
        "metadata": {
            "seed": int(seed),
            "model": "stationary monostatic scalar point target; homogeneous 3-D free space",
            "range_definition": "horizontal separation; target and transceiver at equal depth",
            "receiver_model": "ideal pressure sampler; no receiver filter, clipping or transmit leakage",
            "noise_model": "iid zero-mean real Gaussian; full sampled Nyquist band; not realization-normalized",
            "fractional_delay_method": "evaluate continuous finite LFM at (sample_index - delay_s*fs)/fs",
            "lag_definition": "template start lag k/fs; valid correlation; no filter-length offset",
            "range_estimator": "global maximum absolute signed correlation; integer lag; no detection threshold",
            "omitted_physics": [
                "surface and seabed multipath",
                "refraction",
                "absorption",
                "Doppler",
                "frequency-dependent or extended-target scattering",
                "reverberation",
            ],
            "sound_speed_m_s": _SOUND_SPEED,
            "expected_delay_s": float(delay_s),
            "target_scalar_gain": _TARGET_GAIN,
            "reference_distance_m": _REFERENCE_DISTANCE,
            "two_way_pressure_gain": float(gain),
            "ping_duration_s": _PING_DURATION,
            "ping_f_start_hz": _F_START,
            "ping_f_end_hz": _F_END,
            "ping_taper": "continuous Hann",
            "ping_samples": len(ping),
            "ping_rms_pa": float(np.sqrt(np.mean(ping**2))),
            "record_samples": len(time_s),
            "range_bin_m": _SOUND_SPEED / (2 * fs),
            "range_error_m": float(estimated_range - range_m),
            "search_lag_s": [float(lag_s[0]), float(lag_s[-1])],
            "echo_support_samples": [int(support_indices[0]), int(support_indices[-1] + 1)],
            "echo_mean_square_pa2": echo_power,
            "realized_noise_mean_square_pa2": realized_noise_power,
            "realized_input_snr_db": float(10 * np.log10(echo_power / realized_noise_power)),
            "input_snr_db": float(snr_db),
            "noise_std_pa": float(noise_std),
            "matched_filter_units": "dimensionless; noise-only standard deviation = 1 at each lag",
            "peak_index": peak_index,
            "peak_abs_score": float(abs(matched_filter[peak_index])),
            "input_snr_definition": "10*log10(mean(clean_echo_pa[echo_support]**2)/noise_variance_pa2)",
        },
    }
