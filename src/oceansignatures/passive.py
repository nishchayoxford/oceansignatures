"""Synthetic passive source/channel/receiver; homogeneous 3-D acoustics.

All pressure is Pa. Range is horizontal separation, depths positive down.
These are pedagogical phenomenological sources, not validated animals or ships.
"""

import numpy as np
from scipy import signal as sp_signal
from .metrics import rms, _signal


def _positive(value, name, allow_zero=False):
    if (
        not np.isscalar(value)
        or not np.isfinite(value)
        or (value < 0 if allow_zero else value <= 0)
    ):
        raise ValueError(f"{name} must be finite and {'nonnegative' if allow_zero else 'positive'}")


def _seed(seed):
    if isinstance(seed, (bool, np.bool_)) or not isinstance(seed, (int, np.integer)) or seed < 0:
        raise ValueError("seed must be a nonnegative integer")
    return np.random.default_rng(seed)


def make_source(kind, fs=48000, duration_s=2.0, rms_pa=1.0, seed=2026):
    """Return floor(fs*duration_s) real Pa samples, RMS-normalized over
    the entire unpadded emission record AFTER tapering, at a nominal 1 m.
    fs >= 32 kHz and duration >= 0.1 s are required for these fixed signatures.
    A zero target RMS produces silence. No peak/listening normalization.
    """
    if kind not in ("whistle", "clicks", "vessel"):
        raise ValueError("unknown source kind")
    _positive(fs, "fs")
    _positive(duration_s, "duration_s")
    _positive(rms_pa, "rms_pa", allow_zero=True)
    rng = _seed(seed)
    if fs < 32000 or duration_s < 0.1:
        raise ValueError("source requires fs >= 32000 and duration_s >= 0.1")
    n = int(fs * duration_s)
    t = np.arange(n) / fs
    duration = n / fs
    if kind == "whistle":
        phase = 2 * np.pi * (1600 * t + 1300 * t * t / duration) + 30 * (
            1 - np.cos(2 * np.pi * 5 * t)
        )
        x = np.sin(phase)
    elif kind == "clicks":
        x = np.zeros(n)
        for centre in np.arange(0.05, duration - 0.02, 0.16):
            offset = t - (centre + rng.uniform(-0.003, 0.003))
            support = np.abs(offset) <= 6 * 0.00012
            u = offset[support]
            x[support] += (
                rng.uniform(0.8, 1.2)
                * np.exp(-0.5 * (u / 0.00012) ** 2)
                * np.cos(2 * np.pi * 8000 * u)
            )
    else:
        harmonic = sum(
            np.sin(2 * np.pi * 80 * k * t + rng.uniform(0, 2 * np.pi)) / k**1.2
            for k in range(1, 21)
        )
        frequencies = np.fft.rfftfreq(n, 1 / fs)
        spectrum = np.fft.rfft(rng.normal(size=n)) / np.sqrt(1 + (frequencies / 1000) ** 2)
        spectrum[0] = 0
        broadband = np.fft.irfft(spectrum, n=n)
        x = (harmonic / rms(harmonic) + 0.35 * broadband / rms(broadband)) * (
            0.75 + 0.25 * np.sin(2 * np.pi * 3 * t)
        )
    x *= sp_signal.windows.tukey(n, alpha=min(1.0, 0.04 / duration))
    return x * (rms_pa / rms(x))


def path_parameters(
    range_m, source_depth_m=8.0, receiver_depth_m=12.0, sound_speed=1500.0, surface=True
):
    """Direct and optional ideal pressure-release image paths.

    Gain is signed 1 m / slant distance; delay is slant distance / c.
    range_m is strictly positive horizontal separation; depths may be zero.
    """
    for value, name in ((range_m, "range_m"), (sound_speed, "sound_speed")):
        _positive(value, name)
    for value, name in ((source_depth_m, "source_depth_m"), (receiver_depth_m, "receiver_depth_m")):
        _positive(value, name, allow_zero=True)
    if not isinstance(surface, (bool, np.bool_)):
        raise ValueError("surface must be a boolean")
    geometries = [("direct", source_depth_m - receiver_depth_m, 1.0)]
    if surface:
        geometries.append(("surface", source_depth_m + receiver_depth_m, -1.0))
    paths = []
    for kind, vertical, sign in geometries:
        distance = float(np.hypot(range_m, vertical))
        paths.append(
            {
                "kind": kind,
                "distance_m": distance,
                "delay_s": distance / sound_speed,
                "gain": sign / distance,
            }
        )
    return paths


def propagate(signal, fs, paths):
    """Signed sum with 129-tap Kaiser(beta=8.6) windowed-sinc delays.

    For d=floor(delay*fs), mu=delay*fs-d, h[k]=sinc(k-mu)*Kaiser[k],
    k=-64..64, normalized to unit DC. Linear full convolution, not circular
    FFT shifting; kernel centre is removed, not added as physical latency.
    Retain N+ceil(max_delay*fs)+64 samples, including the numerical tail.
    Up to 64 samples of noncausal interpolation pre-ringing are possible;
    values before t=0 are discarded. Do not infer a physical precursor.
    Accuracy is tested through 0.4*fs, not at Nyquist. No gain correction
    hides frequency-dependent interpolation attenuation.
    """
    x = _signal(signal)
    _positive(fs, "fs")
    paths = list(paths)
    if not paths:
        raise ValueError("at least one path is required")
    for p in paths:
        if not isinstance(p, dict) or not {"delay_s", "gain"} <= p.keys():
            raise ValueError("each path must contain delay_s and gain")
        _positive(p["delay_s"], "delay_s", allow_zero=True)
        if not np.isscalar(p["gain"]) or not np.isfinite(p["gain"]):
            raise ValueError("gain must be finite")
    y = np.zeros(x.size + int(np.ceil(max(p["delay_s"] for p in paths) * fs)) + 64)
    for p in paths:
        samples = p["delay_s"] * fs
        delay = int(np.floor(samples))
        fraction = samples - delay
        if fraction == 0:
            y[delay : delay + x.size] += p["gain"] * x
        else:
            kernel = np.sinc(np.arange(-64, 65) - fraction) * np.kaiser(129, 8.6)
            kernel /= np.sum(kernel)
            delayed = sp_signal.convolve(x, kernel, mode="full", method="direct")
            offset = delay - 64
            left = max(0, offset)
            skip = max(0, -offset)
            y[left : left + len(delayed) - skip] += p["gain"] * delayed[skip:]
    return y


def passive_case(kind, range_m, surface, fs=48000, duration_s=2.0, seed=2026):
    """Executed source → paths → generic receiver → pressure-noise case.

    received_clean_pa is AFTER the causal receiver filter; received_channel_pa
    is BEFORE it. noise_pa is input-referred white noise passed through the
    same filter with independent prehistory. All arrays share the full clock.
    Source RMS is over emission samples; receiver RMS/SNR use the explicit
    half-open analysis window including complete channel and filter tails.
    """
    source = make_source(kind, fs, duration_s, 1.0, seed)
    paths = path_parameters(range_m, surface=surface)
    channel = propagate(source, fs, paths)
    taps = sp_signal.firwin(257, 10000.0, fs=fs, window=("kaiser", 8.6), scale=True)
    clean = sp_signal.convolve(channel, taps, mode="full", method="direct")
    rng = np.random.default_rng(np.random.SeedSequence([int(seed), 1]))
    noise_input = rng.normal(0.0, 0.0005, size=len(clean) + len(taps) - 1)
    noise = sp_signal.convolve(noise_input, taps, mode="valid", method="direct")
    received = clean + noise
    start = max(0, int(np.floor(min(p["delay_s"] for p in paths) * fs)) - 64)
    stop = len(clean)
    window = slice(start, stop)
    clean_rms, noise_rms, received_rms = (
        rms(clean[window]),
        rms(noise[window]),
        rms(received[window]),
    )
    metadata = {
        "kind": kind,
        "range_m": float(range_m),
        "surface": bool(surface),
        "seed": int(seed),
        "paths": paths,
        "pressure_unit": "Pa",
        "reference_pressure_pa": 1e-6,
        "source": {
            "reference_distance_m": 1.0,
            "target_rms_pa": 1.0,
            "rms_convention": "full unpadded emission record after taper",
            "emission_samples": len(source),
            "duration_s": len(source) / fs,
            "requested_duration_s": float(duration_s),
            "depth_m": 8.0,
            "model": "analytic/phenomenological; not species or machinery validation",
        },
        "channel": {
            "sound_speed_m_s": 1500.0,
            "geometry": "homogeneous 3-D point source",
            "range_convention": "horizontal",
            "surface_reflection_coefficient": -1.0,
            "fractional_delay": "129-tap Kaiser beta=8.6 sinc, DC normalized",
            "interpolation_half_support_samples": 64,
            "accuracy_test_band_max_fraction_fs": 0.4,
            "omitted": ["seabed", "refraction", "absorption", "motion/Doppler", "scattering"],
        },
        "receiver": {
            "depth_m": 12.0,
            "model": "assumed generic causal FIR lowpass; not measured hardware",
            "cutoff_hz": 10000.0,
            "numtaps": 257,
            "window": "Kaiser beta=8.6",
            "group_delay_samples": 128,
            "group_delay_s": 128 / fs,
            "fir_coefficients": taps.tolist(),
            "noise_input_std_pa": 0.0005,
            "noise": "Gaussian white input-referred noise; stationary FIR prehistory",
            "noise_rng": "NumPy default_rng(SeedSequence([seed, 1]))",
        },
        "analysis_window_samples": [start, stop],
        "analysis_window_s": [start / fs, stop / fs],
        "analysis_window_convention": "half-open; earliest interpolation support through full receiver tail",
        "measured": {
            "source_rms_pa": rms(source),
            "received_clean_rms_pa": clean_rms,
            "received_rms_pa": received_rms,
            "noise_rms_pa": noise_rms,
            "snr_db": float(20 * np.log10(clean_rms / noise_rms)),
            "received_level_db_re_1uPa": float(20 * np.log10(received_rms / 1e-6)),
        },
        "snr_definition": "20 log10(clean receiver RMS / filtered noise RMS), same analysis window",
        "validation_status": "analytical numerical controls only; no experimental validation",
    }
    return {
        "time_s": np.arange(stop) / fs,
        "source_pa": np.pad(source, (0, stop - len(source))),
        "received_channel_pa": np.pad(channel, (0, stop - len(channel))),
        "received_clean_pa": clean,
        "noise_pa": noise,
        "received_pa": received,
        "fs": fs,
        "metadata": metadata,
    }
