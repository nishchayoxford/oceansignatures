"""Empirical expected acoustic spectra; see docs/NOISE_PREDICTION.md.

One-sided receiver PSDs use Pa²/Hz. Vessel source PSD uses Pa² m²/Hz,
numerically the extrapolated 1 m pressure PSD, not a measured near-field PSD.
"""

from dataclasses import asdict, dataclass

import numpy as np

from .passive import path_parameters


@dataclass(frozen=True)
class SeaWater:
    """Constant path-average properties for Ainslie–McColm absorption.

    NPL's cited comparisons vary properties one at a time, not every corner
    of this guard box. Depth here is water pressure depth, NOT path length.
    """

    temperature_c: float = 10.0
    salinity_ppt: float = 35.0
    depth_m: float = 50.0
    ph: float = 8.0

    def __post_init__(self):
        for name, low, high in (
            ("temperature_c", -6, 35),
            ("salinity_ppt", 5, 50),
            ("ph", 7.7, 8.3),
        ):
            object.__setattr__(
                self, name, _scalar(getattr(self, name), name, low, high, open_bounds=True)
            )
        depth = _scalar(self.depth_m, "depth_m", 0, 7000)
        if depth == 7000:
            raise ValueError("depth_m must be < 7000")
        object.__setattr__(self, "depth_m", depth)


def absorption_db_per_km(frequency_hz, *, water=None, allow_low_frequency_extrapolation=False):
    """Ainslie–McColm α in dB/km; 100 Hz–1 MHz, optionally extrapolate to 20 Hz."""
    water = SeaWater() if water is None else water
    if not isinstance(water, SeaWater):
        raise ValueError("water must be a SeaWater instance")
    f = _frequency(frequency_hz, 100, 1e6, allow_low_frequency_extrapolation) / 1000
    t, s, d, ph = water.temperature_c, water.salinity_ppt, water.depth_m / 1000, water.ph
    f1 = 0.78 * np.sqrt(s / 35) * np.exp(t / 26)
    f2 = 42 * np.exp(t / 17)
    return (
        0.106 * np.exp((ph - 8) / 0.56) * f1 * f * f / (f1 * f1 + f * f)
        + 0.52 * (1 + t / 43) * (s / 35) * np.exp(-d / 6) * f2 * f * f / (f2 * f2 + f * f)
        + 0.00049 * np.exp(-(t / 27 + d / 17)) * f * f
    )


def propagation_loss(
    frequency_hz, *, distance_m, water=None, allow_low_frequency_extrapolation=False
):
    """Monopole 20 log10(r/1 m) + α r/1000, dB, for scalar slant r >= 1 m.

    Absorption spans the full path: source strength is a loss-free extrapolated
    monopole, not a hydrophone measurement made at a physical radius of 1 m.
    """
    r = _scalar(distance_m, "distance_m", 1, 1e7)
    alpha = absorption_db_per_km(
        frequency_hz,
        water=water,
        allow_low_frequency_extrapolation=allow_low_frequency_extrapolation,
    )
    return 20 * np.log10(r) + alpha * r / 1000


def _nonnegative(values, name):
    a = np.asarray(values, dtype=float)
    if not a.size or not np.all(np.isfinite(a)) or np.any(a < 0):
        raise ValueError(f"{name} must be finite and nonnegative")
    return a


def db_to_psd(level_db):
    """dB re 1 µPa²/Hz to Pa²/Hz; source references retain their m² factor."""
    level = np.asarray(level_db, dtype=float)
    if not level.size or not np.all(np.isfinite(level)) or np.any(np.abs(level) > 3000):
        raise ValueError("level_db must be finite within ±3000 dB")
    return 10 ** ((level - 120) / 10)


def psd_to_db(psd):
    """Pa²/Hz to dB re 1 µPa²/Hz; zero maps to -inf, no arbitrary floor."""
    p = _nonnegative(psd, "psd")
    with np.errstate(divide="ignore"):
        return 10 * np.log10(p) + 120


def psd_to_band_power(psd, bandwidth_hz):
    """Flat-band PSD times explicit bandwidth; Pa² (or source Pa² m²)."""
    p = _nonnegative(psd, "psd")
    width = _nonnegative(bandwidth_hz, "bandwidth_hz")
    if np.any(width == 0):
        raise ValueError("bandwidth_hz must be positive")
    return p * width


def band_power_to_psd(band_power, bandwidth_hz):
    """Band mean-square pressure divided by bandwidth, a flat-band approximation."""
    p = _nonnegative(band_power, "band_power")
    width = _nonnegative(bandwidth_hz, "bandwidth_hz")
    if np.any(width == 0):
        raise ValueError("bandwidth_hz must be positive")
    return p / width


def fractional_octave_bandwidth(centre_hz, *, fraction=3):
    """Exact base-10 fractional-octave width; fraction=3 is a decidecade."""
    centre = _nonnegative(centre_hz, "centre_hz")
    n = _scalar(fraction, "fraction", 1, 100)
    if n != int(n) or np.any(centre == 0):
        raise ValueError("positive centres and integer fraction required")
    return centre * (10 ** (0.15 / n) - 10 ** (-0.15 / n))


def _scalar(value, name, low, high, *, open_bounds=False):
    if isinstance(value, (bool, np.bool_)) or not np.isscalar(value):
        raise ValueError(f"{name} must be a finite scalar")
    try:
        value = float(value)
    except (ValueError, TypeError):
        raise ValueError(f"{name} must be a finite scalar") from None
    valid = low < value < high if open_bounds else low <= value <= high
    if not np.isfinite(value) or not valid:
        raise ValueError(f"{name} outside permitted interval {low} to {high}")
    return value


def _frequency(values, low, high, extrapolate=False):
    if not isinstance(extrapolate, (bool, np.bool_)):
        raise ValueError("allow_low_frequency_extrapolation must be boolean")
    f = np.asarray(values, dtype=float)
    lower = 20 if extrapolate else low
    if (
        f.ndim != 1
        or not f.size
        or not np.all(np.isfinite(f))
        or np.any(f < lower)
        or np.any(f > high)
        or np.any(np.diff(f) <= 0)
    ):
        raise ValueError(f"frequency_hz must be finite, increasing, within [{lower}, {high}]")
    return f


def wind_psd(frequency_hz, *, wind_speed_m_s, allow_low_frequency_extrapolation=False):
    """Receiver Pa²/Hz; 100 Hz–100 kHz, optional explicit 20–100 Hz extrapolation.

    Wind 0–20 m/s is a conservative software envelope, not a field-validation claim.
    Calm wind is the nonzero empirical intercept, not silence.
    """
    wind = _scalar(wind_speed_m_s, "wind_speed_m_s", 0, 20)
    f = _frequency(frequency_hz, 100, 100000, allow_low_frequency_extrapolation) / 1000
    level = 50 + 7.5 * np.sqrt(wind) + 20 * np.log10(f) - 40 * np.log10(f + 0.4)
    return 1e-12 * 10 ** (level / 10)


REFERENCE_SPEED_KN = {
    "fishing": 6.4,
    "tug": 3.7,
    "naval": 11.1,
    "recreational": 10.6,
    "government": 8.0,
    "cruise": 17.1,
    "passenger": 9.7,
    "bulker": 13.9,
    "containership": 18.0,
    "vehicle_carrier": 15.8,
    "tanker": 12.4,
    "other": 7.4,
}


def vessel_source_psd(frequency_hz, *, speed_kn, length_m, vessel_class="bulker"):
    """JOMOPANS-ECHO Eqs. 9–11 source PSD, Pa² m²/Hz; 20 Hz–20 kHz.

    Zero speed is rejected: the cruising law cannot predict idle machinery.
    Dredgers are intentionally unsupported (separate operating-state rule).
    Speed (0,60] kn and length (0,500] m are software guards, not validation bounds.
    """
    f = _frequency(frequency_hz, 20, 20000)
    speed = _scalar(speed_kn, "speed_kn", 0, 60)
    length = _scalar(length_m, "length_m", 0, 500)
    if speed == 0 or length == 0:
        raise ValueError("speed_kn and length_m must be positive; idle vessels are unsupported")
    if not isinstance(vessel_class, str) or vessel_class not in REFERENCE_SPEED_KN:
        raise ValueError("unsupported vessel_class")
    vc = REFERENCE_SPEED_KN[vessel_class]
    f1 = 480 / vc
    d = 4 if vessel_class == "cruise" else 3
    level = 191 - 20 * np.log10(f1) - 10 * np.log10((1 - f / f1) ** 2 + d**2)
    if vessel_class in ("bulker", "containership", "vehicle_carrier", "tanker"):
        low = f < 100
        f1 = 600 / vc
        d = 0.8 if vessel_class in ("bulker", "containership") else 1.0
        level[low] = (
            208
            - 40 * np.log10(f1)
            + 10 * np.log10(f[low])
            - 10 * np.log10((1 - (f[low] / f1) ** 2) ** 2 + d**2)
        )
    level += 60 * np.log10(speed / vc) + 20 * np.log10(length / 91.44)
    return 1e-12 * 10 ** (level / 10)


def predict_passage(
    frequency_hz,
    receiver_time_s,
    *,
    speed_kn,
    length_m,
    closest_approach_m,
    wind_speed_m_s,
    vessel_class="bulker",
    source_depth_m=8.0,
    receiver_depth_m=50.0,
    sound_speed_m_s=1500.0,
    closest_approach_emission_time_s=0.0,
    water=None,
    include_wind=True,
    include_vessel=True,
    allow_low_frequency_extrapolation=False,
):
    """Expected quasi-stationary PSD snapshots for a straight, constant-speed track.

    Receiver is fixed at (0,0,zr); source is (v*(te-tCPA), CPA,zs).
    Solve tr=te+r(te)/c for emission time. No Doppler shift or amplitude/time
    Jacobian is applied. No surface/bottom reflection, refraction or flow noise.
    Outputs are a dict of canonical arrays plus JSON-safe model metadata.
    """
    f = _frequency(frequency_hz, 100, 20000, allow_low_frequency_extrapolation)
    tr = np.asarray(receiver_time_s, dtype=float)
    if (
        tr.ndim != 1
        or not tr.size
        or not np.all(np.isfinite(tr))
        or np.any(np.diff(tr) <= 0)
        or np.any(np.abs(tr) > 1e9)
    ):
        raise ValueError("receiver_time_s must be finite, increasing relative seconds within ±1e9")
    cpa = _scalar(closest_approach_m, "closest_approach_m", 1, 1e6)
    zs = _scalar(source_depth_m, "source_depth_m", 0, 6999)
    zr = _scalar(receiver_depth_m, "receiver_depth_m", 0, 6999)
    c = _scalar(sound_speed_m_s, "sound_speed_m_s", 1300, 1700)
    tc = _scalar(closest_approach_emission_time_s, "closest_approach_emission_time_s", -1e9, 1e9)
    for flag in (include_wind, include_vessel):
        if not isinstance(flag, (bool, np.bool_)):
            raise ValueError("component switches must be boolean")
    water = SeaWater() if water is None else water
    source = vessel_source_psd(f, speed_kn=speed_kn, length_m=length_m, vessel_class=vessel_class)
    ambient = wind_psd(
        f,
        wind_speed_m_s=wind_speed_m_s,
        allow_low_frequency_extrapolation=allow_low_frequency_extrapolation,
    )
    alpha = absorption_db_per_km(
        f, water=water, allow_low_frequency_extrapolation=allow_low_frequency_extrapolation
    )
    v = float(speed_kn) * 1852 / 3600
    q = tr - tc
    b2 = cpa * cpa + (zs - zr) ** 2
    te = tc + (c * c * q - np.sqrt(c * c * v * v * q * q + (c * c - v * v) * b2)) / (c * c - v * v)
    horizontal = np.hypot(v * (te - tc), cpa)
    # Reuse the existing validated geometry/delay helper; no duplicate image law.
    paths = [path_parameters(r, zs, zr, c, surface=False)[0] for r in horizontal]
    distance = np.array([p["distance_m"] for p in paths])
    delay = np.array([p["delay_s"] for p in paths])
    loss = 20 * np.log10(distance[:, None]) + distance[:, None] * alpha[None, :] / 1000
    vessel = source[None, :] * 10 ** (-loss / 10) if include_vessel else np.zeros_like(loss)
    wind = np.broadcast_to(ambient if include_wind else np.zeros_like(ambient), loss.shape).copy()
    metadata = {
        "kind": "deterministic expected PSD; not field data or a random waveform",
        "receiver_psd_units": "Pa^2/Hz",
        "source_psd_units": "Pa^2 m^2/Hz",
        "spectral_convention": "one-sided, positive frequency",
        "source_reference": "loss-free monopole strength; nominal 1 m pressure equivalent",
        "clock": "receiver arrival time; retarded emission geometry",
        "doppler": "omitted",
        "quasi_stationary": True,
        "component_sum": "incoherent powers; ambient is already receiver-level, never spread",
        "wind_model": "Stojanovic 2007 Eq. 6 (Coates empirical fit)",
        "vessel_model": "MacGillivray and de Jong 2021 JOMOPANS-ECHO Eqs. 9-11, Table 1",
        "absorption_model": "Ainslie-McColm 1998; NPL domain guidance",
        "propagation": "direct spherical spreading plus full-path absorption",
        "omissions": [
            "Doppler",
            "surface and seabed",
            "refraction",
            "directionality",
            "tonals",
            "flow/self noise",
            "rain",
            "distant shipping",
            "thermal noise",
        ],
        "low_frequency_extrapolation_hz": f[f < 100].tolist(),
        "allow_low_frequency_extrapolation": bool(allow_low_frequency_extrapolation),
        "vessel_model_residual_sd_db": 6.0,
        "uncertainty_note": "population source-model residual, not a receiver confidence interval",
        "water": asdict(water),
        "parameters": {
            "speed_kn": float(speed_kn),
            "length_m": float(length_m),
            "vessel_class": vessel_class,
            "closest_approach_m": cpa,
            "source_depth_m": zs,
            "receiver_depth_m": zr,
            "sound_speed_m_s": c,
            "wind_speed_m_s": float(wind_speed_m_s),
            "closest_approach_emission_time_s": tc,
            "include_wind": bool(include_wind),
            "include_vessel": bool(include_vessel),
        },
    }
    return {
        "frequency_hz": f.copy(),
        "receiver_time_s": tr.copy(),
        "emission_time_s": te,
        "distance_m": distance,
        "travel_time_s": delay,
        "propagation_loss_db": loss,
        "source_psd_pa2_m2_hz": source,
        "wind_psd_pa2_hz": wind,
        "vessel_psd_pa2_hz": vessel,
        "total_psd_pa2_hz": wind + vessel,
        "metadata": metadata,
    }
