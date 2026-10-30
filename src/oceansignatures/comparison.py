"""Compare aligned acoustic spectra, not hydrodynamic self-noise surrogates."""

from datetime import datetime, timedelta

import numpy as np


def _validate_batch(batch):
    axes = {}
    for key in ("frequency_hz", "window_start_s", "window_stop_s", "psd"):
        raw = np.asarray(batch[key])
        if np.iscomplexobj(raw):
            raise ValueError(f"{key} must be real")
        axes[key] = np.asarray(raw, dtype=float)
    frequency = axes["frequency_hz"]
    starts, stops = axes["window_start_s"], axes["window_stop_s"]
    if (
        frequency.ndim != 1
        or frequency.size < 2
        or not np.isfinite(frequency).all()
        or frequency[0] < 0
        or not np.all(np.diff(frequency) > 0)
    ):
        raise ValueError("frequency_hz must be finite, nonnegative and strictly increasing")
    if (
        starts.ndim != 1
        or not starts.size
        or stops.shape != starts.shape
        or not np.isfinite(starts).all()
        or not np.isfinite(stops).all()
        or not np.all(np.diff(starts) > 0)
        or not np.all(stops > starts)
    ):
        raise ValueError("Require finite, ordered window starts and positive durations")
    valid = np.asarray(batch["valid_window"])
    if valid.dtype != np.dtype(bool) or valid.shape != starts.shape:
        raise ValueError("valid_window must be a boolean vector matching the windows")
    if axes["psd"].shape != (starts.size, frequency.size):
        raise ValueError("PSD shape must be (windows, frequencies)")
    if not np.isfinite(axes["psd"][valid]).all() or np.any(axes["psd"][valid] < 0):
        raise ValueError("Valid-window PSD must be finite and nonnegative")


def compare_prediction_observation(prediction, observation):
    """Return observed-minus-predicted residuals on the declared common support."""
    origin = prediction.get("time_origin_utc")
    try:
        origin_time = datetime.fromisoformat(origin)
    except (TypeError, ValueError) as exc:
        raise ValueError("Require an explicit ISO8601 UTC origin") from exc
    if origin_time.utcoffset() != timedelta(0):
        raise ValueError("Require an explicit ISO8601 UTC origin")
    if origin != observation.get("time_origin_utc"):
        raise ValueError("Both spectra must share the exact declared UTC origin")
    if prediction["psd_unit"] != "Pa2/Hz" or observation["psd_unit"] != "Pa2/Hz":
        raise ValueError("Both spectra must have explicit physical units Pa2/Hz")
    if prediction["evidence_kind"] != "predicted_acoustic" or observation["evidence_kind"] not in {
        "synthetic_observation",
        "measured_observation",
    }:
        raise ValueError("Require an acoustic prediction and an acoustic observation")
    for key in ("frequency_hz", "window_start_s", "window_stop_s"):
        if not np.array_equal(prediction[key], observation[key]):
            raise ValueError(f"Spectra must be exactly aligned in {key}; no implicit interpolation")
    _validate_batch(prediction)
    _validate_batch(observation)
    predicted = np.asarray(prediction["psd"], dtype=float)
    observed = np.asarray(observation["psd"], dtype=float)
    compared = np.asarray(prediction["valid_window"]) & np.asarray(observation["valid_window"])
    residual_db = np.full(predicted.shape, np.nan)
    linear_residual = np.full(predicted.shape, np.nan)
    positive = compared[:, None] & (predicted > 0) & (observed > 0)
    residual_db[positive] = 10 * (np.log10(observed[positive]) - np.log10(predicted[positive]))
    residual_db[compared[:, None] & (predicted == 0) & (observed > 0)] = np.inf
    residual_db[compared[:, None] & (predicted > 0) & (observed == 0)] = -np.inf
    linear_residual[compared] = observed[compared] - predicted[compared]
    values = residual_db[positive]
    return {
        "time_origin_utc": prediction["time_origin_utc"],
        "frequency_hz": np.asarray(prediction["frequency_hz"]).copy(),
        "window_start_s": np.asarray(prediction["window_start_s"]).copy(),
        "window_stop_s": np.asarray(prediction["window_stop_s"]).copy(),
        "observed_minus_predicted_db": residual_db,
        "linear_residual_pa2_hz": linear_residual,
        "compared_window": compared,
        "finite_db_comparison": positive,
        "summary": {
            "mean_bias_db": float(np.mean(values)) if values.size else None,
            "rmse_db": float(np.sqrt(np.mean(values**2))) if values.size else None,
            "compared_windows": int(np.count_nonzero(compared)),
            "excluded_windows": int(np.count_nonzero(~compared)),
            "finite_db_pairs": int(np.count_nonzero(positive)),
        },
    }


def _integrated_power(density, valid, selected, df):
    values = density[valid][:, selected]
    peak = np.max(values, axis=1)
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise", under="ignore"):
            scaled = np.divide(
                values, peak[:, None], out=np.zeros_like(values), where=peak[:, None] > 0
            )
            total = np.sum(scaled, axis=1)
            weighted_peak = peak * df
            power = weighted_peak * total
            retry = (peak > 0) & (weighted_peak == 0)
            power[retry] = (peak[retry] * total[retry]) * df
    except FloatingPointError as exc:
        raise ValueError("Integrated band power is not representable as float64") from exc
    if not np.isfinite(power).all() or np.any((peak > 0) & (power == 0)):
        raise ValueError("Integrated band power is not representable as float64")
    result = np.full(len(density), np.nan)
    result[valid] = power
    return result


def _pressure_level(power):
    level = np.full(power.shape, np.nan)
    positive = power > 0
    level[positive] = 10 * (np.log10(power[positive]) + 12)
    level[power == 0] = -np.inf
    return level


def compare_band_levels(prediction, observation, bands_hz):
    """Integrate uniform-grid PSD bins in [low, high) before computing levels."""
    comparison = compare_prediction_observation(prediction, observation)
    frequency = comparison["frequency_hz"]
    df = frequency[1] - frequency[0]
    if not np.allclose(np.diff(frequency), df, rtol=1e-12, atol=0):
        raise ValueError("Band integration requires a uniform frequency grid")
    if not isinstance(bands_hz, dict) or not bands_hz:
        raise ValueError("Specify at least one named frequency band")
    for name, edges in bands_hz.items():
        if not isinstance(name, str) or not name.strip():
            raise ValueError("Each band requires a nonempty name")
        edges = np.asarray(edges, dtype=float)
        if edges.shape != (2,) or not np.isfinite(edges).all():
            raise ValueError("Each band requires two finite frequency edges")
        low, high = edges
        if not frequency[0] <= low < high <= frequency[-1]:
            raise ValueError("Each band must lie within the sampled frequency support")
        if not np.any((frequency >= low) & (frequency < high)):
            raise ValueError("Each band must contain at least one frequency bin centre")
    predicted = np.asarray(prediction["psd"], dtype=float)
    observed = np.asarray(observation["psd"], dtype=float)
    predicted_power, observed_power = [], []
    for low, high in bands_hz.values():
        selected = (frequency >= low) & (frequency < high)
        predicted_power.append(
            _integrated_power(predicted, prediction["valid_window"], selected, df)
        )
        observed_power.append(
            _integrated_power(observed, observation["valid_window"], selected, df)
        )
    predicted_power = np.asarray(predicted_power).T
    observed_power = np.asarray(observed_power).T
    predicted_level = _pressure_level(predicted_power)
    observed_level = _pressure_level(observed_power)
    residual = np.full(predicted_power.shape, np.nan)
    comparable = comparison["compared_window"][:, None] & ~(
        (predicted_power == 0) & (observed_power == 0)
    )
    np.subtract(observed_level, predicted_level, out=residual, where=comparable)
    return {
        "time_origin_utc": comparison["time_origin_utc"],
        "window_start_s": comparison["window_start_s"],
        "window_stop_s": comparison["window_stop_s"],
        "compared_window": comparison["compared_window"],
        "band_names": list(bands_hz),
        "band_edges_hz": np.asarray(list(bands_hz.values())),
        "predicted_power_pa2": predicted_power,
        "observed_power_pa2": observed_power,
        "predicted_level_db_re_1upa": predicted_level,
        "observed_level_db_re_1upa": observed_level,
        "observed_minus_predicted_db": residual,
    }
