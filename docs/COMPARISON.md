# Prediction and observation comparison

`oceansignatures.comparison` compares acoustic pressure spectra on an explicitly shared clock and frequency support. It does not infer calibration or couple local CFD pressure into acoustics.

## Input contract

Both batch dictionaries contain:

- `time_origin_utc`: an explicit ISO8601 UTC string. Origins must be identical; matching relative seconds from different deployments is not alignment.
- `frequency_hz`: finite, nonnegative, strictly increasing frequency centres.
- `window_start_s`, `window_stop_s`: finite relative seconds, ordered starts and positive durations. Corresponding start, stop and frequency arrays must match exactly. There is no interpolation, frequency warping or implicit window averaging.
- `psd`: `[window, frequency]` one-sided pressure PSD in `Pa2/Hz`.
- `psd_unit`: exactly `Pa2/Hz`.
- `valid_window`: Boolean vector. For observations, use the recording pipeline's `qc_ok` mask, not only complete-data coverage.
- `evidence_kind`: `predicted_acoustic` for the prediction; `synthetic_observation` or `measured_observation` for the observation.

The evidence label is a declared provenance class, not authentication of a measurement. A calibrated scalar conversion is not a calibration certificate. Digital-only recordings and hydrodynamic CFD pressure are rejected. Spectral processing conventions and the meaning of the prediction must also be documented: an empirical expected PSD is not an exact finite-record Welch realization.

## Spectral residuals

`compare_prediction_observation(prediction, observation)` returns the common support, the Boolean `compared_window`, linear residual in Pa²/Hz, and

`delta_L = 10 [log10(S_observed) - log10(S_predicted)]`.

Positive values mean more observed power. Mean bias and RMSE are unweighted summaries over finite, strictly positive PSD pairs in included windows. They are not band-integrated errors or uncertainty intervals. `finite_db_comparison` and `finite_db_pairs` state exactly which pairs contribute.

Excluded windows remain present with NaN residuals. True physical zeros are not replaced by a plotting floor: positive/zero gives positive infinity, zero/positive negative infinity, and zero/zero an undefined NaN logarithmic residual. Linear residuals remain defined. No finite pairs means unavailable (`None`) aggregate errors, not zero error.

## Frequency bands

`compare_band_levels(prediction, observation, bands_hz)` requires a uniform frequency grid and a nonempty dictionary of named `(low, high)` edges. Each band must lie within the sampled frequency support and contain at least one centre. The convention is **low-inclusive, high-exclusive**, including the last requested band. This comparison API does not infer that the last supplied frequency is Nyquist. To match the recording pipeline without endpoint ambiguity, select interior bands whose upper edge is below Nyquist, such as 100–500, 500–1000 and 1000–4000 Hz for 8192 Hz audio.

`P_band = delta_f * sum(S[k], low <= f[k] < high)`

`L_band = 10 log10(P_band / 1e-12 Pa²)`

Integrate linear density before taking logarithms. Overlapping named bands are permitted but should not be summed as disjoint total power. No partial-bin weighting or trapezoidal endpoint halving is applied. Peak-scaled arithmetic avoids intermediate summation overflow; unrepresentable nonzero band powers are rejected.

Predicted and observed band powers retain their own validity: a missing observation does not erase an available prediction. Residuals require both masks. Silence, missing data and a flagged passive window remain distinct.

## Verification

Run from the repository root:

```bash
uv run pytest tests/test_comparison.py -q
```

Independent controls include a fourfold PSD ratio, exact linear residuals, mismatched frequency/time/origin rejection, prohibited evidence and digital units, malformed shapes/masks, flagged and all-excluded windows, physical zeros, band integration before dB, and representable versus unrepresentable extreme arithmetic. Synthetic known-discrepancy tests verify this implementation; they do not validate an environmental model against the ocean.
