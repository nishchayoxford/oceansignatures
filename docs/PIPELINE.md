# Recording pipeline

This executable offline path processes ordinary normalized digital audio, not an undocumented proprietary sonar format:

```bash
uv run python -m oceansignatures.pipeline recording.wav build/recording-analysis
```

It writes `summary.json`, `spectrum.csv` and `signal.npz`. Originals are never modified. Identical reruns are permitted; differing existing outputs are not overwritten. Use a fresh output directory for a new analysis.

## What is known and what is not

A decodable WAV is not automatically a calibrated pressure measurement. Without explicit calibration, outputs remain normalized digital amplitude and digital amplitude squared per Hz; the absolute pressure-level field is `null`. A file alone does not establish deployment timing continuity, sensor orientation, source distance, or label truth.

The optional `--calibration-json` file contains `sensitivity_db_v_per_upa`, `gain_db`, and `full_scale_voltage`. This first implementation uses a **declared flat scalar response**, not a frequency-dependent calibration certificate. Its equation is:

`p[Pa] = digital_fraction × V_fullscale × 10^(-gain_dB/20) × 10^(-sensitivity_dB/20) × 10^-6`.

The sensitivity is in dB re 1 V/µPa. The final factor converts µPa to Pa. `full_scale_voltage` is the **ADC input peak voltage corresponding to a normalized digital fraction of +1**; -1 corresponds to its negative. It is not a peak-to-peak or RMS voltage. Gain is removed to obtain the sensor voltage before applying sensitivity. Unknown ADC scaling or gain means this conversion must not be used. The user is responsible for the validity of the calibration over the analysed band. No device calibration or accuracy is implied by the synthetic calibration test.

### Calibration provenance

Both the analysis dictionary and exported `summary.json` retain:

- `calibration`: the actual `sensitivity_db_v_per_upa`, `gain_db`, and `full_scale_voltage` values supplied by the caller, not just a status label;
- `pressure_scale_pa_per_digital_fraction`: the combined conversion factor, such that `signal[Pa] = digital_fraction × pressure_scale_pa_per_digital_fraction`;
- `calibration_assumption`: `flat scalar response over the analysed band`;
- `adc_full_scale_convention`: the explicit ±1/±full-scale peak-voltage convention above.

For example, sensitivity -180 dB re 1 V/µPa, gain 20 dB and full-scale 2 V give a scale of 200 Pa per digital fraction. A 0.2-peak digital sine then has 40 Pa peak pressure. All four provenance fields are `null` when calibration is absent; `calibration_status` remains `unknown`, the PSD remains digital amplitude squared per Hz, and `spl_db_re_1upa` remains `null`. No physical pressure or response assumption is inferred from the recording.

The additive metadata makes new summaries differ from summaries produced before these fields existed. As with any changed analysis, select a fresh output directory; overwrite protection is unchanged.

## Quality controls

The implementation rejects oversized, empty, non-finite and out-of-range normalized recordings, and invalid channel indices. Near-full-scale samples are counted as potential saturation indicators, not proof of actual hardware clipping. The original content hash and explicit channel/sample count are retained. File-size and decoded-sample bounds prevent accidental ingestion of a large deployment in this initial batch tool; deployment-scale streaming is future work.

Finite calibration coefficients alone do not guarantee representable numerical products. Conversion checks that the combined Pa/digital-fraction scale is positive and finite, then checks converted samples. Analysis separately checks each squared signal value and the final Welch density in its original squared units. Overflow, nonfinite products, and nonzero converted samples/squares/rescaled density bins that round to zero are rejected with a stage-specific `ValueError`, without emitting `RuntimeWarning`s or returning nonfinite metrics. For example, `Calibration(-600, -600, 1e200)` can yield finite pressure samples, but analysis of a nonzero ordinary-amplitude tone rejects its unrepresentable pressure squares. A stable RMS calculation alone would not make that PSD valid.

RMS and Welch averaging use a peak-normalized working signal to avoid overflowing intermediate sums when the final products are representable. The waveform and exported PSD retain their original amplitude and squared-amplitude-per-Hz units; this working normalization is not a change of calibration. Silence yields zero RMS and PSD, with a `null` pressure level rather than negative infinity.

## Spectra

The pipeline uses one-sided Welch PSD with a periodic Hann window, 50% overlap, no detrending, and density scaling.[1] Direct waveform RMS and an integrated Welch estimate are distinct reductions: window weighting, partial segments and non-stationarity can make them differ. Never relabel a spectrum-amplitude estimate as a PSD. The spectral CSV header includes the applicable units.

## Independent checks

The tests exercise known-amplitude tones, a hand-derived voltage/gain/sensitivity conversion, calibration provenance in memory and JSON, missing-calibration handling, independent odd/even Hann DFT densities, file hashing, channel selection, potential saturation, invalid inputs, reproducible export and the real subprocess CLI. Extreme finite coefficients test combined-scale, pressure-square and PSD overflow/underflow, rejection without warnings, and stable representable large products. Silence and the existing no-overwrite behavior remain covered. These are software/numerical verification, not validation against a commercial instrument.

## Sources

[1] https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.welch.html
