# Reproducible acoustic prediction → synthetic recording

This offline extension connects the existing empirical acoustic predictor, scalar
calibration, timestamped WAV analysis and comparison APIs. It produces
**predictive approximations + synthetic verification**, never measured ocean
noise or field validation. It does not require hardware, fetch data, bundle NOAA
recordings, change the existing teaching assets, or couple hydrodynamic CFD
pressure into acoustic pressure.

## Run from the repository root

Use the project's existing Python environment (dependencies are already declared
in `pyproject.toml`):

```bash
.venv/bin/python scripts/build_research.py --output outputs/research-run001 \
  --seed 1729 --gain-db 6
.venv/bin/python scripts/build_research.py --output outputs/research-run001 \
  --seed 1729 --gain-db 6   # identical replay: accepted, nothing replaced
.venv/bin/pytest tests/test_research.py tests/test_noise_prediction.py \
  tests/test_observations.py tests/test_comparison.py -q
.venv/bin/ruff check src/oceansignatures/research.py scripts/build_research.py tests/test_research.py
.venv/bin/ruff format --check src/oceansignatures/research.py scripts/build_research.py tests/test_research.py
```

The destination must be **new**: do not pre-create it. A different seed, gain,
producer implementation or software-stack provenance requires another directory.
The entire family is staged next to the destination. Exactly matching existing
files/directories are accepted after regeneration and hashing; differing,
incomplete or extra output is rejected without overwriting it. Symlinks and
`..` components are rejected. This is a single-writer offline contract, not
protection against adversarial concurrent filesystem renames. No publication or
Git action occurs.

The CLI prints strict JSON with the actual file/figure/download counts, manifest
SHA256, destination and computed numerical summary. Determinism is tested both
across different output roots and across separate CLI processes. Bitwise replay
is conditioned on the producer and numerical/rendering software versions, which
are retained; cross-platform font/library identity is not promised.

## One complete research loop

1. **Predict** the source spectrum, direct-path vessel spectrum and ambient wind
   spectrum. Inspect `figures/spectra.png` and `figures/propagation.png`.
2. **Choose a synthetic ensemble** whose frame-local Fourier powers follow those
   received spectra. Inspect `figures/passage.png` for the expected finite-window
   spectrum and time histories.
3. **Record the surrogate** in deterministic, two-channel IEEE float64 WAV chunks,
   apply a declared scalar ADC conversion, and analyse with the existing recording
   kernel. One chunk is deliberately absent and two frames are active-contaminated.
4. **Compare two different questions** in `figures/comparison.png` and
   `figures/residuals.png`: does the software recover a known same-realization
   gain; and how far does one realization differ from its ensemble expectation?
5. **Change one assumption at a time** in `figures/sensitivity.png`. These changes
   are alternatives, not confidence intervals or accuracy bounds.

### Fixed scenario and clocks

| Quantity | Declared value |
|---|---|
| Duration / frame clock | 240 s / 120 half-open 2 s frames |
| Sampling | 8192 samples/s; 16384 samples/frame |
| Spectral kernel | 4096-sample periodic Hann; 2048 overlap; no detrend/padding; arithmetic-mean density |
| Welch segments | Seven complete segments/frame; none crosses a frame boundary |
| Comparison support | 100–4000 Hz at 2 Hz spacing; no low-frequency extrapolation |
| Integrated bands | `[100,500)`, `[500,1000)`, `[1000,4000)` Hz |
| Vessel | 150 m bulker at 12 kn; stationary water assumed |
| Track | 250 m horizontal CPA; source depth 8 m; receiver depth 50 m |
| CPA clock | Emission second 120, **not** receiver arrival second 120 |
| Sound speed | Constant 1500 m/s |
| Absorption water | 10 °C, 35 ppt salinity, pressure depth 50 m, pH 8 |
| Wind | Chosen linear receiver-clock schedule 4→8 m/s over 0–240 s; frozen at each frame midpoint |
| Synthetic UTC | `2000-01-01T00:00:00.000000Z`: arbitrary epoch, not an acquisition date |
| Random generator | NumPy PCG64; seed 1729 by default |
| Injected discrepancy | +6 dB pressure amplitude gain, default |

Retarded geometry solves `arrival = emission + slant_range / c`. Water absorption
depth is a pressure-property parameter, independent of receiver depth. The wind
term is already at the receiver and is not spread or absorbed a second time.
Source units are Pa² m²/Hz; received units are Pa²/Hz. Spectral levels use
`10 log10(PSD) + 120`, with the m² retained in the source reference. Band pressure
levels integrate linear PSD first, then use dB re 1 µPa.

The predictor produces instantaneous expected PSD snapshots. This experiment
freezes each frame's spectrum at its midpoint. The summary quantifies the largest
absolute dB change between that midpoint and the two frame edges over all
supported bins. That is an endpoint stationarity diagnostic, not an exact
continuous-time error bound. The waveform does not implement moving-source
Doppler or a phase-resolved propagation delay; it is a receiver-spectrum surrogate.

## Synthesis and finite-window expectation

A frame has `N=16384` samples and a 0.5 Hz Fourier grid. Evaluate the empirical
receiver PSD `S1[k]` at every frame-grid frequency from 100 through 4000 Hz,
including both endpoint bins; set every other coefficient's power to zero. For
positive interior bins use independent standard normal `a[k]`, `b[k]`:

```text
X[k] = sqrt(N * fs * S1[k] / 4) * (a[k] + i*b[k])
p = irfft(X, n=N)
```

DC and Nyquist are zero. The full arrays of real normal deviates are drawn first,
then the imaginary arrays, with fixed frame/frequency order. Each frame is a
stationary **circular Gaussian** ensemble with the declared expected discrete
Fourier power. There is **no finite-realization RMS normalization**.

Frames are independent and concatenated abruptly, without fades. Individual
frames have bandlimited Fourier support, but their concatenation has seams and
is **not globally bandlimited**. This is intentional and disclosed: every Welch
segment remains wholly inside one frame, so no analysis segment straddles a seam.
Using a different analysis window that crosses seams changes the expected result.
Chunk boundaries, in contrast, do not affect synthesis or reset the analysis
clock.

Raw model PSD is not equated to a finite Hann estimate. Let `S2[m]` be the
corresponding two-sided frame PSD (half the one-sided positive-bin density,
mirrored; DC/Nyquist not halved), and let `w` be the 4096-point periodic Hann
zero-padded to N. Define

```text
K[m] = abs(FFT_N(w)[m])² / (N * sum(w²))
E[P_two_sided[k]] = circular_convolution(S2, K)[k * N/4096]
```

Restore the one-sided factor of two for interior output bins only. This is the
exact mean for this finite Fourier ensemble, up to floating-point arithmetic.
It incorporates support-edge smoothing/leakage, including density outside the
synthesis band. Segment overlap affects variance but not that mean. A direct
covariance/DFT oracle in the tests independently verifies this formula; the full
recording test independently verifies PSDs with SciPy Welch.

`raw_expected_psd_pa2_hz` and `expected_welch_psd_pa2_hz` are exported separately.
`welch_smoothing_max_abs_db` quantifies their maximum difference on the supported
comparison bins. The active signal described below is not part of the empirical
expected PSD; those frames are excluded from every residual statistic.

## Synthetic acquisition, calibration and QC

Each family has eight retained WAV chunks. The intended boundary seconds are
`0, 17.375, 53.125, 80, 82, 109.625, 150.5, 190.875, 216.25, 240`.
The `[80,82)` chunk is absent. Other boundaries are deliberately nonaligned to
2 s frames, and a chunk boundary falls inside the active interval. WAV channel
**1** (zero-based) holds the pressure-derived digital signal. Channel 0 is a
separate 333 Hz, 0.05 digital-amplitude selection distractor, not another sensor
measurement. No channel averaging occurs. Float64 WAV avoids PCM quantization
obscuring the gain-control identity. These files are numerical synthetic
recordings, not gain-normalized listening products.

The same baseline includes a 1200 Hz, 0.3 Pa peak sine over
`[150.25,153.125)` s. It is a controlled active-contamination surrogate, not an
instrument transmission model. The observed family is **the same realization**
multiplied by `10**(gain_db/20)`, including this contamination. Both are converted
to normalized digital fractions with the declared `Calibration`:

```text
sensitivity_db_v_per_upa = -180
instrument gain_db       = 20
full_scale_voltage      = 2 V peak
p[Pa] = digital * 2 * 10**(-(-180+20)/20) * 1e-6
      = digital * 200
```

The injected +6 dB discrepancy is distinct from the declared instrument-chain
20 dB gain. It is a synthetic amplitude multiplier, not a detected hardware fault.
`Calibration` and `calibrate_samples` are reused, not reimplemented. Metadata
retains each coefficient, Pa/digital-fraction scale, flat-response assumption and
ADC ±full-scale **peak** convention. This is not a calibration certificate.

`analyse_chunks` processes the baseline, observed and digital-only control. The
latter uses the same observed WAVs with calibration explicitly absent and has
only digital²/Hz, not pressure or absolute SPL. SHA256 pins are checked by the
recording API. Its internally returned provenance keeps absolute input paths;
export serialization deep-copies that metadata and replaces only chunk paths
with export-relative paths. No private input/root paths are exported.

- Frame index 40, `[80,82)`, has zero coverage and NaN PSD/power, never zeros.
- Indices 75 and 76, `[150,152)` and `[152,154)`, retain complete spectra but are
  active-contaminated and excluded.
- The passive comparison uses **qc_ok**, not merely complete-window `valid`.
- This yields 117 comparable windows and three retained excluded windows.
- Full-scale contacts are checked on the original selected digital channel.
- Each recording analysis allocates 1,967,040 PSD+STFT cells, below its 2,000,000
  cell hard cap. The export retains PSDs/QC/bands, not the large complex STFTs;
  those remain reproducible from the supplied WAVs and processing metadata.

## Three distinct residuals

All comparisons call `compare_prediction_observation` and `compare_band_levels`.
The adapter requires calibrated physical units, maps to `psd_unit="Pa2/Hz"`, uses
`qc_ok`, and passes exact frequency bins, half-open windows and the explicit UTC
origin. The comparison API rejects mismatches; no interpolation, resampling,
origin guessing or unit inference is performed. All bands avoid Nyquist so the
observation and comparison high-exclusive conventions agree.

1. **Matched control:** observed minus the same unscaled realization processed by
   the **same** chunk/calibration/Hann kernel. The expected residual is precisely
   the injected gain, within floating-point tolerance. This verifies a software
   identity, not the empirical source law.
2. **Sampling:** baseline realization minus the ensemble-expected finite-Hann
   PSD. These spectral and band residuals are not forced to zero. Mean and RMSE
   describe this particular seeded sample, not confidence coverage.
3. **Empirical comparison:** observed minus the expected finite-Hann PSD. It
   contains the injected gain **plus** sampling discrepancy. Mean logarithmic
   spectral bias can differ from the gain because logarithm and expectation do
   not commute; the broadband metric integrates powers first.

The synthetic matched control uses the comparison API's `predicted_acoustic`
role solely as an explicitly labelled finite-realization predictor. It must not
be mistaken for the independent empirical prediction. The physical observation
role is `synthetic_observation`; no measured observations are generated.

## One-factor sensitivity, not probabilistic uncertainty

Nineteen cases are exported: baseline plus two alternatives for each factor:

- Vessel source offset −6/+6 dB, leaving ambient unchanged.
- Wind-schedule offset −2/+2 m/s.
- Vessel speed 10/14 kn, affecting both source law and trajectory.
- Horizontal CPA 125/500 m.
- Receiver depth 25/100 m, leaving absorption pressure depth fixed.
- Absorption temperature 5/15 °C, salinity 30/40 ppt, pressure depth 0/100 m,
  and pH 7.8/8.2, each separately.

These summaries use **raw empirical midpoint PSD**, rectangular 2 Hz bin sums
and the time mean of linear band powers before taking a ratio in dB. They do not
use stochastic realizations or imply probabilistic confidence intervals. The
source model's cited population residual scale is not a receiver accuracy bound.
The numerical values remain full precision in NPZ even where a small effect is
rounded for a figure annotation.

## Public Python API

```python
from oceansignatures.research import (
    export_research, comparison_batch, expected_welch_psd,
)

overview = export_research("outputs/research-run001", seed=1729, gain_db=6.0)
batch = comparison_batch(calibrated_analysis, as_prediction=False)
# Only a synthetic matched realization may use as_prediction=True.
expected_psd = expected_welch_psd(full_frame_psd, fs=8192, nperseg=4096)
```

- `export_research(output_directory, *, seed=1729, gain_db=6.0)` returns the
  JSON-safe overview dictionary. Seed is a Python integer in `[0,2**32-1]`;
  gain is a finite real scalar within ±12 dB. Scenario/estimator choices above
  are deliberately fixed. Existing differing/incomplete output raises
  `FileExistsError`; invalid input/path raises `ValueError`.
- `comparison_batch(analysis, *, as_prediction=False)` takes a calibrated
  `analyse_chunks` result, copies its selected 100–4000 Hz support/QC and returns
  the exact comparison contract. It can map declared measurement evidence for
  callers, but this research export uses only synthetic evidence. Unknown
  calibration/units are refused; a measurement cannot become the matched
  synthetic prediction control.
- `expected_welch_psd(frame_psd, *, fs=8192, nperseg=4096)` returns the expected
  one-sided Hann density. The final input axis is a complete even-frame rFFT
  grid, with finite nonnegative densities. Segment length must be even and
  divide the implied frame length. No detrending is supported by this formula.
  Each input spectrum is normalized by its own peak before one-/two-sided
  conversion and FFT convolution, then restored after applying the output
  one-sided weights. This avoids overflow for large finite densities and premature
  halving-to-zero of subnormal inputs; zero spectra remain zero. Nonfinite restored
  density or restoration of a computed positive density to zero raises `ValueError`
  rather than returning NaN or false silence. Representable subnormals are allowed.
  The FFT convolution has absolute roundoff relative to the spectrum's scale;
  arbitrarily weak out-of-support leakage does not have uniform relative accuracy.
  Normalization does not remove this existing floating-point limitation.

## Output schema and numerical downloads

All asset `path` values in the overview are relative to the export directory.
The overview is a small strict-JSON descriptor, not embedded dense spectral data:

```text
schema: "oceansignatures.research.overview.v1"
title: string
kind: "predictive approximations + synthetic verification"
configuration: declared scenario, calibration, estimator, seed, gain and band_names
numericalsummary: actual counts, gain error, sampling statistics, range/delay,
                  band-level extrema, finite-window and midpoint diagnostics
limitations: [string, ...]
figures: [{path, title, caption}, ...]
downloads: [{path, label}, ...]
source_refs: [{label, url}, ...]
```

| Product | Meaning |
|---|---|
| `overview.json` | Description, configuration, computed summary and relative gallery links |
| `metadata.json` | Synthetic declaration, calibration/clock/QC rules, per-chunk SHA256, synthesis/expectation contract, library versions and repository-relative source-code hashes |
| `research.npz` | 120-frame physical/digital PSDs, empirical components, matched/sampling/empirical residuals, expected and actual band powers, geometry, clocks and QC |
| `windows.csv` | 360 window/band rows, QC, retained actual/expected powers and the three residuals; missing/excluded residuals blank |
| `sensitivity.npz` | Case × frame × band linear powers; case × band change of time-mean power |
| `sensitivity.json` | Ordered cases, one changed factor per alternative, band names and integration convention |
| `wav/baseline/*.wav`, `wav/observed/*.wav` | Synthetic two-channel float64 chunks; no absent-chunk zero-fill file |
| `figures/*.png` | Six computed Matplotlib figures |
| `manifest.json` | SHA256 for every exported file except the manifest itself; the CLI separately hashes the manifest |

In `research.npz`, `frequency_hz` has 2049 bins (0–4096 Hz) and is the axis for
baseline/observed/digital/raw-expected/expected-Welch PSD arrays. The
`comparison_frequency_hz` axis has 1951 bins (100–4000 Hz) for source/received
components, loss and spectral residual arrays. Time-dependent arrays have 120
rows. `band_names` explicitly fixes band-column order; never infer column order
from sorted JSON object keys. `source_psd_pa2_m2_hz` is frequency-only.
`receiver_time_s` contains midpoints; `window_start_s`/`window_stop_s` give exact
support; emission seconds and delay retain the retarded-geometry clock.

Missing data remain NaN in NPZ and blank in CSV; strict JSON serializes nonfinite
values as null where needed. Complete contaminated PSD/band values are retained
alongside false `qc_ok`; only comparison residuals are blank/NaN. Neither plot
line interpolation across excluded frames nor zero filling is used. Source-code
hash paths identify repository-relative provenance, not output download links.

## Executed default result

The accepted execution is available in `site/data/research/acoustic/`.
For seed 1729 and injected +6 dB it contains 29 files, including 16 synthetic
WAVs and six figures. Two fresh independent CLI runs reproduced every output
byte, including the manifest. The 120 windows retain one missing frame and two
active frames; 117 windows are compared. Earlier executions remain preserved
as historical evidence rather than being relabelled with the new source identity.

| Computed diagnostic | Value |
|---|---:|
| Same-realization gain bias | 6.000000000000000 dB |
| Maximum spectral gain recovery error | 4.085620730620576e-14 dB |
| Baseline-vs-expected 100–4000 Hz band RMSE | 0.19905314406605468 dB |
| Baseline-vs-expected spectral mean / RMSE | −0.33078428979421526 / 1.7690791731501623 dB |
| Observed-vs-expected spectral mean | 5.669215710205785 dB |
| Largest supported-bin finite-Hann smoothing difference | 2.3706485874557592 dB |
| Largest midpoint-to-frame-edge PSD difference | 0.10605614928233668 dB |
| Minimum sampled slant range | 253.5553489344946 m |
| Sampled travel-time span | 0.16903689928966306–0.5201145973711476 s |

The three sampling band RMSEs are 0.2476561516891248,
0.16616602425248697 and 0.07984742193231777 dB in the explicit `band_names`
order. These values quantify this seeded verification experiment, not ocean-model
accuracy. Research tests include public CLI replay and overwrite-refusal paths.
Fresh independent acoustic acceptance passed 240 declared tests and 621 new
independent probes, including exact silence, weak-bin retention, supported
normalization extremes and representability of the final mean PSD. These are
software checks, not field validation. All 16 synthetic WAVs retain their earlier
bytes; processed roundoff and source provenance are identified separately.
See [the verification report](VERIFICATION.md) for accepted bundle identities
and the scope of the whole-project and browser checks.

## Evidence and limits

Tests exercise the whole WAV→calibration→PSD path against SciPy, a direct
covariance/DFT expected-spectrum oracle, independent range/arrival equations,
source-only sensitivity scaling, exact-support refusal, missing/active QC,
figure line support, CSV/NPZ readback, root-independent deterministic replay,
SHA256 manifests, non-overwrite behavior and the public CLI. These are numerical
and integration tests. No hardware calibration, field-data validation, realistic
waveguide model, hydrodynamic source coupling or source-spectrum accuracy is
established.

The empirical equations and their scope are documented in
[NOISE_PREDICTION.md](NOISE_PREDICTION.md), and the authoritative recording-kernel
conventions in [OBSERVATIONS.md](OBSERVATIONS.md). The export reuses the cited
Stojanovic wind fit, MacGillivray/de Jong JOMOPANS-ECHO source law and
Ainslie–McColm/NPL absorption references from those documents; it makes no new
empirical validity claims.
