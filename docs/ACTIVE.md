# Active echo: a stationary point-target benchmark

This is an analytic, synthetic laboratory example—not a model of a real animal,
vessel, transducer or field detector. It illustrates a known finite ping,
two-way propagation, scalar scattering, additive noise, and matched-filter range
estimation. No classification is performed. No acoustic hardware is used.

## API and reproducibility

```python
from oceansignatures.active import active_case

case = active_case(range_m=150.0, snr_db=-12.0, fs=48000, seed=2026)
```

Inputs have explicit benchmark limits:

- `range_m`: finite real number from 1 to 750 m, inclusive.
- `snr_db`: finite real number from −120 to +120 dB, inclusive.
- `fs`: integer sample rate from 24,000 to 384,000 Hz, inclusive.
- `seed`: nonnegative integer, passed to NumPy's local `default_rng`.

Invalid inputs raise `ValueError` naming the argument. Strings and booleans are
not silently converted. Numerical limits avoid the spreading singularity,
keep the complete echo in the fixed acquisition window, and bound benchmark
resources. A finite taper is not strictly bandlimited; the 24 kHz minimum puts
the nominal 9 kHz chirp endpoint below Nyquist, not all sidelobe energy.

| Return key | Meaning |
|---|---|
| `time_s` | Common pressure-record times, `arange(N)/fs`, starting with transmission at zero |
| `transmitted_pa` | Sampled 1 m reference ping, followed by zeros; not the local transmit leakage at the receiver |
| `received_pa` | Clean scalar echo plus the actual Gaussian noise realization, Pa |
| `received_clean_pa` | Additional diagnostic: full clean echo, Pa |
| `noise_pa` | Additional diagnostic: exact additive noise realization, Pa |
| `lag_s` | Valid-correlation template-start lags, `arange(N-M+1)/fs` |
| `matched_filter` | **Signed** normalized correlation; dimensionless noise-standard-deviation units |
| `estimated_range_m` | Range at the global maximum of `abs(matched_filter)` |
| `true_range_m` | Input horizontal separation |
| `fs` | Sample rate, Hz |
| `metadata` | JSON-compatible computation-derived scalars, support bounds and model assumptions |

All pressure arrays have length `N`. The correlation and lag arrays have length
`N-M+1`. The template length is `M = ceil(0.05*fs)` and
`N = fs + M`; thus the search starts at lag zero and ends at exactly 1 s,
corresponding to 0–750 m. This search interval **does not move around the known
target range**. The last pressure time is `(N-1)/fs`, not `N/fs`.

At 48 kHz the ping has 2400 samples. Metadata reports the actual ping/record
sample counts, ping RMS, chirp parameters, target and propagation gains, sound
speed, expected delay, echo-support sample interval, clean/noise mean squares,
requested and realized SNR, peak index and magnitude, grid spacing, range
error, lag convention, fractional-delay method and omitted physics.

## Signal and propagation

Let `T = 0.05 s`, `f0 = 3000 Hz`, `f1 = 9000 Hz`, and
`b = (f1-f0)/T`. The unscaled real ping is

```text
q(t) = sin²(pi*t/T) * cos(2*pi*(f0*t + b*t²/2)),  0 <= t < T
q(t) = 0,                                       otherwise
s(t) = A*q(t)
A = 1 Pa / sqrt(mean(q(n/fs)²)),  n = 0,...,M-1
```

The continuous Hann taper is evaluated directly, not a library window with an
ambiguous endpoint convention. The sample RMS over the finite **ping window**
is 1 Pa at the assumed 1 m reference; RMS over the padded acquisition record
is smaller. The declared source normalization is synthetic, not a calibrated
source measurement. The nominal time-bandwidth product is 300.

A stationary monostatic transmitter/receiver and point target are at equal
depth in homogeneous 3-D free space. Therefore horizontal separation `R` also
equals the one-way path length. Sound speed is `c = 1500 m/s`:

```text
tau = 2*R/c
g = target_scalar_gain * (reference_distance/R)²
clean_echo(t) = g*s(t-tau)
```

The reference distance is 1 m and the **chosen scalar target gain is 1**.
Each leg has spherical pressure spreading proportional to `1/R`; multiplying
the two legs gives the chosen dimensionless two-way pressure gain. This is a
normalization convention, **not biological or vessel target strength**, a
scattering cross section or an inferred hardware specification. No surface or
seabed multipath, refraction, absorption, Doppler, extended-target scattering,
frequency-dependent scattering or reverberation is included. The receiver is
an ideal pressure sampler with no filter, clipping or transmit leakage.

### Fractional delay and range grid

The implementation evaluates the analytic continuous ping at
`(sample_index - tau*fs)/fs`. It neither rounds `tau` to a sample nor linearly
interpolates an already sampled ping. It is consequently a solution specific
to this known waveform, not a general-purpose bandlimited delay operator.
Support is evaluated in sample coordinates so floating-point subtraction at
an exactly aligned echo endpoint does not add an extra zero sample to the
SNR averaging window.

The range estimator, in contrast, uses **integer correlation lags** without
subsample fitting. At 48 kHz its range-bin spacing is `c/(2*fs) = 0.015625 m`.
A high-SNR correct-main-lobe estimate can have a half-bin quantization error;
this is not an error bound in noise. Noise can select a neighboring carrier
lobe or a completely unrelated record location. The real chirp's signed
correlation has carrier oscillations: maximizing its absolute value is not
Hilbert-envelope processing.

An independent analytic test at `R = 123.4567 m` and 80 dB input SNR returns
123.453125 m at both 48 and 96 kHz, error −0.003575 m. The corresponding
half-bin widths are 0.0078125 and 0.00390625 m. Tests also sample 31 fractional
sample positions at 48 kHz and high SNR. These are numerical checks of the
ideal model, not centimetric field-accuracy claims. The familiar nominal
`c/(2*(f1-f0)) = 0.125 m` bandwidth scale is **not** the sample grid or a
measured two-target resolution; the Hann taper broadens the ambiguity peak.

## Input SNR and seeded noise

Define the echo-support indices by `0 <= (n-tau*fs)/fs < T`. The definition is
independent of the padded record length:

```text
P_echo = mean(clean_echo_pa[echo_support]²)
input_SNR_dB = 10*log10(P_echo / sigma²)
sigma = sqrt(P_echo / 10**(snr_db/10))
noise[n] ~ iid Normal(0, sigma²)
received_pa = received_clean_pa + noise_pa
```

Noise covers the whole acquisition record. It is drawn once from
`default_rng(seed)`; the global RNG is untouched. A realization is **not**
rescaled to have exactly the requested mean square. Therefore the finite
sample `realized_input_snr_db` uses the measured noise mean square on the same
echo support and differs slightly from the requested/ensemble SNR. Both are
reported. These mean squares are not detrended sample variances.

This is **per-sample full-Nyquist-band white-noise SNR**, not SNR inside only
the chirp band, not peak-pulse SNR, not full-record SNR, and not an output
SNR. Noise variance is recomputed for each target range and requested SNR.
Consequently fixed-SNR range comparisons are controlled illustrations, **not
range-performance curves for fixed ambient noise**. Likewise, changing sample
rate at fixed input SNR changes noise power spectral density and sample count;
it is not a fixed physical noise-spectrum experiment.

## Correlation, normalization and estimator

For the unpadded sampled template `s[m]`, length `M`, and received record `y[n]`:

```text
E = sum(s[m]²)
C[k] = sum(y[k+m] * s[m]),  m = 0,...,M-1
z[k] = C[k] / (sigma*sqrt(E)),  k = 0,...,N-M
lag[k] = k/fs
k_hat = first index maximizing abs(z[k])
R_hat = c*lag[k_hat]/2
```

The implementation uses SciPy FFT `correlate(..., mode="valid")`; tests check
selected values against independently computed direct dot products.
`matched_filter` is `z`, not `abs(z)`. Plot its absolute value explicitly if a
magnitude curve is wanted. Lag zero aligns the **start** of the ping with the
start of the received record; there is no causal-filter `M-1` offset to
subtract. The full correlation is not returned and no negative/partial-overlap
lags are searched.

With noise only, each signed `z[k]` has ensemble standard deviation 1, but
neighboring lags are correlated. The score uses the **known simulation noise
standard deviation**, not an independently estimated field noise level.
The maximum of many lag scores is not a standard-normal variable or a
calibrated detection statistic. At an aligned integer delay, its expected
signed signal contribution is `sqrt(M*10**(snr_db/10))`; at −18 dB and 48 kHz
this is approximately 6.16745. That number is neither an observed output SNR
nor a false-alarm-controlled significance level.

The algorithm always reports a range, even if noise wins. There is no
threshold, no target-absent hypothesis test, and no detection probability or
false-alarm claim.

## Executed low-SNR examples

The prescribed six cases below use 48 kHz and seed 2026. These are actual
executed outputs, not labels attached to unrelated arrays. Detailed values
can be regenerated with `uv run python scripts/build_site.py --output build/site`.

| True range (m) | Input SNR (dB) | Expected delay (s) | Estimated range (m) | Error (m) | Peak absolute score |
|---:|---:|---:|---:|---:|---:|
| 75 | −18 | 0.1 | 75.0000 | 0 | 6.42324 |
| 75 | −6 | 0.1 | 75.0000 | 0 | 24.80885 |
| 150 | −18 | 0.2 | 149.9375 | −0.0625 | 6.27642 |
| 150 | −6 | 0.2 | 150.0000 | 0 | 24.10110 |
| 300 | −18 | 0.4 | 300.0000 | 0 | 5.73736 |
| 300 | −6 | 0.4 | 300.0000 | 0 | 24.12297 |

The 150 m, −18 dB example intentionally retains a nearby carrier-lobe error;
it is not corrected using the known range.

### Small repeated-noise illustration—not a detection study

A deterministic verification test runs 24 consecutively seeded cases
(seeds 2026 through 2049), with a target present at 150 m, at each of three
SNRs. The same seed set is reused across SNRs, so comparisons share the same
underlying Gaussian draws scaled to different variances; the SNR groups are
not independent experiments. All trials search the same 0–750 m interval. No trials are rejected and
no seeds are selected for a favorable result. The criterion is absolute
range error at most 0.125 m, chosen as an illustrative nominal bandwidth scale,
not a validated detection threshold or measured two-target resolution.

| Input SNR (dB) | Trials within 0.125 m | RMSE over **all** trials (m) | Largest absolute error (m) |
|---:|---:|---:|---:|
| −30 | 0/24 | 298.716257 | 582.515625 |
| −18 | 23/24 | 56.478607 | 276.687500 |
| −6 | 24/24 | 0 | 0 |

The catastrophic −18 dB outlier is included in RMSE. Reporting only the median
(which is zero there) would conceal it. These small fixed-seed counts are
illustrations of estimator behavior conditional on a target being present,
not accuracy estimates for a population, classification, deployment, or
field detection. The exactly zero high-SNR errors reflect a stationary target
on the chosen grid and ideal matched waveform. This is not evidence of perfect
ranging. No public Monte Carlo API is added; a short ordinary Python loop is
sufficient to reproduce the illustration:

```python
import numpy as np
from oceansignatures.active import active_case

for snr_db in (-30, -18, -6):
    errors = np.array([
        active_case(150, snr_db, seed=seed)["estimated_range_m"] - 150
        for seed in range(2026, 2050)
    ])
    print(snr_db, np.count_nonzero(np.abs(errors) <= 0.125),
          np.sqrt(np.mean(errors**2)), np.max(np.abs(errors)))
```

## Verification and TDD evidence

From the repository root:

```sh
uv run pytest -c pyproject.toml tests/test_active.py -v
uv run ruff check src/oceansignatures/active.py tests/test_active.py
```

The active subsystem has 39 tests. This scoped suite is separate from
full-project integration and browser acceptance.
Checks include analytic delay and waveform evaluation, direct-correlation
agreement, seeded-noise normalization, reproducibility without global RNG
mutation, complete support at domain boundaries, invalid-input validation,
the six required low-SNR cases, fractional-delay checks and the small
repeated-noise illustration.

The public tests independently reconstruct the analytic waveform and selected
correlation values. Echo support is calculated in sample coordinates to avoid
including an extra endpoint in the SNR averaging interval. The site builder
exports every repeated-noise trial, including distant false estimates.
