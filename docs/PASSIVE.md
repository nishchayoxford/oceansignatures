# Passive source, channel and receiver laboratory

This is a standalone, analytical underwater-acoustics teaching model. It is not a species classifier, machinery model, measured hydrophone, CFD calculation, BELLHOP calculation, or experimental validation. All synthetic pressure histories are in Pa. The recorded-reference pathway, if used elsewhere in the project, has a separate calibration status and must not inherit these absolute-pressure units.

## Exact public API

```python
# oceansignatures.passive
make_source(kind, fs=48000, duration_s=2.0, rms_pa=1.0, seed=2026)
path_parameters(range_m, source_depth_m=8., receiver_depth_m=12.,
                sound_speed=1500., surface=True)
propagate(signal, fs, paths)
passive_case(kind, range_m, surface, fs=48000, duration_s=2., seed=2026)

# oceansignatures.metrics
rms(signal)
welch_psd(signal, fs, nperseg=2048)
spectrogram_psd(signal, fs, nperseg=512, noverlap=384)
```

`make_source` returns a 1-D float array with `floor(fs*duration_s)` samples; `kind` is `whistle`, `clicks`, or `vessel`. These fixed-frequency signatures require `fs >= 32000 Hz` and `duration_s >= 0.1 s`. Seeds must be nonnegative integers. No global RNG is used. A zero target RMS gives silence; negative or nonfinite target RMS is rejected.

`path_parameters` returns a list of dictionaries containing `distance_m`, `delay_s`, `gain`, and `kind` (`direct` / `surface`). Horizontal range must be positive; depths are nonnegative, positive downward; sound speed is positive. `surface` must be Boolean. `propagate` requires a nonempty finite real 1-D signal and at least one path with finite gain and nonnegative delay. It does not mutate inputs.

`passive_case` returns `time_s`, `source_pa`, `received_clean_pa`, `received_pa`, `fs`, and `metadata`. Two additional same-clock arrays provide auditability:

- `received_channel_pa`: noise-free channel pressure **before** the receiver filter, zero-padded to the common record.
- `noise_pa`: the actual filtered noise realization, not an inferred residual estimate.
- `received_clean_pa`: noise-free pressure **after** the causal receiver filter, with its entire tail.
- `received_pa = received_clean_pa + noise_pa`.
- `source_pa`: emission-clock source waveform followed by zeros; it is not shifted to the receiver's arrival clock.

All arrays have the same length. `metadata` is JSON-safe and includes actual paths, seed, emission length, source RMS convention, physical omissions, filter coefficients and delay, input-noise standard deviation, analysis sample/time windows, and recomputed measured quantities. Use `metadata['source']['emission_samples']` to recover the unpadded source. Do not interpret the RMS of the zero-padded source as its reference RMS.

## Source phenomenology and reference convention

The source waveform is an amplitude convention at nominal **1 m** in an ideal homogeneous 3-D medium. Its time origin is emission time, not a pre-delayed recording at 1 m. The entire unpadded emission record is tapered first, then scaled so `sqrt(mean(p**2)) == rms_pa`. DC is not subtracted by RMS. For clicks this convention includes the quiet intervals; it is not pulse-active RMS or peak pressure. It is not a claim about natural source levels, equal audibility, or equal peak pressure. Default source RMS is 1 Pa. No per-case received-signal or listening normalization is applied.

All sources have a Tukey gate with `alpha=min(1, 0.04/T)`, approximately 20 ms at each edge, and exactly zero first/last samples. `T=N/fs` is the actual sampled duration.

| Kind | Implemented signature |
|---|---|
| Whistle | Instantaneous frequency `1600 + 2600*t/T + 150*sin(2*pi*5*t)` Hz, analytically integrated phase. Deterministic independent of seed; seed is validated for a uniform interface. |
| Clicks | 8 kHz cosine pulses with Gaussian amplitude standard deviation 0.12 ms, truncated at ±6 standard deviations. Nominal centers `0.05 + 0.16*j` s while less than `T-0.02`; seeded uniform ±3 ms center jitter and 0.8–1.2 pulse weights. |
| Vessel | Twenty 80 Hz harmonics, weights `k**(-1.2)`, seeded uniform phases. Seeded Gaussian broadband shaped by `1/sqrt(1+(f/1000)**2)` in amplitude; its DC bin is zero. Harmonic and broadband components are separately RMS-scaled, mixed with broadband coefficient 0.35, then multiplied by `0.75+0.25*sin(2*pi*3*t)` before the common taper/final RMS scaling. |

The vessel's finite-record Fourier synthesis is a phenomenological colored-noise construction, not a rotating-machine or cavitation solver. Its spectrum has some energy beyond the fractional-delay accuracy band, explicitly assessed below.

## Ideal channel and fractional delay

For horizontal range `R`, source depth `zs`, receiver depth `zr`, and sound speed `c`:

- `rd = sqrt(R**2 + (zs-zr)**2)`;
- `rs = sqrt(R**2 + (zs+zr)**2)`;
- `p_channel(t) = p_source(t-rd/c)/rd - p_source(t-rs/c)/rs` with the surface enabled; omit the second term otherwise.

Distances are in metres, so the gains are the dimensionless ratio `1 m / distance`. The default stationary source/receiver depths are 8/12 m and `c=1500 m/s`. The second term is an image source with ideal pressure-release reflection coefficient −1. There is no seabed, refraction, frequency-dependent absorption, scattering, source motion, Doppler, or spatially distributed source.

Delay is **not rounded** and is not implemented by linear interpolation or a circular FFT shift. For sample delay `d+mu`, with integer `d=floor(delay_s*fs)` and `0<=mu<1`, use

```text
h[k] = sinc(k-mu) * Kaiser_129_beta8.6[k],  k = -64,...,64
h = h / sum(h)
```

Apply direct full linear convolution of the zero-extended input with this kernel, then align the kernel center to physical delay `d`. Integer delays bypass the FIR and are exact. Only unit-DC normalization is applied: there is no frequency-specific correction to conceal attenuation. Output length is `N + ceil(max(delay_s)*fs) + 64`, including all positive-time interpolation tails and one possible spare zero sample.

A finite centered sinc has numerical pre-ringing up to 64 samples before an ideal arrival; this is not a physical precursor. Negative-time values are discarded if a requested delay is shorter than 64 samples. The production 25/100/400 m paths all have longer delays. This method is not an exact causal continuous-time propagator for arbitrary broadband or discontinuous waveforms, and the Nyquist region is deliberately **not** in the accuracy claim.

### Measured delay accuracy, not an all-frequency claim

The retained audit drives the real `propagate` function with impulses and compares the independently evaluated complex response to `exp(-i*2*pi*f*delay)`: 19 fractional delays from 0.05 to 0.95, and 401 frequencies from DC through `0.4*fs`.

- Maximum absolute passband amplitude error: **0.000129722 dB**.
- Maximum absolute passband phase error: **9.12817e-6 rad**.
- Maximum complex relative error: **1.54224e-5** (gate: `<5e-5`).
- Worst amplitude error at `0.49*fs`: **−1.55892 dB**. Near Nyquist the sinc approximation can therefore substantially attenuate pressure; these samples remain in the audit rather than being suppressed.
- For the nine default direct-path source cases, integrated channel amplitude errors relative to exact delay/gain ranged from **−0.000195519 to +0.000000363382 dB**. This is a waveform-specific diagnostic, not a wider-band guarantee.

Separate sinusoidal tests cover 300, 2000, 8000, 12000 and 19200 Hz at 48 kHz with four fractional delays, negative gain, and signed phase. An end-of-record impulse verifies no wraparound or dropped delayed tail. A two-path sinusoid checks signed interference against independently written slant-distance equations, and a source at the pressure-release plane cancels exactly.

## Generic receiver, noise and analysis support

The receiver is an assumed **257-tap causal linear-phase FIR lowpass**, `scipy.signal.firwin`, cutoff 10 kHz, Kaiser beta 8.6, unit DC gain. This is a generic demonstration, not a calibrated instrument or an anti-alias model of acquisition hardware. Its group delay is **128 samples**, retained in the returned physical clock; it is not silently removed. The full convolution adds 256 samples after the channel record.

Input-referred noise is Gaussian white noise with nominal per-sample standard deviation **0.0005 Pa**, generated with `default_rng(SeedSequence([seed, 1]))` and passed through the same FIR. A 256-sample random prehistory makes the recorded noise stationary from its first sample instead of starting with a filter startup ramp. Its realized filtered RMS is measured, not forced to a target. The same seed intentionally supplies common random numbers for comparisons; different cases are not claimed to be independent noise trials. Changing `fs` changes the white-noise spectral density because its per-sample variance, not its Pa²/Hz density, is fixed.

The half-open receiver analysis window is

```text
start = max(0, floor(min_path_delay*fs)-64)
stop  = len(received_pa)
```

It includes interpolation support, all emission samples after propagation, and the complete receiver tail; it excludes most preceding idle time. `analysis_window_samples=[start,stop]` is authoritative; seconds use the same half-open endpoints. The source metric uses its unpadded emission window instead. The receiver metrics are full-window RMS including DC:

```text
snr_db = 20*log10(rms(received_clean_pa[start:stop]) / rms(noise_pa[start:stop]))
received_level_db_re_1uPa = 20*log10(rms(received_pa[start:stop]) / 1e-6)
```

No detector, classifier, detection probability, or field performance is inferred from this SNR. Recomputing full-record RMS including the pre-arrival interval gives a different estimand.

## Spectral metrics

`rms` returns a Python float and uses scaling to avoid overflow while squaring large finite pressures. It includes DC, returns zero for silence, and rejects empty, nonfinite, complex, or non-1-D arrays.

Both PSD functions require uniform sampling supplied by the caller; a numeric `fs` does not certify an irregular acquisition. They reject invalid `fs`, noninteger/too-short segments, and `nperseg` exceeding the input length rather than silently shortening it. Use explicit smaller segment lengths for short records.

- `welch_psd` returns `(frequency_hz, density)`, using periodic Hann windows, per-segment constant detrending, `nperseg//2` overlap, mean averaging, one-sided density, and no padding.
- `spectrogram_psd` returns `(time_s, frequency_hz, density_frequency_by_time)`, using periodic Hann windows, per-frame constant detrending, explicit integer `0<=noverlap<nperseg`, one-sided density, and no boundary extension or padding. Times are centers of half-open frame supports; incomplete final frames are discarded.
- For Pa input, densities are **Pa²/Hz**. Spectral dB, if plotted externally, must use `10*log10(density/1e-12)` and identify any display floor. This is distinct from pressure-level dB and listening normalization.
- Independent DFT tests cover the Hann scaling, odd/even Nyquist endpoint factors and window-weighted detrended Parseval. A PSD integral is not generally the unwindowed full-record RMS squared for a nonstationary click train.

## Verification and retained artifacts

Run the scoped numerical tests:

```sh
uv run pytest -c pyproject.toml tests/test_passive.py tests/test_metrics.py -q
uv run python scripts/build_site.py --output build/site
```

The scoped suite has **52 tests**, including the independent delay, signed-interference, spectral-normalization and deterministic-replay checks described above.

The site catalogue indexes all **18 executed default simulations**: three source kinds × 25/100/400 m × surface off/on, 48 kHz, two-second emission, seed 2026. Each case has downloadable full arrays and metadata. Measured analysis-window SNR spans **12.6084 to 45.3423 dB** for this chosen artificial source/noise scale; it is not a performance envelope.

The delay-accuracy measurements above describe the stated audit grid, not every possible frequency and fractional delay. Reproducibility uses the pinned numerical environment; bit-identical results across arbitrary NumPy/SciPy versions and platforms are not promised.
