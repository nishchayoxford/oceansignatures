# Timestamped offline chunk observations

This additive module leaves `pipeline.py`, `reference.py`, `metrics.py`, their numerical conventions, and version 0.1.0 unchanged. It reads ordinary RIFF/WAVE chunks with **caller-declared** first-sample UTC, channel, evidence kind and optional scalar calibration. It is bounded offline processing, not a demonstration of onboard real-time performance, hardware integration, proprietary formats or deployment acquisition.

## API

```python
from oceansignatures.pipeline import Calibration
from oceansignatures.observations import (
    RecordingChunk, analyse_chunks, frame_spectra, write_observations,
)

chunks = [RecordingChunk(
    path="recording.wav",
    start_utc="2026-01-01T00:00:00.000000Z",  # Caller must know this.
    channel=0,                            # Zero-based; no channel averaging.
    calibration=None,                     # Unknown: digital units only.
    evidence_kind="measurement",          # Required: synthetic OR measurement.
    source="caller-supplied source description",  # Optional, not verified.
    expected_sha256=None,                 # Optional pre-existing lowercase pin.
)]
result = analyse_chunks(
    chunks,
    anchor_utc="2026-01-01T00:00:00.000000Z",
    frame_seconds=2.0,
    nperseg=1024,                         # 2048 also supported at fs=8192.
    noverlap=512,                         # None means floor(nperseg/2).
    window="hann",                       # Periodic Hann or boxcar.
    detrend=False,                        # Or "constant", per inner segment.
    bands_hz={"low": (100.0, 1000.0)},
    active_intervals=[],                  # [(start_UTC, stop_UTC), ...]
)
paths = write_observations(result, "new-observation-products")
```

The timestamp in this example is an API illustration, not an inferred measurement time. Measurement evidence kind is a caller declaration, not authentication. Synthetic and measurement chunks cannot be combined in one result. Unknown and known calibrations cannot be combined; different known scalar calibrations may be applied to their respective chunks, with every coefficient retained.

## Clock and missing-data contract

- Inputs are sorted by declared UTC. UTC must be timezone-aware with zero offset. ISO strings support at most six fractional-second digits; finer strings are rejected rather than truncated. No file creation time, sequence number or audio event is used to invent UTC.
- `anchor_utc` defaults to the earliest declared start. Supply it explicitly when aligning observations with predictions. Integer sample offsets and a common WAV sampling rate are authoritative. There is **no resampling or clock-drift correction**.
- A declared timestamp must agree with an integer sample position within 0.500001 microseconds, accounting only for microsecond serialization rounding. Its signed rounding residual is retained. This is not a tolerance for acquisition jitter, clock drift or arbitrary fractional-sample shifts. Supported WAV rates are 1–384000 Hz.
- Nonoverlapping analysis frames are `[k*N, (k+1)*N)`, with `N=frame_seconds*fs` required to be integral. The output spans the floor-aligned frame containing the first sample through the frame containing the last sample. A supplied anchor may yield negative relative sample/time offsets.
- `window_start_sample` and `window_stop_sample` are exact integer boundaries; `window_start_s` and `window_stop_s` are relative to the anchor. UTC strings are human-readable labels rounded to microseconds, not substitutes for the exact sample clock.
- Chunk boundaries do not reset analysis frames or inner spectral segments. Identical digital samples and scalar conversions produce identical whole/chunk spectra even for nonaligned boundaries.
- Duplicated or overlapping chunk support raises `ValueError`; no deduplication, winner selection or averaging is hidden. `metadata.gaps` records internal missing sample intervals.
- Any incomplete frame (including a frame intersecting an internal gap) has `valid=False` and NaN PSD/STFT/band power. Missing samples are never filled with zeros or interpolated. Entirely missing frames are retained. `partial_window` specifically denotes first/last record-edge truncation; internal gaps are distinguished by coverage and gap intervals.
- `coverage_fraction` uses actual samples divided by the complete frame length, not by the available record length. `valid` means complete support only. `qc_ok = valid & no_near_full_scale_samples & no_active_transmission`.
- Active-transmission intervals are caller-declared half-open UTC intervals. Any positive-duration frame intersection flags contamination; touching an endpoint alone does not. The flag does not subtract a signal or prove source attribution. Complete contaminated frames retain spectra; downstream passive comparisons should use `qc_ok`.
- Near-full-scale contact is `abs(digital)>=0.999` on the selected original digital channel, before calibration. Counts include available samples in invalid frames. This is potential saturation, not proof of hardware clipping.

## Spectra and units

`frame_spectra(samples, fs, *, nperseg=1024, noverlap=None, window="hann", detrend=False)` is the common kernel for later prediction processing. It accepts real finite samples in the caller's amplitude units. It returns `frequency_hz`, `psd[freq]`, complex `stft[segment,freq]`, and inner `segment_start_sample` relative to that frame. It never silently shortens nperseg, pads, or applies an extra whole-frame detrend.

For each complete inner segment, let `y` be the signal after the explicitly selected segment-mean removal (or unchanged when `detrend=False`), `w` the periodic window, and `X=DFT(y*w)`. Under `detrend="constant"`, an **exactly constant** segment uses its shared stored sample as its mean, so true silence remains exact even when a floating-point sum/divide reduction would round to a different value. Constancy is exact sample equality, not approximate equality or an amplitude threshold. Every nonconstant segment retains ordinary arithmetic-mean subtraction; no weak residual is clipped or floored. The raw sample-square representability guard still applies before detrending, including to constant input. The stored complex density transform is

`Z[k] = X[k] / sqrt(fs * sum(w**2)) * sqrt(d[k])`,

where `d=2` at positive-frequency interior bins and `d=1` at DC and the even-length Nyquist bin. Odd transforms have no Nyquist endpoint. Welch density is `mean(abs(Z)**2, axis=segments)`. Thus bin-summed segment PSD equals `sum((y*w)**2)/sum(w**2)`, not necessarily the unwindowed frame mean square. The FFT now uses physical-amplitude samples directly: the existing sample-square guard and frame-length bound keep its intermediates finite, without a frame-wide normalization that can erase a weak segment. The denominator is evaluated as `sqrt(fs)*sqrt(sum(w**2))` to avoid an unnecessary overflowing product. The standalone kernel accepts finite positive `fs` (not just WAV rates); its frequency grid is formed without taking `1/fs`.

**Float64 power guard:** each frequency is scaled independently *after* computing the physical STFT. With `a=abs(Z[:, k])`, `m=max(a)`, and `q=mean((a/m)**2)`, the implementation evaluates `(q*m)*m`, not `mean(abs(frame-normalized_Z)**2)*frame_peak**2`. A zero `m` uses divisor 1. This preserves a representable weak bin alongside an unrelated large bin; it also avoids overflowing an unscaled segment-power sum or `m**2` when the final mean is finite. **Nonzero-bin evidence is retained from the computed physical FFT before density normalization.** Nonfinite PSD, or zero PSD in a bin with any nonzero pre-normalization coefficient, raises `ValueError`. This catches both final-square underflow and earlier division underflow that would otherwise make every coefficient and `m` zero. The decision concerns the final mean PSD at each frequency, not each segment coefficient: a weak coefficient may round to zero while another segment in that same bin supplies representable mean power. Such a mean remains valid. Existing rejection of nonzero sample-square underflow and overflow remains; floating-point overflow/invalid/divide failures in the spectral calculation become `ValueError`, not warnings. Representable subnormal outputs are allowed with float64 rounding precision. No amplitude threshold, clipping, calibration change or spectral-convention change is introduced.

This guard uses the **computed float64 FFT**, not an exact-real DFT oracle. Computed FFT zeros, exactly detrended constants, and inputs annihilated by zero window weights remain valid. Nonzero raw input alone is not sufficient evidence of nonzero spectral power. Cancellation or dynamic range can also put mathematical Fourier content below FFT accuracy; this guard neither recovers that content nor infers it from raw samples. It does prevent *known nonzero computed coefficients* from silently becoming a zero final PSD during normalization. FFT, taper and nonconstant mean-subtraction arithmetic are unchanged.

Default: 2 s frames, 1024-sample periodic Hann, 512-sample overlap, no detrending, one-sided density, arithmetic mean, no padding. Unused final inner-segment samples are explicitly counted in `metadata.processing.unused_tail_samples_per_frame`. The existing `metrics.py` removes segment means and the legacy recording pipeline does not; **do not assume those APIs are interchangeable**. Matched predictions must call this kernel or reproduce its explicit contract.

Result arrays:

| Key | Shape / meaning |
|---|---|
| `frequency_hz` | `[frequency]`, shared for all frames |
| `psd` | `[window, frequency]`, Pa²/Hz or digital²/Hz |
| `stft` | `[window, segment, frequency]`, Pa/√Hz or digital/√Hz |
| `segment_start_sample` | `[window, segment]`, absolute sample offset from anchor; segment stop is start+nperseg; centre is `(start+nperseg/2)/fs` |
| `valid`, `qc_ok`, `active_transmission`, `partial_window` | `[window]` Boolean masks |
| `band_power` | `[window, band]`, Pa² or digital² |
| `band_level_db` | `[window, band]`, dB re (1 µPa)² when calibrated, otherwise dB re 1 digital² |
| `band_spl_db_re_1upa` | `[window, band]` only when calibrated; otherwise Python `None` / JSON `null` |
| `band_names` | Insertion order of `bands_hz`; default single `full` band |

Bands are low-inclusive and high-exclusive, except the Nyquist upper edge is inclusive. Integration is the discrete sum of selected PSD bins multiplied by `df=fs/nperseg`; no trapezoidal endpoint halving or fractional-bin weighting. Empty bands and out-of-range/reversed edges fail. User-defined bands may overlap; therefore summing arbitrary band powers is not necessarily total power.

Calibration **reuses** `pipeline.Calibration` and `pipeline.calibrate_samples`:

`p[Pa] = digital * full_scale_voltage * 10**(-(gain_db+sensitivity_db_v_per_upa)/20) * 1e-6`.

Sensitivity is dB re 1 V/µPa. Digital ±1 corresponds to ADC ±full_scale_voltage **peak**, not RMS or peak-to-peak. Every chunk retains all three coefficients, actual Pa/digital-fraction scale, and flat-response/ADC conventions. Top-level `calibration_status` is `declared scalar calibration` or `unknown`. Unknown calibration retains digital units and no absolute SPL. This declared scalar chain is not a calibration certificate. Overflow and nonzero underflow of conversion, squares and PSD are rejected rather than returned as silence.

True zero-power silence is valid data: power is zero and its logarithmic level is negative infinity in memory/NPZ. JSON cannot represent infinity; it writes `null`, distinguishable from missing data by `valid` and zero `band_power`. CSV uses `-inf` for true silent levels and empty cells for missing/unknown values. This additive module's silent-level convention intentionally differs from the unchanged legacy pipeline's null-only convention.

## Resource, input and export protection

Chunks are inspected and hashed using bounded byte snapshots; decoding uses the exact verified snapshot bytes, not a second unverified path read. Before each contributing decode and again at the end, SHA256 is checked against inspection. An optional predeclared pin makes the first read fail closed too. Inputs are never written. A changed source is rejected; this is not an OS lock preventing another program from writing it.

Resource defaults are also hard ceilings (callers may lower them): 20,000,000 bytes per chunk; 1,000,000 decoded interleaved samples per read; 1,000,000 samples per analysis frame; 10,000 output windows; 2,000,000 total PSD+STFT array cells. At most 1024 chunks, 64 channels, 128 bands and 1024 active intervals. Sparse timestamp spans count toward the output bound. No whole-deployment waveform is concatenated. Bounded chunk bytes are retained temporarily and may be reread/rehashed for successive frames, trading I/O for memory safety; the returned products themselves are bounded in-memory arrays, not an unbounded live stream. The shared standalone kernel also limits frame length and spectral cell count.

`write_observations` stages a complete family before publishing a **new directory**:

- `observations.npz`: all NumPy result arrays, deterministic fixed-timestamp ZIP members, no pickle/object arrays; includes complex STFT and masks.
- `windows.csv`: one row per window/band, clocks, QC flags, band metrics and units.
- `summary.json`: strict finite JSON, metadata, clocks, masks, band metrics, frequencies, array shapes/dtypes and NPZ SHA256. Dense PSD/STFT values are in NPZ, not duplicated in JSON.

An identical existing family is accepted byte-for-byte. Different or incomplete existing output directories are rejected, including pre-created empty directories. Output symlinks, symlink ancestors, `..` paths and targets aliasing source inputs are rejected before publication. Source bytes are checked again for a new publication. Single-writer offline operation is the contract; concurrent adversarial filesystem renames are not supported. No dependencies or original APIs were changed.

## Executed synthetic verification

From the repository root:

```bash
PYTHONPATH=src .venv/bin/python .cache/research-extension/observations/run_synthetic.py
.venv/bin/pytest tests/test_observations.py -q
.venv/bin/ruff check src/oceansignatures/observations.py tests/test_observations.py
.venv/bin/ruff format --check src/oceansignatures/observations.py tests/test_observations.py
```

The retained runner generates explicitly **test-only synthetic** 8192-Hz WAV data with a declared synthetic UTC, 256/1024-Hz tones, stepped amplitudes, four nonaligned chunks, one active interval and a partial last frame. It calls the public API, independently checks expected band powers, reads NPZ back, checks raw hashes and permits only identical reruns. Two executions yielded four windows, `valid=[true,true,true,false]`, `qc_ok=[true,false,true,false]`, **exactly zero** whole/chunk PSD difference, and maximum band-power absolute error `8.673617379884035e-17 Pa²`. Current products are in `products-final/` and the receipt is `execution-final.json`, retained beside the runner. The earlier `products/` and `execution.json` are preserved as the pre-explicit-calibration-status execution, not overwritten. They are numerical/software evidence, not field measurements or experimental validation.

Tests additionally cover independent direct-DFT/Parseval oracles for odd/even periodic Hann/boxcar transforms and both detrend choices; sample-clock alignment; gaps and overlap; active endpoint semantics; calibration and unknown units; SHA changes; bounded decoding and output; extreme numeric inputs; and protected deterministic export/readback.

The OBS-1 correction adds exact DC/Nyquist dynamic-range regressions, including normalized WAV ingestion with differing per-chunk calibrations and 1024-sample segments. It also checks representable subnormal bins, rejection of nonzero PSD underflow, finite means whose unscaled power sum would overflow, and genuine spectral zeros. Correction receipts and before/after ordinary tone/random-recording diagnostics are retained separately under `.cache/research-extension/observation-fix/`. Earlier synthetic bundles and the original failed independent-review evidence are preserved, not regenerated or overwritten. The arithmetic correction can change PSD and band outputs at roundoff level; prior artifact hashes are not expected to identify newly generated corrected products.

The separately authorized constant-detrend correction retains RED/GREEN receipts, before/candidate source bytes and compatibility comparisons under `.cache/research-extension/acoustic-constant-fix/`. Its regressions require exact zero for repeated tiny, ordinary and large supported samples at odd and longer segment lengths, for both windows. Independent complex-DFT controls cover ordinary nonconstant detrending; one-ULP perturbations and mixed constant/weak rows check that exact-constant handling does not suppress real signals. The prior numerical correction and its failed constant-silence verdict remain preserved as historical evidence, not relabelled as accepted.

The further scoped OBS-2-remaining correction retains its normalization design, before/candidate bytes, RED/GREEN receipts and private replay products under `.cache/research-extension/acoustic-normalization-fix/`. Regressions use independent positive-power Decimal oracles for an adjacent-float detrended pair and a window-supported Hann impulse at extreme finite sample rates. Controls preserve a representable mean with a lost weak segment in the same bin, window-annihilated silence and the distinction between computed FFT zeros and normalization loss. The prior failed verdict remains historical evidence; these correction checks are not a new independent acceptance verdict or authorization to promote public products.
