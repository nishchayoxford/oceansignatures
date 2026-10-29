# NOAA reference recording: rights, acquisition, and relative analysis

## Release scope

This release distributes reference-reading and acquisition code, not the NOAA recording itself. The website catalogue has `reference: null`; all displayed experiment results are synthetic. The external source below is an optional example for exercising provenance and signal-quality checks, not a calibration or validation dataset. Rebuilding the website and running the tests require no recording download.

## Recording and rights decision

Source page: <https://www.pmel.noaa.gov/acoustics/multimedia.html>

Original WAV: <https://www.pmel.noaa.gov/acoustics/multimedia/HB-ship-SBNMS.wav>

The source page labels this exact link:

> Humpback whale vocalizations with ship in background in Stellwagen Bank National Marine Sanctuary

The **Credits** text on that same page states:

> Unless otherwise noted (copyrighted material for example), information presented on the PMEL Acoustics Program webpage is considered public information and may be distributed freely.

> If you elect to use materials from this page, please cite NOAA PMEL Acoustics Program as the source, and include the appropriate URL of the page(s) from which the materials were taken.

These passages were read in the complete live HTML, including its sidebar, not inferred from an archive's public-access policy. The extraction service omitted the Credits sidebar; the retained direct source HTML and text contain it. The WAV is directly linked in the Noise Reference Stations section. No clip-specific third-party copyright credit or reuse restriction is marked there. Artist collaborator links are in a separate sidebar section and are not identified as credits on this clip. The source's linked disclaimer was also checked: <https://www.pmel.noaa.gov/acoustics/disclaimer.html>. It concerns external links and non-endorsement, not an additional audio restriction.

**Decision:** redistribution is expressly permitted by the source page's freely-distributable-public-information notice, with NOAA PMEL attribution and URL. No separate Creative Commons license or broader jurisdiction-independent legal conclusion is asserted. Do not apply the project's code license to this audio or claim it as project-created material.

Suggested credit:

> Audio: NOAA PMEL Acoustics Program, https://www.pmel.noaa.gov/acoustics/multimedia.html. Analysis: OceanSignatures. No NOAA affiliation or endorsement is claimed.

## Byte identity and decoding

The acquisition receipt records UTC retrieval at `2026-10-29T12:02:56.632051+00:00`. The HTTP Last-Modified timestamp is **not** the recording date.

| Property | Verified value |
|---|---|
| Bytes | 833,976 |
| SHA256 | `24bf234e0d302ca91fd3e31f3b964185244403f74666569b753ca12080b59750` |
| Magic / container | RIFF/WAVE |
| Encoding | Signed 16-bit PCM, little-endian |
| Channels | 1 |
| Frames / scalar sample count | 416,966 / 416,966 |
| Sample rate | 44,100 Hz |
| Distributed-file duration | 9.45501133786848 s |
| HTTP media type | `audio/x-wav` |

Both libsndfile and ffprobe agree on encoding, channels, sample count, and sample rate. The fetch script was run against the live URL into `.cache/reference/refetch-verification.wav`; `cmp` confirmed identity with the initially acquired file. No complete archive was downloaded.

## Scientific and editorial limitations

This is an **editorial/listening example**, not asserted to be original unprocessed hydrophone data. The source page does not document this clip's gain, filtering, excerpt boundaries, resampling, or playback-speed changes. Those properties are **unknown**, not assumed absent. No speed correction has been applied. Time and frequency axes describe the distributed WAV using its header sample rate, not independently verified original-event timing.

The original recording date, instrument model, hydrophone sensitivity, acquisition gain, and conversion to pressure are unavailable in the inspected source. PCM decoding maps stored integers to digital full-scale fractions. It does **not** calibrate pressure or normalize the recording by its own peak. There are no Pa, underwater SPL, source-level, or species-classifier claims. The animal label is NOAA's description, not an inference from this code.

### Actual QC result — do not conceal this warning

The decoded example has a mean of `-0.19523107919583071`, RMS `0.8564569320584555`, and peak absolute sample `1.0`, all in digital units. Approximately **44.29%** of samples have absolute magnitude at least `0.999`. The analyser therefore flags **possible digital clipping** and a **substantial DC offset relative to digital RMS**. It does not prove which recording, editing, or encoding stage caused these features.

This clip is useful as an authentic public **QC/limitations teaching example**, not as pristine biological truth, a clean spectral template, calibrated level evidence, or a quantitative validation dataset. Show the warning alongside plots. Retain the original waveform; segment-wise constant detrending is used only for density spectra. That detrending cannot undo distortion.

## Reusable Python API

```python
from oceansignatures.reference import load_reference, analyse_reference

recording = load_reference(
    "data/reference/HB-ship-SBNMS.wav",
    expected_sha256="24bf234e0d302ca91fd3e31f3b964185244403f74666569b753ca12080b59750",
)
result = analyse_reference(recording, channel=0)
```

- `load_reference(path, *, expected_sha256, max_bytes=20_000_000)` returns `{samples, fs, metadata}`. `samples` is a float64 array shaped `(frames, channels)`, with no downmixing, gain normalization, or resampling. The digest must be literally 64 lowercase hexadecimal characters. The loader checks byte bounds, pin, RIFF/WAVE magic, RIFF size, decode success, nonempty audio, and finite samples. Metadata contains filename, hash, bytes, format/subtype, channel/frame/sample counts, rate/duration, and explicit uncalibrated digital units. RF64/RIFX and non-WAV formats are deliberately unsupported.
- `analyse_reference(recording, *, channel=0, nperseg=2048, spectrogram_nperseg=1024)` returns `{fs, time_s, waveform, frequency_hz, psd, spectrogram_time_s, spectrogram_frequency_hz, spectrogram_psd, metadata, qc}`. Arrays are NumPy arrays. It selects one channel explicitly and requires at least two frames. The spectrogram density is shaped `(frequency, time)`.
- Density spectra use periodic Hann windows, constant detrending per segment, one-sided density scaling, Welch mean averaging with 50% overlap, and spectrogram overlap of 75%, with no padding. Segment lengths are capped to record length and recorded in metadata. Units are **digital units squared per Hz**. A logarithmic plot must explicitly identify its digital reference or relative maximum; it must never be labelled dB re 1 µPa.
- `qc` includes finite status, mean, RMS, peak, zero fraction, silence flag, near-full-scale threshold/fraction, clipping caveat, and warning strings. A near-full-scale sample flags possible clipping; it is not alone proof of clipping. DC warning threshold is `abs(mean) > 0.1 * RMS`.
- No project-specific passive or active source modules are imported. The analyser is independently reusable.

### Float64 representability domain

Finite decoded DOUBLE-WAV samples are not sufficient for analysis. For the selected channel,
every nonzero sample must have a finite, nonzero float64 square. Silence is supported and
returns zero RMS and zero density. The analyser rejects square overflow or underflow-to-zero
with `ValueError` before spectral analysis; it does not replace an unrepresentable signal with
silence or report successful finite QC.

Mean and RMS use a peak-scaled working array. Welch and the spectrogram use the same working
array, then multiply density by the peak squared to restore **digital units squared per Hz**.
This is numerical scaling, not gain normalization of the recording: the returned waveform is
an exact copy of the selected decoded channel. The original Hann/detrending/overlap/density
conventions are unchanged (ordinary results agree to floating-point rounding).

The supported domain additionally requires finite time/frequency axes, finite spectral
intermediates, and finite restored density. Any nonzero computed working-density bin that
would become zero on restoration is rejected, including very small leakage bins. Subnormal
values are allowed when they remain nonzero; exact-real arithmetic and arbitrary dynamic
range are not promised. Conservatively, `fs * max(effective segment lengths)` must fit in
float64, preventing silent overflow of SciPy's density denominator. These checks can reject
otherwise finite inputs whose final mathematical result might be representable by a different
algorithm. Numerical overflow/division failures become `ValueError`, without RuntimeWarnings.
For example, alternating DOUBLE-WAV samples at ±`1e200` or ±`1e-200` are rejected; ±`1e153`
and ±`1e-150` are supported at 8 Hz with four-sample segments. No physical calibration or
implicit resampling is introduced.

The observed arrays contain 416,966 waveform samples, 1,025 Welch bins, and a 513-by-1,625 spectrogram.

## Fetch and reproduce

From the repository root, after `uv sync --frozen`:

```bash
# Offline test suite; no network by default.
uv run pytest -c pyproject.toml tests/test_reference.py

# Optional acquisition, separate from the site build; review source terms first.
uv run python scripts/fetch_reference.py \
  --destination data/reference/HB-ship-SBNMS.wav

# Offline verification and QC of that local file; no audio writes.
uv run python scripts/fetch_reference.py --verify-only \
  --destination data/reference/HB-ship-SBNMS.wav
```

The script's CLI prints JSON containing metadata, QC, and array shapes; it does not serialize full arrays. `scripts.fetch_reference.fetch_reference(destination, *, max_bytes=20_000_000, timeout_s=20.0, total_timeout_s=60.0)` returns verified file metadata and installs the downloaded bytes atomically only after decoding and pin validation. Its source URL and default hash are fixed in the script. It rejects changed contents rather than updating the pin. The byte cap is enforced on Content-Length and actual bounded reads, and callers cannot increase the maximum. A 20-second socket timeout and a 60-second deadline checked between reads bound ordinary transfer stalls; this is not an OS-level hard real-time deadline and an in-flight read/DNS lookup can outlast a deadline check. Redirected final URLs are rejected. Existing destination bytes survive download/validation failures, and temporary `.part` files are cleaned up.

The default destination is `data/reference/HB-ship-SBNMS.wav`. It is excluded from version control. The script does not fabricate a provenance manifest or infer reuse terms automatically.

## Provenance for any later recording release

Preserve `schema_version`, `id`, `title`, `filename`, `source_page_url`, `download_url`, `retrieved_utc`, `sha256`, `bytes`, `origin`, `attribution`, `rights` (explicit permission, exact notice, source URL, exceptions), `recording` (encoding, channels, frames, sample rate, duration, calibration status), and `processing` (editorial status, known or unknown speed/gain/filtering, project transformations, source-label-only identification). Absence of a processing note is not evidence of original unprocessed data.

## Verification evidence

Offline tests use generated tiny WAV fixtures explicitly labelled as unit-test signals, never as real reference audio. They cover decoding, pin/magic/size/truncation/nonfinite/empty failures, analytical tone energy normalization and selected-channel handling, silent/short input, invalid arguments, fetch byte/deadline/redirect failures, atomic-preservation/cleanup, CLI offline mode, and QC warnings. These tests exercise software behaviour, not the truth of a recording's biological labels or its calibration.
