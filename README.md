# OceanSignatures

**An open, reproducible underwater-acoustics laboratory: source → water → receiver → recorded signature.**

[Open the interactive field guide](https://nishchayoxford.github.io/oceansignatures/)

OceanSignatures connects inspectable mathematical models to computed pressure histories, spectra, spectrograms and active-echo range estimates. It includes an offline recording-analysis pipeline that keeps uncalibrated digital measurements separate from physical pressure.

This is an independent educational project. It does not identify species, simulate a particular commercial instrument, claim manufacturer endorsement, connect to hardware or transmit sound.

## Start with one physical question

1. **Passive:** choose the whistle-like source and direct propagation. Compare 25, 100 and 400 m while holding the source fixed. Pressure decreases and arrival moves later.
2. **Multipath:** at the same range, include the ideal surface reflection. Follow the delayed, sign-reversed copy into frequency-dependent reinforcement and cancellation.
3. **Signatures:** compare whistle-like sweeps, click trains and vessel-like harmonic/broadband noise. Equal reference RMS does not imply equal peak pressure or spectral shape.
4. **Active:** inspect a known chirp and its noisy point-target return. Correlation supplies a candidate delay; two-way travel supplies range. At low SNR, a peak can be wrong.
5. **Evidence:** inspect the repeated-noise results, including every failed estimate. Then read the calibration and validation boundaries before comparing with recordings.

The website explains each step beside the controls. Controls select exact precomputed cases; they do not run a hidden solver or relabel a nearby result.

## Implemented components

| Component | Implementation | Important convention |
|---|---|---|
| Passive sources | Analytic frequency-modulated whistle, seeded Gaussian-envelope clicks, harmonic plus coloured-broadband vessel surrogate | 1 Pa RMS over the unpadded, tapered emission at nominal 1 m; not measured animal/machinery source levels |
| Passive channel | Direct spherical pressure spreading; optional negative image-source surface path; finite windowed-sinc fractional delay | Horizontal range differs from slant path distance; no bottom, absorption, refraction or Doppler |
| Passive receiver | Assumed causal 10 kHz lowpass FIR and filtered seeded pressure noise | Complete filter tail and 128-sample group delay retained; not measured hardware response |
| Active benchmark | 50 ms, 3–9 kHz tapered chirp; stationary scalar point target; outward and return spreading | Chosen target gain, no validated target strength; fixed-SNR range comparisons rescale noise |
| Active estimator | Signed matched correlation, normalized by known simulation noise SD; fixed 0–750 m search | Largest absolute peak always returns a range; no detection threshold or false-alarm guarantee |
| Recording pipeline | Channel selection, original-file hash, QC, waveform and Welch PSD; optional explicit scalar calibration | Unknown calibration stays in digital units; no inferred acquisition timing or proprietary decoder |

The default experiment matrix contains 18 passive cases and 6 active cases. A supplementary active study retains all 72 trials at 150 m: 24 declared seeds at each of −30, −18 and −6 dB input SNR. These are controlled synthetic experiments, not field-performance estimates.

## Reproduce

Requirements: Python 3.11 or later, [uv](https://docs.astral.sh/uv/), and Node.js 22 or later for the frontend tests. The accepted build uses Python 3.11 and the committed `uv.lock`.

```sh
git clone https://github.com/nishchayoxford/oceansignatures.git
cd oceansignatures
uv sync --frozen
uv run pytest -q
node --test tests/test_site.mjs
uv run ruff check src tests scripts
uv run ruff format --check src tests scripts
uv run python scripts/build_site.py --output site
python -m http.server 8000 --bind 127.0.0.1 --directory site
```

Open `http://127.0.0.1:8000/`. The site has no CDN, analytics, online solver or external runtime dependency. After downloading the repository and installing dependencies, building and viewing its synthetic results requires no network. Reference acquisition is separate and never part of the build or test suite.

To regenerate only the data/plot bundle in a separate destination without touching the served site:

```sh
uv run python scripts/build_site.py --output build/site
```

Numerical arrays and metadata remain separate from display reductions. Download full-resolution data for analysis; min/max waveform envelopes preserve transients but are not a replacement for the sampled signal. Renderer/platform changes can affect figure bytes even when numerical results agree; use the pinned environment for reproduction.

## Analyse a recording

```sh
uv run python -m oceansignatures.pipeline recording.wav build/recording-analysis --channel 0
```

This writes `summary.json`, `spectrum.csv` and `signal.npz` without modifying the original recording. Identical reruns are allowed; conflicting existing output is refused. Use `--calibration-json` only when the actual ADC scaling, gain and hydrophone sensitivity are known; see [the calibration convention](docs/PIPELINE.md).

No real recording is bundled with this release. [The reference guide](docs/REFERENCE.md) documents an optional external NOAA example and its substantial signal-quality limitations. An official audio example is not automatically raw, calibrated or suitable validation data.

## Methods and reading route

- [Passive sources, propagation, receiver response and spectral units](docs/PASSIVE.md)
- [Active chirp, two-way delay, correlation, SNR and estimator failures](docs/ACTIVE.md)
- [Recording ingestion, calibration, QC and protected exports](docs/PIPELINE.md)
- [External recording provenance and limitations](docs/REFERENCE.md)
- [Scientific source register](docs/sources.json)

## What would make this a measurement comparison?

The next scientific step is to compare **receiver-domain predictions with receiver-domain observations**, using matched frequency bands, windows and estimator definitions. It requires:

- authorized recordings and documented file/channel/time conventions;
- hydrophone sensitivity and frequency response, gain, ADC full-scale convention and uncertainty;
- source/receiver geometry, environmental sound-speed information, boundaries and operating conditions;
- traceable source/event annotations, including unknown or unreviewed labels;
- independent evaluation grouped by deployment, encounter or original recording—not random overlapping clips from the same event;
- component checks and held-out tests, with noise-only trials if claiming detection performance.

The source recipes, channel and receiver are separate modules so measured responses, better propagation models and justified scattering models can replace the initial assumptions. Those replacements are not implemented merely because the architecture allows them.

## Licence

Project code and original synthetic outputs: [MIT](LICENSE). Third-party recordings and material retain their own terms and attribution requirements; the project licence does not relicense them.
