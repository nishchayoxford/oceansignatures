# Verification report

## Scope and accepted products

This additive extension remains **OceanSignatures 0.1.0**, in the existing
repository and website. The original passive/active experiments and their
scientific outputs are preserved. Environmental acoustic prediction and
calibrated synthetic observations are connected end to end; native CFD remains
a **separate local hydrodynamic-pressure contamination investigation**.

This report records completed numerical, native-evidence, integration and
private-candidate browser checks. It is not itself a deployment receipt, field
validation certificate or manufacturer endorsement. Deployment must use the
same accepted bytes and be verified separately.

| Accepted family | Files, including its manifest | SHA256 of `manifest.json` |
|---|---:|---|
| [Acoustic experiment](../site/data/research/acoustic/manifest.json) | 29 | `ca5de2d069c29aef419c4a8a1810daa93f41de7ad94ccae42bcf5175351c238a` |
| [Native CFD export](../site/data/research/cfd/manifest.json) | 54 | `a697728056aa4f3119a8810a9125767b63e33f1f9e35d0db91fbe8bb4da80768` |
| [Combined research catalogue](../site/data/research/manifest.json) | 85 | `c010360aefb7a71cbf13b3fffbc14b1f65623ce00c2625c3ffe3d77131bfc1be` |

The combined family contains the two complete producer families, a catalogue and
its master manifest. The assembled website has **299 files**, of which **214
nonresearch files** are unchanged by scientific-bundle promotion. Relative to
the original 244-file project baseline, 240 original files retain their bytes;
the four authorized integration edits are the README, site HTML, CI workflow
and ignore rules. New code, tests, documentation and research assets are additive.

## Completed checks

Counts below describe separate, overlapping selections. They must **not** be
summed into a fictitious total number of unique tests.

| Gate | Executed result | Meaning |
|---|---|---|
| Final integrated Python suite after local bundle promotion | **543 passed**, no skips/errors/failures | Collection, execution and JUnit test-identity multisets reconcile; exact project imports pinned |
| Frontend unit tests | **27 passed**, no skips/cancellations/failures | Original and research catalogue/control contracts |
| Ruff lint, format and Git whitespace checks | Passed | Final numerical source and integrated working tree |
| Fresh independent acoustic acceptance | **240 declared tests + 621 new independent tests passed** | Computed float64 FFT, units, scaling, silence, subnormals, representability and compatibility |
| Accepted CFD-runner correction | **75 scoped tests + 23 independent probe groups passed** | Serialized clocks, nested budgets, interruption handling and evidence-preserving re-entry |
| Private-candidate browser acceptance | **1,254 checks passed**: research 430, original laboratory 816, narrow media 8 | Actual candidate content, controls, error handling, media and downloads |
| Independent browser-evidence readback | **599 checks passed** | Separate authentication of the retained browser evidence |

The acoustic reviewer authenticated 59 frozen files, reviewed the complete
correction and latest delta, checked 144 exact-constant cases, 84 independent
nonconstant DFT cases and 84 historical byte-compatibility cases. These subsets
are not additional tests to add to the table. Two fresh public-CLI executions
reproduced **all 29 acoustic files byte-for-byte**, including the manifest;
readback covered 39 scientific arrays and 16 WAV files.

Earlier failed numerical verdicts are retained. The final normalization repair
tracks nonzero FFT evidence before density normalization and assesses
representability of the **final per-frequency mean PSD**. It does not introduce
a weak-signal floor, turn nearly constant input into silence, narrow the finite
positive standalone sample-rate contract, or reject a representable mean merely
because one segment contribution is individually unrepresentable. Extremely
large standalone rates are API-boundary tests, not physically claimed recording
rates. Independent acceptance closed both previously reproduced false-silence
cases without changing the ordinary nonconstant detrending path.

## Numerical demonstration

The declared acoustic case is a **240 s** synthetic passage, with **120 two-second
windows**, **117 compared windows**, one missing window and two active-contaminated
windows. The selected recording rate is 8192 Hz. Processing uses 4096-sample
periodic-Hann segments, 50% overlap and no detrending for this experiment.

- Injected and recovered same-realization discrepancy: **+6 dB**.
- Maximum spectral gain-control error: **4.085620730620576e-14 dB**.
- Baseline-versus-ensemble 100–4000 Hz band RMSE: **0.19905314406605468 dB**.
- Band RMSEs in explicit 100–500, 500–1000, 1000–4000 Hz order:
  **0.2476561516891248**, **0.16616602425248697**, **0.07984742193231777 dB**.

Gain recovery and finite-realization sampling discrepancy answer different
questions; neither establishes empirical source-model accuracy in the ocean.
The one-factor alternatives quantify deterministic sensitivity, not a
probabilistic confidence interval. Equations, axes and exclusions are documented
in [RESEARCH.md](RESEARCH.md), [NOISE_PREDICTION.md](NOISE_PREDICTION.md),
[OBSERVATIONS.md](OBSERVATIONS.md) and [COMPARISON.md](COMPARISON.md).

## Genuine native CFD evidence

Two **OpenFOAM v2206** `pimpleFoam` calculations completed at Re100, each with
**15,360 cells**, reaching **140 s**. The timesteps are **0.01 and 0.005 s**;
probe and force output clocks are matched at 0.02 s. The cases use a 2D laminar
cylinder model, not a three-dimensional turbulent housing-noise model.

All **206 original pinned files per case** were rechecked, including native mesh,
inputs, logs, numerical results and original completion receipts. Solver and
mesh stages actually executed; no synthetic pressure fixture is represented as
native output. The public subset retains 45 original pinned members plus the
completion receipt per case; omitted full-volume and unused-plane members are
explicitly enumerated. Subset verification must not be reported as rechecking
all 206 original pins.

| Late-window result, [100,140) s | Δt 0.01 s | Δt 0.005 s |
|---|---:|---:|
| Mean drag coefficient | 1.36154887 | 1.35837211 |
| Mean-removed lift RMS | 0.24503264 | 0.23905163 |
| Strouhal number from rising signed-lift crossings | 0.16662719 | 0.16681965 |
| Complete rising-to-rising periods | 6 | 6 |
| Lift PSD peak bin | 0.175 Hz | 0.175 Hz |

Native kinematic pressure is converted to Pa using density **once**. Pressure
PSD is in Pa²/Hz, but it is not thereby radiated acoustic PSD. The Strouhal
estimate is from signed lift crossings, not simply the largest pressure bin.
The 0.025 Hz PSD-bin spacing is distinct from the crossing-based estimate.

The wake movie uses **21 actual native snapshots**, 100–140 s, at 1200×500 pixels,
2 fps and **10.5 s**. It uses a common colour scale and no invented intermediate
flow fields. Historical portable-source replays reproduced all 52
non-verification payloads; the two verification/manifest files legitimately
change when moving from full-case verification to the public subset. See
[CFD_RESULTS.md](CFD_RESULTS.md) and [REPRODUCIBILITY.md](REPRODUCIBILITY.md).

## Website acceptance and limitations

The exact assembled candidate was exercised with responsive desktop Chrome at
1440, 390, 320 and 768 px for the research studies, and 1440/390 px for all
24 original laboratory cases. Verification authenticated **85 research HTTP
payloads** and **181 real saved browser downloads**: 48 research and 133 original.
All saved download bytes matched independently supplied SHA256 pins.

All **76 screenshots** were visually inspected. Genuine pointer-started movie
playback reached the end with 21 decoded frames. No unexpected runtime/network
errors or blocking UI findings remained. The only advisory is that dense figure
labels and movie previews are small at 390 px; full-size figures and downloads
were verified. This is desktop responsive testing, not physical-mobile,
touch-screen, cross-browser or screen-reader certification. Fresh keyboard and
mixed-input tests are distinct from the authenticated source-identical prior
headed mouse-only study-select evidence.

Owned QA browsers and local servers were stopped. Parent readback authenticated
**843 evidence payloads plus their seal**, independently rehashed all 181 saved
downloads, and checked the full candidate inventory before local integration.
The earlier served research family is retained separately; promotion does not
overwrite historical source bundles or failed verdicts.

## Reproduction and scientific limits

Use [REPRODUCIBILITY.md](REPRODUCIBILITY.md) for the public acoustic computation,
saved-native reanalysis, original laboratory rebuild and authenticated assembly.
Ordinary testing and saved-native replay do **not** run OpenFOAM. A new native
solve is a separately gated opt-in operation, unnecessary for reproducing these
figures.

The direct propagation model omits shallow-water boundaries, refraction,
directivity, tones and Doppler. The CFD pair supplies two-level timestep
sensitivity, not spatial convergence, extrapolated order or a second-velocity
study; residual startup and short statistical windows remain visible. No
sensor-coupling model joins CFD pressure to the acoustic prediction, and the
PSDs are not summed. No compatible calibrated field recording is bundled; no
field accuracy, biological classification, hardware performance or commercial
integration is claimed.

The official NPL technical-guide index was retrieved and independently confirms
the seawater-absorption guide. Its interactive endpoint failed TLS validation
during access checks; no certificate bypass was used. Optional NOAA acquisition
is outside the build, tests and this accepted bundle.

The [18-minute presentation](PRESENTATION.md) follows the computed results and
these boundaries rather than presenting verification as experimental validation.
