# Native CFD results and portable replay

The accepted export is available at **`site/data/research/cfd/`**. Its retained
local build is `.cache/research-extension/cfd-export/run002/`; public replay does
not require that cache.
Start with `overview.json` (website contract), the three PNGs under `figures/`,
`native_wake.mp4`, and `analysis_summary.json`. These are **native hydrodynamic
CFD**, not acoustic predictions. The original two OpenFOAM cases are unchanged;
this reader never launches a solver or imports a private cache helper.

## Actual completed cases and results

Both original completion receipts and all **206 pinned files per case** were
verified. Both meshes passed native `checkMesh` with 15,360 cells and two active
geometric/solution directions. All three stage receipts have zero exit, no
interruption/timeout, matching log hashes and `End`. Full solver clocks contain
14,000 / 28,000 steps, reaching 140 s within the declared roundoff tolerance.
Each case has 7,000 p/U/force samples at 0.02 s cadence. See [CFD.md](CFD.md) for
the executed geometry, boundaries, initial perturbation, schemes and references.

Primary statistics use **[100,140) s**, 2,000 samples, after the observed startup.
The finer timestep is a separate fresh native simulation, not interpolated data.

| Metric | Δt = 0.01 s | Δt = 0.005 s | Fine versus baseline |
|---|---:|---:|---:|
| Mean Cd | 1.36154887 | 1.35837211 | −0.2333% |
| Mean-removed lift RMS | 0.24503264 | 0.23905163 | −2.4409% |
| Rising signed-lift crossing f (Hz) | 0.16662719 | 0.16681965 | +0.1155% |
| St = fD/U | 0.16662719 | 0.16681965 | +0.1155% |
| Complete rising-to-rising periods | 6 | 6 | — |
| Period standard deviation (s) | 0.00154754 | 0.00162318 | — |
| Lift PSD peak bin (Hz) | 0.175 | 0.175 | Bin width 0.025 Hz |

The primary local pressure fluctuation RMS values are:

| Probe | Position (m), x/y/z | Baseline RMS (Pa) | Fine RMS (Pa) | Change |
|---|---|---:|---:|---:|
| P1 | 0.01 / 0.005 / 0 | 0.0102143431 | 0.0100500838 | −1.6081% |
| P2 | 0.02 / 0.005 / 0 | 0.0132212398 | 0.0132304274 | +0.0695% |
| P3 | 0.04 / 0.005 / 0 | 0.0059147557 | 0.0059511938 | +0.6161% |
| P4 | 0.02 / −0.005 / 0 | 0.0130258429 | 0.0130349606 | +0.0700% |
| P5 | −0.02 / 0 / 0 | 0.0000311048 | 0.0000303931 | −2.2880% |

The comparison checks exact equality of initial fields, mesh, transport/turbulence
properties, spatial discretisation, solution controls and executed generator.
The metadata differ only in timestep. The control dictionaries differ only in
`deltaT` and commensurate output-step counts: physical probe/force/plane cadences
are identical. Probe clocks match bit-for-bit. This is **two-level timestep
sensitivity**, not convergence order, mesh convergence or demonstrated accuracy.

### Adverse stationarity evidence is retained

| Window (s), half-open | Baseline lift RMS | Fine lift RMS |
|---|---:|---:|
| 80–100 | 0.23152606 | 0.22385835 |
| 100–120 | 0.24414899 | 0.23690345 |
| 120–140 | 0.24591119 | 0.24109190 |

The last two blocks grow by **0.7218% / 1.7680%**. The preceding change is
**5.4521% / 5.8274%**. Do not call the whole 80–140 s interval stationary.
The longer [80,140) s sensitivity gives St **0.16640133 / 0.16657334**, nine
complete lift periods and a PSD peak at 0.16666667 Hz with 0.01666667 Hz bins.
At upstream P5 its largest pressure PSD bin is instead 0.01666667 Hz, and its RMS
is more than twice the late-window value: residual drift matters there.
All seven successive 20 s blocks, including startup, remain in the JSON/CSV.

## Definitions, clocks and provenance

- Native p dimensions are verified as `[0 2 -2 0 0 0 0]` (m²/s²).
  **p_Pa = 1000 kg/m³ × p_native**, exactly once. Native p remains in the NPZ.
  Pressure is gauge-referenced at the outlet, not an acoustic receiver pressure.
- Each window removes its own arithmetic mean; RMS is the population root mean
  square of that fluctuation. No linear detrending is applied.
- PSD is the one-sided full-window periodic-Hann periodogram, density-scaled in
  **Pa²/Hz**, with no zero padding. A pressure harmonic is not automatically the
  shedding fundamental.
- Frequency uses rising zero crossings of **signed native Cl**, with linear
  interpolation of crossing times only. It does not use the largest pressure
  peak and does not interpolate/resample the signal. The older private helper
  used crossings after subtracting the finite-window Cl mean, hence its slightly
  different quoted St; this public definition is explicit and reproducible.
- Original p/U and force clocks are both retained. p/U agree exactly. A fixed
  **1e−9 s absolute tolerance, zero relative tolerance** admits native binary64
  accumulation/decimal serialization, not gaps, adaptive steps or resampling.
  Membership is half-open with that tolerance applied only at the boundaries.
  Spectral sample rate is the executed 50 Hz. Crossing locations use the common
  probe clock after force/probe agreement; an independent calculation with the
  original force clock changed f by at most **1.52e−13 Hz**.
- Expected completion digests are explicit API inputs; the CLI defaults are the
  two original digests pinned in the public source. Empty/incomplete maps cannot
  disable checks. Each byte payload is checked against the original expected
  pin and the same retained bytes are parsed and packaged. Original manifests
  are never regenerated or edited to bless changed data.
- Public subset replay checks all **45 retained original pins per case**, plus
  the externally anchored original completion. It does not claim to recheck the
  **161 intentionally omitted files per case**. The complete original run audit
  checked 206; this distinction is visible in `verification.json`.

## Native visualization and archive selection

Each package contains 46 exact native members: the original completion receipt
and 45 pinned inputs/results. They include exact dictionaries, initial fields,
metadata, pre-execution identity map, frozen generator, complete native mesh,
mesh/solver logs, stage receipts, p/U probes, force history and **21 authentic VTP
frames at 100,102,…,140 s**. Actual accumulated-time directory names are preserved.
The other 49 VTP frames and full-volume fields remain local, explicitly omitted
from this public replay subset. No replacement/fake VTK files are created.

VTP arrays declare Float32. The parser respects this dtype and native triangle
connectivity; the numerical extrusion is not misrepresented as a 3D simulation.
Each frame has 107,904 points and 215,040 native triangles, including 75,574
exactly zero-area triangles. Display alone masks 92,160 triangles with twice-area
≤1e−18 m². Native arrays and downloadable bytes remain unchanged. Equal physical
axis scales are used; Gouraud colour shading occurs only inside native triangles,
not a new grid or a fabricated intermediate CFD field.

The movie uses all 21 baseline snapshots with **one common symmetric Uy colour
scale**, 2 fps, 1200×500 pixels and 10.5 s duration. Each 2 s native snapshot is
shown for 0.5 s: nominal 4× time speed, with the final snapshot also held 0.5 s.
No intermediate fields or motion interpolation are generated. `movie.json`
retains exact source times, playback timing and colour bounds.

Archive parts are independent deterministic ZIP files; keep all parts and
`package.json` together. The complete export is **259,714,406 bytes in 54 files**;
the largest file is **29,027,659 bytes**, below the 50 MB limit. There are three
scientific PNG figures and 26 download records in `overview.json`. Every path
there is relative to the export directory. Exact historical native logs/receipts
may contain historical machine paths; ordinary metadata and replay require none.

## Commands and public API

Executed from the repository root (the destination now exists and is protected):

```bash
.venv/bin/python scripts/build_cfd_results.py \
  --baseline .cache/cfd/re100_baseline \
  --sensitivity .cache/cfd/re100_dt005 \
  --output .cache/research-extension/cfd-export/run002
```

Choose a **new** destination for another export. Existing destinations and
symlinked output ancestors are rejected, before input reads or output creation.
A failed validation never creates an accepted results directory. No installs are
needed: only existing numpy/scipy/matplotlib and standard Python are used;
ffmpeg is optional (`--no-movie` disables it).

Portable replay, from the exported folder with its files intact:

```bash
python replay/build_cfd_results.py \
  --baseline native/baseline --sensitivity native/sensitivity \
  --output "$(dirname "$PWD")/replay-new"
```

Use the Python interpreter containing the existing scientific dependencies.
The shell resolves a safe absolute destination; raw `..` path components are
intentionally rejected. No OpenFOAM installation, original cache, sibling
repository or private helper is needed. Replay is the same analysis of saved
native data, not another CFD solve. `analysis-code.zip` contains the same public
reader and a standalone entry point. The source snapshot is bound by the export
manifest.

```python
from oceansignatures.cfd_analysis import (
    BASELINE_COMPLETION_SHA256, read_native_case, analyse_case, compare_cases,
)
case = read_native_case("path/to/native-case-or-packaged-subset",
                        expected_completion_sha256=BASELINE_COMPLETION_SHA256)
summary, arrays = analyse_case(case)
# compare_cases(case, second_case) enforces matched mesh/physics and clocks.
```

Additional public helpers: `parse_probes`, `pressure_pa`, `validate_clock`,
`window_statistics`, `parse_plane`, `package_native_case`, `render_movie`,
`build_export`. `NativeCase` retains verified payload bytes, native histories,
ordered probe locations and evidence reports. Numerical NPZ/CSV outputs retain
raw clocks, dimensional pressure, native pressure, velocity, forces and spectra.

## Verification history and accepted export

The accepted export moves the pressure-history legend outside the scientific
axes. Its geometry regression failed before the correction and passed afterward;
**27 focused CFD-analysis tests** and scoped lint/format checks pass. A fresh
native-data export and visual readback confirm the corrected layout. All **50
other payloads** are byte-identical to `run001`, including every numerical array,
CSV, summary, native archive, other figure and movie. Only that PNG, the bundled
analysis source and its source ZIP change (plus the manifest). Evidence is in
`.cache/research-extension/cfd-legend-readback.json`. The earlier checks below
describe `run001`, not a new integrated acceptance of the current whole project.

- Strict observed RED→GREEN slices covered extraction/conversion, clocks and an
  explicitly **test-only analytical pressure fixture**, fail-closed native
  identity/completion, native connectivity, portable archives, matched analysis,
  public CLI, actual ffmpeg and a documented-replay command correction.
- `pytest tests/test_cfd_analysis.py -q -p no:cacheprovider`: **26 passed**.
- Integrated `pytest -q -p no:cacheprovider`: **360 passed** in the shared tree
  at this verification; scoped Ruff: **all checks passed**.
- Actual native CLI and **two packaged-source CLI replays** completed. The replay
  manifests are byte-identical. All **52 non-verification payloads** match the
  original native export byte-for-byte, including numerical products, original
  ZIP data, figures and MP4. Only `verification.json` changes scope (206 versus
  45 checked original pins); its manifest consequently differs as expected.
- An independent NumPy FFT oracle found maximum PSD absolute errors of
  **1.30e−18 / 1.73e−18 Pa²/Hz**. Independent RMS matched exactly; pressure
  conversion and retained clocks matched exactly. All 412 original pinned
  files were rechecked unchanged after export/replay.
- PNGs and a movie frame were visually inspected; ffprobe confirmed 21 frames,
  2 fps and 10.5 s. No single exported file exceeds 50 MB.

Detailed local verification is retained in
`.cache/research-extension/cfd-export/acceptance.json`, produced by
`verify_export.py` without importing the production analysis functions.
`replay001/` and `replay002/` are verification copies. `run001-initial/` preserves
an earlier export unchanged; it had a rejected raw-`..` replay command in its
instructions. The correction was regression-tested and rebuilt into a fresh
`run001/`, not overwritten onto differing existing data. Both are retained as
historical evidence; the display-corrected **run002** supplies the accepted public
bundle. See [the verification report](VERIFICATION.md) for current acceptance
scope and [the reproduction guide](REPRODUCIBILITY.md) for public-only commands.

Limitations remain: 2D laminar Re100, residual startup, short statistics,
point-interpolated Float32 planes, no spatial/domain convergence, and no acoustic
radiation calculation. Native pressure must not be added to receiver acoustic
PSDs or represented as measured/turbulent broadband hydrophone noise.
