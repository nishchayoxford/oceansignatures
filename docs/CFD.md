# Native laminar cylinder CFD

## Scope and interpretation

This is genuine installed OpenCFD **v2206 `pimpleFoam`**, configured with
`simulationType laminar`, for a stationary circular cylinder in incompressible
water. It is **two-dimensional and spanwise-uniform**, not LES, not turbulent
broadband noise, not a finite hydrophone housing and not an acoustic-radiation
calculation. Native pressure must not be added to acoustic receiver PSDs.

| Quantity | Declared value |
|---|---:|
| Diameter D | 0.01 m |
| Freestream U | 0.01 m/s |
| Kinematic viscosity ν | 1e-6 m²/s |
| Density for pressure/force interpretation ρ | 1000 kg/m³ |
| Reynolds number UD/ν | 100 |
| Circular outer radius | 0.2 m = 20 D |
| Numerical extrusion thickness | 0.001 m, one cell, both end planes `empty` |
| Native cell count | 15,360 hexahedra |
| Cylinder circumference / radial cells | 192 / 80 |
| Radial last/first cell expansion | 50 |
| Baseline timestep / sensitivity timestep | 0.01 / 0.005 s, fixed |
| Declared production endpoint | 140 s, each fresh run from t=0 |
| Probe/force output interval | 0.02 s in both runs |
| Midplane / whole-field intervals | 2 / 10 s |

`cfd/campaign.json` contains the reproducible parameters and proposed analysis
windows. Native cases are confined to ignored `.cache/cfd/` in this repository;
no pre-existing FlowAcoustics CFD cases are reused or modified.

## Boundary conditions, initial condition and numerics

The eight-sector body-fitted blockMesh O-grid uses circular arcs on both cylinder
and outer boundary. The upstream semicircle (x≤0) has fixed freestream U and
zero-gradient p. The downstream semicircle has fixed gauge p=0 and `inletOutlet`
U: zero normal gradient on outward flux, prescribed freestream on reverse flux.
This is a split velocity-inlet/pressure-outlet approximation to an unconfined
external flow, not a circular solid wall or rectangular symmetry channel.
The cylinder has zero no-slip velocity and zero normal pressure gradient.
A 20D radius reduces but does not establish independence from outer-boundary
placement; no domain-size sensitivity has been run.

At t=0 only, the interior velocity is `(U, 0.001 U, 0)`, while the upstream
boundary is `(U,0,0)`. This small physical transverse initial perturbation breaks
reflection symmetry. It is neither random pressure nor a fabricated probe
history. Pressure projection and impulsive no-slip startup generate a large
early transient, deliberately excluded from stationary analysis. The perturbation
is not maintained by a time-dependent forcing.

The scheme is second-order `backward` in time, `Gauss linear` advection and
corrected diffusion, with one PIMPLE outer loop and two pressure correctors
(PISO-style). Pressure uses GAMG and the final solve has zero relative tolerance;
velocity uses smoothSolver. The very low mesh nonorthogonality makes this a
simple laminar benchmark, not a wall-modelled calculation. Near-wall spacing and
wake resolution still need a spatial convergence study for accuracy claims.

Installed solver `pimpleFoam.C`, `createFields.H`, the laminar tutorial's
`turbulenceProperties`, the probes and surfaces case dictionaries, forceCoeffs
configuration and the `inletOutlet` implementation header were inspected before
authoring. The solver explicitly defines `p` as kinematic pressure. Actual
initial and saved fields have dimensions **`[0 2 -2 0 0 0 0]`**, so

`p_hydrodynamic_Pa = rho_kg_m3 * p_native_m2_s2`.

This conversion occurs once, not twice. The gauge datum is the outlet p=0;
fluctuation statistics must state their mean-removal window.

## Primary references and plausibility checks

1. C. H. K. Williamson (1989), *Oblique and parallel modes of vortex shedding in
   the wake of a circular cylinder at low Reynolds numbers*, **JFM 206**, 579–627,
   [doi:10.1017/S0022112089002429](https://doi.org/10.1017/S0022112089002429),
   [primary full text](https://authors.library.caltech.edu/records/wavw6-f7k56/files/WILjfm89a.pdf).
   Equation (6) and Figure 14 give the parallel-shedding fit
   `St = -3.3265/Re + 0.1816 + 1.600e-4 Re`, hence **0.164335 at Re100**.
   The paper distinguishes end-condition-dependent oblique shedding from
   parallel shedding; this one-layer computation only models the latter ideal.
2. R. D. Henderson (1995), *Details of the drag curve near the onset of vortex
   shedding*, **Physics of Fluids 7**, 2102–2104,
   [doi:10.1063/1.868459](https://doi.org/10.1063/1.868459),
   [primary full text](https://authors.library.caltech.edu/records/bqdev-z6q09/files/1_2E868459.pdf).
   Equations (3)–(4), Table I's **vortex-shedding** rows give
   `Cd_viscous=2.5818/Re^0.4369` and
   `Cd_pressure=1.4114-0.2668 Re^0.1648 exp(-0.003375 Re)`.
   Their sum at Re100 is approximately **1.3500**. These are fitted numerical
   reference trends, not experimental truth or an acceptance tolerance.

The native force coefficient reference is `Aref=D*span=1e-5 m²`,
`magUInf=.01 m/s`, with drag along +x and lift along +y. Compare the mean native
Cd and lift-derived `St=fD/U` only after the transient. A plausible frequency or
mean force is not validation of acoustic pressure, pressure RMS, a real housing,
a finite-span cylinder or mesh/time independence.

## Execution, limits and retained evidence

Generator-only operation does not run CFD:

```bash
.venv/bin/python scripts/run_cylinder_cfd.py .cache/cfd/new_case --end-time 140
```

To create and execute a **fresh** baseline case in one command:

```bash
.venv/bin/python scripts/run_cylinder_cfd.py .cache/cfd/new_baseline \
  --end-time 140 --delta-t .01 --execute --timeout-s 650
```

For a fresh timestep sensitivity use `--delta-t .005 --timeout-s 1000`.
The runner refuses existing generation destinations, symlinked case destinations,
unsupported physics/mesh sizes and noncommensurate sample clocks. Native execution
is confined to this repository's `.cache/cfd/`, including nested case directories.
It uses serial execution (one numerical thread), a shared lock, bounded commands
and a cumulative 2400 s native campaign allowance. Do not erase receipts to evade
that allowance. No install, MPI execution or automatic continuation is performed.
The process-group supervisor has real timeout and SIGTERM-cleanup tests for serial
children; it is not a general MPI descendant supervisor.

The generator/runner admission contract also applies to direct Python calls:

- Time and velocity scalars written into native dictionaries round-trip to their
  Python/metadata values. Historical short tokens such as `.01` and `.005`
  (written as `0.01` and `0.005`) are unchanged. The **serialized** timestep must
  divide every sample/output period within an absolute `1e-8` step tolerance,
  with zero relative tolerance and finite step counts, before any directory is
  created. The admitted `delta_t=0.02/3` is no longer rounded to `0.00666667`.
- Under the shared lock, accounting reads **every depth** of `run_record.json`
  beneath the campaign root. Each receipt file must be a nonempty JSON list;
  each `wall_seconds` must be numeric (not boolean), finite and nonnegative,
  and the accumulated charge must remain finite. Malformed/unreadable receipts
  or symlinks anywhere in the campaign tree fail closed. A charge of 2390 s plus
  a new 60 s commitment is rejected equally for flat, nested and mixed layouts;
  a total commitment of exactly 2400 s remains admissible.
- `execute_case` checks freshness **after acquiring the lock and before writing
  identity evidence**. Only initial generated files (`0/{p,U}`, the two constant
  property dictionaries, four system dictionaries, `case_metadata.json` and
  `cylinder.foam`) and their three directories are admissible. Prior logs,
  receipts, identity snapshots, mesh, time/function-output directories and all
  other unexpected entries cause rejection without overwriting existing bytes
  or modification times. It is not a restart or custom-case execution API.
- An `interrupted=true` receipt immediately fails the case, even with zero exit
  and no timeout. No later stage runs and completion cannot say `completed`.
  Logs, run receipts and the failure reason remain preserved.

The measured completed pilot `pilot_re100` reached 6 s (600 steps), with native
stage walls 0.565 s blockMesh, 0.615 s checkMesh and 12.794 s pimpleFoam.
It exercised actual field writes and three VTP snapshots, rather than projecting
throughput from a solver run without output. Its recorded pre-step maximum
Courant number is 1.002521. A final post-solve Courant value is not printed by the
installed solver; the report explicitly does not assert one.

Every attempted native case retains:

- Full `0`, `constant`, `system`, mesh and actual solution directories.
- `log.blockMesh`, `log.checkMesh`, `log.pimpleFoam` and native function outputs.
- `generator_executed.py` and input SHA256 identities before execution.
- `run_record.json`: commands, UTC start/end, elapsed wall time, actual exits,
  timeout/interruption flags and log hashes.
- `completion.json`: achieved versus declared endpoint, completion gates,
  native cell count, logged max Courant, output sizes and SHA256 map.

Mesh admission requires actual `checkMesh -allTopology -allGeometry` output
with **`Mesh OK.` and no `Failed N mesh checks`**. The pilot and baseline mesh
passed: two geometric/solution directions `(1 1 0)`, max skewness
0.007925 and max nonorthogonality 1.71e-6 degrees. A successful mesh check is
separate from flow accuracy. Failed/partial solver output must never be labelled
completed just because some field directories exist.

## Reading the actual producer

The public reader verifies the original native identities, converts pressure
once, and exports numerical results without launching a solver. For retained
local cases, choose a fresh output directory:

```bash
.venv/bin/python scripts/build_cfd_results.py \
  --baseline .cache/cfd/re100_baseline \
  --sensitivity .cache/cfd/re100_dt005 \
  --output .cache/research-extension/cfd-export/new_analysis
```

For replay without the original local cases or OpenFOAM, use the packaged native
subsets and standalone command in [CFD_RESULTS.md](CFD_RESULTS.md). The original
pilot reader is retained locally as historical verification, not a dependency.

Native p/U probes are `postProcessing/probes/0/{p,U}`. Header comments preserve
ordered probe locations (metres):

1. `(0.01, 0.005, 0)`
2. `(0.02, 0.005, 0)`
3. `(0.04, 0.005, 0)`
4. `(0.02,-0.005, 0)`
5. `(-0.02,0,0)`

p has time plus five scalars. U has time plus five parenthesized three-component
vectors. The pilot contains 300 finite rows from .02 through6 s; p/U clocks match
exactly and spacing is uniformly .02 s within decimal serialization roundoff.
Native force columns are named in
`postProcessing/forceCoeffs/0/coefficient.dat`; `Cd` is column1 and `Cl` column4
when indexing from zero including time.

Planes are `postProcessing/wakePlane/<time>/midplane.vtp`, XML VTK PolyData with
ASCII points, triangle connectivity and point-interpolated p/U arrays.
**Native declared array dtype is Float32**, despite ASCII precision17. Each
pilot frame has 107,904 points and 215,040 triangles. Both independent XML/NumPy
parsing and installed `vtkXMLPolyDataReader` read the real output. Maximum |z|
is 5.67e-20 m. The native topo cut contains **75,574 exactly zero-area triangles**;
preserve the original and filter degenerates only in an explicitly labelled
visualization copy. Do not mistake native cutting-plane triangulation count for
volume-cell count or silently replace native point values. Whole fields and
polyMesh are also retained as ASCII; open `cylinder.foam` with an OpenFOAM reader.

## Completion and scientific-analysis gates

The completed local campaign retains its coordinator and status under
`.cache/cfd/`, alongside frozen source/configuration bytes. Both declared cases,
`re100_baseline` and `re100_dt005`, have independently verified completion
receipts. The public generator above can reproduce fresh cases; the private
coordinator is not required for analysis replay. Do not edit a generator during
native execution or rewrite an original completion receipt.

After native completion, verify every pinned file, complete p/U time coverage,
finite arrays and uniform sampling. The initial campaign proposed 80–140 s;
observed amplitude growth motivated the current **[100,140) s** primary window.
The longer [80,140) s sensitivity and every 20 s block are retained. Check lift
amplitudes rather than assuming either window is stationary. Count
like-direction lift zero-crossing periods and report the spread/cycle count,
then cross-check a frequency spectrum using the verified uniform clock. A
60 s window has finite spectral resolution; do not present its peak bin as an
arbitrarily precise frequency. Off-centre wake pressure can contain harmonics;
identify fundamental shedding from lift or transverse velocity rather than
blindly calling the largest pressure peak the shedding frequency.

Compare both cases over matched physical windows and identical probe clocks.
Two timestep values provide a modest sensitivity, not an observed order of
accuracy, a grid-convergence result or domain independence. The public reader,
figures and data are documented in [CFD_RESULTS.md](CFD_RESULTS.md); no data under
the actual native case paths are synthetic fixtures.

## Completed native results

The baseline `re100_baseline` **completed** 14,000 native steps to
140.00000000001765 s (target140 s), with End and zero exit. Its full execution
wall was 275.585 s; output size was 1,497,690,819 bytes. Independent readback
verified all **206 pinned output/input files**, 7,000 finite p/U/force samples,
all solver timestamps and **70 finite native VTP frames**. The `.02 s` p/U clocks
agree bit-for-bit; force output retains more time digits, giving maximum
force/probe timestamp difference 2.86e-11 s. This is fixed-step floating-point
accumulation/serialization, not adaptive sampling. Sample values were not
resampled; common sample-index membership was checked against raw clocks within
1e-9 s for half-open statistic windows.

The originally proposed80–140 s window still contains some amplitude growth:
lift RMS in80–100,100–120,120–140 s is0.23153,0.24415,0.24591 respectively.
Do not describe the whole80–140 s record as perfectly stationary. Over the later
**100–140 s** window, mean **Cd=1.36155**, lift RMS=0.24503, and six complete
signed-lift rising-zero-crossing periods give **f=0.16662719 Hz and
St=0.16662719**. The population period standard deviation is 0.00154754 s.
The 80–140 s record has nine complete periods, St=0.16640133 and
Hann-periodogram peak 0.16666667 Hz with bin width 0.01666667 Hz. This current
definition uses Cl=0, not the finite-window-mean-subtracted crossing threshold
used by the earlier private helper; both are distinguished in the results notes.
The late pressure fluctuation RMS at the five ordered probes is approximately
`[0.0102143, 0.0132212, 0.00591476, 0.0130258, 0.0000311048] Pa`.
These are local hydrodynamic fluctuations. `baseline_analysis.json` preserves
all early and late windows, including adverse startup trends, rather than
selecting only the favourable one.

`baseline_reader_check.json` covers every plane. Native directory names later
contain accumulated time digits; sort them numerically rather than requiring
literal integer directory names. `baseline_wake_140.png` is an inspected
visualization of the actual final transverse-velocity plane; its display copy
masks92,160 triangles with twice-area≤1e-18 m² (including75,574 exactly degenerate
triangles), without altering native files. This figure shows a resolved
alternating wake, with visible downstream mesh coarsening.

The matched Δt=0.005 s case also completed 140 s (28,000 native steps). Both
cases' 412 pinned files have been verified. Its matched late-window mean
Cd=1.35837211 and St=0.16681965 change by −0.2333% and +0.1155%, respectively;
lift RMS changes by −2.4409%. This is sensitivity evidence, not convergence.
Full per-probe results, residual growth and the real native wake movie appear in
[CFD_RESULTS.md](CFD_RESULTS.md). No solver is currently needed for their replay.
Do not rerun generation into an existing case or automatically extend an endpoint.

## Tests

Case generation was developed with observed RED→GREEN tests. The focused suite
covers geometry/units, invalid parameters and clock compatibility, generator-only
CLI, fresh bounded command receipts, timeout and actual SIGTERM child cleanup,
and native output-root rejection. The initial setup verification recorded
**11 passing tests** and **309 integrated tests** in that earlier source tree;
these historical counts are not whole-project acceptance of later additions.
The public analysis has its own numerical, identity and replay tests in
`tests/test_cfd_analysis.py`; see the current verification scope in the results notes.

The subsequent CFD-01–CFD-04 correction has separate retained RED→GREEN
regressions for native scalar serialization, all-depth finite accounting,
interrupted zero-exit stages and immutable direct-call re-entry. Execution tests
load identical generator bytes in **private TEST-ONLY surrogate repositories**;
recording stubs supply synthetic acceptance logs and a real cooperative Python
child establishes the zero-exit SIGTERM receipt. They do not run OpenFOAM.

Correction evidence is under `.cache/research-extension/cfd-runner-fix/`: original
and candidate source snapshots, hashes, diffs, command receipts, collected/JUnit
identities and read-only compatibility checks against the original native cases.
The historical native generator snapshots and completion anchors remain separate
from the corrected source identity: **do not re-pin or relabel old outputs as
produced by this candidate**. Exact historical arguments, including integer
`end_time=140`, reproduce the ten generated input files per production case
without changing the `.01`/`.005` dictionary bytes. No new native solve is needed
for this software correction, and the accepted run002 science is unchanged.
The earlier review's denied realistic accounting rerun remains historically
**unexecuted**; the newly authorized surrogate regressions are distinct evidence.
