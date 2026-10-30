# OceanSignatures: an 18-minute technical presentation

**Title:** From environmental noise to a traceable underwater receiver comparison

**Presenter:** Nishchay Tiwari

**Format:** eleven slides with a short website demonstration; 18 minutes before questions.

This is a suggested talk and rehearsal guide, not a claim of instrument or field
validation. Use the existing computed figures rather than illustrative curves.
The notes below are speaking prompts: explain the mechanism, point to the data,
and pause for the indicated demonstration instead of reading every bullet.

## Timing and scientific story

| Slide | Time | Main claim | Evidence to show |
|---|---|---|---|
| 1 | 0:00–1:00 | A useful signature needs a traceable source-to-receiver chain | Project question and evidence classes |
| 2 | 1:00–2:30 | One workflow can preserve distinct physical branches | Three-branch architecture |
| 3 | 2:30–4:30 | Source level is not received level | Acoustic spectra and propagation |
| 4 | 4:30–6:00 | Predict on the receiver clock | Passage and travel-time history |
| 5 | 6:00–8:00 | Calibration and QC define a comparable observation | Chunk/window workflow and missing/active intervals |
| 6 | 8:00–10:00 | Gain recovery and sampling error answer different questions | Comparison and residual figures |
| 7 | 10:00–11:30 | Sensitivity identifies assumptions worth testing next | One-factor sensitivity figure |
| 8 | 11:30–13:30 | Local flow pressure is a separate sensor-environment problem | Authentic cylinder wake and native logs |
| 9 | 13:30–15:30 | Numerical evidence includes startup and sensitivity | Pressure spectra and force-history comparison |
| 10 | 15:30–17:00 | Results remain inspectable beyond a presentation | Website selector, numerical download, manifest/replay |
| 11 | 17:00–18:00 | The next step is a matched measurement campaign | Specific measurement requirements and conclusion |

For a 15-minute slot, omit the live download on slide 10 and shorten slides
3, 5, 8 and 9 while retaining the distinction between acoustic and hydrodynamic
pressure. For 20 minutes, use the extra time to discuss the estimator on slide 6
and the measurement plan on slide 11. Do not fill the time with more plots.

## 1. The question: what should a receiver record?

**On the slide:** one question, the project name and the three evidence labels:
*acoustic prediction*, *synthetic observation*, *native hydrodynamic CFD*.

**Speaking route:**

“A recorded underwater signature combines a source, propagation, the receiving
system and its local environment. I wanted to build a demonstration in which I
could follow those steps all the way to a quantitative comparison, rather than
stop at an attractive spectrum. OceanSignatures now takes an environmental and
vessel scenario through acoustic prediction, creates a controlled recording,
processes it using explicit calibration and timestamps, and evaluates residuals.
Alongside that, I ran a native flow calculation to examine local pressure around
a simplified sensor housing. Those are related receiver problems, but they are
not automatically the same physical pressure. I will show what the computations
establish and the measurements needed to take the next step.”

**Transition:** a credible combined project begins by keeping the branches separate.

## 2. Architecture: shared workflow, distinct physics

**Visual:** the three branches in the website's Research section. A simple
slide version can use these exact connections:

```text
Vessel + wind → source/ambient PSD → propagation → expected acoustic PSD ─┐
                                                                     ├→ aligned residuals
Synthetic WAV → UTC/channel → calibration → QC → observed acoustic PSD ─┘

Current + cylinder → native transient CFD → local pressure/forces/wake
                                             ↓
                            separate contamination investigation
                            (no established acoustic/sensor coupling)
```

**Speaking route:**

“The upper branch predicts what should reach the receiver. The middle branch
turns a recording into comparable quantities. They meet only when units,
frequency support and time windows agree. The lower branch asks what flow around
the housing could contribute locally. I deliberately did not draw an addition
arrow from CFD pressure into acoustic PSD. A sensor response and coupling model
would be needed before that operation had a physical meaning. The architecture
makes that missing relationship visible instead of hiding it in a total-noise
curve. The original passive and active teaching experiments are still available;
this is an incremental research workflow within the same project.”

**Implementation evidence:** [research](../src/oceansignatures/research.py),
[observations](../src/oceansignatures/observations.py),
[comparison](../src/oceansignatures/comparison.py),
[native reader](../src/oceansignatures/cfd_analysis.py).

## 3. Source strength and propagation

**Visuals:** [source/received spectra](../site/data/research/acoustic/figures/spectra.png)
and [propagation](../site/data/research/acoustic/figures/propagation.png).

**One equation:**

$$S_r(f,t_r)=S_{wind}(f,t_r)+\frac{Q_{vessel}(f)}{r(t_e)^2}
10^{-\alpha(f)r(t_e)/10000},\qquad t_r=t_e+r(t_e)/c.$$

Here $r$ is in metres and $\alpha$ in dB/km; $Q$ is in Pa² m²/Hz and $S$ in
Pa²/Hz. The two components are assumed incoherent.

**Speaking route:**

“The vessel uses the published JOMOPANS-ECHO spectral-density source law, while
the wind contribution uses an empirical ambient-noise fit.[5][2] They enter the
calculation at different places. Vessel source strength still has a distance
reference attached, so it must be propagated. Wind ambient noise is already a
receiver-level contribution: spreading it from the vessel would be a physical
mistake. I keep the two powers separate until the final linear sum. For the
vessel path I use spherical spreading and frequency-dependent absorption; the
absorption convention and environmental assumptions are documented.[3][4]

“On these plots, follow the vessel source spectrum into the received spectrum.
The units change, and the frequency-dependent loss changes the shape. This is a
transparent direct-path baseline, not a shallow-water waveguide solution. The
research model does not include reflections, bathymetry, refraction or Doppler.
That simplification is explicit so a future comparison can tell us which part
of the chain needs to be improved.”

**Methods:** [equations, references and applicability](NOISE_PREDICTION.md).

## 4. A controlled vessel passage

**Visual:** [passage time history and expected spectrogram](../site/data/research/acoustic/figures/passage.png).

**On the slide:** 150 m bulker; 12 kn; 250 m horizontal CPA; source/receiver depths
8/50 m; 240 s; chosen wind schedule 4→8 m/s; comparison band 100–4000 Hz.

**Speaking route:**

“This is a chosen numerical scenario, not an AIS track or weather reconstruction.
The receiver is fixed while the vessel passes, and the ambient wind contribution
changes over the receiver clock. Closest approach is defined at emission time.
A signal arriving at the receiver was emitted earlier, so the geometry is
resolved on that retarded clock rather than evaluated at the same timestamp.

“The displayed spectrum is an expected receiver spectrum. Within each two-second
frame it is frozen at the midpoint. The largest midpoint-to-edge difference on
the supported bins is about 0.106 dB in this case. That gives a useful diagnostic
of the approximation; it is not a universal error bound. The generated recording
is a receiver-spectrum surrogate, not a phase-resolved moving-source waveform.
Retarded geometry therefore does not mean that Doppler has been implemented.”

**Evidence:** [configuration and numerical summary](../site/data/research/acoustic/overview.json).

## 5. From WAV chunks to an observation

**Visual:** [comparison figure with QC intervals](../site/data/research/acoustic/figures/comparison.png)
and this processing line:

`WAV → hash/channel/UTC → scalar calibration → coverage/QC → STFT/Welch → band powers`.

**Speaking route:**

“The recording path is deliberately less tidy than a single array. Files have
explicit first-sample UTCs and a selected channel. Chunk boundaries do not align
with every analysis window. One chunk is absent, and a declared active interval
contaminates two windows. The processing keeps one common clock across files,
so boundaries do not reset the spectral analysis or create invented samples.

“The synthetic calibration declares hydrophone sensitivity, instrument gain and
ADC peak full-scale voltage. In this example the conversion is 200 Pa per digital
fraction. That is a specified synthetic chain, not a calibration certificate.
The digital-only control follows the same files without calibration and retains
digital units; it is not assigned absolute underwater sound levels.

“Out of 120 windows, 117 are suitable for the passive comparison. The absent
window remains missing rather than becoming a zero-noise observation. The two
active-contaminated windows retain their spectra but are excluded from passive
residuals. That distinction matters: missing data, silence and contamination
are different states, and collapsing them would bias a comparison.”

**Methods:** [observation contract](OBSERVATIONS.md), [comparison contract](COMPARISON.md).

## 6. What does the discrepancy actually verify?

**Visual:** [residuals](../site/data/research/acoustic/figures/residuals.png).

**On the slide:**

$$\Delta L=10\log_{10}(S_{observed}/S_{predicted}),\qquad
P_B=\Delta f\sum_{f_k\in B}S(f_k).$$

- Same realization, known amplitude gain: recover **+6 dB**.
- One realization versus finite-Hann ensemble expectation: broadband RMSE
  **approximately 0.1991 dB**.
- Neither number is a field-validation error.

**Speaking route:**

“There are two tests here that should not be confused. First, I multiply the same
pressure realization by a known amplitude gain and process both versions through
the same calibration and spectral kernel. The expected difference is exactly
six decibels, apart from floating-point arithmetic. Recovering it checks the
software chain, residual sign and gain convention.

“Second, I compare one stochastic realization with its ensemble expectation.
That residual is not forced to zero. The approximately 0.1991 dB broadband RMSE
belongs to this seeded verification experiment. It does not demonstrate that an
empirical vessel model predicts the ocean to that accuracy.

“Even the expected spectrum must match the estimator. A finite Hann window
smooths the input spectrum, particularly at the limits of the synthesis band.
The experiment therefore calculates the expected finite-window estimate rather
than simply comparing a Welch estimate with a raw source-law curve. Finally,
band levels integrate linear power before taking logarithms. Averaging spectral
decibels would answer a different question.”

**Evidence:** [finite-window expectation and three residual definitions](RESEARCH.md).

## 7. Sensitivity is a guide to the next experiment

**Visual:** [one-factor sensitivity](../site/data/research/acoustic/figures/sensitivity.png).

**Speaking route:**

“I then change one assumption at a time: source offset, wind, vessel speed,
closest approach, receiver depth and the water properties used for absorption.
The baseline and alternatives use the same declared comparison bands. A source
offset changes only the vessel contribution; changing receiver depth changes
geometry, not the separate pressure-depth parameter in the absorption model.

“These are deterministic alternatives, not confidence intervals. I have not
assigned probability distributions to the inputs or demonstrated their joint
uncertainty. The value of the exercise is diagnostic: it makes clear which
parameters and metadata must be controlled when interpreting a future residual.
A disagreement should not automatically be blamed on the source model when
geometry, calibration or environmental inputs remain uncertain.”

**Evidence:** [ordered sensitivity cases](../site/data/research/acoustic/sensitivity.json).

## 8. The CFD branch: genuine local flow computation

**Visual:** [native wake](../site/data/research/cfd/figures/native_wake.png);
play the [saved-snapshot movie](../site/data/research/cfd/native_wake.mp4) once.

**On the slide:** OpenFOAM v2206; transient `pimpleFoam`; 2D laminar cylinder;
D = 0.01 m; U = 0.01 m/s; Re100; 15,360 cells; 140 s per run.

**Speaking route:**

“The lower branch is a real native OpenFOAM calculation, not a generated pressure
trace. I used a deliberately bounded laminar cylinder case so I could inspect
the mesh, wake, force history and pressure probes and retain the complete
execution evidence. The second run halves the timestep with the mesh and
physical sampling unchanged. The result is a numerical sensitivity pair, not a
second velocity case or a full-scale instrument model.

“The movie shows authentic saved planes using their native connectivity and a
common colour scale. There are no fabricated intermediate flow fields. I retain
which geometrically degenerate triangles were masked for display without
altering the downloaded native arrays.

“The pressure convention matters. OpenFOAM stores kinematic pressure here, in
square metres per square second. Multiplying by the declared water density once
produces pascals. It remains local incompressible hydrodynamic pressure. It is
not a distant propagating acoustic pressure just because it has the same final
unit.”

**Evidence:** [setup](CFD.md), [native provenance and replay](CFD_RESULTS.md).

## 9. Vortex shedding, pressure spectra and numerical limits

**Visuals:** [pressure histories/RMS/PSD](../site/data/research/cfd/figures/hydrodynamic_pressure.png)
and [startup/timestep sensitivity](../site/data/research/cfd/figures/forces_startup_sensitivity.png).

**On the slide:** late interval [100,140) s; St **0.16663 → 0.16682**;
baseline P2 fluctuation RMS **0.01322 Pa**; mean Cd **1.36155 → 1.35837**.

**Speaking route:**

“The signed lift gives a repeatable way to identify shedding periods. I use
rising zero crossings and calculate Strouhal number from that frequency. I do
not assume the largest pressure spectral peak must be the shedding fundamental:
pressure can contain harmonics, and the finite spectral bins limit peak
resolution. In the late interval the lift PSD peak lies at 0.175 Hz on a 0.025 Hz
grid, while the crossing-based frequency is approximately 0.1666 Hz.

“Halving the timestep changes Strouhal by about 0.1155 percent, mean drag by about
minus 0.2333 percent, and lift RMS by about minus 2.4409 percent. Those are
useful observed sensitivities. They are not a convergence order, and agreement
in frequency does not establish that amplitudes are converged.

“The full startup is still shown. The late blocks retain some amplitude growth,
and the longer window changes the statistics. I would rather expose that
limitation than select one favourable interval and call the whole record
stationary. The result supports a local-flow contamination investigation; a
coupled sensor model and a more representative flow regime would be separate
pieces of work.”

**Evidence:** [full statistics and adverse stationarity evidence](CFD_RESULTS.md).

## 10. Demonstrate inspectability, not another calculation

**Live route:**

1. Open the existing site at Research and select the vessel passage.
2. Point to the recovered gain and the separately labelled sampling error.
3. Open the window/band CSV and identify one excluded interval.
4. Select the native cylinder, show the pressure figure and its hydrodynamic label.
5. Point to the native package and the documented saved-data replay command.

**Speaking route:**

“The browser is a viewer for computed products. Selecting a study does not launch
CFD or invent new data. Each displayed result has a downloadable numerical
counterpart. The source files, calibration and processing conventions are
traceable, and output manifests bind the bytes used in the demonstration.

“Reproducibility also has a scope. Replaying the saved native subset reruns the
analysis; it is not a new physical realization or a new solver run. The subset
checks the retained native files and explicitly lists what was omitted. That
separation makes it possible for someone else to inspect the result without
needing my original case directory, while avoiding a claim that a reduced
package contains every original volume field.”

**Offline fallback:** keep the nine original PNG figures and the movie available;
use the downloaded CSV/JSON if the network is unavailable. Do not substitute
illustrative data if a live page fails.

## 11. What I would measure next

**On the slide:** matched deployment data; known instrument response; environmental
and geometry metadata; independent evaluation; sensor-coupling investigation.

**Speaking route:**

“The demonstrated contribution is a working, inspectable chain: empirical
prediction, calibrated synthetic recording processing, quantitative comparison
and genuine native CFD kept in its proper physical role. To turn this into a
field-validation study I would need authorized recordings with calibration and
frequency response, reliable acquisition timing, source and receiver geometry,
environmental information and active-transmission annotations. I would then
compare matched windows and reserve independent deployments or encounters for
evaluation. In parallel, I would investigate how local pressure couples through
the actual mounting and sensor response. That would let the combined workflow
answer a measurement question without confusing software verification with
physical accuracy.”

## Questions to prepare for

- **Why not add CFD PSD to acoustic PSD?** Matching units is necessary but not
  sufficient. No established pressure-to-sensor transfer/coupling model or
  supporting correlation assumptions are available here.
- **Why a low-Reynolds-number 2D case?** It makes a bounded, traceable native
  demonstration and timestep sensitivity feasible. It does not represent
  three-dimensional turbulent broadband housing noise.
- **Why not BELLHOP immediately?** The current direct model is independently
  testable. A waveguide calculation needs appropriate environment/boundary data,
  implementation checks and an experiment that can distinguish its benefit.
- **Does +6 dB mean pressure doubles exactly?** No. The imposed pressure factor
  is $10^{6/20}\approx1.99526$; the power factor is $10^{6/10}\approx3.98107$.
- **Why is a synthetic residual not zero?** A random finite realization differs
  from its expected spectrum. The matched same-realization gain control removes
  that sampling difference; the ensemble comparison does not.
- **What did the software verification cover?** Use the [verification report](VERIFICATION.md)
  for exact accepted counts and identities. Do not quote an earlier integration
  count as final-tree acceptance or describe software tests as physical validation.
- **Is it tied to a commercial product?** No. It uses no proprietary decoder,
  product-specific response, hardware connection or manufacturer endorsement.

## Sources

[2] https://www.mit.edu/~millitsa/resources/pdfs/bwdx.pdf
[3] https://resource.npl.co.uk/acoustics/techguides/seaabsorption
    Official [NPL technical-guide index](https://www.npl.co.uk/research/underwater-acoustics/tech-guides)
    independently confirms this guide; the interactive endpoint could not be
    retrieved with valid TLS during verification. No certificate bypass was used.
[4] https://gorbatschow.github.io/SonarDocs/sound_absorption_sea_ainslie.en
[5] https://publications.tno.nl/publication/34637939/jkMYWh/macgillivray-2021-reference.pdf
