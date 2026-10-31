# Environmental and vessel acoustic prediction

`oceansignatures.noise_prediction` computes deterministic expected, one-sided pressure power spectral densities (PSDs). It does not generate field measurements, waveforms, CFD pressure, or instrument self-noise. This separate extension leaves the existing passive teaching model unchanged.

## Units and empirical equations

### Wind: ambient receiver PSD, not a point source

The wind/surface term is the empirical fit reproduced in Stojanovic (2007), Section II.B, Equation (6), with frequency `g=f_Hz/1000` in kHz and wind speed `U` in m/s:[2]

```text
Lw = 50 + 7.5 sqrt(U) + 20 log10(g) - 40 log10(g + 0.4)
Pw = 10^((Lw - 120)/10) Pa²/Hz
```

`Lw` is interpreted as a power spectral-density level, **dB re 1 µPa²/Hz**. The source text abbreviates the unit as “dB re µPa per Hz” despite defining a power spectral density and using `10 log`; this implementation makes the squared-pressure reference explicit.[2] `Pw` is already ambient noise at the receiver: it is **never** subjected to vessel spreading or absorption. Calm wind retains the fit's nonzero intercept; it is not silence. This is only the wind/surface component, not the full turbulence/shipping/thermal sum in that source.[2]

The paper identifies 100 Hz–100 kHz as the wind-dominated region, not a universal calibrated validity box.[2] The API adopts that conservative frequency envelope and a **software-only** wind guard of 0–20 m/s. No wind-height conversion, site-depth adjustment or universal wind-speed validation is claimed. Setting `allow_low_frequency_extrapolation=True` explicitly permits 20–100 Hz, labels these bins in passage metadata, and does not make them validated wind predictions. Other ambient sources can dominate there.[2]

### Vessel: JOMOPANS-ECHO source PSD

Implemented directly from MacGillivray and de Jong (2021), printed pages 4–5, Equations (9)–(11), and page 7 Table 1. The primary article PDF was retrieved from TNO; the equation page was rendered to check exponents.[5]

Let `f` be the numerical frequency in Hz, `V` the speed in knots, `l` the length in metres, `Vc` the class reference speed, and `l0=300 ft=91.44 m`. The baseline above 100 Hz (and at every frequency for non-cargo classes) is:[5]

```text
f1 = 480 / Vc
D = 4 for cruise, otherwise 3
L0 = 191 - 20 log10(f1) - 10 log10((1 - f/f1)² + D²)
Ls = L0 + 60 log10(V/Vc) + 20 log10(l/91.44)
Qs = 10^((Ls - 120)/10) Pa² m²/Hz
```

For bulker, containership, vehicle carrier and tanker, **replace**, not add to, the baseline below 100 Hz:[5]

```text
fLF = 600 / Vc
DLF = 0.8 for bulker/containership, 1.0 for vehicle_carrier/tanker
L0 = 208 - 40 log10(fLF) + 10 log10(f)
     - 10 log10((1 - (f/fLF)²)² + DLF²)
```

The source spectrum reference is **dB re 1 µPa² m²/Hz**: loss-free monopole strength, numerically the pressure PSD extrapolated to nominal 1 m. It is not a physical near-field hydrophone measurement. The published equations are spectral-density levels; their conversion to decidecade band levels adds `10 log10(0.231 f)`.[5] This implementation returns PSD, never silently the band level.

| API class | Vc (kn) |
|---|---:|
| fishing | 6.4 |
| tug | 3.7 |
| naval | 11.1 |
| recreational | 10.6 |
| government (government/research) | 8.0 |
| cruise | 17.1 |
| passenger | 9.7 |
| bulker | 13.9 |
| containership | 18.0 |
| vehicle_carrier | 15.8 |
| tanker | 12.4 |
| other | 7.4 |

These values are the published Table 1 values; callers supply class explicitly rather than inferring AIS types.[5] Dredgers are deliberately rejected because the publication adds a separate dredging-state speed substitution; this API does not represent that operating state.[5]

The supported range is 20 Hz–20 kHz, matching the reported model/data comparison; below 20 Hz unresolved blade-rate tones are a stated limitation.[5] The reported approximately 6 dB residual standard deviation is a population source-model statistic, **not** a receiver prediction confidence interval.[5] Speed `(0,60] kn` and length `(0,500] m` are numerical/software guards, not empirical validation ranges. Zero speed is rejected: taking the cruising power law to zero is not a model for idle machinery. No current correction is applied; the demonstration assumes stationary water so trajectory speed and model speed agree. No guaranteed source-spectrum continuity is imposed across the published 100 Hz branch switch.

## Absorption and propagation

Ainslie–McColm absorption is evaluated with `g=f_Hz/1000`, temperature `T` °C, salinity `S` ppt, path-average hydrostatic depth `d=depth_m/1000` km and pH. The following coefficients are transcribed from the accessible equation/implementation reference (secondary source); its implementation explicitly converts metres to kilometres before the pressure factors.[4]

```text
f1 = 0.78 sqrt(S/35) exp(T/26)
f2 = 42 exp(T/17)
a = 0.106 exp((pH-8)/0.56) f1 g²/(f1²+g²)
  + 0.52 (1+T/43) (S/35) exp(-d/6) f2 g²/(f2²+g²)
  + 0.00049 exp(-(T/27+d/17)) g²                 [dB/km]
TL = 20 log10(r / 1 m) + a r_m/1000             [dB]
Pr_vessel = Qs / r_m² * 10^(-a r_m/10000)       [Pa²/Hz]
```

The NPL technical guide documents the original Ainslie–McColm reference and a 100 Hz–1 MHz comparison envelope; its approximately 10% agreement with Francois–Garrison is documented for **one-at-a-time** environmental variations, not every joint corner of a rectangular parameter box.[3] Guards are `-6<T<35`, `5<S<50`, `7.7<pH<8.3`, and `0<=depth_m<7000` (zero depth is the source's surface baseline).[3] The depth parameter is not source-to-receiver distance. A constant representative water state is assumed; there is no vertical path integration. The optional 20–100 Hz continuation is explicitly an extrapolation of absorption too, not a supported accuracy claim.

Spherical spreading with absorption is a single unobstructed homogeneous-medium approximation.[2] Absorption spans the complete source-to-receiver distance, not `r-1 m`, because source strength is loss-free extrapolated. Standalone `propagation_loss` accepts scalar slant distance from 1 m to 10,000 km (upper bound is a numerical guard, **not** a useful range of spherical propagation). Near-field source dimensions, waveguide boundaries and bathymetry are not resolved.

## Motion and clocks

The receiver is fixed at `(0,0,zr)`, and the vessel follows `(v*(te-tCPA), CPA,zs)`. Here `v=V*1852/3600 m/s`; depths increase downward. Solve `tr=te+r(te)/c` for emission time at each requested **receiver arrival** time. For `q=tr-tCPA` and `b²=CPA²+(zs-zr)²`, the implemented subsonic root is:

```text
te = tCPA + (c²q - sqrt(c²v²q² + (c²-v²)b²)) / (c²-v²)
r = sqrt(v²(te-tCPA)² + b²)
travel_time = r/c
```

The existing `passive.path_parameters(..., surface=False)` supplies actual slant paths and delays. The existing `reporting.channel_transfer` was inspected: it provides a lossless coherent image response, not the required absorbing direct PSD, so it is not incorrectly reused. An optional image path is **not implemented**: coherent interference and distinct retarded image geometry would need separate treatment. The present model omits surface/bottom reflections, refraction, scattering and source directionality.

These are quasi-stationary **instantaneous expected PSD snapshots**, not STFT/Welch window averages. Retarded geometry does not implement Doppler: frequency warping, moving-source amplitude factors and time-dilation Jacobians are omitted. A comparison must choose windows short enough for geometry to remain approximately stationary, or average the expected spectrum over the exact observation window. The source and wind are assumed incoherent; linear component powers are summed, never dB values. A component can be disabled explicitly; zero wind speed alone does not disable ambient.

## Exact APIs and result contract

```python
SeaWater(temperature_c=10., salinity_ppt=35., depth_m=50., ph=8.)
wind_psd(frequency_hz, *, wind_speed_m_s,
         allow_low_frequency_extrapolation=False)
vessel_source_psd(frequency_hz, *, speed_kn, length_m, vessel_class="bulker")
absorption_db_per_km(frequency_hz, *, water=None,
                     allow_low_frequency_extrapolation=False)
propagation_loss(frequency_hz, *, distance_m, water=None,
                 allow_low_frequency_extrapolation=False)
predict_passage(frequency_hz, receiver_time_s, *, speed_kn, length_m,
                closest_approach_m, wind_speed_m_s, vessel_class="bulker",
                source_depth_m=8., receiver_depth_m=50., sound_speed_m_s=1500.,
                closest_approach_emission_time_s=0., water=None,
                include_wind=True, include_vessel=True,
                allow_low_frequency_extrapolation=False)
db_to_psd(level_db)
psd_to_db(psd)
psd_to_band_power(psd, bandwidth_hz)
band_power_to_psd(band_power, bandwidth_hz)
fractional_octave_bandwidth(centre_hz, *, fraction=3)
```

All frequency inputs are finite, nonempty, strictly increasing 1-D arrays; uniform grids are not required. Frequencies outside the function's envelope are rejected, including DC. Time arrays are finite, increasing relative seconds, not epoch timestamps. Scalar physical inputs are finite and Boolean values are rejected. Positive CPA is required; exact overhead tracks are outside this API. Component switches must be Boolean. Even disabled components' configuration is validated. `water=None` selects the explicitly documented default, independent of geometric receiver depth.

`predict_passage` returns a dictionary:

| Key | Shape / meaning |
|---|---|
| `frequency_hz` | `[F]` Hz |
| `receiver_time_s` | `[T]` receiver arrival seconds |
| `emission_time_s`, `travel_time_s`, `distance_m` | `[T]` retarded geometry |
| `total_psd_pa2_hz` | `[T,F]`, wind + vessel |
| `wind_psd_pa2_hz`, `vessel_psd_pa2_hz` | separate `[T,F]` receiver components |
| `source_psd_pa2_m2_hz` | `[F]` source strength, retained even if disabled |
| `propagation_loss_db` | `[T,F]` vessel path loss |
| `metadata` | JSON-safe parameters, models, clocks, units, omissions, extrapolated bins |

The dB conversion uses a `1e-12` squared-pressure reference; zero linear PSD maps to `-inf` without flooring. Band conversions require explicit positive bandwidth and assume flat PSD inside each band. The exact base-10 decidecade width is `fc*(10^0.05-10^-0.05)`; the published `0.231*fc` is its rounded approximation.[5] Integrate a varying PSD over actual band edges instead of treating a single central value as exact broadband power. Source band power retains m², while receiver band power is Pa²; their dB references must remain distinct.

## Executed example and verification

The reproducible cache-local driver `compute_passage.py` produces `passage.npz`, `passage.json`, `sensitivity.npz` and `sensitivity.json` under `.cache/research-extension/noise/`. The baseline uses a 150 m bulker at 12 kn, 250 m horizontal CPA, 8/50 m source/receiver depths and 6 m/s wind. Frequencies 20–4000 Hz in 2 Hz increments match a subset of positive bins for `fs=8192 Hz`, `nfft=4096`; 20–100 Hz extrapolation is explicitly enabled. Receiver times span −300 to +300 s in 5 s increments. Nothing outside the selected frequency band is synthesized or invented.

Sensitivity products vary wind, speed, CPA, receiver depth, temperature, salinity, absorption depth and pH one at a time; ±6 dB source offsets illustrate the source-model residual scale without claiming probability coverage. The summaries integrate continuous PSD with the trapezoidal rule over 20–4000 Hz; they are not rectangular FFT-bin powers. No random seed or observations are involved.

Tests retain red/green evidence for successive vertical slices: wind benchmark; finite/domain and calm behavior; independent JOMOPANS resonance and speed/length scaling; spectral-reference/bandwidth conversion; independent three-term absorption and range doubling; receiver-clock/slant geometry and component invariance. The test suite verifies zero switches, rejected invalid inputs, and that ambient receiver noise is not spread a second time.

## Sources

[2] https://www.mit.edu/~millitsa/resources/pdfs/bwdx.pdf
[3] https://resource.npl.co.uk/acoustics/techguides/seaabsorption
[4] https://gorbatschow.github.io/SonarDocs/sound_absorption_sea_ainslie.en
[5] https://publications.tno.nl/publication/34637939/jkMYWh/macgillivray-2021-reference.pdf
