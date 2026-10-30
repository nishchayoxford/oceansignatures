"""Reproducible empirical prediction / synthetic-recording verification.

This is not field validation. See docs/RESEARCH.md for the synthesis ensemble,
finite-window interpretation, resource bounds and export schema.
"""

from copy import deepcopy
from dataclasses import asdict, replace
from datetime import datetime, timedelta, timezone
import csv
import hashlib
import io
import json
import os
import platform
from importlib.metadata import version
from pathlib import Path
import tempfile
import zipfile

import numpy as np
from scipy.io import wavfile

from .comparison import compare_band_levels, compare_prediction_observation
from .noise_prediction import SeaWater, predict_passage, wind_psd
from .observations import RecordingChunk, analyse_chunks
from .pipeline import Calibration, calibrate_samples

FS = 8192
FRAME = 16384
NPERSEG = 4096
NFRAMES = 120
ORIGIN = "2000-01-01T00:00:00.000000Z"  # Arbitrary synthetic epoch, not acquisition time.
BANDS = {
    "100–500 Hz": (100.0, 500.0),
    "500–1000 Hz": (500.0, 1000.0),
    "1000–4000 Hz": (1000.0, 4000.0),
}
CALIBRATION = Calibration(-180.0, 20.0, 2.0)


def _utc(seconds):
    anchor = datetime(2000, 1, 1, tzinfo=timezone.utc)
    return (
        (anchor + timedelta(seconds=float(seconds)))
        .isoformat(timespec="microseconds")
        .replace("+00:00", "Z")
    )


def _predict(frequency, times, *, wind_offset_m_s=0.0, source_offset_db=0.0, **changes):
    params = dict(
        speed_kn=12.0,
        length_m=150.0,
        closest_approach_m=250.0,
        source_depth_m=8.0,
        receiver_depth_m=50.0,
        sound_speed_m_s=1500.0,
        closest_approach_emission_time_s=120.0,
        water=SeaWater(),
        wind_speed_m_s=4.0,
        vessel_class="bulker",
    )
    params.update(changes)
    result = predict_passage(frequency, times, **params)
    winds = 4 + 4 * np.asarray(times) / 240 + wind_offset_m_s
    result["vessel_psd_pa2_hz"] *= 10 ** (source_offset_db / 10)
    result["source_psd_pa2_m2_hz"] *= 10 ** (source_offset_db / 10)
    result["wind_psd_pa2_hz"] = np.array([wind_psd(frequency, wind_speed_m_s=w) for w in winds])
    result["total_psd_pa2_hz"] = result["vessel_psd_pa2_hz"] + result["wind_psd_pa2_hz"]
    result["wind_speed_m_s"] = winds
    result["metadata"]["parameters"]["wind_speed_m_s"] = (
        "declared linear 4 to 8 m/s over receiver seconds 0 to 240; held at each frame midpoint"
    )
    result["metadata"]["parameters"]["wind_offset_m_s"] = wind_offset_m_s
    result["metadata"]["parameters"]["source_offset_db"] = source_offset_db
    return result


def _jsonable(value):
    if isinstance(value, np.ndarray):
        return _jsonable(value.tolist())
    if isinstance(value, np.generic):
        return _jsonable(value.item())
    if isinstance(value, dict):
        return {k: _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def _json(path, data):
    path.write_text(json.dumps(_jsonable(data), sort_keys=True, indent=2, allow_nan=False) + "\n")


def _npz(path, arrays):
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for name, values in sorted(arrays.items()):
            buffer = io.BytesIO()
            np.lib.format.write_array(buffer, np.asarray(values), allow_pickle=False)
            entry = zipfile.ZipInfo(name + ".npy", (1980, 1, 1, 0, 0, 0))
            entry.compress_type = zipfile.ZIP_DEFLATED
            entry.external_attr = 0o600 << 16
            archive.writestr(entry, buffer.getvalue())


def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def comparison_batch(analysis, *, as_prediction=False):
    """Adapt calibrated observations on 100–4000 Hz using qc_ok, never valid alone.

    as_prediction=True is reserved for the synthetic matched-realization control,
    not permission to relabel a measurement as an empirical model. No resampling,
    time-origin conversion, interpolation or invented calibration is performed.
    """
    metadata = analysis["metadata"]
    if (
        analysis["unit"] != "Pa"
        or metadata["psd_units"] != "Pa^2/Hz"
        or analysis["calibration_status"] != "declared scalar calibration"
    ):
        raise ValueError("Require physical pressure units and declared calibration")
    kinds = {"synthetic": "synthetic_observation", "measurement": "measured_observation"}
    if metadata["evidence_kind"] not in kinds:
        raise ValueError("Unknown observation evidence kind")
    if as_prediction and metadata["evidence_kind"] != "synthetic":
        raise ValueError("Only a synthetic matched realization can be the prediction control")
    mask = (analysis["frequency_hz"] >= 100) & (analysis["frequency_hz"] <= 4000)
    return dict(
        frequency_hz=analysis["frequency_hz"][mask].copy(),
        window_start_s=analysis["window_start_s"].copy(),
        window_stop_s=analysis["window_stop_s"].copy(),
        psd=analysis["psd"][:, mask].copy(),
        valid_window=analysis["qc_ok"].copy(),
        time_origin_utc=metadata["anchor_utc"],
        psd_unit="Pa2/Hz",
        evidence_kind="predicted_acoustic" if as_prediction else kinds[metadata["evidence_kind"]],
    )


def expected_welch_psd(frame_psd, *, fs=FS, nperseg=NPERSEG):
    """Exact mean Hann periodogram for the finite circular Gaussian ensemble.

    Input last axis is the one-sided PSD on an even frame's full rFFT grid.
    Segment length must divide that frame length; detrending is NOT performed.
    Overlap changes estimator variance, not its mean. DC/Nyquist retain weight 1.
    """
    raw = np.asarray(frame_psd)
    if np.iscomplexobj(raw) or raw.ndim < 1 or raw.shape[-1] < 3:
        raise ValueError("Require a real full one-sided frame PSD")
    p = np.asarray(raw, dtype=float)
    n = 2 * (p.shape[-1] - 1)
    if (
        not np.isfinite(p).all()
        or np.any(p < 0)
        or type(nperseg) is not int
        or not 2 <= nperseg <= n
        or n % nperseg
        or nperseg % 2
        or not np.isscalar(fs)
        or not np.isfinite(fs)
        or fs <= 0
    ):
        raise ValueError("Finite nonnegative PSD, positive fs and even dividing segment required")
    # Normalize each ensemble before halving interior bins or summing in the FFT.
    # Restore density only after the one-sided weights, rounding subnormals once.
    scale = np.max(p, axis=-1, keepdims=True)
    unit = p / np.where(scale == 0, 1.0, scale)
    two_sided = np.concatenate(
        [unit[..., :1], unit[..., 1:-1] / 2, unit[..., -1:], unit[..., -2:0:-1] / 2], axis=-1
    )
    w = 0.5 - 0.5 * np.cos(2 * np.pi * np.arange(nperseg) / nperseg)
    kernel = np.abs(np.fft.fft(w, n=n)) ** 2 / (n * np.sum(w * w))
    smoothed = np.fft.ifft(np.fft.fft(two_sided, axis=-1) * np.fft.fft(kernel), axis=-1).real
    result = np.maximum(smoothed[..., np.arange(nperseg // 2 + 1) * (n // nperseg)], 0)
    result[..., 1:-1] *= 2
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        density = result * scale
    if not np.isfinite(density).all() or np.any((result > 0) & (density == 0)):
        raise ValueError("Expected Welch PSD is not representable as float64")
    return density


def _portable_metadata(analysis, root):
    metadata = deepcopy(analysis["metadata"])
    for chunk in metadata["chunks"]:
        chunk["path"] = Path(chunk["path"]).relative_to(root).as_posix()
    return metadata


def _write_chunks(root, pressure, gain_db):
    scale = float(calibrate_samples([1.0], CALIBRATION)[0])
    edges = np.array([0, 17.375, 53.125, 80, 82, 109.625, 150.5, 190.875, 216.25, 240])
    edges = (edges * FS).astype(int)
    recordings = {}
    for family, gain in (("baseline", 1.0), ("observed", 10 ** (gain_db / 20))):
        directory = root / "wav" / family
        directory.mkdir(parents=True)
        chunks = []
        for i, (a, b) in enumerate(zip(edges[:-1], edges[1:])):
            if a == 80 * FS:  # Omit an entire 2-s chunk, never zero-fill it.
                continue
            selected = gain * pressure[a:b] / scale
            if np.max(np.abs(selected)) >= 0.999:
                raise ValueError("Synthetic digital waveform reaches full scale")
            other = 0.05 * np.sin(2 * np.pi * 333 * np.arange(a, b) / FS)
            path = directory / f"synthetic-{i:02d}.wav"
            # scipy emits deterministic IEEE float64 WAV; no wall-clock PEAK chunk.
            wavfile.write(path, FS, np.column_stack([other, selected]))
            chunks.append(
                RecordingChunk(
                    path,
                    _utc(a / FS),
                    channel=1,
                    calibration=CALIBRATION,
                    evidence_kind="synthetic",
                    source="seeded receiver surrogate; not a recording of the ocean",
                    expected_sha256=_sha(path),
                )
            )
        recordings[family] = chunks
    return recordings


def _analyse(chunks):
    return analyse_chunks(
        chunks,
        anchor_utc=ORIGIN,
        frame_seconds=2.0,
        nperseg=NPERSEG,
        noverlap=2048,
        window="hann",
        detrend=False,
        bands_hz=BANDS,
        active_intervals=[(_utc(150.25), _utc(153.125))],
        max_windows=120,
        max_output_cells=2_000_000,
    )


def _sensitivity(root, frequency, times):
    cases = [{"label": "baseline", "changes": {}}]
    for name, values in (
        ("source_offset_db", (-6.0, 6.0)),
        ("wind_offset_m_s", (-2.0, 2.0)),
        ("speed_kn", (10.0, 14.0)),
        ("closest_approach_m", (125.0, 500.0)),
        ("receiver_depth_m", (25.0, 100.0)),
        ("temperature_c", (5.0, 15.0)),
        ("salinity_ppt", (30.0, 40.0)),
        ("depth_m", (0.0, 100.0)),
        ("ph", (7.8, 8.2)),
    ):
        cases.extend({"label": f"{name} = {value:g}", "changes": {name: value}} for value in values)
    powers = []
    for case in cases:
        changes = case["changes"]
        if set(changes) & set(asdict(SeaWater())):
            changes = {"water": replace(SeaWater(), **changes)}
        prediction = _predict(frequency, times, **changes)
        powers.append(
            np.column_stack(
                [
                    prediction["total_psd_pa2_hz"][:, (frequency >= low) & (frequency < high)].sum(
                        axis=1
                    )
                    * 2
                    for low, high in BANDS.values()
                ]
            )
        )
    powers = np.asarray(powers)
    delta = 10 * np.log10(powers.mean(axis=1) / powers[0].mean(axis=0))
    _npz(
        root / "sensitivity.npz",
        dict(
            receiver_time_s=times,
            band_names=np.array(list(BANDS)),
            band_power_pa2=powers,
            time_mean_power_change_db=delta,
        ),
    )
    _json(
        root / "sensitivity.json",
        dict(
            cases=cases,
            bands_hz=BANDS,
            band_names=list(BANDS),
            interpretation="One-at-a-time deterministic alternatives, not probabilistic confidence intervals.",
            estimator="Raw empirical midpoint PSD; sum bins in [low,high) times 2 Hz. Time mean of linear powers before dB; not smoothed Welch or mean of dB.",
            units="band_power_pa2: Pa2; time_mean_power_change_db: dB relative to baseline",
        ),
    )
    return cases, delta


def _figures(root, a, cases, sensitivity, gain_db):
    import matplotlib as mpl
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    directory = root / "figures"
    directory.mkdir()
    figures = []
    f, t = a["comparison_frequency_hz"], a["receiver_time_s"]
    selected = (a["frequency_hz"] >= 100) & (a["frequency_hz"] <= 4000)

    def db(p):
        return 10 * np.log10(p) + 120

    def figure(rows, height):
        fig = Figure(figsize=(10, height), layout="constrained")
        FigureCanvasAgg(fig)
        return fig, np.atleast_1d(fig.subplots(rows, 1))

    def save(fig, name, title, caption):
        fig.suptitle(title + "\nPredictive approximation / synthetic verification", fontsize=13)
        fig.savefig(directory / name, dpi=140, metadata={"Software": "OceanSignatures research"})
        figures.append(dict(path=f"figures/{name}", title=title, caption=caption))
        fig.clear()

    def exclusions(ax):
        ax.axvspan(80, 82, color="0.6", alpha=0.35)
        ax.axvspan(150, 154, color="#e69f00", alpha=0.18)

    with mpl.rc_context(
        rc={
            **mpl.rcParamsDefault,
            "font.size": 10,
            "axes.grid": True,
            "grid.alpha": 0.2,
            "font.family": "DejaVu Sans",
        }
    ):
        fig, ax = figure(2, 7)
        ax[0].semilogx(f, db(a["source_psd_pa2_m2_hz"]), color="#0072b2")
        ax[0].set(
            ylabel="Source PSD level\n(dB re 1 µPa² m²/Hz)",
            title="150 m bulker, 12 kn — loss-free monopole strength",
        )
        i = int(np.argmin(a["distance_m"]))
        for key, label in (
            ("vessel_psd_pa2_hz", "Vessel at receiver"),
            ("wind_psd_pa2_hz", "Wind at receiver"),
            ("total_psd_pa2_hz", "Incoherent sum"),
        ):
            ax[1].semilogx(f, db(a[key][i]), label=label)
        ax[1].set(
            xlabel="Frequency (Hz)",
            ylabel="Received PSD level\n(dB re 1 µPa²/Hz)",
            title=f"Near-CPA receiver snapshot at {t[i]:g} s",
        )
        ax[1].legend()
        save(
            fig,
            "spectra.png",
            "Source strength is not received pressure",
            "Empirical source spectrum and separate received vessel/wind powers on 100–4000 Hz. Source includes m²; ambient is not spread again. No field measurements.",
        )

        fig, ax = figure(2, 8)
        im = ax[0].pcolormesh(
            t,
            f,
            db(a["expected_welch_psd_pa2_hz"][:, selected]).T,
            shading="nearest",
            cmap="viridis",
        )
        ax[0].set(
            ylabel="Frequency (Hz)", title="Expected finite-Hann PSD; seven segments per frame"
        )
        fig.colorbar(im, ax=ax[0], label="dB re 1 µPa²/Hz")
        for j, name in enumerate(BANDS):
            ax[1].plot(t, db(a["expected_band_power_pa2"][:, j]), label=name)
        ax[1].set(xlabel="Synthetic receiver time (s)", ylabel="Expected band level\n(dB re 1 µPa)")
        ax[1].legend(loc="upper left")
        twin = ax[1].twinx()
        twin.plot(t, a["wind_speed_m_s"], "k--", alpha=0.5)
        twin.set(ylabel="Declared wind (m/s)", ylim=(3, 9))
        save(
            fig,
            "passage.png",
            "A chosen passage with changing wind",
            "240 s, 2 s frames; CPA is at emission second 120. Wind follows a declared 4→8 m/s schedule, not weather data. Spectrogram is the exact ensemble mean for the finite Fourier/Hann surrogate.",
        )

        fig, ax = figure(2, 7)
        ax[0].plot(t, a["distance_m"], label="Slant range")
        ax[0].set(ylabel="Slant range (m)")
        twin = ax[0].twinx()
        twin.plot(t, a["travel_time_s"], "--", color="#d55e00")
        twin.set(ylabel="Travel time (s)")
        for frequency in (100, 1000, 4000):
            ax[1].plot(
                t,
                a["propagation_loss_db"][:, np.flatnonzero(f == frequency)[0]],
                label=f"{frequency} Hz",
            )
        ax[1].set(xlabel="Receiver arrival time (s)", ylabel="Vessel propagation loss (dB)")
        ax[1].legend()
        save(
            fig,
            "propagation.png",
            "Retarded geometry and direct-path loss",
            "Range solves arrival = emission + range/1500. Loss is spherical spreading plus full-path Ainslie–McColm absorption. No waveguide, reflections, refraction or Doppler.",
        )

        fig, ax = figure(3, 9)
        for j, name in enumerate(BANDS):
            ax[j].plot(
                t, db(a["expected_band_power_pa2"][:, j]), label="Empirical → expected Welch"
            )
            ax[j].plot(
                t,
                np.where(a["qc_ok"], db(a["baseline_band_power_pa2"][:, j]), np.nan),
                label="Matched unscaled realization",
                alpha=0.8,
            )
            ax[j].plot(
                t,
                np.where(a["qc_ok"], db(a["observed_band_power_pa2"][:, j]), np.nan),
                label=f"Synthetic observed ({gain_db:+g} dB)",
            )
            active = a["active_transmission"]
            ax[j].plot(
                t[active],
                db(a["observed_band_power_pa2"][active, j]),
                "x",
                color="#d55e00",
                label="Active: excluded",
            )
            exclusions(ax[j])
            ax[j].set(ylabel="Band level\n(dB re 1 µPa)", title=name)
        ax[0].legend(fontsize=8, ncol=2)
        ax[-1].set_xlabel("Synthetic receiver time (s)")
        save(
            fig,
            "comparison.png",
            "Compare like support, then inspect the discrepancy",
            "Half-open bands and identical clocks/processing. Gray: absent 80–82 s chunk; amber: active-contaminated frames 150–154 s. Gaps are retained NaN, not zero. Crosses retain excluded active spectra.",
        )

        fig, ax = figure(3, 9)
        im = ax[0].pcolormesh(
            t, f, a["empirical_residual_db"].T, shading="nearest", cmap="coolwarm", vmin=-4, vmax=14
        )
        ax[0].set(
            ylabel="Frequency (Hz)", title="Observed − expected Welch; display scale −4 to +14 dB"
        )
        fig.colorbar(im, ax=ax[0], label="Residual (dB)", extend="both")
        for j, name in enumerate(BANDS):
            ax[1].plot(t, a["sampling_band_residual_db"][:, j], label=name)
        ax[1].axhline(0, color="k", lw=0.7)
        ax[1].set(
            ylabel="Baseline − expected (dB)",
            title="Finite-realization sampling error, without injected gain",
        )
        ax[1].legend(ncol=3, fontsize=9)
        ax[2].plot(
            t,
            np.mean(a["matched_band_residual_db"] - gain_db, axis=1),
            label="Matched band residual − injected gain",
        )
        ax[2].set(
            xlabel="Synthetic receiver time (s)",
            ylabel="Gain recovery error (dB)",
            ylim=(-1e-11, 1e-11),
        )
        for item in ax:
            exclusions(item)
        save(
            fig,
            "residuals.png",
            "Gain recovery is not model validation",
            "Top: stochastic observed-minus-ensemble spectra (color limits only; arrays are not clipped). Middle: independent band sampling errors. Bottom: same-realization gain verification, excluded intervals shaded. No field-validation inference.",
        )

        fig, ax = figure(1, 9)
        image = ax[0].imshow(sensitivity[1:], aspect="auto", cmap="coolwarm", vmin=-6, vmax=6)
        ax[0].set_yticks(np.arange(len(cases) - 1), [c["label"] for c in cases[1:]])
        ax[0].set_xticks(np.arange(len(BANDS)), list(BANDS))
        for row in range(len(cases) - 1):
            for column in range(len(BANDS)):
                ax[0].text(
                    column,
                    row,
                    f"{sensitivity[row + 1, column]:+.3g}",
                    ha="center",
                    va="center",
                    fontsize=8,
                )
        ax[0].grid(False)
        fig.colorbar(image, ax=ax[0], label="Time-mean band-power change (dB)", extend="both")
        save(
            fig,
            "sensitivity.png",
            "One factor at a time — not confidence intervals",
            "Alternatives change only the named property. Source offsets affect vessel only; speed changes source law and track. Receiver depth changes geometry, not absorption depth. These use raw empirical midpoint PSD bin sums and linear-power time averages, not Welch sampling or probability coverage.",
        )
    return figures


def _compute(root, seed, gain_db):
    times = np.arange(NFRAMES) * 2 + 1.0
    synthesis_frequency = np.fft.rfftfreq(FRAME, 1 / FS)
    supported = (synthesis_frequency >= 100) & (synthesis_frequency <= 4000)
    prediction = _predict(synthesis_frequency[supported], times)
    density = np.zeros((NFRAMES, len(synthesis_frequency)))
    density[:, supported] = prediction["total_psd_pa2_hz"]
    rng = np.random.Generator(np.random.PCG64(seed))
    coefficients = rng.standard_normal(density.shape) + 1j * rng.standard_normal(density.shape)
    coefficients *= np.sqrt(FRAME * FS * density / 4)
    pressure = np.fft.irfft(coefficients, n=FRAME, axis=1).reshape(-1)
    # Controlled active contamination; its entire affected frames are excluded.
    a, b = int(150.25 * FS), int(153.125 * FS)
    pressure[a:b] += 0.3 * np.sin(2 * np.pi * 1200 * np.arange(b - a) / FS)
    chunks = _write_chunks(root, pressure, gain_db)
    baseline = _analyse(chunks["baseline"])
    observed = _analyse(chunks["observed"])
    digital = _analyse([replace(c, calibration=None) for c in chunks["observed"]])
    predicted_batch, observed_batch = (
        comparison_batch(baseline, as_prediction=True),
        comparison_batch(observed),
    )
    matched = compare_prediction_observation(predicted_batch, observed_batch)
    bands = compare_band_levels(predicted_batch, observed_batch, BANDS)
    expected = expected_welch_psd(density)
    selected = (observed["frequency_hz"] >= 100) & (observed["frequency_hz"] <= 4000)
    model_batch = {
        **predicted_batch,
        "psd": expected[:, selected],
        "valid_window": np.ones(NFRAMES, dtype=bool),
    }
    empirical = compare_prediction_observation(model_batch, observed_batch)
    sampling = compare_prediction_observation(model_batch, comparison_batch(baseline))
    empirical_bands = compare_band_levels(model_batch, observed_batch, BANDS)
    sampling_bands = compare_band_levels(model_batch, comparison_batch(baseline), BANDS)
    central = _predict(observed["frequency_hz"][selected], times)
    broadband = selected & (observed["frequency_hz"] < 4000)
    sampling_broadband = 10 * np.log10(
        baseline["psd"][:, broadband].sum(axis=1) / expected[:, broadband].sum(axis=1)
    )
    smoothing = 10 * np.log10(expected[:, selected] / density[:, :: FRAME // NPERSEG][:, selected])
    arrays = dict(
        frequency_hz=observed["frequency_hz"],
        band_names=np.array(list(BANDS)),
        window_start_s=observed["window_start_s"],
        window_stop_s=observed["window_stop_s"],
        baseline_psd_pa2_hz=baseline["psd"],
        observed_psd_pa2_hz=observed["psd"],
        digital_psd=digital["psd"],
        valid=observed["valid"],
        qc_ok=observed["qc_ok"],
        coverage_fraction=observed["coverage_fraction"],
        active_transmission=observed["active_transmission"],
        near_full_scale_samples=observed["near_full_scale_samples"],
        partial_window=observed["partial_window"],
        comparison_frequency_hz=matched["frequency_hz"],
        matched_residual_db=matched["observed_minus_predicted_db"],
        matched_band_residual_db=bands["observed_minus_predicted_db"],
    )
    arrays.update(
        raw_expected_psd_pa2_hz=density[:, :: FRAME // NPERSEG],
        expected_welch_psd_pa2_hz=expected,
        empirical_residual_db=empirical["observed_minus_predicted_db"],
        sampling_residual_db=sampling["observed_minus_predicted_db"],
        sampling_band_residual_db=sampling_bands["observed_minus_predicted_db"],
        empirical_band_residual_db=empirical_bands["observed_minus_predicted_db"],
        expected_band_power_pa2=empirical_bands["predicted_power_pa2"],
        observed_band_power_pa2=observed["band_power"],
        baseline_band_power_pa2=baseline["band_power"],
        receiver_time_s=times,
        **{
            k: central[k]
            for k in (
                "source_psd_pa2_m2_hz",
                "vessel_psd_pa2_hz",
                "wind_psd_pa2_hz",
                "total_psd_pa2_hz",
                "propagation_loss_db",
                "travel_time_s",
                "distance_m",
                "emission_time_s",
                "wind_speed_m_s",
            )
        },
    )
    _npz(root / "research.npz", arrays)
    metadata = dict(
        baseline=_portable_metadata(baseline, root),
        observed=_portable_metadata(observed, root),
        digital_control=_portable_metadata(digital, root),
        prediction=prediction["metadata"],
        synthetic_clock="arbitrary 2000 epoch; not an acquisition date",
        missing_chunk_s=[80.0, 82.0],
        calibration=asdict(CALIBRATION),
    )
    metadata.update(
        schema="oceansignatures.research.metadata.v1",
        provenance={
            "python": platform.python_version(),
            "libraries": {
                name: version(name) for name in ("numpy", "scipy", "matplotlib", "soundfile")
            },
            "source_sha256": {
                f"src/oceansignatures/{name}": _sha(Path(__file__).with_name(name))
                for name in (
                    "research.py",
                    "noise_prediction.py",
                    "observations.py",
                    "comparison.py",
                    "pipeline.py",
                    "passive.py",
                )
            },
        },
        synthesis=dict(
            seed=seed,
            generator="NumPy PCG64; real then imaginary full-frame arrays",
            frame_samples=FRAME,
            frame_count=NFRAMES,
            frequency_step_hz=0.5,
            positive_support_hz=[100.0, 4000.0],
            edges="both synthesis endpoint bins included",
            coefficient="X[k] = sqrt(N*fs*S1[k]/4) * (normal(0,1) + i*normal(0,1)); zero DC/Nyquist",
            transition="Independent periodic Gaussian frames concatenated abruptly, no crossfade; not globally bandlimited. No Welch segment crosses a frame seam.",
            normalization="Ensemble normalization only; no realized RMS rescaling",
            active_contamination=dict(
                interval_s=[150.25, 153.125],
                frequency_hz=1200.0,
                amplitude_pa=0.3,
                treatment="Same baseline contamination is scaled with the injected gain; intersecting frames excluded",
            ),
            gain_db=gain_db,
            gain_interpretation="Synthetic pressure amplitude multiplier before declared ADC conversion; not an instrument fault measurement",
        ),
        expected_estimator=dict(
            method="Circular convolution of two-sided finite-frame density with abs(FFT(zero-padded periodic Hann))^2/(N*sum(w^2)); sample every N/nperseg bins and restore one-sided weights",
            mean="Exact ensemble mean for this finite circular surrogate, including support-edge leakage",
            overlap="Seven 4096-sample segments, step 2048; correlation affects variance, not mean",
            detrend=False,
        ),
        array_units={
            "frequency_hz": "Hz",
            "comparison_frequency_hz": "Hz",
            "receiver_time_s": "s",
            "window_start_s": "s",
            "window_stop_s": "s",
            "emission_time_s": "s",
            "travel_time_s": "s",
            "distance_m": "m",
            "source_psd_pa2_m2_hz": "Pa2 m2/Hz",
            "digital_psd": "digital2/Hz",
            "all_other_psd_arrays": "Pa2/Hz",
            "band_power_arrays": "Pa2",
            "residual_arrays": "dB",
            "propagation_loss_db": "dB",
            "wind_speed_m_s": "m/s",
            "coverage_fraction": "dimensionless",
        },
        serialization="Missing spectra/residuals are NaN in NPZ, blank in CSV, null in JSON. Complete active spectra are retained; use qc_ok for comparisons.",
    )
    _json(root / "metadata.json", metadata)
    with (root / "windows.csv").open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            [
                "start_s",
                "stop_s",
                "valid",
                "qc_ok",
                "coverage_fraction",
                "active_transmission",
                "band",
                "baseline_pa2",
                "observed_pa2",
                "observed_minus_baseline_db",
                "expected_pa2",
                "observed_minus_expected_db",
                "baseline_minus_expected_db",
            ]
        )
        for i in range(NFRAMES):
            for j, name in enumerate(BANDS):
                writer.writerow(
                    [
                        arrays[k][i]
                        for k in (
                            "window_start_s",
                            "window_stop_s",
                            "valid",
                            "qc_ok",
                            "coverage_fraction",
                            "active_transmission",
                        )
                    ]
                    + [
                        name,
                        *[
                            "" if np.isnan(v) else format(v, ".17g")
                            for v in (
                                baseline["band_power"][i, j],
                                observed["band_power"][i, j],
                                bands["observed_minus_predicted_db"][i, j],
                                empirical_bands["predicted_power_pa2"][i, j],
                                empirical_bands["observed_minus_predicted_db"][i, j],
                                sampling_bands["observed_minus_predicted_db"][i, j],
                            )
                        ],
                    ]
                )
    report = dict(
        title="Passage prediction to synthetic calibrated recording",
        kind="predictive approximations + synthetic verification",
        configuration=dict(
            seed=seed,
            gain_db=gain_db,
            sample_rate_hz=FS,
            frame_seconds=2,
            duration_s=240,
            nperseg=NPERSEG,
            noverlap=2048,
            window="hann",
            detrend=False,
            time_origin_utc=ORIGIN,
            bands_hz=BANDS,
            band_names=list(BANDS),
        ),
        numericalsummary=dict(
            compared_windows=matched["summary"]["compared_windows"],
            excluded_windows=matched["summary"]["excluded_windows"],
            gain_bias_db=matched["summary"]["mean_bias_db"],
            gain_max_error_db=float(
                np.nanmax(np.abs(matched["observed_minus_predicted_db"] - gain_db))
            ),
            observation_output_cells=int(observed["psd"].size + observed["stft"].size),
        ),
        limitations=["Synthetic verification, not field validation."],
        figures=[],
        downloads=[],
        source_refs=[],
    )
    report["numericalsummary"].update(
        sampling_broadband_rmse_db=float(
            np.sqrt(np.mean(sampling_broadband[observed["qc_ok"]] ** 2))
        ),
        sampling_broadband_mean_db=float(np.mean(sampling_broadband[observed["qc_ok"]])),
        sampling_spectral_summary=sampling["summary"],
        empirical_spectral_summary=empirical["summary"],
        sampling_band_rmse_db=np.sqrt(
            np.nanmean(sampling_bands["observed_minus_predicted_db"] ** 2, axis=0)
        ),
        sampling_band_mean_db=np.nanmean(sampling_bands["observed_minus_predicted_db"], axis=0),
        welch_smoothing_max_abs_db=float(np.max(np.abs(smoothing))),
    )
    edge_error = max(
        float(
            np.max(
                np.abs(
                    10
                    * np.log10(
                        _predict(central["frequency_hz"], times + shift)["total_psd_pa2_hz"]
                        / central["total_psd_pa2_hz"]
                    )
                )
            )
        )
        for shift in (-1.0, 1.0)
    )
    report["configuration"].update(
        prediction["metadata"]["parameters"],
        water=asdict(SeaWater()),
        calibration=asdict(CALIBRATION),
        selected_channel=1,
    )
    report["numericalsummary"].update(
        quasi_stationary_edge_max_abs_db=edge_error,
        minimum_snapshot_range_m=float(np.min(central["distance_m"])),
        travel_time_min_s=float(np.min(central["travel_time_s"])),
        travel_time_max_s=float(np.max(central["travel_time_s"])),
        expected_band_level_min_db=(
            10 * np.log10(arrays["expected_band_power_pa2"].min(axis=0)) + 120
        ),
        expected_band_level_max_db=(
            10 * np.log10(arrays["expected_band_power_pa2"].max(axis=0)) + 120
        ),
    )
    cases, sensitivity = _sensitivity(root, central["frequency_hz"], times)
    report["figures"] = _figures(root, arrays, cases, sensitivity, gain_db)
    report["limitations"] = [
        "Predictive approximations plus synthetic software verification, never measured or field validated.",
        "Empirical source and wind fits; generic constant water, not a site forecast. No proprietary hardware or NOAA audio.",
        "Direct spherical propagation; no reflections, bathymetry, refraction, directionality, tonals, Doppler, flow/self-noise or CFD coupling.",
        "Independent 2 s circular Gaussian frames, midpoint-frozen PSD; abrupt seams, no crossfade. Welch segments never cross seams. The concatenation is not globally bandlimited.",
        "100–4000 Hz finite-frame synthesis only; no low-frequency extrapolation. Hann leakage is included in expected spectra, including outside synthesis support.",
        "Exact expected Welch is for this discrete Fourier ensemble, not a claim that empirical continuous spectra equal one finite realization.",
        "Same-realization gain control cancels sampling error; comparison against ensemble expectation does not. No stochastic RMS normalization.",
        "Calibration is a declared flat synthetic ADC chain, not a certificate. Digital-only control has no absolute SPL.",
        "Sensitivity alternatives are one-at-a-time deterministic changes, not probabilistic confidence intervals.",
        "Synthetic UTC uses an arbitrary 2000 epoch and is not a real acquisition-date claim.",
    ]
    report["source_refs"] = [
        {
            "label": "Stojanovic 2007, Eq. 6: wind fit",
            "url": "https://www.mit.edu/~millitsa/resources/pdfs/bwdx.pdf",
        },
        {
            "label": "MacGillivray and de Jong 2021, Eqs. 9–11 / Table 1: vessel model",
            "url": "https://publications.tno.nl/publication/34637939/jkMYWh/macgillivray-2021-reference.pdf",
        },
        {
            "label": "NPL: absorption domain guidance",
            "url": "https://resource.npl.co.uk/acoustics/techguides/seaabsorption",
        },
        {
            "label": "Ainslie–McColm equation transcription (secondary)",
            "url": "https://gorbatschow.github.io/SonarDocs/sound_absorption_sea_ainslie.en",
        },
    ]
    report = _jsonable(report)
    report["schema"] = "oceansignatures.research.overview.v1"
    report["downloads"] = [
        dict(path=name, label=label)
        for name, label in (
            ("research.npz", "Prediction, observed, control and residual arrays (NPZ)"),
            ("windows.csv", "Window/band comparison and QC (CSV)"),
            ("metadata.json", "Calibration, estimator and SHA256 source provenance"),
            ("sensitivity.npz", "One-factor sensitivity arrays (NPZ)"),
            ("sensitivity.json", "Sensitivity cases and estimator contract"),
            ("manifest.json", "SHA256 manifest of every exported file except itself"),
        )
    ]
    report["downloads"] += [
        dict(path=p.relative_to(root).as_posix(), label=f"Synthetic {p.parent.name} WAV: {p.name}")
        for p in sorted((root / "wav").rglob("*.wav"))
    ]
    _json(root / "overview.json", report)
    _json(
        root / "manifest.json",
        dict(
            schema="oceansignatures.research.manifest.v1",
            algorithm="SHA256",
            files=_inventory(root)[0],
        ),
    )
    return report


def _safe_output(path):
    if ".." in path.parts:
        raise ValueError("Output path must not contain '..'")
    if any(p.is_symlink() for p in (path, *path.parents)):
        raise ValueError("Output path or ancestor is a symlink")


def _inventory(root):
    files, directories = {}, []
    if not root.is_dir():
        raise FileExistsError("Existing research output is not a directory")
    for path in sorted(root.rglob("*")):
        if path.is_symlink():
            raise ValueError("Export contains a symlink")
        name = path.relative_to(root).as_posix()
        if path.is_file():
            files[name] = _sha(path)
        elif path.is_dir():
            directories.append(name)
        else:
            raise ValueError("Export contains a nonregular file")
    return files, directories


def export_research(output_directory, *, seed=1729, gain_db=6.0):
    """Write a NEW bundle or accept an exactly identical rerun; return overview.

    Fixed 240 s scenario; seed is uint32 and gain_db is finite within ±12 dB.
    Existing differing/incomplete output, extra files and symlinks are rejected.
    All generation is staged beside the destination. Single-writer offline use;
    protection does not promise safety against adversarial concurrent renames.
    """
    if type(seed) is not int or not 0 <= seed <= 2**32 - 1:
        raise ValueError("seed must be an integer in [0, 2**32-1]")
    if (
        isinstance(gain_db, (bool, np.bool_))
        or not np.isscalar(gain_db)
        or not np.isreal(gain_db)
        or not np.isfinite(gain_db)
        or not -12 <= gain_db <= 12
    ):
        raise ValueError("gain_db must be a finite real scalar within ±12 dB")
    gain_db = float(gain_db)
    output = Path(output_directory).absolute()
    _safe_output(output)
    existing = _inventory(output) if output.exists() else None
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".research-stage-", dir=output.parent) as temporary:
        staged = Path(temporary) / "products"
        staged.mkdir()
        report = _compute(staged, seed, gain_db)
        _safe_output(output)
        if existing is not None:
            if _inventory(staged) != existing or _inventory(output) != existing:
                raise FileExistsError("Refusing differing or incomplete existing research output")
            return report
        if output.exists():
            raise FileExistsError("Research output appeared during staging")
        os.rename(staged, output)
    return report
