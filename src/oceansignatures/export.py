"""Deterministic public exports of executed synthetic acoustic experiments.

Builds from core modules and docs/sources.json only; never reads development caches.
The spectral convention is metrics.py's segment-mean removal, deliberately distinct
from the no-detrend audio pipeline. All raw pressure arrays retain full precision.
"""

from pathlib import Path
import io
import json
import zipfile

import numpy as np

from .metrics import welch_psd, spectrogram_psd
from .reporting import envelope, channel_transfer
from .passive import passive_case, propagate
from .active import active_case
from matplotlib.backends.backend_pdf import PdfPages
from scipy import signal as sp_signal
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

STYLE = {
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "figure.titlesize": 16,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.16,
    "savefig.facecolor": "white",
    "pdf.fonttype": 42,
}
BLUE, ORANGE, GREEN, GRAY = "#007c91", "#d25b24", "#657d24", "#727b87"


def _figure(rows=2):
    return plt.subplots(
        rows, 1, figsize=(10, 6.6 if rows == 2 else 5.8), layout="constrained", squeeze=False
    )


def _db(density):
    # Finite display floor only: saved linear densities are unchanged.
    return 10 * np.log10(np.maximum(density, 1e-24) / 1e-12)


def _extrema(ax, t, y, color, label, bins=1200):
    centre, low, high = envelope(t, y, bins)
    ax.vlines(centre, low, high, color=color, linewidth=0.7, label=label)


def _save_figure(output, stem, fig, title, caption, pdf=None, section=None):
    path = Path(output) / "plots" / (stem + ".png")
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160, metadata={"Software": "OceanSignatures " + VERSION})
    if pdf is not None:
        pdf.savefig(fig, dpi=160)
    plt.close(fig)
    result = {
        "src": "plots/" + path.name,
        "title": title,
        "caption": caption,
        "alt": title + ". " + caption,
    }
    if section:
        result["section"] = section
    return result


def render_case(output, entry, pdf=None):
    """Render exact exported arrays, with optional vector/mixed PDF pages."""
    output = Path(output)
    case_id = entry["id"]
    with np.load(output / "data" / (case_id + ".npz"), allow_pickle=False) as archive:
        a = {key: archive[key] for key in archive.files}
    m = json.loads((output / "data" / (case_id + ".json")).read_text())["model_metadata"]
    with plt.rc_context(STYLE):
        if "kind" in entry:
            return _passive_figures(output, entry, a, m, pdf)
        return _active_figures(output, entry, a, m, pdf)


def _passive_figures(output, entry, a, m, pdf):
    case_id = entry["id"]
    results = []
    fig, axes = _figure()
    source_n = m["source"]["emission_samples"]
    start, stop = m["analysis_window_samples"]
    times = a["time_s"]
    fig.suptitle(entry["title"] + " | pressure histories")
    ax, bx = axes[:, 0]
    _extrema(ax, times[:source_n], a["source_pa"][:source_n], BLUE, "Emission min–max")
    ax.set(title="1 m reference source · RMS over unpadded emission", ylabel="Source pressure (Pa)")
    _extrema(bx, times, a["received_pa"], BLUE, "Signal + filtered noise")
    _extrema(bx, times, a["received_clean_pa"], ORANGE, "Clean after receiver")
    bx.axvspan(0, start / FS, color=GRAY, alpha=0.12)
    bx.set(
        title=f"Receiver analysis [{start / FS:.5f}, {stop / FS:.5f}) s",
        xlabel="Time from emission start (s)",
        ylabel="Received pressure (Pa)",
    )
    bx.legend(loc="upper center", bbox_to_anchor=(0.5, -0.20), ncol=2, frameon=False)
    for axis in (ax, bx):
        axis.set_xlim(0, 2.31)
        axis.ticklabel_format(axis="y", style="sci", scilimits=(-2, 3))
    caption = (
        f"Start with the source: 1 Pa RMS over its full unpadded [0, {source_n / FS:g}) s emission, "
        f"including quiet click intervals. Then compare the receiver in Pa (a separate, per-case amplitude axis). "
        f"Each vertical stroke retains a block minimum and maximum; it is not decimated point sampling. "
        f"The receiver analysis is [{start}, {stop}) samples, excluding the shaded idle interval. "
        "Clean means after the 257-tap 10 kHz FIR; its 128-sample group delay is retained. "
        "A centred interpolation kernel may pre-ring numerically; this is not a physical precursor."
    )
    results.append(
        _save_figure(
            output,
            case_id + "-waveform",
            fig,
            "1 · Follow pressure from source to receiver",
            caption,
            pdf,
        )
    )

    fig, axes = _figure()
    fig.suptitle(entry["title"] + " | spectral decomposition")
    ax, bx = axes[:, 0]
    f = a["frequency_hz"] / 1000
    ax.plot(f, _db(a["source_psd_pa2_hz"]), color=BLUE, linewidth=1.1)
    ax.set(title="Source emission · no zero-padded idle time in the estimate")
    for key, color, label, style in (
        ("channel", GRAY, "Before receiver", "--"),
        ("clean", ORANGE, "Clean after FIR", "-"),
        ("noise", GREEN, "Filtered noise", "-"),
        ("received", BLUE, "Signal + noise", "-"),
    ):
        bx.plot(
            f, _db(a[key + "_psd_pa2_hz"]), color=color, linewidth=1, label=label, linestyle=style
        )
    bx.set(
        title="Receiver components · identical receiver analysis window", xlabel="Frequency (kHz)"
    )
    bx.legend(loc="upper center", bbox_to_anchor=(0.5, -0.20), ncol=4, frameon=False, fontsize=9)
    for axis in (ax, bx):
        axis.set(xlim=(0, 12), ylim=(-90, 115), ylabel="PSD (dB re 1 µPa²/Hz)")
        axis.axvline(10, color=GRAY, linewidth=0.8, linestyle=":")
    caption = (
        "Next separate the spectrum into source, pre-filter channel, clean filtered signal, filtered noise and their actual sum. "
        "The sum PSD includes finite-record cross terms; it is not manufactured by adding the two component PSDs. "
        "All cases share frequency/level axes. Welch: 2048-sample periodic Hann, 1024 overlap, segment-mean removal, "
        "one-sided Pa²/Hz; source uses only emission and receiver uses its metadata window. "
        "Spectral dB = 10 log10(PSD / 10⁻¹² Pa²/Hz); displayed 0–12 kHz and −90–115 dB, "
        "with a −120 dB computational display floor only. Exact full-Nyquist linear spectra are downloadable. "
        "Welch window weighting is not the full-record RMS convention."
    )
    results.append(
        _save_figure(
            output,
            case_id + "-psd",
            fig,
            "2 · Separate propagation, receiver and noise",
            caption,
            pdf,
        )
    )

    source_t, f, source_p = spectrogram_psd(a["source_pa"][:source_n], FS)
    receiver_t, _, receiver_p = spectrogram_psd(a["received_pa"][start:stop], FS)
    receiver_t += start / FS
    _npz(
        Path(output) / "data" / (case_id + "-spectrogram.npz"),
        {
            "source_frame_s": source_t,
            "received_frame_s": receiver_t,
            "frequency_hz": f,
            "source_psd_pa2_hz": source_p,
            "received_psd_pa2_hz": receiver_p,
        },
    )
    fig, axes = _figure()
    fig.suptitle(entry["title"] + " | time–frequency structure")
    for axis, frame, density, label in zip(
        axes[:, 0],
        (source_t, receiver_t),
        (source_p, receiver_p),
        ("Source emission", "Received signal + noise"),
    ):
        image = axis.pcolormesh(
            frame,
            f / 1000,
            _db(density),
            vmin=-20,
            vmax=100,
            cmap="magma",
            shading="auto",
            rasterized=True,
        )
        axis.set(title=label, xlim=(0, 2.31), ylim=(0, 12), ylabel="Frequency (kHz)")
        axis.grid(False)
    axes[-1, 0].set_xlabel("Frame-centre time from emission start (s)")
    fig.colorbar(image, ax=list(axes[:, 0]), label="PSD (dB re 1 µPa²/Hz)", pad=0.025)
    caption = (
        "Finally ask when each frequency occurs. These are density spectrograms, not independently normalized pictures: "
        "source and receiver, and every passive case, share −20 to 100 dB re 1 µPa²/Hz. "
        "Values outside this display interval saturate; exact linear matrices are supplied separately. "
        "Frames use 512-sample periodic Hann windows, 384 overlap, per-frame mean removal, no padding; "
        "times are frame centres on the emission clock. Source frames lie inside the unpadded emission; "
        "receiver frames lie inside its analysis window. Blank regions have no complete frame, not measured silence. "
        "The receiver delay, coloration and noise are not removed to improve appearance."
    )
    results.append(
        _save_figure(
            output,
            case_id + "-spectrogram",
            fig,
            "3 · Read the time–frequency signature",
            caption,
            pdf,
        )
    )
    return results


def _active_figures(output, entry, a, m, pdf):
    results = []
    case_id = entry["id"]
    fig, axes = _figure()
    fig.suptitle(entry["title"] + " | ping and echo")
    ax, bx = axes[:, 0]
    ping_n = m["ping_samples"]
    _extrema(ax, a["time_s"][:ping_n], a["transmitted_pa"][:ping_n], BLUE, "Reference ping")
    ax.set(
        title="1 m reference ping · 50 ms, Hann-tapered 3–9 kHz LFM",
        xlabel="Time from transmission (s)",
        ylabel="Reference pressure (Pa)",
        xlim=(0, 0.05),
    )
    _extrema(bx, a["time_s"], a["received_pa"], BLUE, "Echo + noise")
    _extrema(bx, a["time_s"], a["received_clean_pa"], ORANGE, "Clean echo")
    bx.axvspan(m["expected_delay_s"], m["expected_delay_s"] + 0.05, color=ORANGE, alpha=0.10)
    bx.set(
        title="Full acquisition · shaded true echo support",
        xlabel="Time from transmission (s)",
        ylabel="Received pressure (Pa)",
        xlim=(0, 1.05),
    )
    bx.ticklabel_format(axis="y", style="sci", scilimits=(-2, 3))
    bx.legend(loc="upper center", bbox_to_anchor=(0.5, -0.20), ncol=2, frameon=False)
    caption = (
        "Start with the known 1 Pa RMS ping over its unpadded 50 ms support, then its delayed scalar echo. "
        "The panels have different clocks and pressure scales; block min–max strokes retain extrema. "
        f"The analytic two-way delay is {m['expected_delay_s']:g} s and pressure gain is {m['two_way_pressure_gain']:.6g}. "
        "The ideal pressure sampler has no FIR, leakage, reverberation or multipath. "
        "Input SNR divides clean echo-support mean square by expected full-band noise variance; "
        "noise is readjusted for every chosen range/SNR, not held at fixed ambient level. "
        "No normalization is applied to the received waveform."
    )
    results.append(
        _save_figure(
            output,
            case_id + "-history",
            fig,
            "1 · Locate the echo inside a noisy acquisition",
            caption,
            pdf,
        )
    )
    fig, axes = _figure()
    fig.suptitle(entry["title"] + " | matched-filter range estimate")
    range_axis = a["lag_s"] * 1500 / 2
    score = a["matched_filter"]
    peak = m["peak_index"]
    for index, axis in enumerate(axes[:, 0]):
        if index == 0:
            _extrema(axis, range_axis, score, BLUE, "Signed score, block min–max")
        else:
            mask = abs(range_axis - entry["range_m"]) <= 0.65
            axis.plot(
                range_axis[mask], score[mask], color=BLUE, linewidth=1, marker=".", markersize=3
            )
        axis.axvline(
            entry["range_m"], color=ORANGE, linestyle="--", linewidth=1.1, label="True range"
        )
        if index == 0 or abs(range_axis[peak] - entry["range_m"]) <= 0.65:
            axis.plot(
                range_axis[peak],
                score[peak],
                "D",
                color=GREEN,
                markersize=6,
                label="Global |score| winner",
            )
        axis.set(ylabel="Signed score (noise SD)", ylim=(-30, 30))
    axes[0, 0].set(title="Complete 0–750 m search · no truth-based gate", xlim=(0, 750))
    axes[1, 0].set(
        title=f"Truth-centred display only · estimated {range_axis[peak]:.6f} m",
        xlim=(entry["range_m"] - 0.65, entry["range_m"] + 0.65),
        xlabel="Template-start lag × 1500 / 2 (m)",
    )
    handles, labels = axes[0, 0].get_legend_handles_labels()
    axes[1, 0].legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.20),
        ncol=3,
        frameon=False,
        fontsize=9,
    )
    caption = (
        "Now correlate against the entire known ping. The top panel retains min–max score extrema over the full search; "
        "the lower panel plots every local lag sample around truth, solely as a display close-up. "
        "The estimator always selects the global maximum of the absolute signed score over 0–750 m, never a truth-gated maximum. "
        "Score = correlation / (known simulation noise SD × square-root template energy), dimensionless; "
        "noise-only SD is one at each lag, but neighbouring lags are correlated. "
        "The signed oscillations and possible carrier-lobe error are retained. All six cases share ±30 score axes. "
        "This is a range estimator, not a thresholded detector or calibrated probability of detection."
    )
    results.append(
        _save_figure(
            output,
            case_id + "-matched-filter",
            fig,
            "2 · Search globally; inspect the peak locally",
            caption,
            pdf,
        )
    )
    return results


def _verification():
    """Run independent analytic controls now; never quote historical pytest counts."""
    checks = []

    def record(name, error, tolerance, description):
        if not np.isfinite(error) or error > tolerance:
            raise AssertionError(f"{name}: {error} > {tolerance}")
        checks.append(
            {
                "name": name,
                "status": "passed",
                "detail": f"{description} Maximum error {error:.8g}; acceptance ≤ {tolerance:g}.",
            }
        )

    frequency = np.linspace(0, 0.4 * FS, 401)
    impulse = np.zeros(160)
    impulse[0] = 1
    error = 0.0
    for fraction in (0.1, 0.35, 0.7, 0.95):
        delay = (100 + fraction) / FS
        y = propagate(impulse, FS, [{"delay_s": delay, "gain": 1.0}])
        actual = np.exp(-2j * np.pi * frequency[:, None] * np.arange(len(y))[None, :] / FS) @ y
        expected = np.exp(-2j * np.pi * frequency * delay)
        error = max(error, float(np.max(abs(actual - expected))))
    record(
        "Fractional-delay impulse response",
        error,
        5e-5,
        "Four fractional offsets; 401 frequencies from DC to 0.4 fs; complex error against an independent ideal delay.",
    )
    direct = np.sqrt(100**2 + (8 - 12) ** 2)
    reflected = np.sqrt(100**2 + (8 + 12) ** 2)
    expected = (
        np.exp(-2j * np.pi * frequency * direct / 1500) / direct
        - np.exp(-2j * np.pi * frequency * reflected / 1500) / reflected
    )
    record(
        "Continuous two-ray reference",
        float(np.max(abs(channel_transfer(frequency, 100) - expected))),
        1e-12,
        "Independent slant distances and signed phase at 100 m; continuous reference, not a sampled-channel equivalence claim.",
    )
    record(
        "Pressure-release boundary cancellation",
        float(np.max(abs(channel_transfer(frequency, 100, source_depth_m=0)))),
        1e-12,
        "Coincident direct/image path lengths with source depth zero must cancel.",
    )
    x = np.cos(2 * np.pi * 3000 * np.arange(8192) / FS)
    f, density = welch_psd(x, FS)
    record(
        "Welch tone density normalization",
        abs(float(np.sum(density) * (f[1] - f[0])) - 0.5),
        1e-12,
        "A unit-amplitude bin-centred 3 kHz cosine has independently known mean square 0.5 Pa².",
    )
    case = active_case(123.4567, 80, seed=SEED)
    error = abs(case["estimated_range_m"] - 123.4567)
    record(
        "Fractional active-delay range quantization",
        error,
        1500 / (4 * FS),
        "80 dB synthetic control at 123.4567 m; half an integer-lag range bin, not a low-SNR or field error bound.",
    )
    case = active_case(150, -18, seed=SEED)
    template = case["transmitted_pa"][: case["metadata"]["ping_samples"]]
    denominator = case["metadata"]["noise_std_pa"] * np.sqrt(np.sum(template**2))
    lags = [0, 123, 9600, 48000]
    error = max(
        abs(
            np.dot(case["received_pa"][k : k + len(template)], template) / denominator
            - case["matched_filter"][k]
        )
        for k in lags
    )
    record(
        "Signed correlation direct-dot check",
        float(error),
        1e-10,
        "Four valid template-start lags, including both search endpoints, independently evaluated by direct dot products.",
    )
    return {
        "status": "Analytic numerical controls passed",
        "summary": "These controls were executed by this build. They verify ideal-model calculations only; "
        "they are not field validation, a detector evaluation, or an asserted software-suite test count.",
        "checks": checks,
    }


def _robustness(output):
    columns = [
        "seed",
        "snr_db",
        "true_range_m",
        "estimated_range_m",
        "error_m",
        "peak_abs_score",
        "realized_input_snr_db",
        "noise_std_pa",
    ]
    rows = []
    summaries = []
    for snr in (-30, -18, -6):
        for seed in range(2026, 2050):
            case = active_case(150, snr, seed=seed)
            m = case["metadata"]
            rows.append(
                [
                    seed,
                    snr,
                    150,
                    case["estimated_range_m"],
                    m["range_error_m"],
                    m["peak_abs_score"],
                    m["realized_input_snr_db"],
                    m["noise_std_pa"],
                ]
            )
        error = np.array(rows[-24:])[:, 4]
        summaries.append(
            {
                "snr_db": snr,
                "trials": len(error),
                "within_0_125m": int(np.sum(abs(error) <= 0.125)),
                "rmse_all_trials_m": float(np.sqrt(np.mean(error**2))),
                "max_absolute_error_m": float(np.max(abs(error))),
            }
        )
    table = np.array(rows)
    np.savetxt(
        Path(output) / "data/active-robustness.csv",
        table,
        delimiter=",",
        header=",".join(columns),
        comments="",
        fmt="%.17g",
    )
    _npz(Path(output) / "data/active-robustness.npz", dict(zip(columns, table.T)))
    _json(
        Path(output) / "data/active-robustness.json",
        {
            "summaries": summaries,
            "seed_interval_inclusive": [2026, 2049],
            "fs_hz": FS,
            "truth_m": 150,
            "search_m": [0, 750],
            "criterion_m": 0.125,
            "criterion_meaning": "Illustrative nominal c/(2B) scale; not a detection threshold or measured two-target resolution.",
            "trial_policy": "Every seed retained; same seeds across SNRs share scaled Gaussian draws.",
            "interpretation": "Always-return-an-estimate experiment conditional on target presence, not field accuracy.",
        },
    )
    return table, summaries


def _gallery(output, passive, active, pdf):
    results = []
    frequency = np.linspace(50, 12000, 400)
    ranges = np.geomspace(25, 500, 320)
    direct = channel_transfer(frequency[None, :], ranges[:, None], surface=False)
    two_ray = channel_transfer(frequency[None, :], ranges[:, None])
    relative = 20 * np.log10(abs(two_ray / direct))
    _npz(
        Path(output) / "data/interference.npz",
        {
            "frequency_hz": frequency,
            "range_m": ranges,
            "two_ray_relative_to_direct_db": relative,
            "direct_transfer": direct,
            "two_ray_transfer": two_ray,
        },
    )
    fig, axes = _figure(1)
    ax = axes[0, 0]
    fig.suptitle("Why surface interference depends on both range and frequency")
    image = ax.contourf(
        frequency / 1000,
        ranges,
        relative,
        levels=np.linspace(-45, 6, 52),
        cmap="RdBu_r",
        norm=matplotlib.colors.TwoSlopeNorm(vmin=-45, vcenter=0, vmax=6),
        extend="both",
    )
    ax.set(xlabel="Frequency (kHz)", ylabel="Horizontal separation (m)", yscale="log")
    ax.set_yticks([25, 50, 100, 200, 400, 500], labels=["25", "50", "100", "200", "400", "500"])
    ax.grid(False)
    fig.colorbar(image, ax=ax, label="20 log₁₀ |H₂-ray / Hdirect| (dB)")
    caption = (
        "First isolate interference from spherical spreading: this continuous analytic contour divides the coherent "
        "direct-plus-pressure-release image transfer by the direct-path transfer. "
        "Source/receiver depths are 8/12 m, c = 1500 m/s, surface coefficient −1. "
        "Blue bands are cancellation, red bands reinforcement; white is 0 dB relative to direct, not zero pressure. "
        "The logarithmic range axis spans 25–500 m. This is an independently evaluated continuous two-ray reference, "
        "not the discrete propagator, receiver FIR, or a ray-traced ocean model. Values below −45 dB saturate."
    )
    results.append(
        _save_figure(
            output,
            "gallery-interference",
            fig,
            "Interference is a range–frequency structure",
            caption,
            pdf,
            "Channel geometry",
        )
    )

    frequency = np.linspace(0, 12000, 2401)
    hd = channel_transfer(frequency, 100, surface=False)
    ht = channel_transfer(frequency, 100)
    hr = ht - hd
    reference = json.loads((Path(output) / "data" / (passive[0]["id"] + ".json")).read_text())[
        "model_metadata"
    ]
    _, receiver = sp_signal.freqz(reference["receiver"]["fir_coefficients"], worN=frequency, fs=FS)
    _npz(
        Path(output) / "data/transfer.npz",
        {
            "frequency_hz": frequency,
            "direct": hd,
            "surface": hr,
            "coherent_sum": ht,
            "receiver_fir": receiver,
        },
    )
    fig, axes = _figure()
    fig.suptitle("Separate a coherent channel from a causal receiver")
    for h, label, color, style in (
        (hd, "Direct", BLUE, "-"),
        (hr, "Reflected alone", GRAY, "--"),
        (ht, "Coherent sum", ORANGE, "-"),
    ):
        axes[0, 0].plot(
            frequency / 1000,
            20 * np.log10(np.maximum(abs(h), 1e-12)),
            color=color,
            label=label,
            linestyle=style,
        )
    axes[0, 0].set(
        title="100 m horizontal range · phases matter before summation",
        ylabel="Pressure gain (dB re 1)",
        xlim=(0, 12),
        ylim=(-90, -20),
    )
    axes[0, 0].legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3, frameon=False)
    axes[1, 0].plot(frequency / 1000, 20 * np.log10(np.maximum(abs(receiver), 1e-12)), color=GREEN)
    axes[1, 0].set(
        title="Actual 257-tap receiver coefficients · 128-sample group delay retained",
        xlabel="Frequency (kHz)",
        ylabel="FIR gain (dB re 1)",
        xlim=(0, 12),
        ylim=(-120, 5),
    )
    caption = (
        "Read the upper panel as complex pressure addition: reflected and direct magnitudes are nearly coincident at 100 m, "
        "but the reflected path carries a sign reversal and a different delay. Their coherent sum has deep notches; "
        "adding magnitude curves would miss them. Gains reference the nominal 1 m source pressure. "
        "The lower panel is computed from the actual generic receiver coefficients, not a measured instrument: "
        "257-tap Kaiser FIR, 10 kHz cutoff and 128/48000 s group delay. "
        "The continuous channel does not include fractional-delay interpolation error; the receiver is a separate operation."
    )
    results.append(
        _save_figure(
            output,
            "gallery-transfer",
            fig,
            "Direct path, reflection and receiver response",
            caption,
            pdf,
            "Channel → receiver",
        )
    )

    passive_rows = []
    fig, axes = _figure()
    fig.suptitle("Passive range comparison | actual saved case metrics")
    for kind, color in zip(("whistle", "clicks", "vessel"), (BLUE, ORANGE, GREEN)):
        for surface, linestyle in ((False, "--"), (True, "-")):
            entries = [e for e in passive if e["kind"] == kind and e["surface"] == surface]
            meta = [
                json.loads((Path(output) / "data" / (e["id"] + ".json")).read_text())[
                    "model_metadata"
                ]
                for e in entries
            ]
            x = [e["range_m"] for e in entries]
            level = [20 * np.log10(m["measured"]["received_clean_rms_pa"] / 1e-6) for m in meta]
            snr = [m["measured"]["snr_db"] for m in meta]
            label = kind + (" · two paths" if surface else " · direct")
            axes[0, 0].semilogx(x, level, "o", color=color, linestyle=linestyle, label=label)
            axes[1, 0].semilogx(x, snr, "o", color=color, linestyle=linestyle, label=label)
            passive_rows += [
                [("whistle", "clicks", "vessel").index(kind), int(surface), r, lev, s]
                for r, lev, s in zip(x, level, snr)
            ]
    for axis in axes[:, 0]:
        axis.set_xticks([25, 100, 400], labels=["25", "100", "400"])
    axes[0, 0].set(
        title="Clean filtered receiver RMS · complete metadata analysis window",
        ylabel="Level (dB re 1 µPa)",
    )
    axes[1, 0].set(
        title="Signal-to-noise ratio · same window for clean signal and noise",
        xlabel="Horizontal separation (m)",
        ylabel="20 log₁₀(RMSsignal / RMSnoise) (dB)",
    )
    axes[1, 0].legend(
        loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=3, fontsize=9, frameon=False
    )
    np.savetxt(
        Path(output) / "data/passive-comparison.csv",
        passive_rows,
        delimiter=",",
        header="kind_index_whistle_clicks_vessel,surface,range_m,clean_level_db_re_1uPa,snr_db",
        comments="",
        fmt="%.17g",
    )
    caption = (
        "Move from an analytic channel to all 18 saved passive experiments. Markers are computed cases, "
        "lines only connect the three ranges; they are not a dense sampled performance curve. "
        "All sources have 1 Pa unpadded-record RMS, but different spectra and peak pressures. "
        "The receiver and input-referred noise SD remain fixed; RMS and SNR use each complete arrival/tail window. "
        "Solid curves include the ideal surface, dashed curves use only the direct path. "
        "Interference therefore changes broadband energy in a source-dependent way. Seeds are shared, not independent trials; "
        "there is no inference of detectability or field range."
    )
    results.append(
        _save_figure(
            output,
            "gallery-passive-range",
            fig,
            "Range changes level and interference together",
            caption,
            pdf,
            "Actual passive cases",
        )
    )

    fig, axes = _figure()
    fig.suptitle("Active range comparison | controlled SNR, not fixed ambient noise")
    active_rows = []
    for snr, color in ((-18, ORANGE), (-6, BLUE)):
        entries = [e for e in active if e["snr_db"] == snr]
        meta = [
            json.loads((Path(output) / "data" / (e["id"] + ".json")).read_text())["model_metadata"]
            for e in entries
        ]
        x = [e["range_m"] for e in entries]
        noise = [20 * np.log10(m["noise_std_pa"] / 1e-6) for m in meta]
        echo = [20 * np.log10(np.sqrt(m["echo_mean_square_pa2"]) / 1e-6) for m in meta]
        errors = [m["range_error_m"] for m in meta]
        axes[0, 0].plot(x, noise, "o--", color=color, label=f"Noise SD · {snr} dB case")
        if snr == -18:
            axes[0, 0].plot(x, echo, "s-", color=GREEN, label="Clean echo RMS · both SNRs")
        axes[1, 0].plot(x, errors, "o-", color=color, label=f"{snr} dB input SNR")
        active_rows += [[r, snr, n, e, err] for r, n, e, err in zip(x, noise, echo, errors)]
    axes[0, 0].set(
        title="Noise is readjusted at every range to preserve the selected SNR",
        ylabel="Amplitude level (dB re 1 µPa)",
    )
    axes[0, 0].legend(
        loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3, fontsize=9, frameon=False
    )
    axes[1, 0].set(
        title="One fixed seed (2026) · the estimator can pick a neighbouring lobe",
        xlabel="True point-target range (m)",
        ylabel="Estimated − true range (m)",
        ylim=(-0.08, 0.02),
    )
    axes[1, 0].legend(loc="upper center", bbox_to_anchor=(0.5, -0.20), ncol=2, frameon=False)
    for axis in axes[:, 0]:
        axis.set_xticks([75, 150, 300])
    np.savetxt(
        Path(output) / "data/active-comparison.csv",
        active_rows,
        delimiter=",",
        header="range_m,snr_db,noise_sd_db_re_1uPa,echo_rms_db_re_1uPa,error_m",
        comments="",
        fmt="%.17g",
    )
    caption = (
        "These six actual cases do not share a fixed ambient noise level: expected noise variance is recalculated "
        "from echo-support mean square and the selected input SNR. Both echo and noise pressure fall with the chosen "
        "two-way 1/R² gain at fixed SNR. The lower panel reports the actual global-search range error, including "
        "the 150 m, −18 dB neighbouring-lobe error. Shared seed 2026 is a controlled illustration, not a trial average. "
        "Zero errors at these grid-aligned ranges do not establish perfect ranging. The next plate repeats the "
        "experiment and exposes estimates far away from truth."
    )
    results.append(
        _save_figure(
            output,
            "gallery-active-range",
            fig,
            "Why fixed-SNR cases are not fixed-noise range performance",
            caption,
            pdf,
            "Actual active cases",
        )
    )

    table, summaries = _robustness(output)
    fig, axes = _figure()
    fig.suptitle("Repeated-noise range estimation | every trial retained")
    for snr, color, marker in ((-30, GRAY, "x"), (-18, ORANGE, "o"), (-6, BLUE, "+")):
        rows = table[table[:, 1] == snr]
        summary = next(s for s in summaries if s["snr_db"] == snr)
        label = f"{snr} dB: {summary['within_0_125m']}/24 within ±0.125 m"
        for axis in axes[:, 0]:
            axis.scatter(
                rows[:, 0] - 2026, rows[:, 4], color=color, marker=marker, s=38, label=label
            )
    axes[0, 0].set(
        title="All global-search errors · linear axis, no removed failures",
        ylabel="Estimated − true range (m)",
        ylim=(-175, 625),
        xlim=(-0.6, 23.6),
    )
    axes[1, 0].set(
        title="Same 72 errors · symmetric-log scale reveals small errors AND distant outliers",
        ylabel="Error (m), symlog",
        xlabel="Seed offset (seed = 2026 + offset)",
        yscale="symlog",
        ylim=(-175, 625),
        xlim=(-0.6, 23.6),
    )
    axes[1, 0].set_yscale("symlog", linthresh=0.02, linscale=0.8)
    axes[1, 0].axhspan(-0.125, 0.125, color=GREEN, alpha=0.08)
    axes[1, 0].legend(
        loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=3, fontsize=8.7, frameon=False
    )
    for axis in axes[:, 0]:
        axis.axhline(0, color=GRAY, linewidth=0.7)
        axis.set_xticks([0, 4, 8, 12, 16, 20, 23])
    summary_text = "; ".join(
        f"{s['snr_db']} dB: {s['within_0_125m']}/24 within ±0.125 m, all-trial RMSE {s['rmse_all_trials_m']:.3f} m"
        for s in summaries
    )
    caption = (
        "At true range 150 m, run every seed 2026–2049 at each SNR, always searching 0–750 m. "
        "Both panels show every estimate: linear and symmetric-log axes, with no success-only clipping. "
        + summary_text
        + ". The shaded ±0.125 m band is an illustrative c/(2B) scale, not a validated detection threshold. "
        "Matched noise seeds across SNRs share underlying draws. A target is present in every trial; there is no "
        "target-absent experiment or false-alarm calibration. These finite seeded synthetic results are not field accuracy, "
        "population performance, or detection probability. The catastrophic −18 dB outlier is included in RMSE."
    )
    results.append(
        _save_figure(
            output,
            "gallery-active-robustness",
            fig,
            "Do not hide low-SNR failures behind a median",
            caption,
            pdf,
            "Estimator limitations",
        )
    )
    return results


def build_site(output="site"):
    """Build the complete 18+6 case release; preserve frontend files in output.

    Repeated calls overwrite only the generated named data/plot artifacts. No
    internet, cache, reference recording, credentials or browser is required.
    """
    output = Path(output)
    (output / "data").mkdir(parents=True, exist_ok=True)
    (output / "plots").mkdir(parents=True, exist_ok=True)
    verification = _verification()
    passive, active, pages = [], [], []
    with (
        plt.rc_context(STYLE),
        PdfPages(
            output / "data/selected-plates.pdf",
            metadata={
                "Title": "OceanSignatures — selected explanatory plates",
                "Author": "OceanSignatures",
                "Subject": "Synthetic educational models; captions and order in plot-pack.json",
                "CreationDate": None,
                "ModDate": None,
            },
        ) as pdf,
    ):
        for kind in ("whistle", "clicks", "vessel"):
            for distance in (25, 100, 400):
                for surface in (False, True):
                    entry = export_passive_case(output, kind, distance, surface)
                    selected = kind == "whistle" and distance == 100 and surface
                    entry["figures"] = render_case(output, entry, pdf if selected else None)
                    entry["downloads"].append(
                        {
                            "href": f"data/{entry['id']}-spectrogram.npz",
                            "label": "Exact source/receiver spectrogram matrices (NPZ)",
                        }
                    )
                    if selected:
                        pages.extend(entry["figures"])
                    passive.append(entry)
        for distance in (75, 150, 300):
            for snr in (-18, -6):
                entry = export_active_case(output, distance, snr)
                selected = distance == 150 and snr == -18
                entry["figures"] = render_case(output, entry, pdf if selected else None)
                if selected:
                    pages.extend(entry["figures"])
                active.append(entry)
        gallery = _gallery(output, passive, active, pdf)
        pages.extend(gallery)
    _json(
        output / "data/plot-pack.json",
        {
            "format": "PDF; vector text/axes/curves and contour; spectrogram fields and colorbars rasterized at 160 dpi",
            "timestamp_policy": "PDF CreationDate and ModDate omitted",
            "pages": [dict(page_number=i + 1, **page) for i, page in enumerate(pages)],
        },
    )
    with zipfile.ZipFile(
        output / "data/full-waveforms.zip", "w", compression=zipfile.ZIP_STORED
    ) as archive:
        for entry in passive + active:
            name = entry["id"] + ".npz"
            info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
            info.external_attr = 0o644 << 16
            archive.writestr(info, (output / "data" / name).read_bytes())
    registry = json.loads((Path(__file__).resolve().parents[2] / "docs/sources.json").read_text())
    fallback_titles = {
        "1": "SciPy: Welch power spectral density",
        "2": "SciPy: cross-correlation",
        "3": "PAMGuard: simulated sound acquisition",
        "4": "DOSITS: sonar equation",
    }
    sources = [
        {"id": str(s["id"]), "title": s["title"] or fallback_titles[str(s["id"])], "url": s["url"]}
        for s in registry["sources"]
    ]
    catalog = {
        "schema_version": 1,
        "title": "OceanSignatures",
        "version": VERSION,
        "provenance_note": "Executed deterministic synthetic experiments, 48 kHz and seed 2026 unless stated. "
        "All pressure and spectral arrays are saved without display normalization. "
        "No experimental or field validation; no reference recording is bundled. "
        "Rebuild from public core modules, scripts and source registry; no development cache needed.",
        "passive_cases": passive,
        "active_cases": active,
        "reference": None,
        "gallery": gallery,
        "verification": verification,
        "sources": sources,
        "downloads": [
            {
                "href": "data/full-waveforms.zip",
                "label": "All 24 full-resolution waveform NPZs (ZIP)",
            },
            {
                "href": "data/selected-plates.pdf",
                "label": "Selected explanatory plates (10-page vector/mixed PDF)",
            },
            {
                "href": "data/plot-pack.json",
                "label": "Plot-pack page order and full interpretive captions (JSON)",
            },
            {
                "href": "data/active-robustness.csv",
                "label": "All 72 repeated-noise estimates, including outliers (CSV)",
            },
            {
                "href": "data/active-robustness.json",
                "label": "Robustness definitions and all-trial statistics (JSON)",
            },
            {"href": "data/catalog.json", "label": "Complete computed catalogue (JSON)"},
            {
                "href": "data/artifact-manifest.json",
                "label": "Generated asset SHA-256 manifest (JSON)",
            },
        ],
    }
    _json(output / "data/catalog.json", catalog)
    import hashlib
    import scipy

    manifest = {
        "versions": {
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__,
        },
        "files": [],
    }
    for folder in ("data", "plots"):
        for path in sorted((output / folder).iterdir()):
            if path.is_file() and path.name != "artifact-manifest.json":
                manifest["files"].append(
                    {
                        "path": str(path.relative_to(output)),
                        "bytes": path.stat().st_size,
                        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                    }
                )
    _json(output / "data/artifact-manifest.json", manifest)
    return catalog


VERSION = "0.1.0"
SEED = 2026
FS = 48000


def _json(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n")


def _npz(path, arrays):
    """NPZ with fixed ZIP headers, ordered members and no pickled objects."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
        for name, array in sorted(arrays.items()):
            buffer = io.BytesIO()
            np.lib.format.write_array(buffer, np.asarray(array), allow_pickle=False)
            member = zipfile.ZipInfo(name + ".npy", date_time=(1980, 1, 1, 0, 0, 0))
            member.compress_type = zipfile.ZIP_DEFLATED
            member.external_attr = 0o644 << 16
            archive.writestr(member, buffer.getvalue(), compresslevel=6)


def _save_case(output, case_id, case, signals, source_window, receiver_window):
    output = Path(output)
    stem = output / "data" / case_id
    arrays = {k: v for k, v in case.items() if isinstance(v, np.ndarray)}
    arrays["fs_hz"] = np.array(case["fs"])
    columns = []
    names = ["frequency_hz"]
    for name, (key, window) in signals.items():
        f, density = welch_psd(case[key][slice(*window)], case["fs"])
        if not columns:
            columns.append(f)
            arrays["frequency_hz"] = f
        columns.append(density)
        names.append(name + "_pa2_hz")
        arrays[name + "_psd_pa2_hz"] = density
    _npz(stem.with_name(stem.name + ".npz"), arrays)
    np.savetxt(
        stem.with_name(stem.name + ".csv"),
        np.column_stack(columns),
        delimiter=",",
        header=",".join(names),
        comments="",
        fmt="%.17g",
    )
    _json(
        stem.with_name(stem.name + ".json"),
        {
            "schema_version": 1,
            "export_version": VERSION,
            "case_id": case_id,
            "fs_hz": case["fs"],
            "model_metadata": case["metadata"],
            "spectral_method": {
                "window": "periodic Hann",
                "nperseg": 2048,
                "noverlap": 1024,
                "detrend": "constant",
                "scaling": "one-sided density Pa^2/Hz",
                "incomplete_segments": "discarded; no padding",
                "source_window_samples": source_window,
                "receiver_window_samples": receiver_window,
                "note": "Different from the audio pipeline no-detrend convention.",
            },
            "provenance": "Executed synthetic model; no measured or field-validation data.",
            "timestamp_policy": "No wall-clock timestamps; ZIP members fixed to 1980-01-01.",
        },
    )
    return [
        {"href": f"data/{case_id}.{ext}", "label": label}
        for ext, label in (
            ("npz", "Full-resolution pressure, clock and spectrum arrays (NPZ)"),
            ("csv", "One-sided Welch spectrum components (Pa²/Hz CSV)"),
            ("json", "Exact model, support and analysis metadata (JSON)"),
        )
    ]


def export_passive_case(output, kind, range_m, surface):
    """Compute and export one 48 kHz, two-second passive case, seed 2026."""
    case = passive_case(kind, range_m, surface, fs=FS, seed=SEED)
    case_id = f"passive-{kind}-{range_m:.17g}m-surface-{str(surface).lower()}"
    source_window = [0, case["metadata"]["source"]["emission_samples"]]
    window = case["metadata"]["analysis_window_samples"]
    downloads = _save_case(
        output,
        case_id,
        case,
        {
            "source": ("source_pa", source_window),
            "channel": ("received_channel_pa", window),
            "clean": ("received_clean_pa", window),
            "noise": ("noise_pa", window),
            "received": ("received_pa", window),
        },
        source_window,
        window,
    )
    from scipy.io import wavfile

    gain = 0.8 / max(float(np.max(abs(case["source_pa"]))), float(np.max(abs(case["received_pa"]))))
    listening = {
        "digital_fraction_per_pa": gain,
        "physical_playback_calibration": False,
        "encoding": "mono signed PCM16; round(digital_fraction * 32767)",
        "clock": "full shared emission/receiver clock, including zero-padded source tail",
        "interpretation": "Same gain for source and receiver preserves attenuation up to PCM quantization. "
        "Received audio can be very quiet; no loudness normalization. "
        "Playback pressure depends on hardware and volume and is not calibrated.",
    }
    for name in ("source", "received"):
        href = f"data/{case_id}-{name}-same-gain.wav"
        digital = np.rint(case[name + "_pa"] * gain * 32767).astype(np.int16)
        wavfile.write(Path(output) / href, FS, digital)
        listening[name + "_href"] = href
        downloads.append(
            {
                "href": href,
                "label": f"Synthetic {name} listening WAV · same gain, not playback-calibrated",
            }
        )
    metadata_path = Path(output) / "data" / (case_id + ".json")
    metadata = json.loads(metadata_path.read_text())
    metadata["listening"] = listening
    _json(metadata_path, metadata)
    return {
        "id": case_id,
        "kind": kind,
        "range_m": range_m,
        "surface": surface,
        "title": f"{kind.capitalize()} at {range_m:g} m — "
        + ("two paths" if surface else "direct path"),
        "summary": "Synthetic 1 Pa RMS source at a nominal 1 m; 8 m source and 12 m receiver depths. "
        "A generic causal 10 kHz receiver filters signal and noise. No field validation.",
        "metrics": [
            {
                "label": "Analysis-window SNR",
                "value": case["metadata"]["measured"]["snr_db"],
                "unit": "dB",
            },
            {
                "label": "Source emission RMS",
                "value": case["metadata"]["measured"]["source_rms_pa"],
                "unit": "Pa",
            },
        ],
        "figures": [],
        "downloads": downloads,
    }


def export_active_case(output, range_m, snr_db):
    """Compute and export an active case without changing the search or noise."""
    from .active import active_case

    case = active_case(range_m, snr_db, fs=FS, seed=SEED)
    case_id = f"active-{range_m:.17g}m-snr-{snr_db:.17g}db"
    source_window = [0, case["metadata"]["ping_samples"]]
    window = case["metadata"]["echo_support_samples"]
    downloads = _save_case(
        output,
        case_id,
        case,
        {
            "source": ("transmitted_pa", source_window),
            "clean": ("received_clean_pa", window),
            "noise": ("noise_pa", window),
            "received": ("received_pa", window),
        },
        source_window,
        window,
    )
    return {
        "id": case_id,
        "range_m": range_m,
        "snr_db": snr_db,
        "title": f"{range_m:g} m point target · {snr_db:g} dB input SNR",
        "summary": "A finite 50 ms Hann-tapered 3–9 kHz reference ping and scalar target gain 1. "
        "Range is the global absolute signed-correlation maximum, not a detection. "
        "Noise is readjusted for each range/SNR, not fixed ambient noise.",
        "metrics": [
            {"label": "Estimated range", "value": case["estimated_range_m"], "unit": "m"},
            {"label": "Range error", "value": case["metadata"]["range_error_m"], "unit": "m"},
            {
                "label": "Realized echo-support SNR",
                "value": case["metadata"]["realized_input_snr_db"],
                "unit": "dB",
            },
        ],
        "figures": [],
        "downloads": downloads,
    }
