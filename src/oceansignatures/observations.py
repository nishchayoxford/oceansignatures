"""Timestamped offline WAV observations, separate from legacy recording APIs."""

from dataclasses import asdict, dataclass, field
from datetime import datetime, timedelta, timezone
import csv
import hashlib
import json
import os
import tempfile
import zipfile
import io
from itertools import islice
import re
from pathlib import Path

import numpy as np
import soundfile as sf
from scipy.signal import get_window

from .pipeline import Calibration, calibrate_samples


@dataclass(frozen=True)
class RecordingChunk:
    """A declared first-sample UTC and selected channel; no inferred timing."""

    path: Path | str
    start_utc: datetime | str
    channel: int = 0
    calibration: Calibration | None = None
    evidence_kind: str = field(kw_only=True)
    source: str | None = None
    expected_sha256: str | None = None


def _utc(value):
    if isinstance(value, str):
        if re.search(r"[.,]\d{7,}", value):
            raise ValueError("UTC strings support microsecond precision only")
        value = datetime.fromisoformat(value.replace("Z", "+00:00"))
    if not isinstance(value, datetime) or value.utcoffset() != timedelta(0):
        raise ValueError("Timestamp must be timezone-aware UTC")
    return value.astimezone(timezone.utc)


def _iso(value):
    return value.isoformat(timespec="microseconds").replace("+00:00", "Z")


def _spectral_contract(length, fs, nperseg, noverlap, window, detrend):
    if isinstance(fs, (bool, np.bool_)) or not np.isscalar(fs) or not np.isfinite(fs) or fs <= 0:
        raise ValueError("fs must be finite and positive")
    if type(nperseg) is not int or not 2 <= nperseg <= length:
        raise ValueError("nperseg must be an integer in [2, frame samples]; no shortening")
    overlap = nperseg // 2 if noverlap is None else noverlap
    if type(overlap) is not int or not 0 <= overlap < nperseg:
        raise ValueError("noverlap must be an integer in [0, nperseg)")
    if window not in ("hann", "boxcar") or not (detrend is False or detrend == "constant"):
        raise ValueError("Use periodic hann/boxcar and detrend=False or 'constant'")
    return overlap


def frame_spectra(samples, fs, *, nperseg=1024, noverlap=None, window="hann", detrend=False):
    """One-sided density-normalized STFT [segment,freq], Welch=mean(abs(STFT)**2).

    Segment-local detrending only. Periodic windows, no padding, complete
    segments only. This explicit kernel is available for matched predictions.
    """
    raw = np.asarray(samples)
    if np.iscomplexobj(raw) or raw.ndim != 1 or not raw.size or not np.isfinite(raw).all():
        raise ValueError("samples must be a finite nonempty real vector")
    x = np.asarray(raw, dtype=float)
    overlap = _spectral_contract(len(x), fs, nperseg, noverlap, window, detrend)
    nsegments = (len(x) - nperseg) // (nperseg - overlap) + 1
    if len(x) > 1_000_000 or nsegments * (nperseg // 2 + 1) > 2_000_000:
        raise ValueError("Shared spectral work exceeds memory bound")
    starts = np.arange(0, len(x) - nperseg + 1, nperseg - overlap)
    taper = get_window(window, nperseg, fftbins=True)
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise", under="ignore"):
            squared = x * x
            if np.any((x != 0) & (squared == 0)):
                raise ValueError("Signal square is not representable as float64")
            # The raw-square guard and sample bound keep an unscaled FFT finite.
            # A frame-wide normalization can erase a separate weak segment before
            # physical amplitude restoration, bypassing the PSD underflow guard.
            segments = np.array([x[s : s + nperseg] for s in starts])
            if detrend == "constant":
                mean = segments.mean(axis=1, keepdims=True)
                # A reduction can round away from an exactly repeated sample.
                # Use that sample only for exact constants; never threshold a signal.
                constant = np.all(segments == segments[:, :1], axis=1, keepdims=True)
                segments -= np.where(constant, segments[:, :1], mean)
            normalization = np.sqrt(fs) * np.sqrt(np.sum(taper**2))
            stft = np.fft.rfft(segments * taper, axis=1)
            # Keep nonzero-bin evidence before density division can erase it.
            # Reject only an unrepresentable final mean, not a weak segment in
            # a bin whose other segments supply representable power.
            nonzero = np.any(stft != 0, axis=0)
            stft = stft / normalization
            stft[:, 1 : -1 if nperseg % 2 == 0 else None] *= np.sqrt(2)
            # Scale each frequency independently: a strong bin must not erase a weak one.
            amplitude = np.abs(stft)
            peak = amplitude.max(axis=0)
            unit_power = (amplitude / np.where(peak == 0, 1.0, peak)) ** 2
            # Avoid overflowing a sum or peak**2 when the mean power is representable.
            density = (unit_power.mean(axis=0) * peak) * peak
    except (FloatingPointError, OverflowError) as exc:
        raise ValueError("Spectrum or signal square is not representable as float64") from exc
    if not np.isfinite(density).all() or np.any(nonzero & (density == 0)):
        raise ValueError("PSD is not representable as float64")
    return {
        "frequency_hz": np.fft.rfftfreq(nperseg) * fs,
        "psd": density,
        "stft": stft,
        "segment_start_sample": starts,
    }


def _levels(power, reference_power=1.0):
    with np.errstate(divide="ignore", invalid="raise"):
        return 10 * np.log10(power) - 10 * np.log10(reference_power)


def _snapshot(path, limit, expected=None):
    if expected is not None and (
        not isinstance(expected, str) or not re.fullmatch(r"[0-9a-f]{64}", expected)
    ):
        raise ValueError("expected SHA256 must be 64 lowercase hexadecimal characters")
    path = Path(path)
    if path.stat().st_size > limit:
        raise ValueError("WAV exceeds max_chunk_bytes")
    with path.open("rb") as handle:
        payload = handle.read(limit + 1)
    if len(payload) > limit:
        raise ValueError("WAV exceeds max_chunk_bytes")
    digest = hashlib.sha256(payload).hexdigest()
    if expected is not None and digest != expected:
        raise ValueError("SHA256 mismatch; source changed")
    if len(payload) < 12 or payload[:4] != b"RIFF" or payload[8:12] != b"WAVE":
        raise ValueError("Expected RIFF/WAVE")
    if int.from_bytes(payload[4:8], "little") + 8 != len(payload):
        raise ValueError("RIFF size mismatch")
    return payload, digest


def _sample_offset(utc, anchor, fs):
    delta = utc - anchor
    microseconds = (delta.days * 86400 + delta.seconds) * 1000000 + delta.microseconds
    numerator = microseconds * fs
    sample = (numerator + 500000) // 1000000
    residual_us = (numerator - sample * 1000000) / fs
    if abs(residual_us) > 0.500001:
        raise ValueError(
            "Timestamp is not on the anchored sample clock (0.5 us rounding tolerance)"
        )
    return sample, residual_us / 1e6


def analyse_chunks(
    chunks,
    *,
    anchor_utc=None,
    frame_seconds=2.0,
    nperseg=1024,
    noverlap=None,
    window="hann",
    detrend=False,
    bands_hz=None,
    active_intervals=(),
    max_chunk_bytes=20_000_000,
    max_decode_samples=1_000_000,
    max_frame_samples=1_000_000,
    max_windows=10_000,
    max_output_cells=2_000_000,
):
    """Bounded offline WAV decoding; complete anchored windows only have spectra.

    Gaps/partial support yield NaN, never zero fill. Overlap/duplicates fail.
    ``valid`` means complete support; ``qc_ok`` additionally excludes full-scale
    contacts and declared active-transmission intersections. See OBSERVATIONS.md.
    """
    for value, cap in (
        (max_chunk_bytes, 20_000_000),
        (max_decode_samples, 1_000_000),
        (max_frame_samples, 1_000_000),
        (max_windows, 10_000),
        (max_output_cells, 2_000_000),
    ):
        if type(value) is not int or not 1 <= value <= cap:
            raise ValueError(f"Resource bound must be an integer in [1, {cap}]")
    chunks = list(islice(chunks, 1025))
    if not 1 <= len(chunks) <= 1024 or not all(isinstance(c, RecordingChunk) for c in chunks):
        raise ValueError("Supply 1 to 1024 RecordingChunks")
    chunks.sort(key=lambda c: _utc(c.start_utc))
    anchor = _utc(anchor_utc) if anchor_utc is not None else _utc(chunks[0].start_utc)
    if len({c.calibration is not None for c in chunks}) != 1:
        raise ValueError("Cannot mix known and unknown calibration")
    calibrated = chunks[0].calibration is not None
    if any(c.evidence_kind not in ("synthetic", "measurement") for c in chunks):
        raise ValueError("Declare evidence_kind='synthetic' or 'measurement'")
    if len({c.evidence_kind for c in chunks}) != 1:
        raise ValueError("Do not combine synthetic and measurement evidence")
    records = []
    for chunk in chunks:
        if chunk.calibration is not None and not isinstance(chunk.calibration, Calibration):
            raise ValueError("calibration must be an existing Calibration instance or None")
        payload, digest = _snapshot(chunk.path, max_chunk_bytes, chunk.expected_sha256)
        info = sf.info(io.BytesIO(payload))
        if info.frames < 1 or not 1 <= info.channels <= min(64, max_decode_samples):
            raise ValueError("WAV must be nonempty with supported channel count")
        if not 1 <= info.samplerate <= 384000:
            raise ValueError("Unsupported sample rate")
        if records and info.samplerate != records[0]["sample_rate_hz"]:
            raise ValueError("All chunks must have the same sample rate; no resampling")
        fs = info.samplerate
        if type(chunk.channel) is not int or not 0 <= chunk.channel < info.channels:
            raise ValueError("Select a valid zero-based integer channel")
        start, residual = _sample_offset(_utc(chunk.start_utc), anchor, fs)
        records.append(
            {
                "path": str(Path(chunk.path).absolute()),
                "sha256": digest,
                "bytes": len(payload),
                "channel": chunk.channel,
                "evidence_kind": chunk.evidence_kind,
                "source": chunk.source,
                "start_utc": _iso(_utc(chunk.start_utc)),
                "start_sample": start,
                "stop_sample": start + info.frames,
                "timestamp_rounding_residual_s": residual,
                "sample_rate_hz": fs,
                "channels": info.channels,
                "subtype": info.subtype,
                "calibration": asdict(chunk.calibration) if calibrated else None,
                "pressure_scale_pa_per_digital_fraction": float(
                    calibrate_samples([1.0], chunk.calibration)[0]
                )
                if calibrated
                else None,
                "adc_full_scale_convention": "digital +/-1 = +/-full_scale_voltage V (peak, not RMS or peak-to-peak)"
                if calibrated
                else None,
                "calibration_assumption": "flat scalar response over the analysed band"
                if calibrated
                else None,
            }
        )
    gaps = []
    for a, b in zip(records[:-1], records[1:]):
        if b["start_sample"] < a["stop_sample"]:
            raise ValueError("Chunk overlap or duplicate timing is not permitted")
        if b["start_sample"] > a["stop_sample"]:
            gaps.append({"start_sample": a["stop_sample"], "stop_sample": b["start_sample"]})
    if not np.isscalar(frame_seconds) or not np.isfinite(frame_seconds) or frame_seconds <= 0:
        raise ValueError("frame_seconds must be finite and positive")
    requested = frame_seconds * fs
    if requested > max_frame_samples or abs(requested - round(requested)) > 1e-8:
        raise ValueError("Frame must be an integer number of samples within max_frame_samples")
    frame = round(requested)
    overlap = _spectral_contract(frame, fs, nperseg, noverlap, window, detrend)
    end = records[-1]["stop_sample"]
    first = records[0]["start_sample"] // frame * frame
    count_windows = (end - first + frame - 1) // frame
    nsegments = (frame - nperseg) // (nperseg - overlap) + 1
    if (
        count_windows > max_windows
        or count_windows * (nperseg // 2 + 1) * (nsegments + 1) > max_output_cells
    ):
        raise ValueError("Output exceeds max_windows or max_output_cells")
    starts = np.arange(first, end, frame, dtype=np.int64)
    segment_starts = np.arange(nsegments) * (nperseg - overlap)
    frequency = np.fft.rfftfreq(nperseg, 1 / fs)
    bands = {"full": (0, fs / 2)} if bands_hz is None else bands_hz
    if not isinstance(bands, dict) or not 1 <= len(bands) <= 128:
        raise ValueError("bands_hz must contain 1 to 128 named frequency pairs")
    masks = []
    for name, edges in bands.items():
        if not isinstance(name, str) or not name or len(edges) != 2:
            raise ValueError("Each band needs a nonempty name and two edges")
        low, high = edges
        if not np.isfinite(edges).all() or not 0 <= low < high <= fs / 2:
            raise ValueError("Band edges must satisfy 0 <= low < high <= Nyquist")
        mask = (frequency >= low) & (frequency < high if high < fs / 2 else frequency <= high)
        if not mask.any():
            raise ValueError("Band contains no frequency bins")
        masks.append(mask)
    intervals = list(islice(active_intervals, 1025))
    if len(intervals) > 1024:
        raise ValueError("At most 1024 active intervals")
    intervals = [
        ((_utc(a) - anchor).total_seconds(), (_utc(b) - anchor).total_seconds())
        for a, b in intervals
    ]
    if any(b <= a for a, b in intervals):
        raise ValueError("Active interval stop must follow start")
    density = np.full((len(starts), len(frequency)), np.nan)
    stft = np.full((len(starts), nsegments, len(frequency)), complex(np.nan, np.nan))
    coverage, near, active = [], [], []
    for i, start in enumerate(starts):
        x = np.full(frame, np.nan)
        count, contacts = 0, 0
        for chunk, record in zip(chunks, records):
            left, right = (
                max(start, record["start_sample"]),
                min(start + frame, record["stop_sample"]),
            )
            if left >= right:
                continue
            payload, _ = _snapshot(chunk.path, max_chunk_bytes, record["sha256"])
            with sf.SoundFile(io.BytesIO(payload)) as stream:
                stream.seek(int(left - record["start_sample"]))
                block_size = max_decode_samples // stream.channels
                for offset in range(int(left), int(right), block_size):
                    take = min(block_size, int(right) - offset)
                    digital = stream.read(take, dtype="float64", always_2d=True)[:, chunk.channel]
                    if (
                        len(digital) != take
                        or not np.isfinite(digital).all()
                        or np.any(np.abs(digital) > 1)
                    ):
                        raise ValueError("Expected finite normalized digital samples in [-1, 1]")
                    x[offset - start : offset - start + take] = (
                        calibrate_samples(digital, chunk.calibration) if calibrated else digital
                    )
                    contacts += np.count_nonzero(np.abs(digital) >= 0.999)
            count += right - left
        coverage.append(count / frame)
        near.append(contacts)
        active.append(any(a < (start + frame) / fs and b > start / fs for a, b in intervals))
        if count == frame:
            spectrum = frame_spectra(
                x, fs, nperseg=nperseg, noverlap=overlap, window=window, detrend=detrend
            )
            density[i], stft[i] = spectrum["psd"], spectrum["stft"]
    # Detect changes even after a chunk's last decoded frame. Never rewrite input.
    for record in records:
        _snapshot(record["path"], max_chunk_bytes, record["sha256"])
    with np.errstate(over="raise", invalid="raise", under="ignore"):
        try:
            power = np.column_stack(
                [(density[:, mask] * (fs / nperseg)).sum(axis=1) for mask in masks]
            )
        except FloatingPointError as exc:
            raise ValueError("Band power is not representable as float64") from exc
    valid = np.array(coverage) == 1
    return {
        "unit": "Pa" if calibrated else "normalized digital amplitude",
        "calibration_status": "declared scalar calibration" if calibrated else "unknown",
        "sample_rate_hz": fs,
        "window_start_sample": starts,
        "window_stop_sample": starts + frame,
        "window_start_s": starts / fs,
        "window_stop_s": (starts + frame) / fs,
        "window_start_utc": np.array(
            [_iso(anchor + timedelta(seconds=int(s) / fs)) for s in starts]
        ),
        "window_stop_utc": np.array(
            [_iso(anchor + timedelta(seconds=int(s + frame) / fs)) for s in starts]
        ),
        "frequency_hz": frequency,
        "psd": density,
        "stft": stft,
        "valid": valid,
        "coverage_fraction": np.array(coverage),
        "partial_window": (starts < records[0]["start_sample"]) | (starts + frame > end),
        "near_full_scale_samples": np.array(near),
        "active_transmission": np.array(active),
        "qc_ok": valid & (np.array(near) == 0) & ~np.array(active),
        "segment_start_sample": starts[:, None] + segment_starts[None, :],
        "band_names": list(bands),
        "band_power": power,
        "band_level_db": _levels(power, 1e-12 if calibrated else 1.0),
        "band_spl_db_re_1upa": _levels(power, 1e-12) if calibrated else None,
        "metadata": {
            "schema": "oceansignatures.observations.v1",
            "evidence_kind": chunks[0].evidence_kind,
            "chunks": records,
            "anchor_utc": _iso(anchor),
            "gaps": gaps,
            "clock": "integer sample offsets / fs; UTC labels rounded to microseconds; half-open windows",
            "timestamp_rounding_tolerance_s": 0.500001e-6,
            "processing": {
                "frame_samples": frame,
                "nperseg": nperseg,
                "noverlap": overlap,
                "window": window,
                "window_periodic": True,
                "detrend": detrend,
                "scaling": "density",
                "onesided": True,
                "average": "mean",
                "padding": False,
                "unused_tail_samples_per_frame": frame - (int(segment_starts[-1]) + nperseg),
            },
            "psd_units": "Pa^2/Hz" if calibrated else "digital^2/Hz",
            "stft_units": "Pa/sqrt(Hz)" if calibrated else "digital/sqrt(Hz)",
            "band_power_units": "Pa^2" if calibrated else "digital^2",
            "band_level_reference": "1 uPa^2" if calibrated else "1 digital^2",
            "bands_hz": {name: list(edges) for name, edges in bands.items()},
            "band_edge_policy": "low inclusive; high exclusive except Nyquist inclusive; sum PSD bins * df",
            "active_intervals_relative_s": intervals,
            "near_full_scale_threshold": 0.999,
            "clipping_caveat": "potential full-scale contact, not proof of hardware clipping",
            "execution": "bounded offline chunk decoding; no onboard realtime performance claim",
        },
    }


def _json_safe(value):
    if isinstance(value, np.ndarray):
        return _json_safe(value.tolist())
    if isinstance(value, np.generic):
        return _json_safe(value.item())
    if isinstance(value, dict):
        return {key: _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def _safe_output(path):
    if ".." in path.parts:
        raise ValueError("Output path must not contain '..'")
    for item in (path, *path.parents):
        if item.is_symlink():
            raise ValueError("Output path or ancestor is a symlink")


def write_observations(analysis, output_directory):
    """Publish deterministic NPZ/CSV/JSON in a fresh directory, or verify replay.

    Never replace existing differing/incomplete products. Symlinks are rejected.
    The complete family is staged before directory promotion; concurrent external
    writers/renames are outside this offline single-writer contract.
    """
    output = Path(output_directory).absolute()
    _safe_output(output)
    names = ("observations.npz", "windows.csv", "summary.json")
    sources = {Path(record["path"]).resolve() for record in analysis["metadata"]["chunks"]}
    for name in names:
        target = output / name
        _safe_output(target)
        if target.resolve() in sources or any(
            target.exists() and target.samefile(source) for source in sources
        ):
            raise ValueError("Output would replace an immutable source")
    archive = io.BytesIO()
    arrays = {key: value for key, value in analysis.items() if isinstance(value, np.ndarray)}
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as writer:
        for key in sorted(arrays):
            body = io.BytesIO()
            np.lib.format.write_array(body, arrays[key], allow_pickle=False)
            info = zipfile.ZipInfo(key + ".npy", date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o600 << 16
            writer.writestr(info, body.getvalue())
    summary = {key: value for key, value in analysis.items() if key not in ("psd", "stft")}
    summary["array_archive"] = {
        "filename": "observations.npz",
        "sha256": hashlib.sha256(archive.getvalue()).hexdigest(),
        "arrays": {
            key: {"shape": list(value.shape), "dtype": str(value.dtype)}
            for key, value in sorted(arrays.items())
        },
    }
    summary["serialization"] = {
        "nonfinite": "JSON null for missing/unknown values and negative infinity (true zero-power level); distinguish with valid and band_power. NPZ retains NaN/-inf; CSV uses blank for missing and -inf for silence.",
    }
    table = io.StringIO(newline="")
    clock_keys = (
        "window_start_utc",
        "window_stop_utc",
        "window_start_sample",
        "window_stop_sample",
        "window_start_s",
        "window_stop_s",
    )
    qc_keys = (
        "valid",
        "qc_ok",
        "coverage_fraction",
        "partial_window",
        "near_full_scale_samples",
        "active_transmission",
    )
    writer = csv.writer(table, lineterminator="\n")
    writer.writerow(
        (
            *clock_keys,
            *qc_keys,
            "band",
            "low_hz",
            "high_hz",
            "band_power",
            "band_level_db",
            "band_spl_db_re_1upa",
            "power_unit",
        )
    )

    def cell(value):
        if value is None or isinstance(value, (float, np.floating)) and np.isnan(value):
            return ""
        if isinstance(value, (float, np.floating)):
            return format(value, ".17g")
        return value

    for i in range(len(analysis["valid"])):
        for j, name in enumerate(analysis["band_names"]):
            spl = analysis["band_spl_db_re_1upa"]
            row = [analysis[key][i] for key in (*clock_keys, *qc_keys)]
            row += [
                name,
                *analysis["metadata"]["bands_hz"][name],
                analysis["band_power"][i, j],
                analysis["band_level_db"][i, j],
                None if spl is None else spl[i, j],
                analysis["metadata"]["band_power_units"],
            ]
            writer.writerow([cell(value) for value in row])
    bodies = {
        "observations.npz": archive.getvalue(),
        "windows.csv": table.getvalue().encode(),
        "summary.json": (
            json.dumps(_json_safe(summary), sort_keys=True, indent=2, allow_nan=False) + "\n"
        ).encode(),
    }
    if output.exists():
        if not output.is_dir() or any(
            not (output / name).is_file() or (output / name).read_bytes() != body
            for name, body in bodies.items()
        ):
            raise FileExistsError("Refusing differing or incomplete existing observation output")
        return [output / name for name in names]
    for record in analysis["metadata"]["chunks"]:
        _snapshot(record["path"], 20_000_000, record["sha256"])
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".observation-stage-", dir=output.parent) as temporary:
        staged = Path(temporary) / "products"
        staged.mkdir()
        for name, body in bodies.items():
            with (staged / name).open("xb") as handle:
                handle.write(body)
        _safe_output(output)
        if output.exists():
            raise FileExistsError("Output appeared during staging")
        os.rename(staged, output)
    return [output / name for name in names]
