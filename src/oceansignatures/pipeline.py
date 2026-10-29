"""Calibration-aware offline ingestion; no proprietary instrument formats."""

from dataclasses import dataclass
import argparse
import hashlib
import io
import json
from pathlib import Path

import numpy as np
import soundfile as sf
from scipy.signal import welch


@dataclass(frozen=True)
class Calibration:
    """Flat response assumption; sensitivity is dB re 1 V per micropascal."""

    sensitivity_db_v_per_upa: float
    gain_db: float
    full_scale_voltage: float

    def __post_init__(self):
        values = (self.sensitivity_db_v_per_upa, self.gain_db, self.full_scale_voltage)
        if not np.all(np.isfinite(values)) or self.full_scale_voltage <= 0:
            raise ValueError("Calibration must be finite with positive full-scale voltage")
        if abs(self.sensitivity_db_v_per_upa) > 600 or abs(self.gain_db) > 600:
            raise ValueError("Calibration exponent exceeds supported numeric range")


def _pressure_scale(calibration: Calibration) -> float:
    # Combine the dB factors before multiplication: intermediate voltages can
    # overflow even when the final Pa/digital-fraction scale is representable.
    exponent = -(float(calibration.gain_db) + float(calibration.sensitivity_db_v_per_upa)) / 20 - 6
    scale = float(calibration.full_scale_voltage) * 10.0**exponent
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError(
            "Combined calibration pressure scale is not representable as positive float64"
        )
    return scale


def calibrate_samples(samples, calibration: Calibration) -> np.ndarray:
    """Convert normalized ADC fractions to Pa using the declared scalar chain."""
    digital = np.asarray(samples, dtype=float)
    scale = _pressure_scale(calibration)
    error = "Calibrated pressure is not representable as float64"
    try:
        with np.errstate(over="raise", invalid="raise", under="ignore"):
            pressure = digital * scale
    except FloatingPointError as exc:
        raise ValueError(error) from exc
    if not np.all(np.isfinite(pressure)) or np.any((digital != 0) & (pressure == 0)):
        raise ValueError(error)
    return pressure


def analyse_recording(path, calibration=None, channel=0, max_bytes=20_000_000):
    """Analyse a normalized digital audio channel; calibration is never inferred.

    The optional scalar calibration assumes a flat sensor response over the
    retained band. It is a declared conversion, not a calibration certificate.
    """
    path = Path(path)
    if path.stat().st_size > max_bytes:
        raise ValueError("Recording exceeds the configured file-size bound")
    info = sf.info(path)
    if info.frames < 2 or info.frames * info.channels > 10_000_000:
        raise ValueError("Recording must contain 2 to 10000000 decoded samples")
    if (
        isinstance(channel, bool)
        or not isinstance(channel, int)
        or not 0 <= channel < info.channels
    ):
        raise ValueError("Select a valid zero-based integer channel")
    samples, fs = sf.read(path, dtype="float64", always_2d=True)
    signal = samples[:, channel]
    if not np.all(np.isfinite(signal)) or np.any(np.abs(signal) > 1.0):
        raise ValueError("Only finite normalized digital samples in [-1, 1] are supported")
    near_full_scale = int(np.count_nonzero(np.abs(signal) >= 1 - 1 / 32768))
    if calibration is not None:
        signal = calibrate_samples(signal, calibration)
    square_error = "Signal pressure-square/amplitude-square is not representable as float64"
    try:
        with np.errstate(over="raise", invalid="raise", under="ignore"):
            squared = signal**2
    except FloatingPointError as exc:
        raise ValueError(square_error) from exc
    if not np.all(np.isfinite(squared)) or np.any((signal != 0) & (squared == 0)):
        raise ValueError(square_error)
    # Bound the FFT and averaging intermediates, then restore the original
    # squared units. Finite individual squares need not have a finite sum.
    normalizer = float(np.max(np.abs(signal))) or 1.0
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise", under="ignore"):
            normalized = signal / normalizer
            frequency, normalized_density = welch(
                normalized,
                fs=fs,
                window="hann",
                nperseg=min(2048, signal.size),
                detrend=False,
                scaling="density",
                return_onesided=True,
            )
            density = normalized_density * normalizer**2
            rms_value = float(normalizer * np.sqrt(np.mean(normalized**2)))
    except FloatingPointError as exc:
        raise ValueError("PSD is not representable as float64") from exc
    if not np.all(np.isfinite(density)) or np.any((normalized_density != 0) & (density == 0)):
        raise ValueError("PSD is not representable as float64")
    return {
        "sample_rate_hz": fs,
        "sample_count": signal.size,
        "channel": channel,
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "calibration_status": "declared scalar calibration" if calibration else "unknown",
        "calibration": (
            {
                "sensitivity_db_v_per_upa": float(calibration.sensitivity_db_v_per_upa),
                "gain_db": float(calibration.gain_db),
                "full_scale_voltage": float(calibration.full_scale_voltage),
            }
            if calibration is not None
            else None
        ),
        "pressure_scale_pa_per_digital_fraction": (
            _pressure_scale(calibration) if calibration is not None else None
        ),
        "calibration_assumption": (
            "flat scalar response over the analysed band" if calibration is not None else None
        ),
        "adc_full_scale_convention": (
            "digital fractions +/-1 correspond to +/-full_scale_voltage V at the ADC input "
            "(peak, not RMS or peak-to-peak)"
            if calibration is not None
            else None
        ),
        "unit": "Pa" if calibration else "normalized digital amplitude",
        "rms": rms_value,
        "spl_db_re_1upa": (
            float(20 * np.log10(rms_value / 1e-6))
            if calibration is not None and rms_value > 0
            else None
        ),
        "time_s": np.arange(signal.size) / fs,
        "signal": signal,
        "frequency_hz": frequency,
        "psd": density,
        "psd_method": "Welch; periodic Hann; 50% overlap; no detrending; one-sided density",
        "qc": {
            "near_full_scale_samples": near_full_scale,
            "clipping_interpretation": "potential full-scale contact, not proof of hardware clipping",
            "timing_continuity": "not established from a single audio file",
            "finite_samples": True,
        },
    }


def write_analysis(analysis, output_directory):
    """Write reproducible products; refuse to replace differing prior outputs."""
    output = Path(output_directory)
    array_keys = ("time_s", "signal", "frequency_hz", "psd")
    summary = {key: value for key, value in analysis.items() if key not in array_keys}
    stream = io.StringIO()
    np.savetxt(
        stream,
        np.column_stack((analysis["frequency_hz"], analysis["psd"])),
        delimiter=",",
        header="frequency_hz,psd_"
        + ("Pa2_per_Hz" if analysis["unit"] == "Pa" else "digital2_per_Hz"),
        comments="",
    )
    archive = io.BytesIO()
    np.savez_compressed(archive, **{key: analysis[key] for key in array_keys})
    bodies = {
        "summary.json": (
            json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n"
        ).encode(),
        "spectrum.csv": stream.getvalue().encode(),
        "signal.npz": archive.getvalue(),
    }
    for name, body in bodies.items():
        path = output / name
        if path.exists() and path.read_bytes() != body:
            raise FileExistsError(f"Refusing to overwrite differing output: {name}")
    output.mkdir(parents=True, exist_ok=True)
    for name, body in bodies.items():
        path = output / name
        if not path.exists():
            with path.open("xb") as handle:
                handle.write(body)
    return [output / name for name in bodies]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("recording", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--channel", type=int, default=0)
    parser.add_argument("--calibration-json", type=Path)
    args = parser.parse_args(argv)
    calibration = None
    if args.calibration_json is not None:
        calibration = Calibration(**json.loads(args.calibration_json.read_text()))
    analysis = analyse_recording(args.recording, calibration, args.channel)
    paths = write_analysis(analysis, args.output)
    print(
        json.dumps(
            {
                "calibration_status": analysis["calibration_status"],
                "outputs": [str(path) for path in paths],
            }
        )
    )


if __name__ == "__main__":
    main()
