"""Read hash-bound OpenFOAM cylinder evidence; hydrodynamic, never acoustic.

Only numpy/scipy/matplotlib and the Python standard library are required.
Native inputs are read-only. No solver is imported or launched.
"""

from __future__ import annotations

import io
import re
import json
import hashlib
from pathlib import Path, PurePosixPath
from dataclasses import dataclass
import xml.etree.ElementTree as ET
import zipfile
import numpy as np
from scipy.signal import periodogram


def _json_bytes(value):
    return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode()


def _fresh_directory(path):
    path = Path(path)
    if ".." in path.parts or any(p.is_symlink() for p in (path, *path.parents)):
        raise ValueError("unsafe output destination")
    if path.exists():
        raise FileExistsError("output already exists; choose a NEW directory")
    return path


def _zip_bytes(items):
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
        for name, raw in sorted(items):
            info = zipfile.ZipInfo(_safe_name(name), date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            archive.writestr(info, raw, compresslevel=6)
    return buffer.getvalue()


def package_native_case(case, destination):
    """Copy verified bytes into deterministic independent zip parts (<45 MB each).

    Preserve original completion unchanged, including pins of intentionally omitted
    volume fields and unused snapshots. These remain in the original local case.
    """
    dest = _fresh_directory(destination)
    groups, group, total = [], [], 0
    for name, raw in sorted(case.payloads.items()):
        if total + len(raw) > 120_000_000 and group:
            groups.append(group)
            group, total = [], 0
        group.append((name, raw))
        total += len(raw)
    if group:
        groups.append(group)
    compressed = []
    while groups:
        group = groups.pop(0)
        raw = _zip_bytes(group)
        if len(raw) >= 45_000_000:
            if len(group) < 2:
                raise ValueError("single archive member exceeds public download limit")
            cut = len(group) // 2
            groups[0:0] = [group[:cut], group[cut:]]
        else:
            compressed.append((group, raw))
    members, archives = {}, {}
    dest.mkdir(parents=True)
    for index, (group, raw) in enumerate(compressed, 1):
        name = f"original-data-{index:02d}.zip"
        (dest / name).write_bytes(raw)
        archives[name] = sha256(raw)
        for member, _ in group:
            members[member] = name
    manifest = {
        "schema": "oceansignatures-native-subset-v1",
        "completion_sha256": sha256(case.payloads["completion.json"]),
        "archives": archives,
        "members": members,
        "omitted_files": sorted(set(case.completion["output_hashes"]) - set(case.payloads)),
        "omission_policy": "Full volume fields and unused 2..98 s VTP frames excluded from public subset; original local cases unchanged.",
        "contains": "Exact native dictionaries, metadata, mesh, mesh checks, solver logs, run receipts, probes, forces, and 21 original VTP frames at 100..140 s.",
    }
    (dest / "package.json").write_bytes(_json_bytes(manifest))
    return manifest


def _source_reader(root):
    """Return read-once byte accessor and available members; never extract archives."""
    if any(p.is_symlink() for p in (root, *root.parents)):
        raise ValueError("symlinked native input is not accepted")

    def direct(name):
        path = root / _safe_name(name)
        if any(p.is_symlink() for p in (path, *path.parents)):
            raise ValueError("symlinked native input is not accepted")
        return path.read_bytes()

    if not (root / "package.json").exists():
        return direct, None
    package = json.loads(direct("package.json"))
    if package.get("schema") != "oceansignatures-native-subset-v1":
        raise ValueError("unsupported packaged native subset")
    archives = {}
    for name, digest in package["archives"].items():
        archives[name] = zipfile.ZipFile(io.BytesIO(_verify(direct(name), digest, name)))
    members = package["members"]
    actual = {}
    for part, z in archives.items():
        for name in z.namelist():
            _safe_name(name)
            if name in actual:
                raise ValueError("duplicate packaged native member")
            actual[name] = part
    if actual != members:
        raise ValueError("package member coverage mismatch")

    def read(name):
        if name not in members:
            raise ValueError("missing packaged native input: " + name)
        return archives[members[name]].read(name)

    return read, set(members)


def parse_plane(raw: bytes, *, expected_time_s):
    """Read native ASCII Float32 VTP with native triangle connectivity, never regrid."""
    root = ET.fromstring(raw)
    piece = root.find("./PolyData/Piece")
    if root.get("type") != "PolyData" or piece is None:
        raise ValueError("expected native PolyData")
    data = {}
    for a in root.iter("DataArray"):
        name = a.get("Name")
        if name in data or a.get("format") != "ascii":
            raise ValueError("duplicate or non-ASCII native array")
        values = _finite(np.fromstring(a.text or "", sep=" "))
        if name in ("Points", "p", "U", "TimeValue"):
            if a.get("type") != "Float32":
                raise ValueError("unexpected native VTP dtype")
            values = values.astype(np.float32).astype(float)
        elif name in ("connectivity", "offsets"):
            if not np.equal(values, np.floor(values)).all():
                raise ValueError("non-integer native connectivity")
            values = values.astype(np.int64)
        data[name] = values
    points = data["Points"].reshape(-1, 3)
    tris = data["connectivity"].reshape(-1, 3)
    p, u = data["p"], data["U"].reshape(-1, 3)
    if (
        len(points) != int(piece.attrib["NumberOfPoints"])
        or len(tris) != int(piece.attrib["NumberOfPolys"])
        or p.shape != (len(points),)
        or u.shape != points.shape
    ):
        raise ValueError("native plane shape mismatch")
    if (
        not np.array_equal(data["offsets"], np.arange(1, len(tris) + 1) * 3)
        or tris.min() < 0
        or tris.max() >= len(points)
    ):
        raise ValueError("invalid native triangle connectivity")
    if np.max(np.abs(points[:, 2])) > 1e-12:
        raise ValueError("native midplane is off z=0")
    # TimeValue is Float32, distinct from the binary64 native directory clock.
    time_tolerance = float(np.spacing(np.float32(expected_time_s))) / 2 + CLOCK_ATOL_S
    if (
        data["TimeValue"].shape != (1,)
        or abs(data["TimeValue"][0] - expected_time_s) > time_tolerance
    ):
        raise ValueError("native VTP time mismatch")
    a, b, c = points[tris[:, 0]], points[tris[:, 1]], points[tris[:, 2]]
    area2 = np.linalg.norm(np.cross(b - a, c - a), axis=1)
    mask = area2 <= 1e-18
    return {
        "points_m": points,
        "triangles": tris,
        "pressure_native_m2_s2": p,
        "velocity_m_s": u,
        "display_mask": mask,
        "declared_dtype": "Float32",
        "report": {
            "time_s": expected_time_s,
            "points": len(points),
            "triangles": len(tris),
            "zero_area_triangles": int((area2 == 0).sum()),
            "display_masked_triangles": int(mask.sum()),
            "display_mask_rule": "twice triangle area <= 1e-18 m2; native arrays unchanged",
            "vtk_time_tolerance_s": time_tolerance,
            "max_abs_z_m": float(np.max(np.abs(points[:, 2]))),
        },
    }


BASELINE_COMPLETION_SHA256 = "1190826790986ac808016f7a68d1e9af3871ce28375262c97329152c1ca0f5ea"
SENSITIVITY_COMPLETION_SHA256 = "3eec7f859fa357e54311319847aadf9c6e2d0a7bd2dd77eaa7c6cfae04b58680"
P_PATH = "postProcessing/probes/0/p"
U_PATH = "postProcessing/probes/0/U"
FORCE_PATH = "postProcessing/forceCoeffs/0/coefficient.dat"
INITIAL_FILES = {
    "0/p",
    "0/U",
    "constant/transportProperties",
    "constant/turbulenceProperties",
    "system/blockMeshDict",
    "system/controlDict",
    "system/fvSchemes",
    "system/fvSolution",
    "case_metadata.json",
    "generator_executed.py",
    "cylinder.foam",
}
MESH_FILES = {
    "constant/polyMesh/" + s for s in ("points", "faces", "owner", "neighbour", "boundary")
}
ESSENTIAL_FILES = (
    INITIAL_FILES
    | MESH_FILES
    | {
        P_PATH,
        U_PATH,
        FORCE_PATH,
        "input_hashes.json",
        "run_record.json",
        "log.blockMesh",
        "log.checkMesh",
        "log.pimpleFoam",
    }
)


def sha256(raw):
    return hashlib.sha256(raw).hexdigest()


def _safe_name(name):
    if (
        not isinstance(name, str)
        or not name
        or PurePosixPath(name).is_absolute()
        or ".." in name.split("/")
        or "\\" in name
    ):
        raise ValueError("unsafe archive identity path")
    return name


def _pins(pins, required):
    if not isinstance(pins, dict) or not pins or not required.issubset(pins):
        raise ValueError("empty or missing essential identity pins")
    for name, digest in pins.items():
        _safe_name(name)
        if not isinstance(digest, str) or not re.fullmatch("[0-9a-f]{64}", digest):
            raise ValueError("malformed identity pin")


def _verify(raw, digest, name):
    if (
        not isinstance(digest, str)
        or not re.fullmatch("[0-9a-f]{64}", digest)
        or sha256(raw) != digest
    ):
        raise ValueError(f"identity mismatch: {name}")
    return raw


@dataclass
class NativeCase:
    metadata: dict
    completion: dict
    payloads: dict[str, bytes]
    frames: dict[float, str]
    time_s: object
    pressure_native_m2_s2: object
    pressure_Pa: object
    velocity_m_s: object
    force_time_s: object
    drag: object
    lift: object
    locations_m: object
    clock: dict
    report: dict


def _number(text, key):
    match = re.search(r"\b" + re.escape(key) + r"\s+([^;]+);", text)
    if not match:
        raise ValueError(f"missing executed dictionary value {key}")
    return float(match[1])


def read_native_case(source, *, expected_completion_sha256):
    """Verify a complete original case, retaining only replay essentials and 21 planes.

    The completion digest is an EXTERNAL trust anchor, never regenerated by this reader.
    Every original pin is checked. Parsing and export consume the same retained bytes.
    """
    read, available = _source_reader(Path(source))
    completion_raw = _verify(read("completion.json"), expected_completion_sha256, "completion.json")
    c = json.loads(completion_raw)
    if c.get("status") != "completed" or c.get("failure_reason") is not None:
        raise ValueError("incomplete native case")
    pins = c.get("output_hashes")
    _pins(pins, ESSENTIAL_FILES)
    meta_raw = _verify(read("case_metadata.json"), pins["case_metadata.json"], "case_metadata.json")
    meta = json.loads(meta_raw)
    if (
        meta.get("schema") != "oceansignatures-native-cylinder-v1"
        or meta.get("adaptive_timestep") is not False
    ):
        raise ValueError("unsupported native metadata contract")
    end, dt = meta["end_time_s"], meta["delta_t_s"]
    if (
        abs(end - 140) > CLOCK_ATOL_S
        or c["requested_end_time_s"] != end
        or abs(c["actual_end_time_s"] - end) > CLOCK_ATOL_S
        or c["solved_steps"] != round(end / dt)
    ):
        raise ValueError("incomplete/inconsistent native endpoint or solved steps")
    plane_names = sorted(
        (
            p
            for p in pins
            if p.startswith("postProcessing/wakePlane/") and p.endswith("/midplane.vtp")
        ),
        key=lambda p: float(PurePosixPath(p).parent.name),
    )
    plane_times = np.array([float(PurePosixPath(p).parent.name) for p in plane_names])
    validate_clock(plane_times, meta["plane_interval_s"], end)
    frames = {
        float(t): p
        for t, p in zip(plane_times, plane_names)
        if 100 - CLOCK_ATOL_S <= t <= 140 + CLOCK_ATOL_S
    }
    if len(frames) != 21:
        raise ValueError("missing required plane identity pins")
    keep = ESSENTIAL_FILES | set(frames.values())
    if available is not None and available != keep | {"completion.json"}:
        raise ValueError("packaged essential input/plane coverage mismatch")
    checked = set(pins) if available is None else keep
    payloads = {"completion.json": completion_raw, "case_metadata.json": meta_raw}
    for name in sorted(checked):
        digest = pins[name]
        raw = meta_raw if name == "case_metadata.json" else _verify(read(name), digest, name)
        if name in keep:
            payloads[name] = raw
    initial = json.loads(payloads["input_hashes.json"])
    _pins(initial, INITIAL_FILES)
    for name, digest in initial.items():
        if name not in payloads or pins.get(name) != digest:
            raise ValueError("executed input identity disagrees with completion")
    receipts = json.loads(payloads["run_record.json"])
    if [r["stage"] for r in receipts] != ["blockMesh", "checkMesh", "pimpleFoam"]:
        raise ValueError("incomplete native stage receipts")
    for r in receipts:
        if (
            r.get("returncode") != 0
            or r.get("timed_out") is not False
            or r.get("interrupted") is not False
        ):
            raise ValueError("incomplete or failed native stage")
        raw = payloads["log." + r["stage"]]
        _verify(raw, r["log_sha256"], "stage log")
        if not re.search(rb"(?m)^End\s*$", raw):
            raise ValueError("native stage log has no End")
    mesh = payloads["log.checkMesh"].decode()
    cells = re.search(r"\bcells:\s*(\d+)", mesh)
    if (
        "Mesh OK." not in mesh
        or re.search(r"Failed\s+\d+\s+mesh checks", mesh)
        or not cells
        or int(cells[1]) != meta["cells"]
        or c["native_cells"] != meta["cells"]
    ):
        raise ValueError("native mesh not verified OK")
    if (
        "2 geometric (non-empty/wedge) directions (1 1 0)" not in mesh
        or "2 solution (non-empty) directions (1 1 0)" not in mesh
    ):
        raise ValueError("mesh is not the declared 2D model")
    solver = payloads["log.pimpleFoam"].decode()
    solver_times = np.array([float(x) for x in re.findall(r"^Time = (\S+)", solver, re.M)])
    solver_clock = validate_clock(solver_times, dt, end)
    for name, dimensions in [("0/p", KINEMATIC_DIMENSIONS), ("0/U", [0, 1, -1, 0, 0, 0, 0])]:
        match = re.search(rb"dimensions\s*\[([^]]+)\]", payloads[name])
        if not match or [float(x) for x in match[1].split()] != dimensions:
            raise ValueError("native unit dimensions do not match contract")
    if meta["pressure_dimensions"] != KINEMATIC_DIMENSIONS:
        raise ValueError("metadata pressure dimensions disagree")
    control = payloads["system/controlDict"].decode()
    for key, value in [
        ("endTime", end),
        ("deltaT", dt),
        ("rhoInf", meta["density_kg_m3"]),
        ("magUInf", meta["velocity_m_s"]),
        ("lRef", meta["diameter_m"]),
        ("Aref", meta["diameter_m"] * meta["span_m"]),
    ]:
        if not np.isclose(_number(control, key), value, rtol=1e-12, atol=0):
            raise ValueError("executed configuration disagrees with metadata: " + key)
    if not re.search(r"adjustTimeStep\s+no\s*;", control) or not re.search(
        r"simulationType\s+laminar\s*;", payloads["constant/turbulenceProperties"].decode()
    ):
        raise ValueError("unsupported physics or timestep policy")
    for obj, interval in [
        ("probes", meta["probe_interval_s"]),
        ("forceCoeffs", meta["probe_interval_s"]),
        ("wakePlane", meta["plane_interval_s"]),
    ]:
        body = (
            control.split(obj + " {")[-1]
            if obj + " {" in control
            else re.split(r"\b" + obj + r"\s*\{", control)[-1]
        )
        if not np.isclose(_number(body, "writeInterval") * dt, interval, rtol=0, atol=1e-12):
            raise ValueError("executed sampling interval disagrees with metadata")
    time, p, xyz = parse_probes(payloads[P_PATH])
    utime, u, uxyz = parse_probes(payloads[U_PATH], vector=True)
    if (
        not np.array_equal(time, utime)
        or not np.array_equal(xyz, uxyz)
        or not np.array_equal(xyz, meta["probe_locations_m"])
    ):
        raise ValueError("probe clocks or locations disagree")
    clock = validate_clock(time, meta["probe_interval_s"], end)
    force_text = payloads[FORCE_PATH].decode()
    header = re.search(r"^# Time\s+(.+)$", force_text, re.M)
    if not header:
        raise ValueError("missing native force columns")
    columns = ["Time"] + header[1].split()
    force = _finite(np.loadtxt(io.StringIO(force_text), ndmin=2))
    if force.shape[1] != len(columns) or "Cd" not in columns or "Cl" not in columns:
        raise ValueError("native force columns disagree")
    force_clock = validate_clock(force[:, 0], meta["probe_interval_s"], end)
    if not np.allclose(force[:, 0], time, rtol=0, atol=CLOCK_ATOL_S):
        raise ValueError("force/probe clock mismatch")
    report = {
        "completion_sha256": expected_completion_sha256,
        "verified_original_pins": len(checked),
        "original_pins": len(pins),
        "omitted_original_pins": len(set(pins) - checked),
        "retained_pins": len(keep),
        "native_cells": meta["cells"],
        "mesh_status": "Mesh OK.",
        "solver_clock": solver_clock,
        "force_clock": force_clock,
        "max_force_probe_clock_difference_s": float(np.max(np.abs(force[:, 0] - time))),
        "selected_plane_times_s": list(frames),
        "original_plane_count": len(plane_names),
        "input_identity": "original completion and pre-execution pins verified",
        "reconstruction": "not_required: serial native probes and direct sampled planes",
    }
    return NativeCase(
        meta,
        c,
        payloads,
        frames,
        time,
        p,
        pressure_pa(p, meta["density_kg_m3"], meta["pressure_dimensions"]),
        u,
        force[:, 0],
        force[:, columns.index("Cd")],
        force[:, columns.index("Cl")],
        xyz,
        clock,
        report,
    )


def validate_clock(time, interval_s, endpoint_s):
    """Check native samples against a fixed cadence without changing any samples."""
    time = _finite(np.asarray(time))
    if interval_s <= 0 or endpoint_s <= 0:
        raise ValueError("invalid clock interval/endpoint")
    n = round(endpoint_s / interval_s)
    expected = np.arange(1, n + 1) * interval_s
    if time.shape != expected.shape or not np.allclose(time, expected, rtol=0, atol=CLOCK_ATOL_S):
        raise ValueError("native clock has missing/irregular samples or incomplete endpoint")
    return {
        "samples": n,
        "interval_s": interval_s,
        "tolerance_s": CLOCK_ATOL_S,
        "max_accumulation_serialization_error_s": float(np.max(np.abs(time - expected))),
        "first_s": float(time[0]),
        "last_s": float(time[-1]),
        "resampled": False,
    }


def window_mask(time, window):
    """Half-open [start,end), with only <=1e-9 s clock-boundary roundoff admitted."""
    start, end = window
    return (time >= start - CLOCK_ATOL_S) & (time < end - CLOCK_ATOL_S)


def window_statistics(time, pressure_Pa, lift, drag, window, interval_s, diameter_m, velocity_m_s):
    """Demean within window; population RMS and one-sided periodic-Hann PSD.

    Lift frequency uses rising crossings of zero lift (not the pressure peak).
    Linear crossing interpolation estimates event times, NOT signal resampling.
    """
    mask = window_mask(np.asarray(time), window)
    t, p, cl, cd = (np.asarray(x)[mask] for x in (time, pressure_Pa, lift, drag))
    if len(t) < 4:
        raise ValueError("analysis window has insufficient samples")
    for a in (t, p, cl, cd):
        _finite(a)
    fluct = p - p.mean(axis=0)
    freq, psd = periodogram(
        fluct, fs=1 / interval_s, window="hann", detrend=False, scaling="density", axis=0
    )
    _, lift_psd = periodogram(
        cl - cl.mean(), fs=1 / interval_s, window="hann", detrend=False, scaling="density"
    )
    ix = np.flatnonzero((cl[:-1] <= 0) & (cl[1:] > 0))
    crossings = t[ix] - cl[ix] * (t[ix + 1] - t[ix]) / (cl[ix + 1] - cl[ix])
    periods = np.diff(crossings)
    f = float(1 / periods.mean()) if len(periods) else None
    return (
        {
            "window_s": list(window),
            "membership": "[start,end)",
            "samples": len(t),
            "mean_removal": "separate arithmetic mean within this window; no linear detrend",
            "pressure_mean_Pa": p.mean(axis=0).tolist(),
            "pressure_rms_Pa": np.sqrt(np.mean(fluct**2, axis=0)).tolist(),
            "mean_Cd": float(cd.mean()),
            "mean_Cl": float(cl.mean()),
            "lift_rms": float(np.std(cl)),
            "lift_frequency_Hz": f,
            "strouhal": f * diameter_m / velocity_m_s if f is not None else None,
            "rising_lift_crossings_s": crossings.tolist(),
            "complete_lift_periods": len(periods),
            "period_mean_s": float(periods.mean()) if len(periods) else None,
            "period_std_s": float(periods.std()) if len(periods) else None,
            "period_range_s": [float(periods.min()), float(periods.max())]
            if len(periods)
            else None,
            "lift_psd_peak_Hz": float(freq[1 + np.argmax(lift_psd[1:])]),
            "pressure_psd_peak_Hz": freq[1 + np.argmax(psd[1:], axis=0)].tolist(),
            "psd_bin_width_Hz": float(freq[1] - freq[0]),
            "psd_units": "Pa2/Hz",
            "pressure_units": "Pa",
        },
        freq,
        psd,
    )


KINEMATIC_DIMENSIONS = [0, 2, -2, 0, 0, 0, 0]
CLOCK_ATOL_S = 1e-9  # explicit fixed-step accumulation/ASCII serialization bound


def analyse_case(case):
    """Deterministic reduction of retained native histories (no source reopens)."""
    meta = case.metadata

    def window(bounds):
        return window_statistics(
            case.time_s,
            case.pressure_Pa,
            case.lift,
            case.drag,
            bounds,
            meta["probe_interval_s"],
            meta["diameter_m"],
            meta["velocity_m_s"],
        )

    primary, f, psd = window((100, 140))
    sensitivity, f80, psd80 = window((80, 140))
    blocks = [window((i, i + 20))[0] for i in range(0, 140, 20)]
    summary = {
        "primary": primary,
        "sensitivity_80_140": sensitivity,
        "blocks_20_s": blocks,
        "clock": case.clock,
        "configuration": meta,
        "pressure_interpretation": "native local hydrodynamic pressure; NOT acoustic pressure",
        "pressure_conversion": "p_Pa = density_kg_m3 * native_p_m2_s2 (once)",
        "frequency_method": "rising zero crossings of signed native lift; linear event-time interpolation, not sample resampling",
        "psd_method": "one-sided periodogram, periodic Hann window, density scaling, full selected record, no zero padding",
        "startup_policy": "Primary 100–140 s follows observed growth; 80–140 s and all 20 s blocks retained. Stationarity is not assumed.",
        "block_lift_rms_change_80_100_to_100_120_percent": 100
        * (blocks[5]["lift_rms"] / blocks[4]["lift_rms"] - 1),
        "block_lift_rms_change_100_120_to_120_140_percent": 100
        * (blocks[6]["lift_rms"] / blocks[5]["lift_rms"] - 1),
    }
    arrays = {
        "time_s": case.time_s,
        "native_pressure_m2_s2": case.pressure_native_m2_s2,
        "pressure_Pa": case.pressure_Pa,
        "velocity_m_s": case.velocity_m_s,
        "probe_locations_m": case.locations_m,
        "force_time_s": case.force_time_s,
        "Cd": case.drag,
        "Cl": case.lift,
        "frequency_Hz": f,
        "pressure_psd_Pa2_Hz": psd,
        "sensitivity_frequency_Hz": f80,
        "sensitivity_pressure_psd_Pa2_Hz": psd80,
    }
    return summary, arrays


def compare_cases(baseline, second):
    """Require identical mesh/physics and physical sampling clocks; only dt varies."""
    a, b = baseline.metadata, second.metadata
    if {k: v for k, v in a.items() if k != "delta_t_s"} != {
        k: v for k, v in b.items() if k != "delta_t_s"
    } or a["delta_t_s"] == b["delta_t_s"]:
        raise ValueError("cases are not matched except timestep")
    matched_files = (INITIAL_FILES | MESH_FILES) - {"case_metadata.json", "system/controlDict"}
    for name in matched_files:
        if baseline.payloads[name] != second.payloads[name]:
            raise ValueError("cases have unmatched mesh/physics bytes: " + name)

    def normalise(case):
        text = case.payloads["system/controlDict"].decode()
        text = re.sub(r"\bdeltaT\s+[^;]+;", "deltaT VARIABLE;", text)
        return re.sub(
            r"writeControl\s+timeStep;\s*writeInterval\s+(\d+);",
            lambda m: (
                "writeControl runTime; writeInterval "
                + format(int(m[1]) * case.metadata["delta_t_s"], ".12g")
                + ";"
            ),
            text,
        )

    if normalise(baseline) != normalise(second):
        raise ValueError("cases have unmatched executed controls besides dt/output step counts")
    if not np.allclose(baseline.time_s, second.time_s, rtol=0, atol=CLOCK_ATOL_S):
        raise ValueError("cases have unmatched probe clocks")
    x, y = analyse_case(baseline)[0]["primary"], analyse_case(second)[0]["primary"]
    metrics = ["mean_Cd", "lift_rms", "lift_frequency_Hz", "strouhal", "pressure_rms_Pa"]
    delta = {
        k: ((np.asarray(y[k]) - np.asarray(x[k])) / np.asarray(x[k]) * 100).tolist()
        for k in metrics
    }
    return {
        "interpretation": "two-level timestep sensitivity; not convergence order",
        "window_s": [100, 140],
        "baseline_delta_t_s": a["delta_t_s"],
        "second_delta_t_s": b["delta_t_s"],
        "matched_mesh_physics_files": sorted(matched_files),
        "control_equivalence": "Only deltaT and commensurate output-step counts differ; physical output intervals identical.",
        "max_probe_clock_difference_s": float(np.max(np.abs(baseline.time_s - second.time_s))),
        "changes_percent": delta,
        "change_definition": "100*(second-baseline)/baseline",
        "baseline": x,
        "second": y,
    }


def _finite(values):
    if not np.isfinite(values).all():
        raise ValueError("non-finite native values")
    return values


def parse_probes(raw: bytes, *, vector: bool = False):
    """Return original clock, native scalar/vector values and ordered locations."""
    text = raw.decode("utf-8")
    headers = re.findall(r"^# Probe (\d+) \(([^)]+)\)", text, re.M)
    if not headers or [int(i) for i, _ in headers] != list(range(len(headers))):
        raise ValueError("missing or unordered probe locations")
    xyz = _finite(np.array([[float(x) for x in s.split()] for _, s in headers]))
    data = _finite(np.loadtxt(io.StringIO(text.replace("(", " ").replace(")", " ")), ndmin=2))
    components = 3 if vector else 1
    if xyz.shape != (len(headers), 3) or data.shape[1] != 1 + len(headers) * components:
        raise ValueError("probe shape does not match locations")
    values = data[:, 1:]
    if vector:
        values = values.reshape(-1, len(headers), 3)
    return data[:, 0], values, xyz


def pressure_pa(native_p, density_kg_m3: float, dimensions):
    """Convert kinematic m²/s² to Pa exactly once; reject dimensional pressure."""
    if list(dimensions) != KINEMATIC_DIMENSIONS:
        raise ValueError("expected kinematic pressure dimensions [0 2 -2 0 0 0 0]")
    if not np.isfinite(density_kg_m3) or density_kg_m3 <= 0:
        raise ValueError("positive finite density required")
    return _finite(np.asarray(native_p, dtype=float) * density_kg_m3)


def _write_numerical(case, summary, arrays, label, out):
    data = out / "data"
    data.mkdir(exist_ok=True)
    items = []
    for name, array in sorted(arrays.items()):
        buf = io.BytesIO()
        np.save(buf, array, allow_pickle=False)
        items.append((name + ".npy", buf.getvalue()))
    (data / f"{label}.npz").write_bytes(_zip_bytes(items))
    n = len(case.locations_m)
    table = np.column_stack(
        [case.time_s, case.pressure_Pa, case.force_time_s, case.drag, case.lift]
    )
    np.savetxt(
        data / f"{label}_history.csv",
        table,
        delimiter=",",
        fmt="%.17g",
        comments="",
        header=",".join(
            ["probe_time_s"] + [f"p{i + 1}_Pa" for i in range(n)] + ["force_time_s", "Cd", "Cl"]
        ),
    )
    np.savetxt(
        data / f"{label}_pressure_psd.csv",
        np.column_stack([arrays["frequency_Hz"], arrays["pressure_psd_Pa2_Hz"]]),
        delimiter=",",
        fmt="%.17g",
        comments="",
        header=",".join(["frequency_Hz"] + [f"p{i + 1}_Pa2_Hz" for i in range(n)]),
    )
    rows = [summary["primary"], summary["sensitivity_80_140"]] + summary["blocks_20_s"]
    table = [
        [*s["window_s"], s["samples"], s["mean_Cd"], s["lift_rms"], *s["pressure_rms_Pa"]]
        for s in rows
    ]
    np.savetxt(
        data / f"{label}_windows.csv",
        table,
        delimiter=",",
        fmt="%.17g",
        comments="",
        header=",".join(
            ["start_s", "end_exclusive_s", "samples", "mean_Cd", "lift_rms"]
            + [f"p{i + 1}_rms_Pa" for i in range(n)]
        ),
    )


def _plotting():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )
    return plt


def _wake_axes(ax, plane, case, bound, *, wide=False):
    from matplotlib.tri import Triangulation
    from matplotlib.patches import Circle

    xyz, tri = plane["points_m"], plane["triangles"]
    native_tri = Triangulation(xyz[:, 0], xyz[:, 1], tri, mask=plane["display_mask"])
    artist = ax.tripcolor(
        native_tri,
        plane["velocity_m_s"][:, 1],
        shading="gouraud",
        cmap="RdBu_r",
        vmin=-bound,
        vmax=bound,
        rasterized=True,
    )
    ax.add_patch(
        Circle(
            (0, 0),
            case.metadata["diameter_m"] / 2,
            facecolor="#313844",
            edgecolor="black",
            zorder=3,
        )
    )
    for i, (x, y, _) in enumerate(case.locations_m, 1):
        ax.plot(x, y, "o", ms=4, color="black", mec="white", mew=0.6, zorder=4)
        if not wide:
            ax.annotate(
                str(i), (x, y), xytext=(4, 5), textcoords="offset points", fontsize=9, zorder=5
            )
    ax.set_aspect("equal", adjustable="box")
    if wide:
        r = case.metadata["outer_radius_m"]
        ax.set_xlim(-r, r)
        ax.set_ylim(-r, r)
    else:
        ax.set_xlim(-0.03, 0.12)
        ax.set_ylim(-0.04, 0.04)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    return artist


def _render_figures(cases, summaries, arrays, planes, out):
    plt = _plotting()
    figures = out / "figures"
    figures.mkdir()
    baseline = cases["baseline"]
    values = [np.max(np.abs(p["velocity_m_s"][:, 1])) for p in planes]
    bound = float(max(values))
    if bound == 0:
        bound = 1e-12
    fig, ax = plt.subplots(
        1, 2, figsize=(12, 5), gridspec_kw={"width_ratios": [1, 1.7]}, layout="constrained"
    )
    artist = _wake_axes(ax[0], planes[-1], baseline, bound, wide=True)
    _wake_axes(ax[1], planes[-1], baseline, bound)
    ax[0].set_title("Native circular domain")
    ax[1].set_title("Wake + probe locations, t = 140 s")
    fig.colorbar(artist, ax=ax, label="Native transverse velocity Uy (m/s)", shrink=0.7)
    fig.suptitle("2D laminar cylinder wake — Re = 100")
    fig.savefig(figures / "native_wake.png", dpi=150)
    plt.close(fig)
    fig, ax = plt.subplots(3, 1, figsize=(10, 10), layout="constrained")
    primary = summaries["baseline"]["primary"]
    mask = window_mask(baseline.time_s, (100, 140))
    means = np.array(primary["pressure_mean_Pa"])
    n = len(baseline.locations_m)
    colors = plt.get_cmap("tab10")(np.arange(n))
    for i in range(n):
        ax[0].plot(
            baseline.time_s[mask],
            baseline.pressure_Pa[mask, i] - means[i],
            color=colors[i],
            lw=1,
            label=f"P{i + 1}",
        )
        f, psd = arrays["baseline"]["frequency_Hz"], arrays["baseline"]["pressure_psd_Pa2_Hz"][:, i]
        ax[2].semilogy(f[1:], psd[1:], color=colors[i], lw=1, label=f"P{i + 1}")
    ax[0].set(
        xlabel="Native time (s)",
        ylabel="p − window mean (Pa)",
        title="Local hydrodynamic pressure, 100–140 s",
    )
    ax[0].legend(loc="upper left", bbox_to_anchor=(1.01, 1), borderaxespad=0)
    x = np.arange(n)
    ax[1].bar(x - 0.18, primary["pressure_rms_Pa"], width=0.36, label="100–140 s", color="#166e94")
    ax[1].bar(
        x + 0.18,
        summaries["baseline"]["sensitivity_80_140"]["pressure_rms_Pa"],
        width=0.36,
        label="80–140 s",
        color="#e59b41",
    )
    ax[1].set(
        xticks=x,
        xticklabels=[f"P{i + 1}" for i in range(n)],
        ylabel="Pressure RMS (Pa)",
        yscale="log",
    )
    ax[1].legend()
    ax[2].set(
        xlim=(0, 1),
        xlabel="Frequency (Hz)",
        ylabel="One-sided PSD (Pa²/Hz)",
        title=f"Hann periodogram; bin width {primary['psd_bin_width_Hz']:.3f} Hz (no zero padding)",
    )
    if primary["lift_frequency_Hz"] is not None:
        ax[2].axvline(
            primary["lift_frequency_Hz"], color="black", ls="--", label="Lift-derived shedding f"
        )
    ax[2].legend(ncol=3)
    for a in ax:
        a.grid(alpha=0.2)
    fig.savefig(figures / "hydrodynamic_pressure.png", dpi=150)
    plt.close(fig)
    fig, ax = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
    for label, case in cases.items():
        title = f"Δt = {case.metadata['delta_t_s']:g} s"
        s = summaries[label]
        ax[0, 0].plot(case.force_time_s, case.drag, lw=0.8, label=title)
        ax[0, 1].plot(case.force_time_s, case.lift, lw=0.8, label=title)
        ax[1, 0].plot(
            [np.mean(b["window_s"]) for b in s["blocks_20_s"]],
            [b["lift_rms"] for b in s["blocks_20_s"]],
            "o-",
            label=title,
        )
        late = window_mask(case.force_time_s, (100, 140))
        ax[1, 1].plot(case.force_time_s[late], case.drag[late], lw=1, label=title)
    ax[0, 0].set(
        title="Drag startup (full range, including impulse)",
        ylabel="Cd",
        yscale="symlog",
        ylim=(-100, 10),
    )
    ax[0, 1].set(title="Lift startup: amplitude continues growing", ylabel="Cl")
    ax[1, 0].set(
        title="20 s block RMS — adverse trend retained",
        ylabel="Mean-removed lift RMS",
        xlabel="Block centre (s)",
    )
    ax[1, 1].set(title="Late drag: matched 100–140 s", ylabel="Cd", xlabel="Native time (s)")
    for a in ax.flat:
        a.grid(alpha=0.2)
        a.legend()
    for a in ax[0]:
        a.axvspan(100, 140, color="grey", alpha=0.12)
        a.set_xlabel("Native time (s)")
    fig.savefig(figures / "forces_startup_sensitivity.png", dpi=150)
    plt.close(fig)
    return bound


def render_movie(case, planes, bound, out):
    """21 native snapshots, 2 fps, no intermediate fields or temporal interpolation."""
    import shutil
    import subprocess

    if not shutil.which("ffmpeg"):
        return {"status": "not_created", "reason": "optional ffmpeg executable absent"}
    out = Path(out)
    frames = out / "movie_frames"
    frames.mkdir()
    plt = _plotting()
    times = [p["report"]["time_s"] for p in planes]
    if len(planes) != 21 or not np.allclose(
        times, np.arange(100, 141, 2), rtol=0, atol=CLOCK_ATOL_S
    ):
        raise ValueError("movie requires the 21 native 100..140 s snapshots")
    for index, plane in enumerate(planes):
        fig, ax = plt.subplots(figsize=(12, 5), layout="constrained")
        artist = _wake_axes(ax, plane, case, bound)
        fig.colorbar(artist, ax=ax, label="Native Uy (m/s)", shrink=0.8)
        prefix = "TEST-ONLY fixture" if case.metadata.get("test_only") else "Native Re100 wake"
        ax.set_title(
            f"{prefix} | t = {times[index]:.0f} s | 2 s native snapshots, 2 fps playback (4×)"
        )
        fig.savefig(frames / f"frame_{index:03d}.png", dpi=100)
        plt.close(fig)
    command = [
        "ffmpeg",
        "-v",
        "error",
        "-nostdin",
        "-n",
        "-framerate",
        "2",
        "-i",
        str(frames / "frame_%03d.png"),
        "-c:v",
        "libx264",
        "-threads",
        "1",
        "-crf",
        "20",
        "-pix_fmt",
        "yuv420p",
        "-fflags",
        "+bitexact",
        "-flags:v",
        "+bitexact",
        "-map_metadata",
        "-1",
        "-movflags",
        "+faststart",
        str(out / "native_wake.mp4"),
    ]
    subprocess.run(command, check=True, capture_output=True, timeout=120)
    record = {
        "status": "created",
        "path": "native_wake.mp4",
        "source_snapshot_count": len(planes),
        "source_times_s": times,
        "playback_fps": 2,
        "playback_duration_s": len(planes) / 2,
        "physical_time_span_s": times[-1] - times[0],
        "interpolated_fields": False,
        "common_colour_bounds_m_s": [-bound, bound],
        "playback_note": "Each native 2 s snapshot shown for 0.5 s (4x nominal time speed); last snapshot held 0.5 s. No in-between fields.",
    }
    (out / "movie.json").write_bytes(_json_bytes(record))
    return record


def build_export(
    baseline,
    destination,
    *,
    sensitivity=None,
    baseline_sha256=BASELINE_COMPLETION_SHA256,
    sensitivity_sha256=SENSITIVITY_COMPLETION_SHA256,
    movie=False,
):
    """Validate before creating output; publish a new immutable directory atomically."""
    import shutil
    import tempfile

    out = _fresh_directory(destination)
    cases = {"baseline": read_native_case(baseline, expected_completion_sha256=baseline_sha256)}
    if sensitivity is not None:
        cases["sensitivity"] = read_native_case(
            sensitivity, expected_completion_sha256=sensitivity_sha256
        )
    comparison = (
        compare_cases(cases["baseline"], cases["sensitivity"]) if sensitivity is not None else None
    )
    summaries, arrays = {}, {}
    for label, case in cases.items():
        summaries[label], arrays[label] = analyse_case(case)
    # Validate all selected native frames before creating any result directory.
    all_planes = {
        label: [
            parse_plane(case.payloads[name], expected_time_s=t) for t, name in case.frames.items()
        ]
        for label, case in cases.items()
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix="." + out.name + "-", dir=out.parent))
    try:
        packages = {}
        for label, case in cases.items():
            _write_numerical(case, summaries[label], arrays[label], label, staging)
            packages[label] = package_native_case(case, staging / "native" / label)
        summary = {
            "schema": "oceansignatures-native-analysis-v1",
            "cases": summaries,
            "timestep_sensitivity": comparison,
        }
        (staging / "analysis_summary.json").write_bytes(_json_bytes(summary))
        bound = _render_figures(cases, summaries, arrays, all_planes["baseline"], staging)
        (staging / "verification.json").write_bytes(
            _json_bytes(
                {
                    label: {**case.report, "planes": [p["report"] for p in all_planes[label]]}
                    for label, case in cases.items()
                }
            )
        )
        # Bundle the actual analysis implementation, not private cache helpers.
        replay = staging / "replay"
        replay.mkdir()
        module_raw = Path(__file__).read_bytes()
        cli_raw = b'#!/usr/bin/env python3\nfrom cfd_analysis import main\nif __name__ == "__main__":\n    raise SystemExit(main())\n'
        (replay / "cfd_analysis.py").write_bytes(module_raw)
        (replay / "build_cfd_results.py").write_bytes(cli_raw)
        command = "python replay/build_cfd_results.py --baseline native/baseline"
        if sensitivity is not None:
            command += " --sensitivity native/sensitivity"
        command += ' --output "$(dirname "$PWD")/replay-new"'
        readme = (
            "Native hydrodynamic CFD analysis replay\n\n" + command + "\n\n"
            "Run from this export folder. Requires Python, numpy, scipy, matplotlib; ffmpeg is optional.\n"
            "No OpenFOAM, original cache, sibling repository or native solver is required.\n"
            "Keep native/*/package.json and every original-data-*.zip together. Do not extract or edit them.\n"
            "The native completion receipts retain every original expected pin. Each included member\n"
            "is verified against those pins. Full volume fields and unused VTP frames stay local.\n"
            "Replay means the same analysis of the same saved data, not another simulation.\n"
            "Only the exact historical native logs/receipts retain machine paths; no ordinary output metadata needs them.\n"
        )
        (staging / "REPLAY.txt").write_text(readme)
        (staging / "analysis-code.zip").write_bytes(
            _zip_bytes(
                [
                    ("cfd_analysis.py", module_raw),
                    ("build_cfd_results.py", cli_raw),
                    ("REPLAY.txt", readme.encode()),
                ]
            )
        )
        figures = [
            {
                "path": "figures/native_wake.png",
                "title": "Native cylinder wake and probes",
                "caption": "Actual t=140 s Float32 sampled plane; native triangle connectivity, equal physical axis scales. Only degenerate/tiny triangles masked for display (twice area ≤1e-18 m²). Common Uy colour bounds across 100–140 s; linear display shading inside native cells, no regridding.",
            },
            {
                "path": "figures/hydrodynamic_pressure.png",
                "title": "Hydrodynamic pressure histories, RMS and PSD",
                "caption": "p_native × density once gives Pa. Each statistic removes its own window mean. Primary [100,140) s, sensitivity [80,140) s; all native samples retained. PSD in Pa²/Hz, not an acoustic receiver spectrum. Fundamental frequency is inferred from signed lift, not a pressure harmonic.",
            },
            {
                "path": "figures/forces_startup_sensitivity.png",
                "title": "Startup and timestep sensitivity",
                "caption": "Full native drag/lift startup and 20 s block RMS retain amplitude growth. Late-window comparison uses identical physical sampling and mesh/physics; two timesteps do not establish convergence order or mesh independence.",
            },
        ]
        downloads = [
            {"path": "analysis_summary.json", "label": "Full numerical summary"},
            {"path": "verification.json", "label": "Identity, clocks and plane checks"},
            {"path": "analysis-code.zip", "label": "Portable analysis Python source"},
            {"path": "REPLAY.txt", "label": "Replay instructions"},
        ]
        for p in sorted((staging / "data").iterdir()):
            downloads.append({"path": p.relative_to(staging).as_posix(), "label": p.name})
        for label, package in packages.items():
            downloads.append(
                {
                    "path": f"native/{label}/package.json",
                    "label": f"{label}: archive index and omissions",
                }
            )
            for name in package["archives"]:
                downloads.append(
                    {
                        "path": f"native/{label}/{name}",
                        "label": f"{label}: exact original native data ({name})",
                    }
                )
        limitations = [
            "2D spanwise-uniform laminar Re100 cylinder; not LES, turbulent broadband noise or a finite housing.",
            "Incompressible local hydrodynamic pressure is not propagating acoustic pressure. Never add these PSDs to acoustic predictions.",
            "100–140 s primary follows observed startup; 80–140 s and 20 s blocks retain residual growth. A short late record does not prove stationarity.",
            "No mesh or domain-size independence demonstrated; two timestep values give sensitivity, not convergence order.",
            "Native planes are Float32 point-interpolated fields. Degenerate cutting-plane triangles are retained in downloads and masked only in display.",
            "Full volume fields and earlier VTP frames are intentionally omitted from the public subset; unchanged originals remain local.",
        ]
        kind = (
            "test-only analytic fixture; not native CFD"
            if cases["baseline"].metadata.get("test_only")
            else "native hydrodynamic CFD"
        )
        overview = {
            "title": "Native laminar cylinder: pressure and vortex shedding",
            "kind": kind,
            "configuration": cases["baseline"].metadata,
            "metrics": {
                "baseline": summaries["baseline"]["primary"],
                "timestep_sensitivity": comparison,
            },
            "limitations": limitations,
            "figures": figures,
            "downloads": downloads,
            "wake_colour_bounds_m_s": [-bound, bound],
        }
        if movie:
            record = render_movie(cases["baseline"], all_planes["baseline"], bound, staging)
            overview["movie"] = record
            if record["status"] == "created":
                overview["moviepath"] = record["path"]
                downloads.extend(
                    [
                        {
                            "path": record["path"],
                            "label": "Native wake movie: 21 authentic snapshots",
                        },
                        {"path": "movie.json", "label": "Movie source times and playback contract"},
                    ]
                )
        (staging / "overview.json").write_bytes(_json_bytes(overview))
        manifest = {
            p.relative_to(staging).as_posix(): sha256(p.read_bytes())
            for p in sorted(staging.rglob("*"))
            if p.is_file()
        }
        if any(p.stat().st_size >= 50_000_000 for p in staging.rglob("*") if p.is_file()):
            raise ValueError("public file exceeds 50 MB")
        (staging / "manifest.json").write_bytes(
            _json_bytes({"schema": "sha256-files-v1", "files": manifest})
        )
        if out.exists():
            raise FileExistsError("destination appeared during build; refusing overwrite")
        staging.rename(out)
        return overview
    except BaseException:
        shutil.rmtree(staging)
        raise


def main(argv=None):
    import argparse

    parser = argparse.ArgumentParser(
        description="Verify, analyse and package native hydrodynamic CFD; never run a solver."
    )
    parser.add_argument("--baseline", required=True)
    parser.add_argument("--sensitivity")
    parser.add_argument("--baseline-sha256", default=BASELINE_COMPLETION_SHA256)
    parser.add_argument("--sensitivity-sha256", default=SENSITIVITY_COMPLETION_SHA256)
    parser.add_argument("--output", required=True)
    parser.add_argument("--no-movie", action="store_true")
    args = parser.parse_args(argv)
    try:
        overview = build_export(
            args.baseline,
            args.output,
            sensitivity=args.sensitivity,
            baseline_sha256=args.baseline_sha256,
            sensitivity_sha256=args.sensitivity_sha256,
            movie=not args.no_movie,
        )
    except (ValueError, OSError, KeyError) as error:
        parser.exit(2, f"CFD export rejected: {error}\n")
    print(
        json.dumps(
            {
                "kind": overview["kind"],
                "figures": len(overview["figures"]),
                "downloads": len(overview["downloads"]),
                "baseline_metrics": overview["metrics"]["baseline"],
            },
            indent=2,
        )
    )
    return 0
