"""Assemble authenticated computed bundles into the existing static field guide.

Single-writer offline assembly. No calculations, solver calls, invented results,
network access, archive extraction or changes to the original experiment bundle.
"""

import hashlib
import json
import math
from pathlib import Path, PurePosixPath
import re
import tempfile


SCHEMAS = {
    "acoustic": "oceansignatures.research.manifest.v1",
    "cfd": "sha256-files-v1",
}
KINDS = {
    "acoustic": "predictive approximations + synthetic verification",
    "cfd": "native hydrodynamic CFD",
}


def _sha(raw):
    return hashlib.sha256(raw).hexdigest()


def _json(raw):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError("Duplicate JSON identity")
            result[key] = value
        return result

    def bad_constant(value):
        raise ValueError("Nonfinite JSON constant")

    return json.loads(raw, object_pairs_hook=unique, parse_constant=bad_constant)


def _encode(value):
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()


def _root(value):
    path = Path(value).absolute()
    if ".." in path.parts or any(p.is_symlink() for p in (path, *path.parents)):
        raise ValueError("Unsafe source or destination path")
    return path


def _name(value):
    if (
        not isinstance(value, str)
        or not re.fullmatch(r"[A-Za-z0-9_./-]+", value)
        or PurePosixPath(value).is_absolute()
        or any(part in {"", ".", ".."} for part in value.split("/"))
    ):
        raise ValueError("Unsafe bundle member path")
    return value


def _digest(value):
    if not isinstance(value, str) or re.fullmatch(r"[0-9a-f]{64}", value) is None:
        raise ValueError("Expected SHA256 must be an exact lowercase hexadecimal digest")
    return value


def _inventory(root):
    if not root.is_dir():
        raise ValueError("Bundle must be a directory")
    files = {}
    for path in sorted(root.rglob("*")):
        if path.is_symlink() or not (path.is_file() or path.is_dir()):
            raise ValueError("Bundle contains a symlink or nonregular member")
        if path.is_file():
            files[_name(path.relative_to(root).as_posix())] = path
    return files


def _bundle(root, expected, family):
    paths = _inventory(root)
    if "manifest.json" not in paths:
        raise ValueError("Missing bundle manifest")
    raw = paths["manifest.json"].read_bytes()
    if _sha(raw) != _digest(expected):
        raise ValueError("Bundle manifest identity mismatch")
    manifest = _json(raw)
    pins = manifest.get("files")
    if manifest.get("schema") != SCHEMAS[family] or not isinstance(pins, dict) or not pins:
        raise ValueError("Unsupported or empty bundle manifest")
    if len(pins) > 1000 or set(paths) != set(pins) | {"manifest.json"}:
        raise ValueError("Bundle membership mismatch")
    total = 0
    payload = {"manifest.json": raw}
    for name, digest in pins.items():
        _name(name)
        _digest(digest)
        size = paths[name].stat().st_size
        total += size
        if size >= 50_000_000 or total > 512_000_000:
            raise ValueError("Bundle exceeds bounded publication size")
        content = paths[name].read_bytes()
        if _sha(content) != digest:
            raise ValueError("Bundle payload identity mismatch")
        payload[name] = content
    overview = _json(payload["overview.json"])
    if overview.get("kind") != KINDS[family]:
        raise ValueError("Unsupported scientific evidence kind")
    return overview, payload


def _text(value):
    if not isinstance(value, str) or not value.strip():
        raise ValueError("Nonempty display text required")
    return value


def _list(value):
    if not isinstance(value, list) or not value:
        raise ValueError("Nonempty display list required")
    return value


def _metric(label, value, unit):
    if unit in {"windows", "cells"} and (
        type(value) is not int or not (int(unit == "cells") <= value <= 2**53 - 1)
    ):
        raise ValueError("Invalid integer count or loss of browser integer precision")
    if type(value) not in (int, float) or not math.isfinite(value):
        raise ValueError("Finite numeric metric required")
    return {"label": label, "value": value, "unit": unit}


def _study(overview, payload, family):
    def asset(name):
        if _name(name) not in payload:
            raise ValueError("Display asset is absent from authenticated bundle")
        return f"data/research/{family}/{name}"

    study = {
        "id": "vessel-passage" if family == "acoustic" else "native-cylinder",
        "title": _text(overview["title"]),
        "figures": [
            {
                "src": asset(f["path"]),
                "title": _text(f["title"]),
                "caption": _text(f["caption"]),
                "alt": _text(f["title"]) + ". " + _text(f["caption"]),
            }
            for f in _list(overview["figures"])
        ],
        "downloads": [
            {"href": asset(d["path"]), "label": _text(d["label"])}
            for d in _list(overview["downloads"])
        ],
        "limitations": [_text(x) for x in _list(overview["limitations"])],
    }
    if family == "acoustic":
        s = overview["numericalsummary"]
        study.update(
            evidence_kind="predicted_acoustic_with_synthetic_observation",
            summary="Empirical vessel and wind prediction, calibrated synthetic WAV chunks, quality flags and matched spectral comparisons. The gain control and finite-realization sampling errors are separate diagnostics, not field validation.",
            metrics=[
                _metric("Compared windows", s["compared_windows"], "windows"),
                _metric("Excluded windows", s["excluded_windows"], "windows"),
                _metric("Recovered known gain", s["gain_bias_db"], "dB"),
                _metric("Maximum gain-control error", s["gain_max_error_db"], "dB"),
                _metric("Broadband sampling RMSE", s["sampling_broadband_rmse_db"], "dB"),
                _metric(
                    "Observed − expected spectral mean",
                    s["empirical_spectral_summary"]["mean_bias_db"],
                    "dB",
                ),
            ],
        )
    else:
        b = overview["metrics"]["baseline"]
        pair = overview["metrics"]["timestep_sensitivity"]
        study.update(
            evidence_kind="native_hydrodynamic",
            summary="Two completed native OpenFOAM runs with matched mesh and physical sampling. Late-window pressure, signed lift and timestep sensitivity quantify local hydrodynamic fluctuations, not radiated acoustic pressure.",
            metrics=[
                _metric("Native volume cells", overview["configuration"]["cells"], "cells"),
                _metric("Baseline Strouhal", b["strouhal"], "dimensionless"),
                _metric("Halved-timestep Strouhal", pair["second"]["strouhal"], "dimensionless"),
                _metric("Strouhal timestep change", pair["changes_percent"]["strouhal"], "%"),
                _metric("P2 local pressure RMS", b["pressure_rms_Pa"][1], "Pa"),
                _metric("Baseline mean drag coefficient", b["mean_Cd"], "dimensionless"),
            ],
        )
        if "moviepath" in overview:
            study["video"] = {
                "src": asset(overview["moviepath"]),
                "poster": study["figures"][0]["src"],
                "title": "Native cylinder wake — authentic saved snapshots",
                "caption": _text(overview["movie"]["playback_note"]),
            }
    return study


def assemble_research_site(acoustic, cfd, output, *, acoustic_sha256, cfd_sha256):
    """Copy only authenticated bundles; differing existing destinations are refused."""
    out, a, c = (_root(p) for p in (output, acoustic, cfd))
    for source in (a, c):
        if out.is_relative_to(source) or source.is_relative_to(out):
            raise ValueError("Output must be separate from source bundles")
    expected = {"acoustic": acoustic_sha256, "cfd": cfd_sha256}
    studies, content = [], {}
    for family, source in (("acoustic", a), ("cfd", c)):
        overview, payload = _bundle(source, expected[family], family)
        studies.append(_study(overview, payload, family))
        content.update({f"{family}/{name}": raw for name, raw in payload.items()})
    catalog = {
        "schema_version": 1,
        "project_version": "0.1.0",
        "source_manifests": expected,
        "studies": studies,
    }
    content["catalog.json"] = _encode(catalog)
    content["manifest.json"] = _encode(
        {
            "schema": "oceansignatures.research-site.v1",
            "files": {name: _sha(raw) for name, raw in sorted(content.items())},
        }
    )
    if out.exists():
        existing = _inventory(out)
        if set(existing) != set(content) or any(
            existing[n].read_bytes() != raw for n, raw in content.items()
        ):
            raise FileExistsError("Refusing differing or incomplete research site output")
        return catalog
    out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".research-site-stage-", dir=out.parent) as temporary:
        staged = Path(temporary) / "products"
        staged.mkdir()
        for name, raw in content.items():
            target = staged / name
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(raw)
        _root(out)
        if out.exists():
            raise FileExistsError("Destination appeared during assembly")
        staged.rename(out)
    return catalog
