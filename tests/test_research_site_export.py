"""Test-only packaging fixtures, never public scientific results."""

import hashlib
import importlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import pytest


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def bundle(root, *, native=False):
    root.mkdir(parents=True)
    overview = {
        "title": "TEST ONLY packaging fixture",
        "kind": "native hydrodynamic CFD"
        if native
        else "predictive approximations + synthetic verification",
        "figures": [
            {"path": "figure.png", "title": "TEST ONLY", "caption": "Not a scientific result."}
        ],
        "downloads": [{"path": "values.csv", "label": "TEST ONLY data"}],
        "limitations": ["Synthetic packaging fixture; never promote this fixture."],
    }
    if native:
        overview["configuration"] = {"cells": 3}
        overview["metrics"] = {
            "baseline": {"strouhal": 0.2, "mean_Cd": 1.3, "pressure_rms_Pa": [0.001, 0.002]},
            "timestep_sensitivity": {
                "second": {"strouhal": 0.201},
                "changes_percent": {"strouhal": 0.5},
            },
        }
    else:
        overview["numericalsummary"] = {
            "compared_windows": 2,
            "excluded_windows": 1,
            "gain_bias_db": 6.0,
            "gain_max_error_db": 1e-14,
            "sampling_broadband_rmse_db": 0.2,
            "empirical_spectral_summary": {"mean_bias_db": 5.7},
        }
    payloads = {
        "overview.json": json.dumps(overview).encode(),
        "figure.png": b"TEST ONLY nonimage payload",
        "values.csv": b"test_only\n1\n",
    }
    for name, raw in payloads.items():
        (root / name).write_bytes(raw)
    manifest = {
        "schema": "sha256-files-v1" if native else "oceansignatures.research.manifest.v1",
        "files": {name: sha(raw) for name, raw in payloads.items()},
    }
    raw = json.dumps(manifest).encode()
    (root / "manifest.json").write_bytes(raw)
    return sha(raw)


def api():
    assert importlib.util.find_spec("oceansignatures.research_site") is not None, (
        "research site assembler missing"
    )
    return importlib.import_module("oceansignatures.research_site")


def test_real_payload_assembly_is_additive_and_preserves_identity(tmp_path):
    acoustic, native = tmp_path / "acoustic", tmp_path / "native"
    apin, cpin = bundle(acoustic), bundle(native, native=True)
    old = tmp_path / "site/data/catalog.json"
    old.parent.mkdir(parents=True)
    old.write_bytes(b"original experiment sentinel")
    output = old.parent / "research"
    catalog = api().assemble_research_site(
        acoustic, native, output, acoustic_sha256=apin, cfd_sha256=cpin
    )
    assert old.read_bytes() == b"original experiment sentinel"
    assert catalog["project_version"] == "0.1.0"
    assert [s["evidence_kind"] for s in catalog["studies"]] == [
        "predicted_acoustic_with_synthetic_observation",
        "native_hydrodynamic",
    ]
    assert catalog["source_manifests"] == {"acoustic": apin, "cfd": cpin}
    assert catalog["studies"][0]["metrics"][2]["value"] == 6.0
    assert catalog["studies"][1]["metrics"][1]["value"] == 0.2
    for source, prefix in ((acoustic, "acoustic"), (native, "cfd")):
        for p in source.iterdir():
            assert (output / prefix / p.name).read_bytes() == p.read_bytes()
    assert catalog["studies"][0]["figures"][0]["src"] == "data/research/acoustic/figure.png"
    for study in catalog["studies"]:
        assert study["figures"][0]["alt"]
    assert json.loads((output / "catalog.json").read_text()) == catalog
    pins = json.loads((output / "manifest.json").read_text())["files"]
    assert pins == {
        p.relative_to(output).as_posix(): sha(p.read_bytes())
        for p in output.rglob("*")
        if p.is_file() and p != output / "manifest.json"
    }


def test_public_assembly_cli_and_browser_schema(tmp_path):
    script = Path(__file__).resolve().parents[1] / "scripts/build_research_site.py"
    assert script.is_file(), "public research assembly CLI missing"
    acoustic, native, output = tmp_path / "acoustic", tmp_path / "native", tmp_path / "output"
    apin, cpin = bundle(acoustic), bundle(native, native=True)
    result = subprocess.run(
        [
            sys.executable,
            str(script),
            "--acoustic",
            str(acoustic),
            "--cfd",
            str(native),
            "--acoustic-sha256",
            apin,
            "--cfd-sha256",
            cpin,
            "--output",
            str(output),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout)["studies"] == 2
    adapter = script.parents[1] / "site/research.mjs"
    node = subprocess.run(
        [
            "node",
            "--input-type=module",
            "-e",
            f"import fs from 'node:fs'; import {{validateResearch}} from {json.dumps(adapter.as_uri())}; validateResearch(JSON.parse(fs.readFileSync(process.argv[1], 'utf8')));",
            str(output / "catalog.json"),
        ],
        capture_output=True,
        text=True,
    )
    assert node.returncode == 0, node.stderr


@pytest.mark.parametrize("damage", ["wrong_anchor", "corrupt", "unlisted", "symlink", "missing"])
def test_authenticated_bundle_rejection_is_before_destination_creation(tmp_path, damage):
    # Verification-only coverage of identity gates used by the first vertical slice.
    acoustic, native = tmp_path / "acoustic", tmp_path / "native"
    apin, cpin = bundle(acoustic), bundle(native, native=True)
    if damage == "wrong_anchor":
        apin = "a" * 64
    elif damage == "corrupt":
        (acoustic / "values.csv").write_bytes(b"changed finite value\n2\n")
    elif damage == "unlisted":
        (acoustic / "extra.txt").write_bytes(b"unlisted")
    elif damage == "symlink":
        (acoustic / "values.csv").unlink()
        (acoustic / "values.csv").symlink_to(native / "values.csv")
    else:
        (acoustic / "figure.png").unlink()
    output = tmp_path / "absent-parent/output"
    with pytest.raises(ValueError):
        api().assemble_research_site(
            acoustic, native, output, acoustic_sha256=apin, cfd_sha256=cpin
        )
    assert not output.parent.exists()


def test_identical_repeat_and_conflicting_output_preservation(tmp_path):
    acoustic, native, output = tmp_path / "acoustic", tmp_path / "native", tmp_path / "output"
    apin, cpin = bundle(acoustic), bundle(native, native=True)
    kwargs = dict(acoustic_sha256=apin, cfd_sha256=cpin)
    first = api().assemble_research_site(acoustic, native, output, **kwargs)
    assert api().assemble_research_site(acoustic, native, output, **kwargs) == first
    victim = output / "acoustic/values.csv"
    victim.write_bytes(b"preserve the existing differing data")
    with pytest.raises(FileExistsError):
        api().assemble_research_site(acoustic, native, output, **kwargs)
    assert victim.read_bytes() == b"preserve the existing differing data"


def test_destination_symlink_rejected_before_missing_source_reads(tmp_path):
    external = tmp_path / "external"
    external.mkdir()
    link = tmp_path / "link"
    link.symlink_to(external, target_is_directory=True)
    with pytest.raises(ValueError, match="Unsafe"):
        api().assemble_research_site(
            "missing", "missing", link / "output", acoustic_sha256="a" * 64, cfd_sha256="b" * 64
        )
    assert not (external / "output").exists()


def revise_overview(root, change):
    path = root / "overview.json"
    data = json.loads(path.read_bytes())
    change(data)
    path.write_text(json.dumps(data))
    manifest_path = root / "manifest.json"
    manifest = json.loads(manifest_path.read_bytes())
    manifest["files"]["overview.json"] = sha(path.read_bytes())
    manifest_path.write_text(json.dumps(manifest))
    return sha(manifest_path.read_bytes())


@pytest.mark.parametrize(
    "field,value",
    [("compared_windows", -1), ("excluded_windows", 1.5), ("cells", 0), ("cells", 2**53 + 1)],
)
def test_counts_must_have_their_physical_integer_meaning(tmp_path, field, value):
    acoustic, native, output = (
        tmp_path / "acoustic",
        tmp_path / "native",
        tmp_path / "absent/output",
    )
    apin, cpin = bundle(acoustic), bundle(native, native=True)
    if field == "cells":
        cpin = revise_overview(native, lambda d: d["configuration"].update(cells=value))
    else:
        apin = revise_overview(acoustic, lambda d: d["numericalsummary"].update({field: value}))
    with pytest.raises(ValueError, match="count"):
        api().assemble_research_site(
            acoustic, native, output, acoustic_sha256=apin, cfd_sha256=cpin
        )
    assert not output.parent.exists()


@pytest.mark.parametrize("damage", ["evidence", "asset", "metric", "limitations"])
def test_display_contract_rejection_before_writes(tmp_path, damage):
    # Verification-only coverage; changed expected pins belong to TEST ONLY fixtures.
    acoustic, native, output = (
        tmp_path / "acoustic",
        tmp_path / "native",
        tmp_path / "absent/output",
    )
    apin, cpin = bundle(acoustic), bundle(native, native=True)

    def mutate(d):
        if damage == "evidence":
            d["kind"] = "field measurement"
        elif damage == "asset":
            d["figures"][0]["path"] = "not-in-the-manifest.png"
        elif damage == "metric":
            d["numericalsummary"]["gain_bias_db"] = True
        else:
            d["limitations"] = []

    apin = revise_overview(acoustic, mutate)
    with pytest.raises(ValueError):
        api().assemble_research_site(
            acoustic, native, output, acoustic_sha256=apin, cfd_sha256=cpin
        )
    assert not output.parent.exists()


def test_staging_failure_never_promotes_partial_results(tmp_path, monkeypatch):
    acoustic, native, output = tmp_path / "acoustic", tmp_path / "native", tmp_path / "output"
    apin, cpin = bundle(acoustic), bundle(native, native=True)
    original = Path.write_bytes

    def failure(path, data):
        if path.name == "values.csv" and any(
            p.name.startswith(".research-site-stage-") for p in path.parents
        ):
            raise OSError("test-only simulated staging failure")
        return original(path, data)

    monkeypatch.setattr(Path, "write_bytes", failure)
    with pytest.raises(OSError, match="staging failure"):
        api().assemble_research_site(
            acoustic, native, output, acoustic_sha256=apin, cfd_sha256=cpin
        )
    assert not output.exists()
    assert not list(tmp_path.glob(".research-site-stage-*"))
    assert sha((acoustic / "manifest.json").read_bytes()) == apin


def test_ci_includes_research_frontend_and_real_catalog_gate():
    workflow = (Path(__file__).resolve().parents[1] / ".github/workflows/pages.yml").read_text()
    commands = [
        line.strip()[7:] for line in workflow.splitlines() if line.strip().startswith("- run: ")
    ]
    assert any("node --test" in c and "tests/test_research_site.mjs" in c for c in commands), (
        "CI omits the additive research frontend tests"
    )
    assert "validateResearch" in workflow and "site/data/research/manifest.json" in workflow
