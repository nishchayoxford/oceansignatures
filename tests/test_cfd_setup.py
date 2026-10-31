"""Case-generation tests; no native solver is run by pytest."""

import importlib.util
import json
import sys

import pytest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def module():
    path = ROOT / "scripts/run_cylinder_cfd.py"
    assert path.exists(), "Native cylinder case generator is not implemented"
    spec = importlib.util.spec_from_file_location("run_cylinder_cfd", path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_generate_body_fitted_2d_water_case(tmp_path):
    mod = module()
    case = tmp_path / "case"
    metadata = mod.generate_case(case, end_time=4.0)
    assert metadata["cells"] == 15360
    assert metadata["reynolds"] == pytest.approx(100.0, abs=1e-12)
    assert metadata["pressure_dimensions"] == [0, 2, -2, 0, 0, 0, 0]
    assert metadata["spanwise_model"] == "2D spanwise-uniform laminar"
    assert len(metadata["probe_locations_m"]) == 5
    mesh = (case / "system/blockMeshDict").read_text()
    assert mesh.count("hex (") == 8
    assert mesh.count("arc ") == 32
    assert "type empty" in mesh
    assert "[0 2 -2 0 0 0 0]" in (case / "0/p").read_text()
    assert "simulationType laminar;" in (case / "constant/turbulenceProperties").read_text()
    control = (case / "system/controlDict").read_text()
    assert "adjustTimeStep no;" in control
    assert "endTime 4;" in control
    assert "fields (p U);" in control
    assert "format ascii;" in control
    assert json.loads((case / "case_metadata.json").read_text()) == metadata


@pytest.mark.parametrize(
    "kwargs",
    [
        {"delta_t": 0.02 / 3},
        {"end_time": 140.000000000002, "volume_interval": 10.000000000002},
        {"velocity": 0.010000001},
    ],
)
def test_cfd01_native_scalars_round_trip_and_match_metadata(tmp_path, kwargs):
    import re

    case = tmp_path / "TEST_ONLY"
    metadata = module().generate_case(case, **kwargs)
    control = (case / "system/controlDict").read_text()

    def scalar(name):
        return float(re.search(r"\b" + name + r"\s+([^;]+);", control)[1])

    assert scalar("deltaT") == metadata["delta_t_s"]
    assert scalar("endTime") == metadata["end_time_s"]
    assert scalar("writeInterval") == metadata["volume_interval_s"]
    assert scalar("magUInf") == metadata["velocity_m_s"]
    intervals = [
        int(s) for s in re.findall(r"writeControl timeStep; writeInterval (\d+);", control)
    ]
    assert [n * scalar("deltaT") for n in intervals] == [
        metadata["probe_interval_s"],
        metadata["probe_interval_s"],
        metadata["plane_interval_s"],
    ]
    assert intervals[0] * scalar("deltaT") == pytest.approx(0.02, rel=0, abs=1e-15)
    assert round(metadata["end_time_s"] / scalar("deltaT")) * scalar("deltaT") == pytest.approx(
        metadata["end_time_s"], rel=0, abs=1e-10
    )


@pytest.mark.parametrize("kwargs", [{"delta_t": 0.02 / (3 + 1e-9)}, {"end_time": 1e308}])
def test_cfd01_clock_validation_is_finite_and_not_relative_to_step_count(tmp_path, kwargs):
    case = tmp_path / "TEST_ONLY_rejected"
    with pytest.raises(ValueError, match="clocks"):
        module().generate_case(case, **kwargs)
    assert not case.exists()


def test_bounded_command_records_real_timeout_and_nonzero_exit(tmp_path):
    mod = module()
    assert hasattr(mod, "run_command"), "Bounded subprocess runner is missing"
    result = mod.run_command(
        tmp_path, "timeout", [sys.executable, "-c", "import time; time.sleep(5)"], 0.1
    )
    assert result["timed_out"] is True
    assert result["returncode"] != 0
    assert result["wall_seconds"] < 3
    result = mod.run_command(
        tmp_path, "failure", [sys.executable, "-c", "print('witness'); raise SystemExit(7)"], 3
    )
    assert result["returncode"] == 7
    assert result["timed_out"] is False
    assert len(result["log_sha256"]) == 64
    assert (tmp_path / "log.failure").read_text().strip() == "witness"
    assert len(json.loads((tmp_path / "run_record.json").read_text())) == 2


@pytest.mark.parametrize(
    "kwargs",
    [
        {"velocity": float("nan")},
        {"delta_t": 0},
        {"end_time": -1},
        {"delta_t": 0.003},
        {"radial_cells": 2},
        {"sector_cells": 2},
    ],
)
def test_invalid_physics_and_cadence_rejected_before_write(tmp_path, kwargs):
    mod = module()
    case = tmp_path / "bad"
    with pytest.raises(ValueError):
        mod.generate_case(case, **kwargs)
    assert not case.exists()


def test_cli_generates_fresh_case_without_running_solver(tmp_path):
    import subprocess

    case = tmp_path / "cli"
    result = subprocess.run(
        [sys.executable, str(ROOT / "scripts/run_cylinder_cfd.py"), str(case), "--end-time", "2"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert (case / "system/controlDict").exists()
    assert not (case / "log.pimpleFoam").exists()


def test_native_execution_rejects_outside_cache_before_commands(tmp_path):
    mod = module()
    assert hasattr(mod, "execute_case"), "Native execution entry point is missing"
    with pytest.raises(ValueError, match=".cache/cfd"):
        mod.execute_case(tmp_path / "unowned", timeout_s=60)
    assert not (tmp_path / "unowned").exists()


def surrogate_runner(tmp_path):
    """Load identical generator bytes with an isolated __file__-derived CFD root."""
    import hashlib

    source = ROOT / "scripts/run_cylinder_cfd.py"
    copied = tmp_path / "TEST_ONLY_repo/scripts/run_cylinder_cfd.py"
    copied.parent.mkdir(parents=True)
    copied.write_bytes(source.read_bytes())
    assert (
        hashlib.sha256(copied.read_bytes()).digest() == hashlib.sha256(source.read_bytes()).digest()
    )
    spec = importlib.util.spec_from_file_location("test_only_cfd_runner", copied)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    root = copied.parents[1] / ".cache/cfd"
    root.mkdir(parents=True)
    bashrc = copied.parents[1] / "TEST_ONLY_bashrc"
    bashrc.write_text("# Recording stubs only: never source or invoke OpenFOAM.\n")
    return mod, root, bashrc


def recording_stages(mod, monkeypatch, *, interrupted_stage=None):
    """Synthetic acceptance logs, not native science or measured stage timings."""
    calls = []

    def record(case, stage, argv, timeout_s):
        calls.append(stage)
        metadata = json.loads((case / "case_metadata.json").read_text())
        text = "TEST ONLY recording stub\n"
        if stage == "checkMesh":
            text += f"cells: {metadata['cells']}\nMesh OK.\n"
        if stage == "pimpleFoam":
            text += f"Time = {metadata['end_time_s']}\nEnd\n"
        (case / ("log." + stage)).write_text(text)
        receipt = dict(
            stage=stage,
            command=argv,
            timeout_seconds=timeout_s,
            wall_seconds=0.001,
            returncode=0,
            timed_out=False,
            interrupted=stage == interrupted_stage,
            test_only=True,
        )
        path = case / "run_record.json"
        receipts = json.loads(path.read_text()) if path.exists() else []
        path.write_text(json.dumps(receipts + [receipt]))
        return receipt

    monkeypatch.setattr(mod, "run_command", record)
    return calls


def case_identity(case):
    return {
        str(p.relative_to(case)): (p.read_bytes(), p.stat().st_mtime_ns)
        for p in case.rglob("*")
        if p.is_file()
    }


@pytest.mark.parametrize("layout", ["flat", "nested", "mixed"])
def test_cfd02_every_admitted_depth_charges_realistic_prior_cases(tmp_path, monkeypatch, layout):
    mod, root, bashrc = surrogate_runner(tmp_path)
    calls = recording_stages(mod, monkeypatch)
    for i in range(2):
        parent = root if layout == "flat" or (layout == "mixed" and i == 0) else root / "group/deep"
        old = parent / f"old-{i}"
        old.mkdir(parents=True)
        receipts = [
            dict(
                stage=s,
                wall_seconds=w,
                timeout_seconds=t,
                returncode=0,
                timed_out=False,
                interrupted=False,
                test_only=True,
            )
            for s, w, t in [
                ("blockMesh", 0.5, 30),
                ("checkMesh", 0.5, 30),
                ("pimpleFoam", 1194, 1199),
            ]
        ]
        assert sum(r["wall_seconds"] for r in receipts) == 1195
        assert all(r["wall_seconds"] < r["timeout_seconds"] for r in receipts)
        assert sum(r["wall_seconds"] for r in receipts) < 1200 <= 1800
        (old / "run_record.json").write_text(json.dumps(receipts))
    case = root / ("new" if layout == "flat" else "group/deep/new")
    mod.generate_case(case)
    before = case_identity(root)
    with pytest.raises(ValueError, match="2400 s campaign cap: charged=2390"):
        mod.execute_case(case, timeout_s=60, foam_bashrc=bashrc)
    assert calls == []
    after = case_identity(root)
    # Lock creation is permitted, but no existing receipt or case byte/mtime changes.
    after.pop(".native.lock", None)
    assert after == before


@pytest.mark.parametrize(
    "records",
    [
        [{"wall_seconds": float("nan")}],
        [{"wall_seconds": float("inf")}],
        [{"wall_seconds": -1}],
        [{"wall_seconds": True}],
        [{"wall_seconds": "1"}],
        [{}],
        {},
        [],
        [{"wall_seconds": 1e308}, {"wall_seconds": 1e308}],
    ],
)
def test_cfd02_invalid_or_nonfinite_accounting_fails_closed(tmp_path, monkeypatch, records):
    mod, root, bashrc = surrogate_runner(tmp_path)
    calls = recording_stages(mod, monkeypatch)
    old = root / "old"
    old.mkdir()
    (old / "run_record.json").write_text(json.dumps(records))
    case = root / "new"
    mod.generate_case(case)
    with pytest.raises(ValueError, match="campaign accounting"):
        mod.execute_case(case, timeout_s=60, foam_bashrc=bashrc)
    assert calls == []
    assert not (case / "generator_executed.py").exists()


@pytest.mark.parametrize("symlink_kind", ["directory", "receipt"])
def test_cfd02_symlinked_history_is_not_silently_omitted(tmp_path, monkeypatch, symlink_kind):
    mod, root, bashrc = surrogate_runner(tmp_path)
    calls = recording_stages(mod, monkeypatch)
    external = tmp_path / "TEST_ONLY_history"
    external.mkdir()
    (external / "run_record.json").write_text('[{"wall_seconds": 2341}]')
    if symlink_kind == "directory":
        (root / "old").symlink_to(external, target_is_directory=True)
    else:
        (root / "old").mkdir()
        (root / "old/run_record.json").symlink_to(external / "run_record.json")
    case = root / "new"
    mod.generate_case(case)
    with pytest.raises(ValueError, match="campaign accounting"):
        mod.execute_case(case, timeout_s=60, foam_bashrc=bashrc)
    assert calls == []


@pytest.mark.parametrize("charge", [0.0, 2340.0])
def test_cfd02_fresh_nested_case_and_exact_budget_boundary_remain_allowed(
    tmp_path, monkeypatch, charge
):
    mod, root, bashrc = surrogate_runner(tmp_path)
    calls = recording_stages(mod, monkeypatch)
    for name in ("old-a", "group/old-b"):
        old = root / name
        old.mkdir(parents=True)
        (old / "run_record.json").write_text(json.dumps([{"wall_seconds": charge / 2}]))
    case = root / "group/deeper/new"
    mod.generate_case(case)
    result = mod.execute_case(case, timeout_s=60, foam_bashrc=bashrc)
    assert result["status"] == "completed"
    assert calls == ["blockMesh", "checkMesh", "pimpleFoam"]


@pytest.mark.parametrize("stage", ["blockMesh", "checkMesh", "pimpleFoam"])
def test_cfd03_zero_exit_interruption_stops_at_each_stage(tmp_path, monkeypatch, stage):
    mod, root, bashrc = surrogate_runner(tmp_path)
    calls = recording_stages(mod, monkeypatch, interrupted_stage=stage)
    case = root / "new"
    mod.generate_case(case)
    with pytest.raises(RuntimeError, match=stage + " interrupted"):
        mod.execute_case(case, timeout_s=60, foam_bashrc=bashrc)
    stages = ["blockMesh", "checkMesh", "pimpleFoam"]
    assert calls == stages[: stages.index(stage) + 1]
    completion = json.loads((case / "completion.json").read_text())
    assert completion["status"] != "completed"
    assert stage + " interrupted" in completion["failure_reason"]
    receipts = json.loads((case / "run_record.json").read_text())
    assert receipts[-1]["interrupted"] is True
    assert receipts[-1]["returncode"] == 0 and receipts[-1]["timed_out"] is False
    for later in stages[stages.index(stage) + 1 :]:
        assert not (case / ("log." + later)).exists()


def test_cfd03_real_cooperative_python_child_can_exit_zero_on_sigterm(tmp_path):
    import os
    import signal
    import subprocess
    import time

    mod, root, _ = surrogate_runner(tmp_path)
    case = root / "TEST_ONLY_signal"
    case.mkdir()
    pid_file = case / "child.pid"
    child_code = (
        "import os,sys,signal,time; from pathlib import Path; "
        "signal.signal(signal.SIGTERM, lambda s,f: sys.exit(0)); "
        + f"Path({str(pid_file)!r}).write_text(str(os.getpid())); time.sleep(10)"
    )
    wrapper_code = (
        f"import runpy; m=runpy.run_path({mod.__file__!r}); "
        f"m['run_command']({str(case)!r}, 'cooperative', "
        f"{[sys.executable, '-c', child_code]!r}, 5)"
    )
    wrapper = subprocess.Popen([sys.executable, "-c", wrapper_code])
    try:
        for _ in range(200):
            if pid_file.exists():
                break
            time.sleep(0.01)
        assert pid_file.exists(), "cooperative child did not install its handler"
        wrapper.send_signal(signal.SIGTERM)
        assert wrapper.wait(timeout=3) == 0
        receipt = json.loads((case / "run_record.json").read_text())[0]
        assert receipt["interrupted"] is True
        assert receipt["returncode"] == 0
        assert receipt["timed_out"] is False
        with pytest.raises(ProcessLookupError):
            os.kill(int(pid_file.read_text()), 0)
    finally:
        if wrapper.poll() is None:
            wrapper.kill()
            wrapper.wait()
        if pid_file.exists():
            try:
                os.kill(int(pid_file.read_text()), signal.SIGKILL)
            except ProcessLookupError:
                pass


PRIOR_EXECUTION_ARTIFACTS = [
    "generator_executed.py",
    "input_hashes.json",
    "completion.json",
    "run_record.json",
    "log.blockMesh",
    "log.checkMesh",
    "log.pimpleFoam",
    "log.other",
    "constant/polyMesh/points",
    "postProcessing/",
    "1/",
    "processor0/",
    "unexpected.dat",
]


def add_prior_artifact(case, artifact):
    path = case / artifact
    if artifact.endswith("/"):
        path.mkdir(parents=True)
    else:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            '[{"wall_seconds": 0.5}]'
            if artifact == "run_record.json"
            else "TEST ONLY: immutable historical evidence\n"
        )


def tree_identity(case):
    files = case_identity(case)
    directories = {
        str(p.relative_to(case)): p.stat().st_mtime_ns
        for p in [case, *case.rglob("*")]
        if p.is_dir()
    }
    return files, directories


@pytest.mark.parametrize("artifact", [*PRIOR_EXECUTION_ARTIFACTS, "all"])
def test_cfd04_each_prior_artifact_rejects_before_any_evidence_write(
    tmp_path, monkeypatch, artifact
):
    mod, root, bashrc = surrogate_runner(tmp_path)
    calls = recording_stages(mod, monkeypatch)
    case = root / "old"
    mod.generate_case(case)
    for name in PRIOR_EXECUTION_ARTIFACTS if artifact == "all" else [artifact]:
        add_prior_artifact(case, name)
    before = tree_identity(case)
    with pytest.raises(FileExistsError, match="prior execution artifacts"):
        mod.execute_case(case, timeout_s=60, foam_bashrc=bashrc)
    assert calls == []
    assert tree_identity(case) == before


def test_cfd04_freshness_is_rechecked_inside_lock(tmp_path, monkeypatch):
    import fcntl

    mod, root, bashrc = surrogate_runner(tmp_path)
    calls = recording_stages(mod, monkeypatch)
    case = root / "new"
    mod.generate_case(case)
    flock = fcntl.flock
    observed = []

    def acquire(lock, flags):
        flock(lock, flags)
        add_prior_artifact(case, "generator_executed.py")
        observed.append(tree_identity(case))

    monkeypatch.setattr(fcntl, "flock", acquire)
    with pytest.raises(FileExistsError, match="prior execution artifacts"):
        mod.execute_case(case, timeout_s=60, foam_bashrc=bashrc)
    assert calls == []
    assert len(observed) == 1
    assert tree_identity(case) == observed[0]


def test_cfd04_completed_surrogate_cannot_be_reentered(tmp_path, monkeypatch):
    mod, root, bashrc = surrogate_runner(tmp_path)
    calls = recording_stages(mod, monkeypatch)
    case = root / "fresh"
    mod.generate_case(case)
    assert mod.execute_case(case, timeout_s=60, foam_bashrc=bashrc)["status"] == "completed"
    assert calls == ["blockMesh", "checkMesh", "pimpleFoam"]
    before = tree_identity(case)
    calls.clear()
    with pytest.raises(FileExistsError, match="prior execution artifacts"):
        mod.execute_case(case, timeout_s=60, foam_bashrc=bashrc)
    assert calls == []
    assert tree_identity(case) == before


def test_sigterm_wrapper_cleans_serial_child_and_records_interruption(tmp_path):
    import os
    import signal
    import subprocess
    import time

    pid_file = tmp_path / "child.pid"
    child_code = (
        "import os,time; from pathlib import Path; Path("
        + repr(str(pid_file))
        + ").write_text(str(os.getpid())); time.sleep(30)"
    )
    code = (
        "import runpy; m=runpy.run_path("
        + repr(str(ROOT / "scripts/run_cylinder_cfd.py"))
        + "); m['run_command']("
        + repr(str(tmp_path))
        + ", 'signal', "
        + repr([sys.executable, "-c", child_code])
        + ", 10)"
    )
    wrapper = subprocess.Popen([sys.executable, "-c", code])
    try:
        for _ in range(100):
            if pid_file.exists():
                break
            time.sleep(0.01)
        assert pid_file.exists()
        child_pid = int(pid_file.read_text())
        wrapper.send_signal(signal.SIGTERM)
        wrapper.wait(timeout=3)
        record = tmp_path / "run_record.json"
        assert record.exists(), "SIGTERM discarded receipt and orphaned serial child"
        assert json.loads(record.read_text())[0]["interrupted"] is True
        with pytest.raises(ProcessLookupError):
            os.kill(child_pid, 0)
    finally:
        if wrapper.poll() is None:
            wrapper.kill()
            wrapper.wait()
        if pid_file.exists():
            try:
                os.kill(int(pid_file.read_text()), signal.SIGKILL)
            except ProcessLookupError:
                pass
