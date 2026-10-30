"""Analytical fixtures are TEST-ONLY, never native CFD evidence."""

from pathlib import Path
import importlib
import importlib.util
import hashlib
import json
import numpy as np
import pytest


def write_test_only_native_fixture(root):
    """Tiny native-FORMAT fixture, not a CFD solution or mesh-validation proof."""
    root.mkdir()
    metadata = dict(
        schema="oceansignatures-native-cylinder-v1",
        test_only=True,
        diameter_m=0.01,
        outer_radius_m=0.2,
        span_m=0.001,
        velocity_m_s=0.01,
        viscosity_m2_s=1e-6,
        density_kg_m3=1000,
        cells=3,
        end_time_s=140,
        delta_t_s=0.01,
        adaptive_timestep=False,
        probe_interval_s=0.02,
        plane_interval_s=2,
        volume_interval_s=10,
        pressure_dimensions=[0, 2, -2, 0, 0, 0, 0],
        probe_locations_m=[[0.01, 0.005, 0]],
    )
    control = (
        "application pimpleFoam; startTime 0; endTime 140; deltaT 0.01; adjustTimeStep no; "
        "rhoInf 1000; magUInf 0.01; lRef 0.01; Aref 1e-05; "
        "probes { writeControl timeStep; writeInterval 2; } "
        "forceCoeffs { writeControl timeStep; writeInterval 2; } "
        "wakePlane { writeControl timeStep; writeInterval 200; }"
    )
    files = {
        "case_metadata.json": json.dumps(metadata).encode(),
        "cylinder.foam": b"",
        "generator_executed.py": b"# TEST ONLY; no solver or real mesh\n",
        "system/controlDict": control.encode(),
        "system/blockMeshDict": b"// TEST ONLY",
        "system/fvSchemes": b"// TEST ONLY",
        "system/fvSolution": b"// TEST ONLY",
        "constant/transportProperties": b"nu [0 2 -1 0 0 0 0] 1e-06;",
        "constant/turbulenceProperties": b"simulationType laminar;",
        "0/p": b"dimensions [0 2 -2 0 0 0 0];",
        "0/U": b"dimensions [0 1 -1 0 0 0 0];",
    }
    initial = {k: hashlib.sha256(v).hexdigest() for k, v in files.items()}
    files["input_hashes.json"] = json.dumps(initial).encode()
    for name in ("points", "faces", "owner", "neighbour", "boundary"):
        files["constant/polyMesh/" + name] = b"// TEST ONLY mesh-format placeholder"
    files["log.blockMesh"] = b"TEST ONLY\nEnd\n"
    files["log.checkMesh"] = (
        b"TEST ONLY\n cells: 3\nMesh has 2 geometric (non-empty/wedge) directions (1 1 0)\n"
        b"Mesh has 2 solution (non-empty) directions (1 1 0)\nMesh OK.\nEnd\n"
    )
    files["log.pimpleFoam"] = (
        "TEST ONLY\n" + "\n".join(f"Time = {i * 0.01:.17g}" for i in range(1, 14001)) + "\nEnd\n"
    ).encode()
    files["run_record.json"] = json.dumps(
        [
            dict(
                stage=s,
                returncode=0,
                timed_out=False,
                interrupted=False,
                log_sha256=hashlib.sha256(files["log." + s]).hexdigest(),
            )
            for s in ("blockMesh", "checkMesh", "pimpleFoam")
        ]
    ).encode()
    t = np.arange(1, 7001) * 0.02
    p = (7 + 2 * np.sin(2 * np.pi * t / 5)) / 1000
    files["postProcessing/probes/0/p"] = (
        "# Probe 0 (0.01 0.005 0)\n" + "\n".join(f"{a:.17g} {b:.17g}" for a, b in zip(t, p))
    ).encode()
    files["postProcessing/probes/0/U"] = (
        "# Probe 0 (0.01 0.005 0)\n" + "\n".join(f"{a:.17g} (0.01 0 0)" for a in t)
    ).encode()
    files["postProcessing/forceCoeffs/0/coefficient.dat"] = (
        "# Time Cd Cd(f) Cd(r) Cl\n"
        + "\n".join(f"{a:.17g} 1.3 0.65 0.65 {np.sin(2 * np.pi * a / 5):.17g}" for a in t)
    ).encode()
    for time in range(2, 141, 2):
        files[f"postProcessing/wakePlane/{time}/midplane.vtp"] = (
            f"<VTKFile type='PolyData'><PolyData><FieldData><DataArray type='Float32' Name='TimeValue' format='ascii'>{time}</DataArray></FieldData><Piece NumberOfPoints='3' NumberOfPolys='2'><Points><DataArray type='Float32' Name='Points' format='ascii'>.01 0 0 .02 0 0 .01 .01 0</DataArray></Points><Polys><DataArray type='Int32' Name='connectivity' format='ascii'>0 1 2 0 0 0</DataArray><DataArray type='Int32' Name='offsets' format='ascii'>3 6</DataArray></Polys><PointData><DataArray type='Float32' Name='p' format='ascii'>.001 .002 .003</DataArray><DataArray type='Float32' Name='U' format='ascii'>.01 0 0 .01 .001 0 .01 -.001 0</DataArray></PointData></Piece></PolyData></VTKFile>"
        ).encode()
    completion = dict(
        status="completed",
        failure_reason=None,
        native_cells=3,
        requested_end_time_s=140,
        actual_end_time_s=140,
        solved_steps=14000,
        output_hashes={k: hashlib.sha256(v).hexdigest() for k, v in files.items()},
    )
    files["completion.json"] = json.dumps(completion).encode()
    for name, raw in files.items():
        dest = root / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(raw)
    return hashlib.sha256(files["completion.json"]).hexdigest()


def test_native_identity_complete_clocks_and_same_verified_bytes(tmp_path):
    m = api()
    assert hasattr(m, "read_native_case"), "verified native case reader missing"
    root = tmp_path / "TEST_ONLY"
    pin = write_test_only_native_fixture(root)
    case = m.read_native_case(root, expected_completion_sha256=pin)
    assert case.pressure_Pa.shape == (7000, 1)
    assert case.clock["resampled"] is False
    assert len(case.frames) == 21
    assert case.report["verified_original_pins"] == len(
        json.loads((root / "completion.json").read_bytes())["output_hashes"]
    )
    # A finite corruption is an identity failure, not merely a finite-value check.
    p = root / "postProcessing/probes/0/p"
    p.write_bytes(p.read_bytes().replace(b"0.007", b"0.008", 1))
    with pytest.raises(ValueError, match="identity"):
        m.read_native_case(root, expected_completion_sha256=pin)
    np.testing.assert_array_equal(
        case.pressure_Pa[0],
        [1000 * float(case.payloads["postProcessing/probes/0/p"].splitlines()[1].split()[1])],
    )


@pytest.mark.parametrize(
    "damage",
    ["empty", "missing_probe", "missing_plane", "missing_mesh", "incomplete", "wrong_anchor"],
)
def test_native_identity_fails_closed(tmp_path, damage):
    m = api()
    assert hasattr(m, "read_native_case"), "verified native case reader missing"
    root = tmp_path / "TEST_ONLY"
    pin = write_test_only_native_fixture(root)
    completion_path = root / "completion.json"
    c = json.loads(completion_path.read_bytes())
    if damage == "empty":
        c["output_hashes"] = {}
    if damage == "missing_probe":
        del c["output_hashes"]["postProcessing/probes/0/p"]
    if damage == "missing_mesh":
        del c["output_hashes"]["constant/polyMesh/points"]
    if damage == "missing_plane":
        del c["output_hashes"]["postProcessing/wakePlane/100/midplane.vtp"]
    if damage == "incomplete":
        c["status"] = "running"
    completion_path.write_text(json.dumps(c))
    pin = (
        hashlib.sha256(completion_path.read_bytes()).hexdigest()
        if damage != "wrong_anchor"
        else "a" * 64
    )
    with pytest.raises(ValueError):
        m.read_native_case(root, expected_completion_sha256=pin)


@pytest.mark.parametrize(
    "target",
    [
        "postProcessing/forceCoeffs/0/coefficient.dat",
        "postProcessing/wakePlane/100/midplane.vtp",
        "system/fvSchemes",
        "constant/polyMesh/points",
    ],
)
def test_other_finite_identity_corruption_is_rejected(tmp_path, target):
    m = api()
    root = tmp_path / "TEST_ONLY"
    pin = write_test_only_native_fixture(root)
    p = root / target
    p.write_bytes(p.read_bytes() + b"\n0\n")
    with pytest.raises(ValueError, match="identity"):
        m.read_native_case(root, expected_completion_sha256=pin)


@pytest.mark.parametrize(
    "damage",
    [
        "empty_input_pins",
        "missing_input_pin",
        "mesh_failed",
        "stage_timeout",
        "stage_pending",
        "bad_dimensions",
        "missing_solver_step",
    ],
)
def test_reanchored_test_only_fixture_semantic_gates(tmp_path, damage):
    # Verification-only coverage: these existing gates must remain fail-closed.
    m = api()
    root = tmp_path / "TEST_ONLY"
    write_test_only_native_fixture(root)
    if damage in ("empty_input_pins", "missing_input_pin"):
        path = root / "input_hashes.json"
        pins = json.loads(path.read_bytes())
        if damage == "empty_input_pins":
            pins = {}
        else:
            del pins["0/p"]
        path.write_text(json.dumps(pins))
    if damage == "mesh_failed":
        path = root / "log.checkMesh"
        path.write_bytes(path.read_bytes() + b"Failed 1 mesh checks\n")
    if damage == "missing_solver_step":
        path = root / "log.pimpleFoam"
        path.write_bytes(path.read_bytes().replace(b"Time = 1\n", b""))
    if damage in ("stage_timeout", "stage_pending"):
        path = root / "run_record.json"
        records = json.loads(path.read_bytes())
        if damage == "stage_timeout":
            records[-1]["timed_out"] = True
        else:
            records.pop()
        path.write_text(json.dumps(records))
    if damage == "bad_dimensions":
        path = root / "0/p"
        path.write_bytes(b"dimensions [1 -1 -2 0 0 0 0];")
        input_path = root / "input_hashes.json"
        pins = json.loads(input_path.read_bytes())
        pins["0/p"] = hashlib.sha256(path.read_bytes()).hexdigest()
        input_path.write_text(json.dumps(pins))
    records_path = root / "run_record.json"
    records = json.loads(records_path.read_bytes())
    for r in records:
        r["log_sha256"] = hashlib.sha256((root / ("log." + r["stage"])).read_bytes()).hexdigest()
    records_path.write_text(json.dumps(records))
    path = root / "completion.json"
    completion = json.loads(path.read_bytes())
    for name in completion["output_hashes"]:
        completion["output_hashes"][name] = hashlib.sha256((root / name).read_bytes()).hexdigest()
    path.write_text(json.dumps(completion))
    with pytest.raises(ValueError):
        m.read_native_case(
            root, expected_completion_sha256=hashlib.sha256(path.read_bytes()).hexdigest()
        )


def test_real_ffmpeg_animation_uses_each_test_only_snapshot_once(tmp_path):
    import shutil
    import subprocess

    if not shutil.which("ffmpeg") or not shutil.which("ffprobe"):
        pytest.skip("optional ffmpeg/ffprobe absent")
    m = api()
    assert hasattr(m, "render_movie"), "native snapshot movie missing"
    root = tmp_path / "TEST_ONLY"
    pin = write_test_only_native_fixture(root)
    case = m.read_native_case(root, expected_completion_sha256=pin)
    planes = [m.parse_plane(case.payloads[p], expected_time_s=t) for t, p in case.frames.items()]
    out = tmp_path / "movie"
    out.mkdir()
    record = m.render_movie(case, planes, 0.002, out)
    assert record["source_snapshot_count"] == 21
    assert len(list((out / "movie_frames").glob("*.png"))) == 21
    result = subprocess.run(
        [
            "ffprobe",
            "-v",
            "error",
            "-count_frames",
            "-select_streams",
            "v:0",
            "-show_entries",
            "stream=nb_read_frames",
            "-of",
            "json",
            str(out / record["path"]),
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    assert int(json.loads(result.stdout)["streams"][0]["nb_read_frames"]) == 21
    assert record["interpolated_fields"] is False


def test_public_cli_export_and_deterministic_portable_replay(tmp_path):
    import subprocess
    import sys

    script = Path(__file__).resolve().parents[1] / "scripts/build_cfd_results.py"
    assert script.exists(), "public CLI missing"
    root = tmp_path / "TEST_ONLY"
    pin = write_test_only_native_fixture(root)
    output = tmp_path / "first"
    command = [
        sys.executable,
        str(script),
        "--baseline",
        str(root),
        "--baseline-sha256",
        pin,
        "--output",
        str(output),
        "--no-movie",
    ]
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    overview = json.loads((output / "overview.json").read_text())
    assert overview["kind"] == "test-only analytic fixture; not native CFD"
    assert len(overview["figures"]) == 3
    for record in overview["figures"] + overview["downloads"]:
        assert (output / record["path"]).is_file()
        assert not Path(record["path"]).is_absolute()
    replay = tmp_path / "replay"
    result = subprocess.run(
        [
            sys.executable,
            str(output / "replay/build_cfd_results.py"),
            "--baseline",
            str(output / "native/baseline"),
            "--baseline-sha256",
            pin,
            "--output",
            str(replay),
            "--no-movie",
        ],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert result.returncode == 0, result.stderr
    for p in (output / "data").glob("*"):
        assert p.read_bytes() == (replay / "data" / p.name).read_bytes(), p.name
    assert (output / "analysis_summary.json").read_bytes() == (
        replay / "analysis_summary.json"
    ).read_bytes()
    assert (output / "manifest.json").is_file()
    # The documented replay command must obey the same no-raw-dotdot write guard.
    assert "--output ../" not in (output / "REPLAY.txt").read_text()
    original = (output / "overview.json").read_bytes()
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode != 0
    assert (output / "overview.json").read_bytes() == original


def test_public_cli_rejects_incomplete_and_symlink_before_creation(tmp_path):
    import subprocess
    import sys

    script = Path(__file__).resolve().parents[1] / "scripts/build_cfd_results.py"
    assert script.exists(), "public CLI missing"
    root = tmp_path / "TEST_ONLY"
    pin = write_test_only_native_fixture(root)
    p = root / "completion.json"
    c = json.loads(p.read_bytes())
    c["status"] = "running"
    p.write_text(json.dumps(c))
    pin = hashlib.sha256(p.read_bytes()).hexdigest()
    out = tmp_path / "rejected"
    result = subprocess.run(
        [
            sys.executable,
            str(script),
            "--baseline",
            str(root),
            "--baseline-sha256",
            pin,
            "--output",
            str(out),
            "--no-movie",
        ],
        capture_output=True,
    )
    assert result.returncode != 0 and not out.exists()
    external = tmp_path / "external"
    external.mkdir()
    link = tmp_path / "linked"
    link.symlink_to(external, target_is_directory=True)
    result = subprocess.run(
        [
            sys.executable,
            str(script),
            "--baseline",
            "missing",
            "--output",
            str(link / "result"),
            "--no-movie",
        ],
        capture_output=True,
    )
    assert result.returncode != 0 and not (external / "result").exists()


def test_analysis_blocks_and_matched_two_level_sensitivity(tmp_path):
    m = api()
    assert hasattr(m, "analyse_case") and hasattr(m, "compare_cases"), "analysis/comparison missing"
    root = tmp_path / "TEST_ONLY"
    pin = write_test_only_native_fixture(root)
    case = m.read_native_case(root, expected_completion_sha256=pin)
    summary, arrays = m.analyse_case(case)
    assert summary["primary"]["window_s"] == [100, 140]
    assert summary["sensitivity_80_140"]["samples"] == 3000
    assert [b["window_s"] for b in summary["blocks_20_s"]] == [
        [i, i + 20] for i in range(0, 140, 20)
    ]
    assert arrays["pressure_Pa"].shape == (7000, 1)
    import copy

    second = copy.copy(case)
    second.metadata = dict(case.metadata, delta_t_s=0.005)
    second.payloads = dict(case.payloads)
    second.payloads["system/controlDict"] = (
        case.payloads["system/controlDict"]
        .replace(b"deltaT 0.01;", b"deltaT 0.005;")
        .replace(b"writeInterval 2;", b"writeInterval 4;")
        .replace(b"writeInterval 200;", b"writeInterval 400;")
    )
    pair = m.compare_cases(case, second)
    assert pair["interpretation"] == "two-level timestep sensitivity; not convergence order"
    assert pair["changes_percent"]["mean_Cd"] == 0
    second.payloads["constant/polyMesh/points"] += b" changed"
    with pytest.raises(ValueError, match="matched"):
        m.compare_cases(case, second)


def test_packaged_subset_replays_without_original_or_cache(tmp_path):
    m = api()
    assert hasattr(m, "package_native_case"), "portable archive missing"
    root = tmp_path / "TEST_ONLY"
    pin = write_test_only_native_fixture(root)
    case = m.read_native_case(root, expected_completion_sha256=pin)
    package = tmp_path / "portable"
    manifest = m.package_native_case(case, package)
    root.rename(tmp_path / "unreachable_source")
    replay = m.read_native_case(package, expected_completion_sha256=pin)
    np.testing.assert_array_equal(replay.pressure_Pa, case.pressure_Pa)
    assert replay.frames == case.frames
    assert replay.report["verified_original_pins"] == 45
    assert replay.report["omitted_original_pins"] > 0
    assert manifest["omitted_files"]
    assert all((package / p).stat().st_size < 50_000_000 for p in manifest["archives"])
    with pytest.raises(FileExistsError):
        m.package_native_case(case, package)
    # Original expected pins survive, rather than rehashing and accepting altered p.
    meta = json.loads((package / "package.json").read_text())
    del meta["members"]["postProcessing/probes/0/p"]
    (package / "package.json").write_text(json.dumps(meta))
    with pytest.raises(ValueError):
        m.read_native_case(package, expected_completion_sha256=pin)


def test_native_vtp_connectivity_and_display_only_degeneracy(tmp_path):
    m = api()
    assert hasattr(m, "parse_plane"), "native VTP reader missing"
    root = tmp_path / "TEST_ONLY"
    write_test_only_native_fixture(root)
    raw = (root / "postProcessing/wakePlane/100/midplane.vtp").read_bytes()
    plane = m.parse_plane(raw, expected_time_s=100)
    np.testing.assert_array_equal(plane["triangles"], [[0, 1, 2], [0, 0, 0]])
    assert plane["display_mask"].tolist() == [False, True]
    assert plane["declared_dtype"] == "Float32"
    np.testing.assert_array_equal(
        plane["pressure_native_m2_s2"],
        np.array([0.001, 0.002, 0.003], dtype=np.float32).astype(float),
    )
    assert plane["report"]["zero_area_triangles"] == 1
    with pytest.raises(ValueError, match="time"):
        m.parse_plane(raw, expected_time_s=102)
    with pytest.raises(ValueError):
        m.parse_plane(raw.replace(b"0 1 2 0 0 0", b"0 1 9 0 0 0"), expected_time_s=100)


def test_analytic_fixture_only_clocks_rms_psd_and_lift_frequency():
    m = api()
    assert hasattr(m, "window_statistics"), "clock/statistics feature missing"
    # TEST-ONLY manufactured signal: stronger pressure second harmonic is not shedding.
    t = np.arange(1, 7001) * 0.02
    lifted = np.sin(2 * np.pi * t / 5)
    p = (7 + 2 * np.sin(4 * np.pi * t / 5))[:, None]
    c = m.validate_clock(t + 3e-11, 0.02, 140)
    assert c["resampled"] is False
    with pytest.raises(ValueError, match="clock"):
        m.validate_clock(t + np.where(t > 100, 1e-5, 0), 0.02, 140)
    with pytest.raises(ValueError, match="clock"):
        m.validate_clock(t[:-1], 0.02, 140)
    stats, frequency, psd = m.window_statistics(
        t + 3e-11, p, lifted, 1.3 + 0.1 * lifted, (100, 140), 0.02, 0.01, 0.01
    )
    assert stats["samples"] == 2000
    assert stats["pressure_rms_Pa"][0] == pytest.approx(np.sqrt(2), rel=1e-12)
    assert stats["lift_frequency_Hz"] == pytest.approx(0.2, rel=1e-10)
    assert stats["strouhal"] == pytest.approx(0.2, rel=1e-10)
    assert stats["pressure_psd_peak_Hz"] == [0.4]
    assert stats["psd_bin_width_Hz"] == pytest.approx(0.025)
    assert np.sum(psd[:, 0]) * 0.025 == pytest.approx(2, rel=1e-12)
    assert frequency[np.argmax(psd[:, 0])] == 0.4


def api():
    assert importlib.util.find_spec("oceansignatures.cfd_analysis") is not None, (
        "public reader missing"
    )
    return importlib.import_module("oceansignatures.cfd_analysis")


def test_extract_native_probe_bytes_and_convert_kinematic_once():
    m = api()
    raw = b"# Probe 0 (0.01 0.005 0)\n# Probe 1 (-0.02 0 0)\n0.02 0.001 -0.002\n0.04 0.003 0.004\n"
    t, p, xyz = m.parse_probes(raw)
    np.testing.assert_array_equal(t, [0.02, 0.04])
    np.testing.assert_array_equal(xyz, [[0.01, 0.005, 0], [-0.02, 0, 0]])
    np.testing.assert_array_equal(m.pressure_pa(p, 1000, [0, 2, -2, 0, 0, 0, 0]), [[1, -2], [3, 4]])
    with pytest.raises(ValueError, match="kinematic"):
        m.pressure_pa(p, 1000, [1, -1, -2, 0, 0, 0, 0])
    _, u, _ = m.parse_probes(
        raw.replace(b"0.001 -0.002", b"(1 2 3) (4 5 6)").replace(
            b"0.003 0.004", b"(7 8 9) (10 11 12)"
        ),
        vector=True,
    )
    assert u.shape == (2, 2, 3)
    with pytest.raises(ValueError):
        m.parse_probes(raw.replace(b"0.001", b"nan"))


def test_pressure_history_legend_stays_outside_scientific_axes(tmp_path, monkeypatch):
    m = api()
    root = tmp_path / "TEST_ONLY"
    pin = write_test_only_native_fixture(root)
    case = m.read_native_case(root, expected_completion_sha256=pin)
    summary, arrays = m.analyse_case(case)
    plane = m.parse_plane(case.payloads[case.frames[max(case.frames)]], expected_time_s=140)
    from matplotlib.figure import Figure

    original = Figure.savefig
    checked = []

    def inspect(figure, path, **kwargs):
        original(figure, path, **kwargs)
        if Path(path).name == "hydrodynamic_pressure.png":
            figure.set_dpi(kwargs["dpi"])
            figure.canvas.draw()
            renderer = figure.canvas.get_renderer()
            ax = figure.axes[0]
            legend_box = ax.get_legend().get_window_extent(renderer)
            assert not legend_box.overlaps(ax.get_window_extent(renderer)), (
                "pressure legend obscures the scientific time-history region"
            )
            assert figure.bbox.contains(legend_box.x0, legend_box.y0)
            assert figure.bbox.contains(legend_box.x1, legend_box.y1)
            checked.append(True)

    monkeypatch.setattr(Figure, "savefig", inspect)
    m._render_figures(
        {"baseline": case}, {"baseline": summary}, {"baseline": arrays}, [plane], tmp_path
    )
    assert checked == [True]
