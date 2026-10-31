#!/usr/bin/env python3
"""Generate a small, explicit OpenCFD v2206 2-D laminar cylinder case.

Native execution is opt-in; generated cases are never synthetic flow evidence.
"""

import hashlib
import json
import math
import os
import signal
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path


def run_command(case, stage, argv, timeout_s):
    """Run one serial command with a hard wall bound and a durable receipt.

    Not an MPI supervisor: native solver invocations deliberately remain serial.
    """
    case = Path(case)
    log = case / ("log." + stage)
    if log.exists():
        raise FileExistsError(log)
    started = datetime.now(timezone.utc).isoformat()
    clock = time.monotonic()
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    timed_out = False
    interrupted = False

    def on_signal(_signum, _frame):
        nonlocal interrupted
        interrupted = True

    previous = {sig: signal.signal(sig, on_signal) for sig in (signal.SIGTERM, signal.SIGINT)}
    try:
        with log.open("wb") as output:
            proc = subprocess.Popen(
                argv,
                stdout=output,
                stderr=subprocess.STDOUT,
                start_new_session=True,
                env=env,
                cwd=case,
            )
            while proc.poll() is None and not interrupted:
                remaining = timeout_s - (time.monotonic() - clock)
                if remaining <= 0:
                    timed_out = True
                    break
                try:
                    proc.wait(timeout=min(0.1, remaining))
                except subprocess.TimeoutExpired:
                    pass
            if proc.poll() is None:
                os.killpg(proc.pid, signal.SIGTERM)
                try:
                    proc.wait(timeout=1)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid, signal.SIGKILL)
                    proc.wait()
    finally:
        for sig, handler in previous.items():
            signal.signal(sig, handler)
    receipt = dict(
        stage=stage,
        command=argv,
        started_utc=started,
        finished_utc=datetime.now(timezone.utc).isoformat(),
        timeout_seconds=timeout_s,
        timed_out=timed_out,
        interrupted=interrupted,
        returncode=proc.returncode,
        wall_seconds=time.monotonic() - clock,
        log_sha256=hashlib.sha256(log.read_bytes()).hexdigest(),
        numerical_threads=1,
        mpi_ranks=0,
    )
    path = case / "run_record.json"
    records = json.loads(path.read_text()) if path.exists() else []
    records.append(receipt)
    path.write_text(json.dumps(records, indent=2) + "\n")
    return receipt


def generate_case(
    case,
    *,
    velocity=0.01,
    end_time=140.0,
    delta_t=0.01,
    radial_cells=80,
    sector_cells=24,
    volume_interval=10.0,
):
    """Write a fresh body-fitted O-grid case, SI units throughout."""
    values = (velocity, end_time, delta_t, volume_interval)
    if not all(math.isfinite(v) and v > 0 for v in values):
        raise ValueError("Physical and time parameters must be finite and positive")
    if not 0.008 <= velocity <= 0.015 or not 0.001 <= delta_t <= 0.02:
        raise ValueError("This bounded laminar demonstration supports Re80–150 and dt=0.001–0.02 s")
    if not (
        isinstance(radial_cells, int)
        and isinstance(sector_cells, int)
        and 60 <= radial_cells <= 120
        and 20 <= sector_cells <= 32
        and 8 * radial_cells * sector_cells <= 30000
    ):
        raise ValueError("Mesh must remain a resolved 10k–30k cell O-grid")

    def scalar(value):
        # Keep historical short tokens when exact; otherwise round-trip the float.
        short = f"{value:g}"
        return short if float(short) == value else repr(float(value))

    native_dt = float(scalar(delta_t))
    for period in (0.02, 2.0, end_time, volume_interval):
        steps = float(scalar(period)) / native_dt
        if not math.isfinite(steps) or not math.isclose(
            steps, round(steps), rel_tol=0, abs_tol=1e-8
        ):
            raise ValueError(
                "Fixed dt must exactly divide finite serialized sample and output clocks"
            )
    case = Path(case)
    if any(p.is_symlink() for p in (case, *case.parents)):
        raise ValueError("Case destination must not traverse symlinks")
    case.mkdir(parents=True)
    for sub in ("0", "constant", "system"):
        (case / sub).mkdir()

    def put(name, body, cls="dictionary"):
        header = (
            "FoamFile\n{\n version 2.0; format ascii; class "
            + cls
            + "; object "
            + Path(name).name
            + ";\n}\n"
        )
        (case / name).write_text(header + body + "\n")

    diameter, radius, outer, span, nu, rho = 0.01, 0.005, 0.2, 0.001, 1e-6, 1000.0
    angles = [math.radians(-90 + 45 * i) for i in range(8)]
    vertices = []
    for z in (-span / 2, span / 2):
        for r in (radius, outer):
            vertices.extend((r * math.cos(a), r * math.sin(a), z) for a in angles)

    def vec(v):
        return "(" + " ".join(f"{x:.17g}" for x in v) + ")"

    blocks, arcs = [], []
    patches = {"cylinder": [], "inlet": [], "outlet": [], "frontAndBack": []}
    for i in range(8):
        j = (i + 1) % 8
        a, b, c, d = i, i + 8, j + 8, j
        blocks.append(
            f"hex ({a} {b} {c} {d} {a + 16} {b + 16} {c + 16} {d + 16}) "
            f"({radial_cells} {sector_cells} 1) simpleGrading (50 1 1)"
        )
        patches["cylinder"].append((d, a, a + 16, d + 16))
        patches["outlet" if i < 4 else "inlet"].append((b, c, c + 16, b + 16))
        patches["frontAndBack"].extend([(d, c, b, a), (a + 16, b + 16, c + 16, d + 16)])
        mid = angles[i] + math.pi / 8
        for layer, z in ((0, -span / 2), (16, span / 2)):
            for off, r in ((0, radius), (8, outer)):
                arcs.append(
                    f"arc {i + off + layer} {j + off + layer} "
                    + vec((r * math.cos(mid), r * math.sin(mid), z))
                )
    boundary = []
    for name, faces in patches.items():
        typ = "wall" if name == "cylinder" else "empty" if name == "frontAndBack" else "patch"
        boundary.append(
            name + " { type " + typ + "; faces (\n" + "\n".join(vec(f) for f in faces) + "\n); }"
        )
    put(
        "system/blockMeshDict",
        "convertToMeters 1;\nvertices (\n"
        + "\n".join(vec(v) for v in vertices)
        + "\n);\nblocks (\n"
        + "\n".join(blocks)
        + "\n);\nedges (\n"
        + "\n".join(arcs)
        + "\n);\nboundary (\n"
        + "\n".join(boundary)
        + "\n);\nmergePatchPairs ();",
    )
    put(
        "constant/transportProperties", f"transportModel Newtonian;\nnu [0 2 -1 0 0 0 0] {nu:.17g};"
    )
    put("constant/turbulenceProperties", "simulationType laminar;")
    put(
        "0/p",
        "dimensions [0 2 -2 0 0 0 0];\ninternalField uniform 0;\n"
        "boundaryField { cylinder { type zeroGradient; } inlet { type zeroGradient; } "
        "outlet { type fixedValue; value uniform 0; } frontAndBack { type empty; } }",
        "volScalarField",
    )
    put(
        "0/U",
        f"dimensions [0 1 -1 0 0 0 0];\ninternalField uniform ({velocity:.17g} {velocity * 0.001:.17g} 0);\n"
        "boundaryField { cylinder { type fixedValue; value uniform (0 0 0); } "
        f"inlet {{ type fixedValue; value uniform ({velocity:.17g} 0 0); }} "
        f"outlet {{ type inletOutlet; inletValue uniform ({velocity:.17g} 0 0); value uniform ({velocity:.17g} 0 0); }} "
        "frontAndBack { type empty; } }",
        "volVectorField",
    )
    put(
        "system/fvSchemes",
        """
    ddtSchemes { default backward; }
    gradSchemes { default Gauss linear; }
    divSchemes { default none; div(phi,U) Gauss linear;
        div((nuEff*dev2(T(grad(U))))) Gauss linear; }
    laplacianSchemes { default Gauss linear corrected; }
    interpolationSchemes { default linear; }
    snGradSchemes { default corrected; }
    fluxRequired { default no; p; }
    """,
    )
    put(
        "system/fvSolution",
        """
    solvers {
        p { solver GAMG; tolerance 1e-9; relTol 0.05; smoother GaussSeidel; }
        pFinal { $p; relTol 0; }
        U { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-10; relTol 0; }
        UFinal { $U; }
    }
    PIMPLE { nOuterCorrectors 1; nCorrectors 2; nNonOrthogonalCorrectors 0;
        momentumPredictor yes; pRefCell 0; pRefValue 0; }
    """,
    )
    probes = [
        [0.01, 0.005, 0.0],
        [0.02, 0.005, 0.0],
        [0.04, 0.005, 0.0],
        [0.02, -0.005, 0.0],
        [-0.02, 0.0, 0.0],
    ]
    probe_interval = round(0.02 / delta_t)
    plane_interval = round(2.0 / delta_t)
    put(
        "system/controlDict",
        f"""
    application pimpleFoam;
    startFrom startTime; startTime 0;
    stopAt endTime; endTime {scalar(end_time)}; deltaT {scalar(delta_t)};
    writeControl runTime; writeInterval {scalar(volume_interval)};
    purgeWrite 0; writeFormat ascii; writePrecision 12; writeCompression off;
    timeFormat general; timePrecision 12; runTimeModifiable no;
    adjustTimeStep no;
    functions {{
        probes {{
            type probes; libs ("libsampling.so");
            writeControl timeStep; writeInterval {probe_interval};
            writePrecision 12; fixedLocations true; interpolationScheme cell;
            fields (p U); probeLocations ({" ".join(vec(p) for p in probes)});
        }}
        forceCoeffs {{
            type forceCoeffs; libs ("libforces.so");
            patches (cylinder); p p; U U; rho rhoInf; rhoInf {rho:g};
            CofR (0 0 0); liftDir (0 1 0); dragDir (1 0 0); pitchAxis (0 0 1);
            magUInf {scalar(velocity)}; lRef {diameter:g}; Aref {diameter * span:g};
            writeControl timeStep; writeInterval {probe_interval};
            writePrecision 12; log false;
        }}
        wakePlane {{
            type surfaces; libs ("libsampling.so");
            writeControl timeStep; writeInterval {plane_interval};
            surfaceFormat vtk; interpolationScheme cellPoint; fields (p U);
            formatOptions {{ vtk {{ format ascii; precision 17; }} }}
            surfaces (midplane {{
                type cuttingPlane; planeType pointAndNormal;
                pointAndNormalDict {{ basePoint (0 0 0); normalVector (0 0 1); }}
                interpolate true; triangulate true;
                isoMethod topo; regularise false; snap false;
            }});
        }}
    }}
    """,
    )
    metadata = dict(
        schema="oceansignatures-native-cylinder-v1",
        diameter_m=diameter,
        outer_radius_m=outer,
        span_m=span,
        velocity_m_s=velocity,
        viscosity_m2_s=nu,
        density_kg_m3=rho,
        reynolds=velocity * diameter / nu,
        cells=8 * radial_cells * sector_cells,
        radial_cells=radial_cells,
        circumferential_cells=8 * sector_cells,
        radial_expansion_ratio=50,
        spanwise_model="2D spanwise-uniform laminar",
        pressure_dimensions=[0, 2, -2, 0, 0, 0, 0],
        pressure_units="m2/s2 (kinematic); multiply once by density for Pa",
        end_time_s=end_time,
        delta_t_s=delta_t,
        adaptive_timestep=False,
        probe_interval_s=probe_interval * delta_t,
        plane_interval_s=plane_interval * delta_t,
        volume_interval_s=volume_interval,
        probe_locations_m=probes,
        initial_transverse_velocity_m_s=0.001 * velocity,
        perturbation="Uniform interior Uy=0.001*U at t=0 only; inlet Uy=0. Not a time-series forcing.",
        interpretation="Native incompressible hydrodynamic pressure, NOT propagating acoustic pressure.",
    )
    (case / "case_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (case / "cylinder.foam").touch()
    return metadata


def execute_case(case, *, timeout_s=900, foam_bashrc=None):
    """Execute serial native stages, bounded by the shared 2400 s campaign cap."""
    import fcntl
    import re
    import shutil

    root = Path(__file__).resolve().parents[1] / ".cache/cfd"
    case = Path(case).absolute()
    if any(p.is_symlink() for p in (case, *case.parents)) or not case.resolve().is_relative_to(
        root.resolve()
    ):
        raise ValueError("Native runs must be fresh cases beneath this repository's .cache/cfd")
    if not math.isfinite(timeout_s) or not 0 < timeout_s <= 1800:
        raise ValueError("Native invocation wall bound must be in (0, 1800] seconds")
    foam_bashrc = Path(foam_bashrc or Path.home() / "OpenFOAM/OpenFOAM-v2206/etc/bashrc")
    if not foam_bashrc.is_file():
        raise FileNotFoundError(foam_bashrc)
    metadata = json.loads((case / "case_metadata.json").read_text())
    with (root / ".native.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        # Only the generator's initial inventory is admissible. Check under the
        # campaign lock, before copying identities or entering completion writes.
        initial_files = {
            "0/p",
            "0/U",
            "constant/transportProperties",
            "constant/turbulenceProperties",
            "system/blockMeshDict",
            "system/fvSchemes",
            "system/fvSolution",
            "system/controlDict",
            "case_metadata.json",
            "cylinder.foam",
        }
        for path in case.rglob("*"):
            name = path.relative_to(case).as_posix()
            if path.is_symlink() or not (
                (name in initial_files and path.is_file())
                or (name in {"0", "constant", "system"} and path.is_dir())
            ):
                raise FileExistsError(
                    f"Case contains prior execution artifacts or unexpected entries: {path}"
                )
        charged = 0.0
        for path in sorted(root.rglob("*")):
            if path.is_symlink():
                raise ValueError(f"Invalid campaign accounting: symlink {path}")
            if path.name != "run_record.json":
                continue
            try:
                records = json.loads(path.read_text())
                if not isinstance(records, list) or not records:
                    raise ValueError("expected nonempty receipt list")
                for receipt in records:
                    wall = receipt["wall_seconds"]
                    if type(wall) not in (int, float) or not math.isfinite(wall) or wall < 0:
                        raise ValueError("wall_seconds must be finite, numeric and nonnegative")
                    charged += wall
                    if not math.isfinite(charged):
                        raise ValueError("nonfinite cumulative wall_seconds")
            except (OSError, ValueError, TypeError, KeyError, OverflowError) as exc:
                raise ValueError(f"Invalid campaign accounting in {path}: {exc}") from exc
        if charged + timeout_s > 2400:
            raise ValueError(f"2400 s campaign cap: charged={charged}, requested={timeout_s}")
        if sum(p.stat().st_size for p in root.rglob("*") if p.is_file()) > 4_500_000_000:
            raise ValueError("Campaign is approaching its 5 GB output cap")
        snapshot = case / "generator_executed.py"
        shutil.copyfile(__file__, snapshot)
        inputs = {
            str(p.relative_to(case)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in case.rglob("*")
            if p.is_file()
        }
        (case / "input_hashes.json").write_text(json.dumps(inputs, indent=2) + "\n")
        clock = time.monotonic()
        status = "failed"
        reason = None
        mesh_cells = None
        try:
            for stage, opts in [
                ("blockMesh", []),
                ("checkMesh", ["-allTopology", "-allGeometry"]),
                ("pimpleFoam", []),
            ]:
                remaining = timeout_s - (time.monotonic() - clock)
                if remaining <= 0:
                    raise RuntimeError("Native case wall budget exhausted")
                # bashrc forwards positional parameters, so clear them during source.
                shell = 'rc=$1; shift; argv=("$@"); set --; source "$rc" >/dev/null 2>&1 || exit $?; exec "${argv[@]}"'
                argv = [
                    "bash",
                    "-c",
                    shell,
                    "native",
                    str(foam_bashrc),
                    stage,
                    "-case",
                    str(case),
                    *opts,
                ]
                receipt = run_command(
                    case, stage, argv, min(remaining, 30) if stage != "pimpleFoam" else remaining
                )
                if receipt["interrupted"]:
                    raise RuntimeError(f"{stage} interrupted; see preserved log")
                text = (case / ("log." + stage)).read_text()
                if receipt["returncode"] != 0 or receipt["timed_out"]:
                    raise RuntimeError(f"{stage} failed/timeout; see preserved log")
                if stage == "checkMesh":
                    if "Mesh OK." not in text or re.search(r"Failed\s+\d+\s+mesh checks", text):
                        raise RuntimeError("Full checkMesh acceptance failed")
                    cell_match = re.search(r"cells:\s+(\d+)", text)
                    if cell_match is None:
                        raise RuntimeError("Full mesh log contains no native cell count")
                    mesh_cells = int(cell_match.group(1))
                    if mesh_cells != metadata["cells"]:
                        raise RuntimeError("Native cell count does not match requested mesh")
                if stage == "pimpleFoam":
                    times = [float(t) for t in re.findall(r"^Time = ([0-9.eE+-]+)", text, re.M)]
                    if (
                        not re.search(r"^End\s*$", text, re.M)
                        or not times
                        or not math.isclose(times[-1], metadata["end_time_s"], abs_tol=1e-8)
                    ):
                        raise RuntimeError("Native completion/end-time gate failed")
            status = "completed"
        except Exception as exc:
            reason = str(exc)
        log = (case / "log.pimpleFoam").read_text() if (case / "log.pimpleFoam").exists() else ""
        times = [float(t) for t in re.findall(r"^Time = ([0-9.eE+-]+)", log, re.M)]
        co = [float(c) for c in re.findall(r"Courant Number mean:.*?max: ([0-9.eE+-]+)", log)]
        hashes = {
            str(p.relative_to(case)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(case.rglob("*"))
            if p.is_file() and p.name != "completion.json"
        }
        report = dict(
            status=status,
            failure_reason=reason,
            native_cells=mesh_cells,
            requested_end_time_s=metadata["end_time_s"],
            actual_end_time_s=times[-1] if times else None,
            solved_steps=len(times),
            logged_max_courant=max(co) if co else None,
            courant_note="pimpleFoam reports pre-step Courant; no final post-solve Courant claim",
            wall_seconds=time.monotonic() - clock,
            output_hashes=hashes,
            output_bytes=sum(p.stat().st_size for p in case.rglob("*") if p.is_file()),
        )
        (case / "completion.json").write_text(json.dumps(report, indent=2) + "\n")
        if status != "completed":
            raise RuntimeError(reason)
        return {k: v for k, v in report.items() if k != "output_hashes"}


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path)
    parser.add_argument("--velocity", type=float, default=0.01)
    parser.add_argument("--end-time", type=float, default=140.0)
    parser.add_argument("--delta-t", type=float, default=0.01)
    parser.add_argument("--radial-cells", type=int, default=80)
    parser.add_argument("--sector-cells", type=int, default=24)
    parser.add_argument("--volume-interval", type=float, default=10.0)
    parser.add_argument("--execute", action="store_true", help="Run native stages (serial only)")
    parser.add_argument("--timeout-s", type=float, default=900)
    parser.add_argument("--foam-bashrc", type=Path)
    args = parser.parse_args()
    result = generate_case(
        args.case,
        velocity=args.velocity,
        end_time=args.end_time,
        delta_t=args.delta_t,
        radial_cells=args.radial_cells,
        sector_cells=args.sector_cells,
        volume_interval=args.volume_interval,
    )
    print(json.dumps(result, indent=2))
    if args.execute:
        print(
            json.dumps(
                execute_case(args.case, timeout_s=args.timeout_s, foam_bashrc=args.foam_bashrc),
                indent=2,
            )
        )
