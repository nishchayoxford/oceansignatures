#!/usr/bin/env python3
"""Build a new empirical-prediction / synthetic-verification bundle offline."""

import argparse
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from oceansignatures.research import export_research  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="New directory; an exactly identical replay is accepted",
    )
    parser.add_argument("--seed", type=int, default=1729)
    parser.add_argument(
        "--gain-db",
        type=float,
        default=6.0,
        help="Synthetic pressure-amplitude discrepancy in dB (default +6)",
    )
    args = parser.parse_args()
    try:
        overview = export_research(args.output, seed=args.seed, gain_db=args.gain_db)
    except (ValueError, OSError) as exc:
        parser.exit(2, f"Research export refused: {exc}\n")
    manifest = args.output / "manifest.json"
    report = dict(
        export_directory=str(args.output),
        files_total=sum(p.is_file() for p in args.output.rglob("*")),
        figures=len(overview["figures"]),
        downloads=len(overview["downloads"]),
        manifest_sha256=hashlib.sha256(manifest.read_bytes()).hexdigest(),
        numericalsummary=overview["numericalsummary"],
    )
    print(json.dumps(report, sort_keys=True, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
