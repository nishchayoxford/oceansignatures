#!/usr/bin/env python3
"""Rebuild public numerical assets without a browser or development cache."""

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from oceansignatures.export import build_site  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("site"),
        help="Generated data/ and plots/ destination; existing frontend is preserved",
    )
    args = parser.parse_args()
    catalog = build_site(args.output)
    files = [
        p for folder in ("data", "plots") for p in (args.output / folder).rglob("*") if p.is_file()
    ]
    print(
        json.dumps(
            {
                "output": str(args.output),
                "passive_cases": len(catalog["passive_cases"]),
                "active_cases": len(catalog["active_cases"]),
                "gallery": len(catalog["gallery"]),
                "png_figures": sum(p.suffix == ".png" for p in files),
                "generated_files": len(files),
                "bytes": sum(p.stat().st_size for p in files),
                "analytic_checks": len(catalog["verification"]["checks"]),
                "verification": catalog["verification"]["status"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
