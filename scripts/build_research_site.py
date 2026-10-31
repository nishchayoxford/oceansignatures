#!/usr/bin/env python3
"""Assemble verified acoustic/CFD bundles; never compute or run a solver."""

import argparse
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from oceansignatures.research_site import assemble_research_site


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acoustic", type=Path, required=True)
    parser.add_argument("--cfd", type=Path, required=True)
    parser.add_argument(
        "--acoustic-sha256", required=True, help="Accepted acoustic manifest digest"
    )
    parser.add_argument("--cfd-sha256", required=True, help="Accepted CFD manifest digest")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    catalog = assemble_research_site(
        args.acoustic,
        args.cfd,
        args.output,
        acoustic_sha256=args.acoustic_sha256,
        cfd_sha256=args.cfd_sha256,
    )
    manifest = (args.output / "manifest.json").read_bytes()
    print(
        json.dumps(
            {
                "project_version": catalog["project_version"],
                "studies": len(catalog["studies"]),
                "files": len(json.loads(manifest)["files"]) + 1,
                "manifest_sha256": hashlib.sha256(manifest).hexdigest(),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
