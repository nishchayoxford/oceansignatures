#!/usr/bin/env python3
"""Public read-only native CFD analysis/export entry point; never launches CFD."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from oceansignatures.cfd_analysis import main

if __name__ == "__main__":
    raise SystemExit(main())
