#!/usr/bin/env python3
"""Fetch one explicitly redistributable NOAA PMEL example, pinned by SHA256.

No archive download; no automatic re-pinning when remote contents change.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import tempfile
import time
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from oceansignatures.reference import MAX_BYTES, analyse_reference, load_reference  # noqa: E402

SOURCE_URL = "https://www.pmel.noaa.gov/acoustics/multimedia/HB-ship-SBNMS.wav"
EXPECTED_SHA256 = "24bf234e0d302ca91fd3e31f3b964185244403f74666569b753ca12080b59750"
DEFAULT_DESTINATION = ROOT / "data" / "reference" / "HB-ship-SBNMS.wav"


def fetch_reference(
    destination=DEFAULT_DESTINATION, *, max_bytes=MAX_BYTES, timeout_s=20.0, total_timeout_s=60.0
):
    """Fetch the pinned example, decode-check it, and atomically install bytes."""
    if type(max_bytes) is not int or not 0 < max_bytes <= MAX_BYTES:
        raise ValueError(f"max_bytes must be an integer in [1, {MAX_BYTES}]")
    for value, upper in ((timeout_s, 20), (total_timeout_s, 60)):
        if not isinstance(value, (int, float)) or not 0 < value <= upper:
            raise ValueError(f"timeout must be finite, positive, and <= {upper} seconds")
    destination = Path(destination)
    deadline = time.monotonic() + total_timeout_s
    request = urllib.request.Request(SOURCE_URL, headers={"User-Agent": "OceanSignatures/0.1"})
    with urllib.request.urlopen(request, timeout=timeout_s) as response:
        if response.geturl() != SOURCE_URL:
            raise ValueError("Unexpected redirect; review source before re-pinning")
        length = response.headers.get("Content-Length")
        if length is not None and int(length) > max_bytes:
            raise ValueError("Download exceeds size limit")
        payload = bytearray()
        while True:
            if time.monotonic() >= deadline:
                raise ValueError("Transfer deadline exceeded")
            # HTTPResponse.read1 returns after at most one underlying read;
            # read(n) can otherwise wait indefinitely on a slow byte trickle.
            chunk = response.read1(min(65536, max_bytes + 1 - len(payload)))
            if time.monotonic() >= deadline:
                raise ValueError("Transfer deadline exceeded")
            if not chunk:
                break
            payload.extend(chunk)
            if len(payload) > max_bytes:
                raise ValueError("Download exceeds size limit")
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=destination.parent, suffix=".part", delete=False) as f:
            temporary = Path(f.name)
            f.write(payload)
        result = load_reference(temporary, expected_sha256=EXPECTED_SHA256, max_bytes=max_bytes)
        temporary.replace(destination)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    result["metadata"]["filename"] = destination.name
    return result["metadata"]


def main(argv=None):
    """Fetch (unless verify-only), verify the pin, and print actual QC metadata."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--destination", type=Path, default=DEFAULT_DESTINATION)
    parser.add_argument("--verify-only", action="store_true", help="No network or audio writes")
    args = parser.parse_args(argv)
    try:
        if not args.verify_only:
            fetch_reference(args.destination)
        recording = load_reference(args.destination, expected_sha256=EXPECTED_SHA256)
        analysis = analyse_reference(recording)
        print(
            json.dumps(
                {
                    "metadata": analysis["metadata"],
                    "qc": analysis["qc"],
                    "array_shapes": {
                        key: list(value.shape)
                        for key, value in analysis.items()
                        if hasattr(value, "shape")
                    },
                },
                indent=2,
                allow_nan=False,
            )
        )
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"reference error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
