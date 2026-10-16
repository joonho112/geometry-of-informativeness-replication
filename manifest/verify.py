#!/usr/bin/env python3
"""Check the supplied input files against their recorded checksums."""
import csv
import hashlib
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]


def main():
    ledger = ROOT / "manifest/pipeline_manifest.csv"
    with ledger.open(newline="") as stream:
        records = list(csv.DictReader(stream))
    failures = []
    for item in records:
        path = ROOT / item["path"]
        if not path.is_file():
            failures.append(f"Missing: {item['path']}")
        elif (
            path.stat().st_size != int(item["bytes"])
            or hashlib.sha256(path.read_bytes()).hexdigest() != item["sha256"]
        ):
            failures.append(f"Changed: {item['path']}")
    known = {x["path"] for x in records}
    for folder in ("analysis/inputs", "reference"):
        for path in (ROOT / folder).rglob("*"):
            if (
                path.is_file()
                and path.name != ".DS_Store"
                and path.relative_to(ROOT).as_posix() not in known
            ):
                failures.append(f"Unlisted input: {path.relative_to(ROOT)}")
    if failures:
        print("\n".join(failures), file=sys.stderr)
        return 1
    print(f"PASS: {len(records)} input and reference files match their SHA-256 checksums.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
