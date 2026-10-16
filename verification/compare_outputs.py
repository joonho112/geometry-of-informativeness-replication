#!/usr/bin/env python3
"""Compare recalculated CSV files with the results used in the paper."""
import argparse
import csv
import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]


def equivalent(expected, actual, atol, rtol):
    if expected == actual:
        return True
    try:
        left, right = float(expected), float(actual)
    except ValueError:
        return False
    if math.isnan(left) and math.isnan(right):
        return True
    return math.isclose(left, right, rel_tol=rtol, abs_tol=atol)


def compare(item):
    actual = ROOT / item["actual"]
    if not actual.is_file():
        return {
            "file": item["actual"],
            "passed": False,
            "error": "missing generated output",
        }
    with (ROOT / item["expected"]).open(newline="") as stream:
        expected_rows = list(csv.reader(stream))
    with actual.open(newline="") as stream:
        actual_rows = list(csv.reader(stream))
    errors = []
    cells = 0
    if len(expected_rows) != len(actual_rows):
        errors.append(f"row count {len(actual_rows)} != {len(expected_rows)}")
    for row, (left, right) in enumerate(zip(expected_rows, actual_rows), 1):
        if len(left) != len(right):
            errors.append(f"row {row}: column count mismatch")
            continue
        for column, (a, b) in enumerate(zip(left, right), 1):
            cells += 1
            equal = (
                a == b
                if row == 1
                else equivalent(
                    a, b, item["absolute_tolerance"], item["relative_tolerance"]
                )
            )
            if not equal and len(errors) < 10:
                errors.append(f"row {row}, column {column}: {b!r} != {a!r}")
    return {
        "file": item["actual"],
        "passed": not errors,
        "cells": cells,
        "errors": errors,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--group",
        choices=["simulation", "posterior", "sensitivity", "firststage", "all"],
        required=True,
    )
    args = parser.parse_args()
    items = json.loads((ROOT / "verification/output-comparisons.json").read_text())
    results = [
        compare(x) for x in items if args.group == "all" or x["group"] == args.group
    ]
    out = ROOT / "verification/generated"
    out.mkdir(parents=True, exist_ok=True)
    passed = all(x["passed"] for x in results)
    (out / f"comparison-{args.group}.json").write_text(
        json.dumps({"passed": passed, "results": results}, indent=2) + "\n"
    )
    for result in results:
        print(("PASS" if result["passed"] else "FAIL") + ": " + result["file"])
        if not result["passed"]:
            print(result.get("error", result.get("errors")))
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
