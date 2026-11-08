"""Run one fixed 31-point PB family, keeping every actual solve in the ledger."""

from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path

from mechanism_radial_pb import LEDGER, MATRIX, OUT as FIRST_RESULT, ROOT, solve


OUT = ROOT / "data/processed/salt_transport/22_pb_main_20261028.csv"
FIELDS = ["id", "a_nm", "sigma_C_m2", "c_mM", "c_mol_m3", "cells", "psi_axis",
          "psi_wall_cell", "k_plus", "k_minus", "D_s_local_m2_s",
          "max_fv_residual", "gauss_relative", "newton_iterations", "status", "note"]
LEDGER_FIELDS = ["call_number", "utc_time", "point_id", "cells", "status",
                 "command", "exit_code", "reason"]


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("family_index", type=int, choices=range(12))
    args = parser.parse_args()
    with MATRIX.open(encoding="utf-8-sig") as stream:
        matrix = json.load(stream)
    assert len(matrix) == 372
    family = matrix[31 * args.family_index:31 * (args.family_index + 1)]
    assert len(family) == 31 and len({(p["a_nm"], p["sigma_C_m2"]) for p in family}) == 1
    ledger = read_csv(LEDGER)
    results = read_csv(OUT)
    known = {r["id"] for r in results}
    if args.family_index == 0 and family[0]["id"] not in known:
        if not any(r["point_id"] == family[0]["id"] and r["status"] == "PASS" for r in ledger):
            raise RuntimeError("First point must already have a passing ledger entry")
        with FIRST_RESULT.open(encoding="utf-8") as stream:
            first = json.load(stream)
        results.append(first)
        known.add(first["id"])
        write_csv(OUT, results, FIELDS)
    for point in family:
        point_id = point["id"]
        if point_id in known:
            continue
        if any(r["point_id"] == point_id and r["status"] == "PASS" for r in ledger):
            raise RuntimeError(f"Passing ledger entry without stored result: {point_id}")
        if len(ledger) >= 444:
            raise RuntimeError("PB solve cap reached")
        entry = {
            "call_number": len(ledger) + 1,
            "utc_time": datetime.now(timezone.utc).isoformat(),
            "point_id": point_id, "cells": point["cells"], "status": "STARTED",
            "command": f"python src/mechanism_pb_main_batch.py {args.family_index}",
            "exit_code": "", "reason": "",
        }
        ledger.append(entry)
        write_csv(LEDGER, ledger, LEDGER_FIELDS)
        try:
            result = solve(point)
        except Exception as error:
            entry["status"] = "FAILED"
            entry["exit_code"] = "1"
            entry["reason"] = str(error)
            write_csv(LEDGER, ledger, LEDGER_FIELDS)
            raise
        results.append(result)
        known.add(point_id)
        write_csv(OUT, results, FIELDS)
        entry["status"] = "PASS"
        entry["exit_code"] = "0"
        write_csv(LEDGER, ledger, LEDGER_FIELDS)
    completed = [r for r in results if r["id"] in {p["id"] for p in family}]
    print(f"family={args.family_index} completed={len(completed)}/31 actual_calls={len(ledger)} "
          f"max_gauss={max(float(r['gauss_relative']) for r in completed):.3g}")


if __name__ == "__main__":
    main()
