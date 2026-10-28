"""Run the 48 fixed radial PB grid refinements with counted attempts."""

from __future__ import annotations

import json
from datetime import datetime, timezone

from mechanism_pb_main_batch import FIELDS, LEDGER_FIELDS, OUT as MAIN_OUT, read_csv, write_csv
from mechanism_radial_pb import LEDGER, ROOT, solve


MATRIX = ROOT / "原始资料/赵珈绮/22_机制预测与实验验证方案_v1.0.0/matrices/pb_refinement.json"
OUT = ROOT / "data/processed/salt_transport/22_pb_refinement_20261028.csv"


def main() -> None:
    with MATRIX.open(encoding="utf-8-sig") as stream:
        matrix = json.load(stream)
    assert len(matrix) == 48
    main_rows = {r["id"]: r for r in read_csv(MAIN_OUT)}
    assert len(main_rows) == 372
    ledger = read_csv(LEDGER)
    results = read_csv(OUT)
    known = {r["id"] for r in results}
    for entry in matrix:
        parent = main_rows[entry["parent_id"]]
        assert (float(entry["a_nm"]) == float(parent["a_nm"])
                and float(entry["sigma_C_m2"]) == float(parent["sigma_C_m2"])
                and float(entry["c_mM"]) == float(parent["c_mM"]))
        point = dict(entry, c_mol_m3=float(entry["c_mM"]))  # 1 mM = 1 mol/m^3
        if point["id"] in known:
            continue
        if any(r["point_id"] == point["id"] and r["status"] == "PASS" for r in ledger):
            raise RuntimeError(f"Passing ledger entry without stored result: {point['id']}")
        if len(ledger) >= 444:
            raise RuntimeError("PB solve cap reached")
        record = {
            "call_number": len(ledger) + 1,
            "utc_time": datetime.now(timezone.utc).isoformat(),
            "point_id": point["id"], "cells": point["cells"], "status": "STARTED",
            "command": "python src/mechanism_pb_refinement.py", "exit_code": "", "reason": "",
        }
        ledger.append(record)
        write_csv(LEDGER, ledger, LEDGER_FIELDS)
        try:
            result = solve(point)
        except Exception as error:
            record["status"] = "FAILED"
            record["exit_code"] = "1"
            record["reason"] = str(error)
            write_csv(LEDGER, ledger, LEDGER_FIELDS)
            raise
        results.append(result)
        known.add(point["id"])
        write_csv(OUT, results, FIELDS)
        record["status"] = "PASS"
        record["exit_code"] = "0"
        write_csv(LEDGER, ledger, LEDGER_FIELDS)
    print(f"refinement={len(results)}/48 actual_calls={len(ledger)} "
          f"max_gauss={max(float(r['gauss_relative']) for r in results):.3g}")


if __name__ == "__main__":
    main()
