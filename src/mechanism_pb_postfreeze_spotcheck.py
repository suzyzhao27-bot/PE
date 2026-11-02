"""Spend four fixed reserve calls on new T12 concentration midpoints after freeze."""

from __future__ import annotations

import math
from datetime import datetime, timezone

from mechanism_pb_main_batch import FIELDS, LEDGER_FIELDS, OUT as MAIN, read_csv, write_csv
from mechanism_radial_pb import LEDGER as FROZEN_LEDGER, ROOT, solve


OUT = ROOT / "data/processed/salt_transport/22_pb_postfreeze_spotcheck_20261029.csv"
POST_LEDGER = ROOT / "data/processed/salt_transport/22_pb_postfreeze_call_ledger_20261029.csv"
COMMAND = "python src/mechanism_pb_postfreeze_spotcheck.py"


def main() -> None:
    main_rows = read_csv(MAIN)
    frozen_ledger = read_csv(FROZEN_LEDGER)
    post_ledger = read_csv(POST_LEDGER)
    results = read_csv(OUT)
    if len(main_rows) != 372 or len(frozen_ledger) != 437 or len(frozen_ledger) + len(post_ledger) > 444:
        raise ArithmeticError("Unexpected main PB table or call budget")
    anchors = sorted(float(r["c_mol_m3"]) for r in main_rows
                     if float(r["a_nm"]) == 25 and float(r["sigma_C_m2"]) == -.008
                     and 1 < float(r["c_mol_m3"]) < 10)
    if len(anchors) != 5 or not (1.16 < anchors[0] < 1.17 and 9.34 < anchors[-1] < 9.35):
        raise ArithmeticError("Concentration anchors changed")
    # The two dominant 16-node intervals span anchors[0:3] and anchors[2:5].
    midpoints = [(index + 1, math.sqrt(left * right))
                 for index, (left, right) in enumerate(zip(anchors, anchors[1:]))]
    known_main = {float(r["c_mol_m3"]) for r in main_rows
                  if float(r["a_nm"]) == 25 and float(r["sigma_C_m2"]) == -.008}
    if len(midpoints) != 4 or any(c in known_main for _, c in midpoints):
        raise ArithmeticError("Spotcheck nodes are not new")
    planned_ids = {f"PB_postfreeze_a25_s08_q{index}" for index, _ in midpoints}
    if any(r["id"] not in planned_ids for r in results):
        raise ArithmeticError("Unexpected post-freeze result")
    for index, c in midpoints:
        point_id = f"PB_postfreeze_a25_s08_q{index}"
        if any(r["id"] == point_id for r in results):
            continue
        if any(r["point_id"] == point_id for r in post_ledger):
            raise RuntimeError(f"Prior unaccepted call exists: {point_id}")
        if len(frozen_ledger) + len(post_ledger) >= 444:
            raise RuntimeError("PB solve cap reached")
        point = {"id": point_id, "a_nm": 25, "sigma_C_m2": -.008,
                 "c_mM": c, "c_mol_m3": c, "cells": 96}
        record = {"call_number": len(frozen_ledger) + len(post_ledger) + 1,
                  "utc_time": datetime.now(timezone.utc).isoformat(),
                  "point_id": point_id, "cells": 96, "status": "STARTED",
                  "command": COMMAND, "exit_code": "",
                  "reason": "Post-freeze T12 independent midpoint in dominant 31/16 interval"}
        post_ledger.append(record)
        write_csv(POST_LEDGER, post_ledger, LEDGER_FIELDS)
        try:
            result = solve(point)
        except Exception as error:
            record["status"] = "FAILED"
            record["exit_code"] = "1"
            record["reason"] += f"; {error}"
            write_csv(POST_LEDGER, post_ledger, LEDGER_FIELDS)
            raise
        result["note"] = "Post-freeze T12 interpolation spotcheck; no material validation"
        results.append(result)
        write_csv(OUT, results, FIELDS)
        record["status"] = "PASS"
        record["exit_code"] = "0"
        write_csv(POST_LEDGER, post_ledger, LEDGER_FIELDS)
    if len(results) != 4:
        raise ArithmeticError("Incomplete post-freeze spotcheck")
    total_calls = len(frozen_ledger) + len(post_ledger)
    print(f"new_points={len(results)} actual_calls={total_calls} "
          f"remaining={444-total_calls} max_gauss={max(float(r['gauss_relative']) for r in results):.6g}")


if __name__ == "__main__":
    main()
