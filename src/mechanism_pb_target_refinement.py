"""Use the fixed technical reserve on exact 1/10 mM integration endpoints."""

from __future__ import annotations

from datetime import datetime, timezone

from mechanism_pb_main_batch import FIELDS, LEDGER_FIELDS, read_csv, write_csv
from mechanism_pb_integral_check import OUT as CHECK_OUT
from mechanism_pb_refinement import OUT as REFINEMENT_OUT
from mechanism_radial_pb import LEDGER, ROOT, solve


OUT = ROOT / "data/processed/salt_transport/22_pb_target_refinement_20261028.csv"


def main() -> None:
    risk = [r for r in read_csv(CHECK_OUT) if r["status"] == "NUMERIC_RISK_OVER_2_PERCENT"]
    assert len(risk) == 12 and {float(r["feed_mM"]) for r in risk} == {1.0, 10.0}
    assert len(read_csv(REFINEMENT_OUT)) == 48
    ledger = read_csv(LEDGER)
    results = read_csv(OUT)
    known = {r["id"] for r in results}
    for item in risk:
        a = int(float(item["a_nm"]))
        sigma = float(item["sigma_C_m2"])
        feed = int(float(item["feed_mM"]))
        point = {"id": f"PB_extra_a{a}_s{int(-sigma * 1000):02d}_feed{feed}",
                 "a_nm": a, "sigma_C_m2": sigma, "c_mM": feed,
                 "c_mol_m3": float(feed), "cells": 96}
        if point["id"] in known:
            continue
        if any(r["point_id"] == point["id"] and r["status"] == "PASS" for r in ledger):
            raise RuntimeError(f"Passing ledger entry without stored result: {point['id']}")
        if len(ledger) >= 444 or len(ledger) - 420 >= 24:
            raise RuntimeError("PB technical reserve or total cap reached")
        record = {"call_number": len(ledger) + 1,
                  "utc_time": datetime.now(timezone.utc).isoformat(),
                  "point_id": point["id"], "cells": 96, "status": "STARTED",
                  "command": "python src/mechanism_pb_target_refinement.py",
                  "exit_code": "", "reason": f"T12 initial 31/16 difference {item['relative_difference']}"}
        ledger.append(record)
        write_csv(LEDGER, ledger, LEDGER_FIELDS)
        try:
            result = solve(point)
        except Exception as error:
            record["status"] = "FAILED"
            record["exit_code"] = "1"
            record["reason"] += f"; {error}"
            write_csv(LEDGER, ledger, LEDGER_FIELDS)
            raise
        results.append(result)
        known.add(point["id"])
        write_csv(OUT, results, FIELDS)
        record["status"] = "PASS"
        record["exit_code"] = "0"
        write_csv(LEDGER, ledger, LEDGER_FIELDS)
    print(f"supplement={len(results)}/12 actual_calls={len(ledger)} reserve_remaining={444-len(ledger)}")


if __name__ == "__main__":
    main()
