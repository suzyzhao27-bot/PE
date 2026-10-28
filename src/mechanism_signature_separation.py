"""Compare W/C/WC treatment responses in the already frozen 54 pairs."""

from __future__ import annotations

import csv
import json
import math
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FREEZE = ROOT / "data/processed/salt_transport/freeze_predictions.json"
OUT = ROOT / "data/processed/salt_transport/22_W_C_WC_处理响应差距_20261028.csv"
FIELDS = ("h_nm", "salt_feed_mM", "external", "time_min", "mechanism_a", "mechanism_b",
          "salt_response_a_uM", "salt_response_b_uM", "salt_gap_uM",
          "urea_response_a_uM", "urea_response_b_uM", "urea_gap_uM",
          "model_difference_type", "numeric_status", "experimental_identifiability")


def states(group: dict[str, object]) -> dict[tuple[str, str, int], float]:
    return {(r["state"], r["probe"], r["time_min"]): r["receiver_mM"]
            for r in group["pool_states"]}


def main() -> None:
    frozen = json.loads(FREEZE.read_text(encoding="utf-8"))
    if frozen["status"] != "FROZEN_CODEX_TECHNICAL_STUDENT_UNREVIEWED":
        raise ArithmeticError("Unexpected freeze status")
    groups = {}
    for pair in frozen["prediction_groups"]:
        key = (pair["h_nm"], pair["salt_feed_mM"], pair["external"])
        groups.setdefault(key, {})[pair["mechanism"]] = pair
    if len(groups) != 18 or any(set(value) != {"W", "C", "WC"} for value in groups.values()):
        raise ArithmeticError("Incomplete W/C/WC strata")
    output = []
    for (h, feed, external), mechanisms in sorted(groups.items()):
        pool = {name: states(pair) for name, pair in mechanisms.items()}
        for minute in (60, 120, 240):
            for probe in ("salt", "urea"):
                e_values = [pool[name][("E", probe, minute)] for name in ("W", "C", "WC")]
                if max(e_values) - min(e_values) > 1e-12:
                    raise ArithmeticError(f"E state differs across mechanisms: {h} {feed} {minute}")
            response = {}
            for name in ("W", "C", "WC"):
                response[name] = {probe: 1000 * (pool[name][("E", probe, minute)]
                                                  - pool[name][("B", probe, minute)])
                                  for probe in ("salt", "urea")}
            for a, b in combinations(("W", "C", "WC"), 2):
                s_gap = abs(response[a]["salt"] - response[b]["salt"])
                u_gap = abs(response[a]["urea"] - response[b]["urea"])
                salt_diff = not math.isclose(s_gap, 0.0, abs_tol=1e-9)
                urea_diff = not math.isclose(u_gap, 0.0, abs_tol=1e-9)
                kind = ("BOTH_PROBES_DIFFER" if salt_diff and urea_diff else
                        "SALT_ONLY_DIFFERENCE" if salt_diff else
                        "UREA_ONLY_DIFFERENCE" if urea_diff else "MODEL_RESPONSES_EQUAL")
                status = ("NUMERIC_UNRESOLVED_T12" if feed == 10 and ("C" in (a, b) or "WC" in (a, b))
                          else "MODEL_ONLY_TARGET_CHECK")
                output.append({"h_nm": h, "salt_feed_mM": feed, "external": external,
                               "time_min": minute, "mechanism_a": a, "mechanism_b": b,
                               "salt_response_a_uM": response[a]["salt"],
                               "salt_response_b_uM": response[b]["salt"], "salt_gap_uM": s_gap,
                               "urea_response_a_uM": response[a]["urea"],
                               "urea_response_b_uM": response[b]["urea"], "urea_gap_uM": u_gap,
                               "model_difference_type": kind, "numeric_status": status,
                               "experimental_identifiability": "UNDETERMINED_ACTUAL_LOQ_NOISE_FLOW_NOT_MEASURED"})
    if len(output) != 162:
        raise ArithmeticError(f"Expected 162 pairwise comparisons, got {len(output)}")
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(output)
    print(f"comparisons={len(output)} strata={len(groups)} "
          f"numeric_unresolved={sum(r['numeric_status'] == 'NUMERIC_UNRESOLVED_T12' for r in output)}")


if __name__ == "__main__":
    main()
