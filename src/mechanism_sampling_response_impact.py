"""Check B/E treatment responses against omitted-sampling counterfactuals."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/processed/salt_transport"
MATRIX = ROOT / "原始资料/赵珈绮/22_机制预测与实验验证方案_v1.0.0/matrices/paired_predictions.json"
INPUTS = (DATA / "22_中性曲线取样省略影响_20261029.csv",
          DATA / "22_带电盐曲线取样省略影响_20261029.csv")
OUT = DATA / "22_处理响应对取样事件敏感性_20261029.csv"
FIELDS = ("pair_id", "h_nm", "feed_mM", "mechanism", "external", "probe", "time_min",
          "sampled_E_minus_B_uM", "no_withdrawal_E_minus_B_uM",
          "response_difference_uM", "direction_status", "numeric_status", "interpretation")


def main() -> None:
    pairs = json.loads(MATRIX.read_text(encoding="utf-8-sig"))
    source = []
    for path in INPUTS:
        with path.open(encoding="utf-8-sig", newline="") as stream:
            source.extend(csv.DictReader(stream))
    if len(pairs) != 54 or len(source) != 648:
        raise ArithmeticError("Unexpected fixed response inputs")
    indexed = {(r["pair_id"], r["state"], r["probe"], int(r["time_min"])): r for r in source}
    if len(indexed) != 648:
        raise ArithmeticError("Duplicate sampled/no-withdrawal state")
    rows = []
    for pair in pairs:
        for probe in ("salt", "urea"):
            for minute in (60, 120, 240):
                b, e = (indexed[pair["pair_id"], state, probe, minute] for state in ("B", "E"))
                sampled = float(e["frozen_sampled_receiver_uM"]) - float(b["frozen_sampled_receiver_uM"])
                no_withdrawal = float(e["no_withdrawal_receiver_uM"]) - float(b["no_withdrawal_receiver_uM"])
                # 1e-9 uM only suppresses floating-point sign noise near exact zero.
                significant = abs(sampled) > 1e-9 and abs(no_withdrawal) > 1e-9
                direction = ("REVERSED" if significant and sampled * no_withdrawal < 0 else
                             "SAME_DIRECTION" if significant else "ZERO_OR_NUMERIC_NOISE")
                risk = probe == "salt" and int(pair["c_feed_mM"]) == 10 and pair["mechanism"] in ("C", "WC")
                rows.append({"pair_id": pair["pair_id"], "h_nm": pair["h_nm"],
                             "feed_mM": pair["c_feed_mM"], "mechanism": pair["mechanism"],
                             "external": pair["external"], "probe": probe, "time_min": minute,
                             "sampled_E_minus_B_uM": sampled,
                             "no_withdrawal_E_minus_B_uM": no_withdrawal,
                             "response_difference_uM": sampled - no_withdrawal,
                             "direction_status": direction,
                             "numeric_status": "NUMERIC_UNRESOLVED_T12" if risk else "MODEL_ONLY_TARGET_CHECK",
                             "interpretation": "NO_WITHDRAWAL_COUNTERFACTUAL_NOT_PROTOCOL"})
    if len(rows) != 324:
        raise ArithmeticError("Expected 54 pairs x two probes x three times")
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    later = [r for r in rows if r["time_min"] > 60]
    print(f"responses={len(rows)} max_abs_response_difference_uM="
          f"{max(abs(r['response_difference_uM']) for r in later):.8g} "
          f"direction_reversals={sum(r['direction_status'] == 'REVERSED' for r in rows)} "
          f"T12_rows={sum(r['numeric_status'] == 'NUMERIC_UNRESOLVED_T12' for r in rows)}")


if __name__ == "__main__":
    main()
