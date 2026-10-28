"""Propagate four new post-freeze PB points through the twelve charged 10 mM curves."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from mechanism_transport import PoreTable, run_curve


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/processed/salt_transport"
MATRIX = ROOT / "原始资料/赵珈绮/22_机制预测与实验验证方案_v1.0.0/matrices/paired_predictions.json"
OUT = DATA / "22_pb_postfreeze_10mM完整时序影响_20261029.csv"
FIELDS = ("pair_id", "h_nm", "mechanism", "external", "time_min",
          "frozen_receiver_mM", "augmented_receiver_mM", "signed_change_uM",
          "relative_change_to_frozen", "augmented_mass_relative_error",
          "numeric_status")


def read(name: str) -> list[dict[str, str]]:
    with (DATA / name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    pairs = json.loads(MATRIX.read_text(encoding="utf-8-sig"))
    spots, frozen_curves = read("22_pb_postfreeze_spotcheck_20261029.csv"), read("22_paired_illustrations_20261028.csv")
    post_ledger = read("22_pb_postfreeze_call_ledger_20261029.csv")
    if len(pairs) != 54 or len(spots) != 4 or len(post_ledger) != 4 or len(frozen_curves) != 864:
        raise ArithmeticError("Unexpected source size")
    if any(r["status"] != "PASS" for r in post_ledger):
        raise ArithmeticError("Unaccepted new PB call")
    table = PoreTable()
    table.nodes = sorted(table.nodes + [(float(r["c_mol_m3"]), float(r["D_s_local_m2_s"])) for r in spots])
    frozen = {(r["pair_id"], int(float(r["time_min"]))): r for r in frozen_curves
              if r["state"] == "B" and r["probe"] == "salt" and float(r["c_feed_mM"]) == 10
              and r["mechanism"] in ("C", "WC")}
    if len(frozen) != 48:
        raise ArithmeticError("Unexpected frozen charged curves")
    output = []
    for pair in pairs:
        if pair["c_feed_mM"] != 10 or pair["mechanism"] not in ("C", "WC"):
            continue
        if pair["sigma_B_C_m2"] != -.008 or pair["sigma_E_C_m2"] != 0:
            raise ArithmeticError("Unexpected charge-state assumptions")
        curve = run_curve(10.0, "salt", float(pair["h_nm"]), float(pair["b_B"]),
                          -.008, pair["external"], table)
        for item in curve:
            minute = int(item["time_min"])
            old = float(frozen[pair["pair_id"], minute]["receiver_mol_m3"])
            new = item["receiver_mol_m3"]
            output.append({
                "pair_id": pair["pair_id"], "h_nm": pair["h_nm"],
                "mechanism": pair["mechanism"], "external": pair["external"],
                "time_min": minute, "frozen_receiver_mM": old,
                "augmented_receiver_mM": new,
                "signed_change_uM": 1000 * (new - old),
                "relative_change_to_frozen": abs(new - old) / old if old else 0.0,
                "augmented_mass_relative_error": item["mass_relative_error"],
                "numeric_status": "POSTFREEZE_SPOTCHECK_T12_UNRESOLVED",
            })
    if len(output) != 48 or len({r["pair_id"] for r in output}) != 12:
        raise ArithmeticError("Expected 12 pairs times four states")
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(output)
    later = [r for r in output if r["time_min"]]
    print(f"pairs=12 states=48 max_absolute_change_uM={max(abs(r['signed_change_uM']) for r in later):.8g} "
          f"max_relative={max(r['relative_change_to_frozen'] for r in later):.8g} "
          f"over_2pct={sum(r['relative_change_to_frozen'] > .02 for r in later)}")


if __name__ == "__main__":
    main()
