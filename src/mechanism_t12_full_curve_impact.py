"""Compare whole 10 mM charged trajectories using existing 31/16 PB nodes."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from mechanism_transport import PoreTable, run_curve


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/processed/salt_transport"
MATRIX = ROOT / "原始资料/赵珈绮/22_机制预测与实验验证方案_v1.0.0/matrices/paired_predictions.json"
OUT = DATA / "22_T12_10mM完整时序差_20261029.csv"
FIELDS = ("pair_id", "h_nm", "mechanism", "external", "time_min", "receiver_31_mM",
          "receiver_16_mM", "absolute_difference_uM", "relative_to_31_receiver",
          "mass_error_31", "mass_error_16", "numeric_status")


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    pairs = json.loads(MATRIX.read_text(encoding="utf-8-sig"))
    main_pb = rows(DATA / "22_pb_main_20261028.csv")
    extras = rows(DATA / "22_pb_target_refinement_20261028.csv")
    frozen_curves = rows(DATA / "22_paired_illustrations_20261028.csv")
    if (len(pairs), len(main_pb), len(frozen_curves)) != (54, 372, 864):
        raise ArithmeticError("Unexpected source size")
    family = sorted((float(r["c_mol_m3"]), float(r["D_s_local_m2_s"])) for r in main_pb
                    if float(r["a_nm"]) == 25 and float(r["sigma_C_m2"]) == -.008)
    endpoints = sorted((float(r["c_mol_m3"]), float(r["D_s_local_m2_s"])) for r in extras
                       if float(r["a_nm"]) == 25 and float(r["sigma_C_m2"]) == -.008)
    fine = PoreTable()
    if len(family) != 31 or [x for x, _ in endpoints] != [1.0, 10.0] or fine.nodes != sorted([(0.0, 0.0)] + family + endpoints):
        raise ArithmeticError("Fine PB nodes differ from original run")
    coarse = PoreTable()
    coarse.nodes = sorted([(0.0, 0.0)] + family[::2] + endpoints)
    if len(coarse.nodes) != 19:
        raise ArithmeticError("Unexpected coarse PB table")
    frozen = {(r["pair_id"], int(float(r["time_min"]))): r for r in frozen_curves
              if r["state"] == "B" and r["probe"] == "salt" and float(r["c_feed_mM"]) == 10
              and r["mechanism"] in ("C", "WC")}
    if len(frozen) != 48:
        raise ArithmeticError("Unexpected frozen charged states")

    output = []
    for pair in pairs:
        if pair["c_feed_mM"] != 10 or pair["mechanism"] not in ("C", "WC"):
            continue
        b = float(pair["b_B"])
        sigma = float(pair["sigma_B_C_m2"])
        if sigma != -.008:
            raise ArithmeticError("Unexpected charge")
        args = (10.0, "salt", float(pair["h_nm"]), b, sigma, pair["external"])
        fine_curve = run_curve(*args, fine)
        coarse_curve = run_curve(*args, coarse)
        for a, c in zip(fine_curve, coarse_curve):
            minute = int(a["time_min"])
            original = frozen[pair["pair_id"], minute]
            r31, r16 = a["receiver_mol_m3"], c["receiver_mol_m3"]
            if abs(r31 - float(original["receiver_mol_m3"])) > 1e-10 or abs(a["donor_mol_m3"] - float(original["donor_mol_m3"])) > 1e-10:
                raise ArithmeticError("Replayed fine trajectory differs from frozen source")
            output.append({
                "pair_id": pair["pair_id"], "h_nm": pair["h_nm"],
                "mechanism": pair["mechanism"], "external": pair["external"],
                "time_min": minute, "receiver_31_mM": r31, "receiver_16_mM": r16,
                "absolute_difference_uM": abs(r31 - r16) * 1000,
                "relative_to_31_receiver": abs(r31 - r16) / r31 if r31 else 0.0,
                "mass_error_31": a["mass_relative_error"],
                "mass_error_16": c["mass_relative_error"],
                "numeric_status": "NUMERIC_UNRESOLVED_T12",
            })
    if len(output) != 48 or len({r["pair_id"] for r in output}) != 12:
        raise ArithmeticError("Expected 12 charged pairs times four states")
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(output)
    nonzero = [r for r in output if r["time_min"]]
    print(f"pairs=12 states=48 max_gap_uM={max(r['absolute_difference_uM'] for r in nonzero):.8g} "
          f"max_relative={max(r['relative_to_31_receiver'] for r in nonzero):.8g} "
          f"over_2pct={sum(r['relative_to_31_receiver'] > .02 for r in nonzero)}")


if __name__ == "__main__":
    main()
