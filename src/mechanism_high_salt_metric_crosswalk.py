"""Reconcile fixed high-salt C-pair membrane gain with system initial flux."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/processed/salt_transport"
MATRIX = ROOT / "原始资料/赵珈绮/22_机制预测与实验验证方案_v1.0.0/matrices/paired_predictions.json"
EFFECTS = DATA / "22_paired_effects_20261028.csv"
FLUX = DATA / "22_高盐初始通量荷电中性对照_20261029.csv"
OUT = DATA / "22_高盐C情景膜内与装置通量口径对账_20261029.csv"
FIELDS = ("pair_id", "h_nm", "external", "sigma_B_mC_m2", "b_B", "b_E",
          "frozen_E_s_membrane_E_over_B", "membrane_B_over_E_ratio",
          "system_B_over_E_initial_flux_ratio", "membrane_excess_percent",
          "system_excess_percent", "numeric_status", "interpretation")


def read(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    pairs = json.loads(MATRIX.read_text(encoding="utf-8-sig"))
    effects = {r["pair_id"]: r for r in read(EFFECTS)}
    flux = {(int(r["h_nm"]), int(r["sigma_mC_m2"]), r["external"]): r for r in read(FLUX)}
    if len(pairs) != 54 or len(effects) != 54 or len(flux) != 18:
        raise ArithmeticError("Fixed pair or flux source changed")
    rows = []
    for pair in pairs:
        if pair["mechanism"] != "C" or int(pair["c_feed_mM"]) != 600:
            continue
        if (float(pair["b_B"]), float(pair["b_E"]), float(pair["sigma_B_C_m2"]),
                float(pair["sigma_E_C_m2"])) != (0.01, 0.01, -0.008, 0.0):
            raise ArithmeticError("C-pair baseline changed")
        effect = effects[pair["pair_id"]]
        source = flux[int(pair["h_nm"]), -8, pair["external"]]
        membrane = 1 / float(effect["E_s_initial"])
        system = float(source["charged_to_neutral_J_ratio"])
        if membrane <= 1 or system <= 1:
            raise ArithmeticError("Unexpected high-salt direction")
        rows.append({"pair_id": pair["pair_id"], "h_nm": pair["h_nm"],
                     "external": pair["external"], "sigma_B_mC_m2": -8,
                     "b_B": pair["b_B"], "b_E": pair["b_E"],
                     "frozen_E_s_membrane_E_over_B": effect["E_s_initial"],
                     "membrane_B_over_E_ratio": membrane,
                     "system_B_over_E_initial_flux_ratio": system,
                     "membrane_excess_percent": 100 * (membrane - 1),
                     "system_excess_percent": 100 * (system - 1),
                     "numeric_status": effect["numeric_status"],
                     "interpretation": "DIFFERENT_METRICS_DO_NOT_SUBSTITUTE"})
    if len(rows) != 6:
        raise ArithmeticError("Expected six high-salt C pairs")
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"pairs={len(rows)} membrane_B_over_E={min(r['membrane_B_over_E_ratio'] for r in rows):.8g}.."
          f"{max(r['membrane_B_over_E_ratio'] for r in rows):.8g} "
          f"system_B_over_E={min(r['system_B_over_E_initial_flux_ratio'] for r in rows):.8g}.."
          f"{max(r['system_B_over_E_initial_flux_ratio'] for r in rows):.8g}")


if __name__ == "__main__":
    main()
