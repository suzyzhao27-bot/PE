"""Calculate frozen initial-state membrane E_s, E_u, and Pi for 54 pairs."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

from mechanism_neutral_baseline import A_M, D_SALT, TAU, membrane_permeability
from mechanism_transport import PoreTable, external_resistance, salt_flux


ROOT = Path(__file__).resolve().parents[1]
MATRIX = ROOT / "原始资料/赵珈绮/22_机制预测与实验验证方案_v1.0.0/matrices/paired_predictions.json"
OUT = ROOT / "data/processed/salt_transport/22_paired_effects_20261028.csv"
FIELDS = ["pair_id", "h_nm", "c_feed_mM", "mechanism", "external",
          "P_mem_s_B_m_s", "P_mem_s_E_m_s", "P_mem_u_B_m_s", "P_mem_u_E_m_s",
          "E_s_initial", "E_u_initial", "Pi_initial", "F_common_initial", "numeric_status"]


def salt_membrane_p(feed: float, h_nm: float, b: float, sigma: float,
                    external: str, table: PoreTable) -> float:
    if sigma == 0:
        return membrane_permeability(h_nm, b, "salt_zero_charge")
    flux = salt_flux(feed, 0.0, h_nm, b, sigma, external, table)
    epsilon_f = b * TAU
    side = external_resistance("salt", external) / 2 + math.pi * A_M / (4 * epsilon_f * D_SALT)
    interface_difference = feed - 2 * flux * side
    if interface_difference <= 0:
        raise ArithmeticError("Nonpositive membrane concentration difference")
    return flux / interface_difference


def main() -> None:
    with MATRIX.open(encoding="utf-8-sig") as stream:
        pairs = json.load(stream)
    table = PoreTable()
    rows = []
    for pair in pairs:
        h = float(pair["h_nm"])
        feed = float(pair["c_feed_mM"])
        external = pair["external"]
        p_s = {}
        p_u = {}
        for state in ("B", "E"):
            b = float(pair[f"b_{state}"])
            sigma = float(pair[f"sigma_{state}_C_m2"])
            p_s[state] = salt_membrane_p(feed, h, b, sigma, external, table)
            p_u[state] = membrane_permeability(h, b, "urea")
        e_s = p_s["E"] / p_s["B"]
        e_u = p_u["E"] / p_u["B"]
        pi = e_s / e_u
        # W/WC have a large planned salt gain. C is a charge-only diagnostic,
        # including cases with E_s near or below one; leave this index undefined.
        f_common = math.log(e_u) / math.log(e_s) if pair["mechanism"] != "C" and e_s > 1 else "NA"
        rows.append({
            "pair_id": pair["pair_id"], "h_nm": h, "c_feed_mM": feed,
            "mechanism": pair["mechanism"], "external": external,
            "P_mem_s_B_m_s": p_s["B"], "P_mem_s_E_m_s": p_s["E"],
            "P_mem_u_B_m_s": p_u["B"], "P_mem_u_E_m_s": p_u["E"],
            "E_s_initial": e_s, "E_u_initial": e_u, "Pi_initial": pi,
            "F_common_initial": f_common,
            "numeric_status": ("NUMERIC_UNRESOLVED_T12" if feed == 10 and
                               float(pair["sigma_B_C_m2"]) != 0
                               else "MODEL_ONLY_TARGET_CHECK"),
        })
    if len(rows) != 54:
        raise ArithmeticError("Missing fixed pairs")
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"effects={len(rows)} Pi_min={min(r['Pi_initial'] for r in rows):.6g} "
          f"Pi_max={max(r['Pi_initial'] for r in rows):.6g}")


if __name__ == "__main__":
    main()
