"""Re-express frozen paired salt effects with the pore-entry resistance included.

This is algebraic postprocessing of saved membrane permeabilities. It does not
solve PB or integrate reservoir curves, and leaves the frozen source untouched.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

from mechanism_neutral_baseline import A_M, D_SALT, TAU


ROOT = Path(__file__).resolve().parents[1]
PAIRS = ROOT / "原始资料/赵珈绮/22_机制预测与实验验证方案_v1.0.0/matrices/paired_predictions.json"
SOURCE = ROOT / "data/processed/salt_transport/22_paired_effects_20261028.csv"
OUT = ROOT / "data/processed/salt_transport/23_成对效应膜阻口径勘误_20261030.csv"
FIELDS = [
    "pair_id", "h_nm", "c_feed_mM", "mechanism", "external",
    "P_mem_s_B_original_m_s", "P_mem_s_E_original_m_s",
    "P_mem_s_B_entry_included_m_s", "P_mem_s_E_entry_included_m_s",
    "E_s_original", "E_s_entry_included", "E_u_original",
    "Pi_original", "Pi_entry_included", "F_common_entry_included",
    "numeric_status", "correction_scope",
]


def main() -> None:
    with PAIRS.open(encoding="utf-8-sig") as stream:
        pairs = {row["pair_id"]: row for row in json.load(stream)}
    with SOURCE.open(encoding="utf-8-sig", newline="") as stream:
        original = list(csv.DictReader(stream))
    if len(pairs) != 54 or len(original) != 54 or {r["pair_id"] for r in original} != set(pairs):
        raise ValueError("Frozen pair identity mismatch")

    output = []
    changed = 0
    for old in original:
        pair = pairs[old["pair_id"]]
        for key in ("h_nm", "c_feed_mM"):
            if not math.isclose(float(old[key]), float(pair[key]), rel_tol=0, abs_tol=1e-12):
                raise ValueError(f"Pair {old['pair_id']} {key} mismatch")
        if any(old[key] != pair[key] for key in ("mechanism", "external")):
            raise ValueError(f"Pair {old['pair_id']} label mismatch")

        corrected = {}
        for state in ("B", "E"):
            p = float(old[f"P_mem_s_{state}_m_s"])
            if p <= 0:
                raise ValueError("Nonpositive saved permeability")
            if float(pair[f"sigma_{state}_C_m2"]) == 0:
                corrected[state] = p  # Neutral source already includes both pore entries.
            else:
                b = float(pair[f"b_{state}"])
                if b <= 0:
                    raise ValueError("Nonpositive b")
                two_entries = math.pi * A_M / (2 * b * TAU * D_SALT)
                corrected[state] = 1 / (1 / p + two_entries)
                changed += 1

        e_s = corrected["E"] / corrected["B"]
        e_u = float(old["E_u_initial"])
        pi = e_s / e_u
        if not math.isclose(float(old["E_s_initial"]),
                            float(old["P_mem_s_E_m_s"]) / float(old["P_mem_s_B_m_s"]),
                            rel_tol=1e-12):
            raise ValueError("Original E_s ratio mismatch")
        f_common = (math.log(e_u) / math.log(e_s)
                    if old["mechanism"] != "C" and e_s > 1 else "NA")
        output.append({
            "pair_id": old["pair_id"], "h_nm": old["h_nm"],
            "c_feed_mM": old["c_feed_mM"], "mechanism": old["mechanism"],
            "external": old["external"],
            "P_mem_s_B_original_m_s": old["P_mem_s_B_m_s"],
            "P_mem_s_E_original_m_s": old["P_mem_s_E_m_s"],
            "P_mem_s_B_entry_included_m_s": corrected["B"],
            "P_mem_s_E_entry_included_m_s": corrected["E"],
            "E_s_original": old["E_s_initial"], "E_s_entry_included": e_s,
            "E_u_original": old["E_u_initial"], "Pi_original": old["Pi_initial"],
            "Pi_entry_included": pi, "F_common_entry_included": f_common,
            "numeric_status": old["numeric_status"],
            "correction_scope": "INITIAL_MEMBRANE_PERMEABILITY_ONLY",
        })

    if changed != 36:
        raise ValueError(f"Expected 36 charged states, got {changed}")
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, FIELDS)
        writer.writeheader()
        writer.writerows(output)
    print(f"pairs={len(output)} charged_states={changed} output={OUT}")


if __name__ == "__main__":
    main()
