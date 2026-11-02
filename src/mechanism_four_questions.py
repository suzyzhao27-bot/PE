"""Extract numerical answers to the four fixed pre-experiment questions."""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/processed/salt_transport"
BOOK2 = ROOT / "原始资料/赵珈绮/20_学生Stage2_B1成果_v1.0.0/WP1/book2_long.csv"
OUT = DATA / "22_four_science_questions_20261028.json"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def one(items: list[dict[str, str]], **conditions: object) -> dict[str, str]:
    match = [r for r in items if all((float(r[k]) == float(v) if isinstance(v, (int, float))
                                    else r[k] == v) for k, v in conditions.items())]
    if len(match) != 1:
        raise ArithmeticError(f"Expected one row for {conditions}; found {len(match)}")
    return match[0]


def numerical_range(items: list[dict[str, str]], field: str) -> list[float]:
    values = [float(r[field]) for r in items]
    return [min(values), max(values)]


def main() -> None:
    neutral = rows(DATA / "22_neutral_baseline_20261028.csv")
    effects = rows(DATA / "22_paired_effects_20261028.csv")
    curves = rows(DATA / "22_paired_illustrations_20261028.csv")
    detection = rows(DATA / "22_G4_模型检测需求_20261028.csv")
    legacy = rows(DATA / "22_legacy_urea_conditions_20261028.csv")
    book2 = rows(BOOK2)
    if (len(neutral), len(effects), len(curves), len(detection), len(legacy), len(book2)) != (
            288, 54, 864, 108, 164, 41):
        raise ArithmeticError("Unexpected source length")

    fixed = {}
    for probe in ("salt_zero_charge", "urea"):
        for external in ("zero", "reference_2600_neutral"):
            p = {}
            for h in (65, 700, 2600):
                row = one(neutral, probe=probe, external=external, h_nm=float(h), b=.01,
                          salt_feed_mM=1.0, time_min=0)
                p[str(h)] = float(row["permeability_total_m_s"])
            fixed[f"{probe}|{external}"] = {"P_total_m_s_by_h_nm": p,
                                             "P_2600_over_P_65": p["2600"] / p["65"]}
    old = {}
    for h in (65, 700, 2600):
        increments = {}
        for state in ("B", "R"):
            sample = f"{state}{h}"
            initial = one(book2, sample_id=sample, time_min=0.0)
            at100 = one(book2, sample_id=sample, time_min=100.0)
            increments[state] = float(at100["kappa_uScm"]) - float(initial["kappa_uScm"])
        old[str(h)] = {"increment_100_uScm": increments,
                       "R_over_B_increment_100": increments["R"] / increments["B"]}

    charge = {}
    for feed in (1, 10, 600):
        sub = [r for r in effects if r["mechanism"] == "C" and float(r["c_feed_mM"]) == feed]
        charge[str(feed)] = {"E_s_initial_range": numerical_range(sub, "E_s_initial"),
                            "E_u_initial_range": numerical_range(sub, "E_u_initial"),
                            "Pi_initial_range": numerical_range(sub, "Pi_initial"),
                            "numeric_statuses": sorted({r["numeric_status"] for r in sub})}

    urea = {}
    for h in (65, 700, 2600):
        variants = []
        for scenario in ("C1_T85", "C2_T85"):
            row = one(legacy, scenario=scenario, sample_id=f"R{h}", source_time_min=100.0)
            variants.append(row)
        urea[str(h)] = {
            "R100_W_receiver_mM_C1_C2": [float(r["urea_receiver_W_mol_m3"]) for r in variants],
            "B100_C_anchor_receiver_mM_C1_C2": [float(r["urea_receiver_C_anchor_mol_m3"]) for r in variants],
            "condition": variants[0]["status"],
        }

    signatures = {}
    for mechanism in ("W", "C", "WC"):
        sub = [r for r in effects if r["mechanism"] == mechanism and float(r["h_nm"]) == 700
               and float(r["c_feed_mM"]) == 1 and r["external"] == "zero"]
        if len(sub) != 1:
            raise ArithmeticError(mechanism)
        pid = sub[0]["pair_id"]
        sig = {"pair_id": pid, "initial_E_s": float(sub[0]["E_s_initial"]),
               "initial_E_u": float(sub[0]["E_u_initial"]),
               "initial_Pi": float(sub[0]["Pi_initial"])}
        for probe in ("salt", "urea"):
            for state in ("B", "E"):
                row = one(curves, pair_id=pid, probe=probe, state=state, time_min=60.0)
                sig[f"{probe}_{state}_receiver_60_mM"] = float(row["receiver_mol_m3"])
                det = one(detection, pair_id=pid, state=state)
                sig[f"{probe}_{state}_LOQ_ceiling_60_uM"] = float(
                    det[f"{probe}_LOQ_at_most_uM"])
        signatures[mechanism] = sig

    result = {"source_status": "MODEL_ONLY_AND_CONDITIONAL_OLD_DATA",
              "question_1_fixed_structure": {"neutral_P": fixed, "old_Book2_100min": old,
                                             "comparison_boundary": "Old conductivity increments are not membrane permeability or paired treatment effects."},
              "question_2_charge_only": charge,
              "question_3_legacy_urea_targets": urea,
              "question_4_paired_signatures_700nm_1mM_zero_external": signatures,
              "boundaries": ["W/C/WC use prescribed b and sigma, not fitted material parameters.",
                             "10 mM charged cases remain T12 numeric unresolved.",
                             "Actual LOQ, flow and shared geometry/external scaling are unmeasured."]}
    with OUT.open("w", encoding="utf-8") as stream:
        json.dump(result, stream, ensure_ascii=False, indent=2)
    print(f"old_h={len(old)} charge_feeds={len(charge)} urea_h={len(urea)} "
          f"signatures={len(signatures)}")


if __name__ == "__main__":
    main()
