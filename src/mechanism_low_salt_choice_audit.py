"""Summarize fixed low-salt alternatives without choosing one before calibration."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/processed/salt_transport"
MODEL = DATA / "22_G4_模型检测需求_20261028.csv"
INTERFACE = DATA / "22_模型与实验条件接口_20261029.csv"
OUT = DATA / "22_低盐选择模型检测接口_20261029.csv"
FIELDS = ("low_salt_mM", "linked_core_model_conditions", "linked_comparisons",
          "illustrated_B_E_states", "min_salt_receiver_60_uM", "min_salt_pair_id",
          "illustrative_salt_LOQ_at_most_uM_for_all_states", "min_urea_receiver_60_uM",
          "illustrative_urea_LOQ_at_most_uM_for_all_states",
          "T12_unresolved_states", "T12_unresolved_linked_comparisons",
          "actual_salt_LOQ", "decision_status")


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    model, interface = rows(MODEL), rows(INTERFACE)
    if len(model) != 108 or len(interface) != 18:
        raise ArithmeticError("Unexpected source size")
    output = []
    for feed in (1, 10):
        flag = f"link_if_low_{feed}mM"
        links = [r for r in interface if r[flag] == "1" and int(r["salt_feed_mM"]) == feed]
        cases = [r for r in model if int(r["salt_feed_mM"]) == feed
                 and int(r["h_nm"]) in (65, 2600)]
        if len(links) != 4 or len(cases) != 24:
            raise ArithmeticError("Unexpected linked low-salt matrix")
        min_salt = min(cases, key=lambda r: float(r["salt_receiver_60_uM"]))
        min_urea = min(float(r["urea_receiver_60_uM"]) for r in cases)
        unresolved = sum(r["numeric_status"] == "NUMERIC_UNRESOLVED_T12" for r in cases)
        linked_comparisons = sum(int(r["model_comparison_rows"]) for r in links)
        if (feed == 1 and (unresolved, linked_comparisons) != (0, 36)) or (feed == 10 and (unresolved, linked_comparisons) != (8, 36)):
            raise ArithmeticError("Unexpected low-salt numeric risk")
        output.append({
            "low_salt_mM": feed,
            "linked_core_model_conditions": len(links),
            "linked_comparisons": linked_comparisons,
            "illustrated_B_E_states": len(cases),
            "min_salt_receiver_60_uM": min_salt["salt_receiver_60_uM"],
            "min_salt_pair_id": min_salt["pair_id"],
            "illustrative_salt_LOQ_at_most_uM_for_all_states": float(min_salt["salt_receiver_60_uM"]) / 5,
            "min_urea_receiver_60_uM": min_urea,
            "illustrative_urea_LOQ_at_most_uM_for_all_states": min_urea / 5,
            "T12_unresolved_states": unresolved,
            "T12_unresolved_linked_comparisons": 0 if feed == 1 else linked_comparisons,
            "actual_salt_LOQ": "NOT_MEASURED",
            "decision_status": "WAIT_FOR_GLOBAL_PRETEST_LOQ_DECISION",
        })
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(output)
    for row in output:
        print(f"feed={row['low_salt_mM']}mM min_salt_60_uM={float(row['min_salt_receiver_60_uM']):.8g} "
              f"illustrative_LOQ_ceiling_uM={row['illustrative_salt_LOQ_at_most_uM_for_all_states']:.8g} "
              f"unresolved_states={row['T12_unresolved_states']}")


if __name__ == "__main__":
    main()
