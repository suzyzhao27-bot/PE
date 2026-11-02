"""Measure actual new-PB-point effect on the frozen 25 nm charged integral."""

from __future__ import annotations

import csv
from pathlib import Path

from mechanism_transport import PoreTable


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/processed/salt_transport"
SPOTS = DATA / "22_pb_postfreeze_spotcheck_20261029.csv"
FROZEN_LEDGER = DATA / "22_pb_call_ledger_20261028.csv"
POST_LEDGER = DATA / "22_pb_postfreeze_call_ledger_20261029.csv"
OUT = DATA / "22_pb_postfreeze_10mM积分独立补点差_20261029.csv"
FIELDS = ("lower_mM", "new_node_mM", "upper_mM", "old_integral_mol_m_s",
          "augmented_integral_mol_m_s", "signed_change_mol_m_s",
          "total_old_0_10_mol_m_s", "total_augmented_0_10_mol_m_s",
          "total_relative_change", "numeric_status")


def read(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    spots, frozen_ledger, post_ledger = read(SPOTS), read(FROZEN_LEDGER), read(POST_LEDGER)
    if len(spots) != 4 or len(frozen_ledger) != 437 or len(post_ledger) != 4:
        raise ArithmeticError("Unexpected post-freeze PB evidence")
    if any(r["status"] != "PASS" for r in post_ledger) or {r["point_id"] for r in post_ledger} != {r["id"] for r in spots} or [int(r["call_number"]) for r in post_ledger] != list(range(438, 442)):
        raise ArithmeticError("New PB rows lack passing call records")
    original = PoreTable()
    refined = PoreTable()
    nodes = [(float(r["c_mol_m3"]), float(r["D_s_local_m2_s"])) for r in spots]
    if any(not (1 < c < 10) for c, _ in nodes):
        raise ArithmeticError("New node outside target interval")
    refined.nodes = sorted(original.nodes + nodes)
    old_total, new_total = original.integral(0, 10), refined.integral(0, 10)
    output = []
    for c, _ in sorted(nodes):
        lower = max(x for x, _ in original.nodes if x < c)
        upper = min(x for x, _ in original.nodes if x > c)
        if sum(lower < other < upper for other, _ in nodes) != 1:
            raise ArithmeticError("More than one new PB point per old interval")
        old_piece, new_piece = original.integral(lower, upper), refined.integral(lower, upper)
        output.append({
            "lower_mM": lower, "new_node_mM": c, "upper_mM": upper,
            "old_integral_mol_m_s": old_piece,
            "augmented_integral_mol_m_s": new_piece,
            "signed_change_mol_m_s": new_piece - old_piece,
            "total_old_0_10_mol_m_s": old_total,
            "total_augmented_0_10_mol_m_s": new_total,
            "total_relative_change": abs(new_total - old_total) / old_total,
            "numeric_status": "INDEPENDENT_SPOTCHECK_ONLY_T12_NOT_FULLY_CLEARED",
        })
    if abs(sum(r["signed_change_mol_m_s"] for r in output) - (new_total - old_total)) > 1e-10 * abs(new_total - old_total):
        raise ArithmeticError("Four-interval sum does not reconcile")
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(output)
    print(f"new_nodes=4 old_integral={old_total:.10g} augmented_integral={new_total:.10g} "
          f"relative_change={abs(new_total-old_total)/old_total:.8g}")


if __name__ == "__main__":
    main()
