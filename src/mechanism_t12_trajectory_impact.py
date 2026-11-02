"""Compare 31/16 PB concentration tables at saved 10 mM charged pool states."""

from __future__ import annotations

import csv
import math
from pathlib import Path

from mechanism_transport import PoreTable, salt_flux


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/processed/salt_transport"
OUT = DATA / "22_T12_10mM池状态通量差_20261029.csv"
FIELDS = ("pair_id", "h_nm", "mechanism", "external", "time_min", "donor_mM",
          "receiver_mM", "J_31_mol_m2_s", "J_16_mol_m2_s", "absolute_difference_mol_m2_s",
          "relative_difference", "numeric_status")


def rows(name: str) -> list[dict[str, str]]:
    with (DATA / name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    main_pb = rows("22_pb_main_20261028.csv")
    extras = rows("22_pb_target_refinement_20261028.csv")
    curves = rows("22_paired_illustrations_20261028.csv")
    ledger = rows("22_pb_call_ledger_20261028.csv")
    if (len(main_pb), len(curves), len(ledger)) != (372, 864, 437):
        raise ArithmeticError("Frozen input count changed")
    family = sorted((float(r["c_mol_m3"]), float(r["D_s_local_m2_s"])) for r in main_pb
                    if float(r["a_nm"]) == 25 and float(r["sigma_C_m2"]) == -0.008)
    endpoints = sorted((float(r["c_mol_m3"]), float(r["D_s_local_m2_s"])) for r in extras
                       if float(r["a_nm"]) == 25 and float(r["sigma_C_m2"]) == -0.008)
    if len(family) != 31 or [c for c, _ in endpoints] != [1.0, 10.0]:
        raise ArithmeticError("Unexpected 25 nm / -8 mC/m² PB nodes")
    fine = PoreTable()
    if fine.nodes != sorted([(0.0, 0.0)] + family + endpoints):
        raise ArithmeticError("Fine table differs from actual transport input")
    coarse = PoreTable()
    coarse.nodes = sorted([(0.0, 0.0)] + family[::2] + endpoints)
    if len(fine.nodes) != 34 or len(coarse.nodes) != 19:
        raise ArithmeticError("Unexpected table size")

    output = []
    for row in curves:
        if not (float(row["c_feed_mM"]) == 10 and row["mechanism"] in ("C", "WC")
                and row["state"] == "B" and row["probe"] == "salt"):
            continue
        donor = float(row["donor_mol_m3"])
        receiver = float(row["receiver_mol_m3"])
        h = float(row["h_nm"])
        b = .01 if row["mechanism"] == "C" else .0001
        args = (donor, receiver, h, b, -.008, row["external"])
        j31, j16 = salt_flux(*args, fine), salt_flux(*args, coarse)
        difference = abs(j31 - j16)
        relative = difference / abs(j31) if j31 else (0.0 if j16 == 0 else math.inf)
        if not all(math.isfinite(x) and x >= 0 for x in (j31, j16, difference, relative)):
            raise ArithmeticError("Invalid flux comparison")
        output.append({"pair_id": row["pair_id"], "h_nm": h,
                       "mechanism": row["mechanism"], "external": row["external"],
                       "time_min": row["time_min"], "donor_mM": donor,
                       "receiver_mM": receiver, "J_31_mol_m2_s": j31,
                       "J_16_mol_m2_s": j16, "absolute_difference_mol_m2_s": difference,
                       "relative_difference": relative,
                       "numeric_status": "NUMERIC_UNRESOLVED_T12"})
    keys = {(r["pair_id"], r["time_min"]) for r in output}
    if len(output) != 48 or len(keys) != 48 or len({r["pair_id"] for r in output}) != 12:
        raise ArithmeticError("Expected 12 pairs x 4 saved times")
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(output)
    print(f"states={len(output)} max_relative={max(r['relative_difference'] for r in output):.6g} "
          f"over_2pct={sum(r['relative_difference'] > .02 for r in output)} pb_calls=437")


if __name__ == "__main__":
    main()
