"""Check whether four new PB points change 10 mM subwindow direction labels."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/processed/salt_transport"
OUT = DATA / "22_T12_10mM补点后分窗方向复核_20261029.csv"
FIELDS = ("pair_id", "h_nm", "mechanism", "external",
          "old_increment_60_120_uM", "new_increment_60_120_uM",
          "old_increment_120_240_uM", "new_increment_120_240_uM",
          "old_opposite_sign", "new_opposite_sign", "direction_label_stable",
          "numeric_status")


def read(name: str) -> list[dict[str, str]]:
    with (DATA / name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    original = read("22_W_C_WC_处理响应分窗方向_20261029.csv")
    curves = read("22_pb_postfreeze_10mM完整时序影响_20261029.csv")
    frozen = read("22_paired_illustrations_20261028.csv")
    if len(original) != 108 or len(curves) != 48 or len(frozen) != 864:
        raise ArithmeticError("Unexpected source size")
    old = {r["pair_id"]: r for r in original if r["probe"] == "salt"
           and int(float(r["salt_feed_mM"])) == 10 and r["mechanism"] in ("C", "WC")}
    e = {(r["pair_id"], int(float(r["time_min"]))): float(r["receiver_mol_m3"])
         for r in frozen if r["state"] == "E" and r["probe"] == "salt"
         and float(r["c_feed_mM"]) == 10 and r["mechanism"] in ("C", "WC")}
    b = {(r["pair_id"], int(r["time_min"])): float(r["augmented_receiver_mM"])
         for r in curves}
    if len(old) != 12 or len(e) != 48 or len(b) != 48:
        raise ArithmeticError("Expected twelve charged 10 mM series")
    output = []
    for pair_id, row in sorted(old.items()):
        response = {t: 1000 * (e[pair_id, t] - b[pair_id, t]) for t in (60, 120, 240)}
        new_first = response[120] - response[60]
        new_second = response[240] - response[120]
        old_first = float(row["increment_60_120_uM"])
        old_second = float(row["increment_120_240_uM"])
        new_flag = int(new_first * new_second < 0)
        old_flag = int(row["opposite_increment_sign"])
        output.append({
            "pair_id": pair_id, "h_nm": row["h_nm"],
            "mechanism": row["mechanism"], "external": row["external"],
            "old_increment_60_120_uM": old_first,
            "new_increment_60_120_uM": new_first,
            "old_increment_120_240_uM": old_second,
            "new_increment_120_240_uM": new_second,
            "old_opposite_sign": old_flag, "new_opposite_sign": new_flag,
            "direction_label_stable": int(old_flag == new_flag and
                                          old_first * new_first >= 0 and
                                          old_second * new_second >= 0),
            "numeric_status": "POSTFREEZE_SPOTCHECK_T12_UNRESOLVED",
        })
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(output)
    print(f"series={len(output)} old_reversals={sum(r['old_opposite_sign'] for r in output)} "
          f"new_reversals={sum(r['new_opposite_sign'] for r in output)} "
          f"direction_stable={sum(r['direction_label_stable'] for r in output)}")


if __name__ == "__main__":
    main()
