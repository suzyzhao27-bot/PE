"""Audit fixed 60–120 and 120–240 min treatment-response increments."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/processed/salt_transport"
SOURCE = DATA / "22_paired_illustrations_20261028.csv"
OUT = DATA / "22_W_C_WC_处理响应分窗方向_20261029.csv"
FIELDS = ("pair_id", "h_nm", "salt_feed_mM", "mechanism", "external", "probe",
          "response_60_uM", "response_120_uM", "response_240_uM",
          "increment_60_120_uM", "increment_120_240_uM",
          "opposite_increment_sign", "numeric_status")


def main() -> None:
    with SOURCE.open(encoding="utf-8-sig", newline="") as stream:
        source = list(csv.DictReader(stream))
    if len(source) != 864:
        raise ArithmeticError("Unexpected frozen curve size")
    grouped = defaultdict(dict)
    metadata = {}
    for row in source:
        key = row["pair_id"], row["probe"]
        state_time = row["state"], int(float(row["time_min"]))
        if state_time in grouped[key]:
            raise ArithmeticError("Duplicate frozen state")
        grouped[key][state_time] = float(row["receiver_mol_m3"])
        metadata[row["pair_id"]] = row
    if len(grouped) != 108:
        raise ArithmeticError("Expected 54 pairs times two probes")
    output = []
    for (pair_id, probe), states in sorted(grouped.items()):
        if set(states) != {(state, minute) for state in ("B", "E") for minute in (0, 60, 120, 240)}:
            raise ArithmeticError("Incomplete B/E time states")
        response = {minute: 1000 * (states["E", minute] - states["B", minute])
                    for minute in (60, 120, 240)}
        first = response[120] - response[60]
        second = response[240] - response[120]
        info = metadata[pair_id]
        output.append({
            "pair_id": pair_id, "h_nm": info["h_nm"],
            "salt_feed_mM": info["c_feed_mM"],
            "mechanism": info["mechanism"], "external": info["external"],
            "probe": probe, "response_60_uM": response[60],
            "response_120_uM": response[120], "response_240_uM": response[240],
            "increment_60_120_uM": first, "increment_120_240_uM": second,
            "opposite_increment_sign": int(first * second < 0),
            "numeric_status": ("NUMERIC_UNRESOLVED_T12" if probe == "salt"
                               and int(float(info["c_feed_mM"])) == 10
                               and info["mechanism"] in ("C", "WC") else
                               "MODEL_ONLY_TARGET_CHECK"),
        })
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(output)
    reversals = [r for r in output if r["opposite_increment_sign"]]
    print(f"series={len(output)} opposite_subwindow_increments={len(reversals)} "
          f"numeric_unresolved_reversals={sum(r['numeric_status']=='NUMERIC_UNRESOLVED_T12' for r in reversals)}")


if __name__ == "__main__":
    main()
