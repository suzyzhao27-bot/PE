"""Quantify finite-pool approach to equal-reservoir equilibrium in frozen curves."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/processed/salt_transport"
SOURCE = DATA / "22_paired_illustrations_20261028.csv"
OUT = DATA / "22_双池条件曲线平衡接近程度_20261029.csv"
FIELDS = ("pair_id", "h_nm", "salt_feed_mM", "mechanism", "external", "state", "probe",
          "equilibrium_concentration_mM", "receiver_fraction_60", "receiver_fraction_120",
          "receiver_fraction_240", "remaining_gradient_fraction_240",
          "within_1pct_of_equilibrium_240", "numeric_status")


def main() -> None:
    with SOURCE.open(encoding="utf-8-sig", newline="") as stream:
        source = list(csv.DictReader(stream))
    if len(source) != 864:
        raise ArithmeticError("Unexpected frozen curve size")
    grouped = defaultdict(dict)
    metadata = {}
    for row in source:
        key = row["pair_id"], row["state"], row["probe"]
        minute = int(float(row["time_min"]))
        if minute in grouped[key]:
            raise ArithmeticError("Duplicate time state")
        feed = float(row["c_feed_mM"]) if row["probe"] == "salt" else 10.0
        donor, receiver = float(row["donor_mol_m3"]), float(row["receiver_mol_m3"])
        if abs(donor + receiver - feed) > 1e-8 * feed:
            raise ArithmeticError("Equal-volume pool concentration sum changed")
        grouped[key][minute] = receiver
        metadata[key] = row
    if len(grouped) != 216:
        raise ArithmeticError("Expected 54 pairs x B/E x salt/urea")
    output = []
    for key, times in sorted(grouped.items()):
        if set(times) != {0, 60, 120, 240} or abs(times[0]) > 1e-12:
            raise ArithmeticError("Incomplete curve time states")
        pair_id, state, probe = key
        row = metadata[key]
        feed = float(row["c_feed_mM"]) if probe == "salt" else 10.0
        fraction = {minute: 2 * times[minute] / feed for minute in (60, 120, 240)}
        if any(not (-1e-12 <= x <= 1 + 1e-9) for x in fraction.values()):
            raise ArithmeticError("Receiver exceeds equal-pool equilibrium")
        output.append({
            "pair_id": pair_id, "h_nm": row["h_nm"],
            "salt_feed_mM": row["c_feed_mM"], "mechanism": row["mechanism"],
            "external": row["external"], "state": state, "probe": probe,
            "equilibrium_concentration_mM": feed / 2,
            "receiver_fraction_60": fraction[60],
            "receiver_fraction_120": fraction[120],
            "receiver_fraction_240": fraction[240],
            "remaining_gradient_fraction_240": 1 - fraction[240],
            "within_1pct_of_equilibrium_240": int(fraction[240] >= .99),
            "numeric_status": row["numeric_status"],
        })
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(output)
    print(f"series={len(output)} within_1pct_240={sum(r['within_1pct_of_equilibrium_240'] for r in output)} "
          f"within_1pct_60={sum(r['receiver_fraction_60']>=.99 for r in output)}")


if __name__ == "__main__":
    main()
