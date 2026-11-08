"""Conditional urea targets from archived Book2 salt curves under shared resistance scaling."""

from __future__ import annotations

import csv
import math
from pathlib import Path

from mechanism_neutral_baseline import D_SALT, D_UREA, hindrance_urea


ROOT = Path(__file__).resolve().parents[1]
BOOK = ROOT / "原始资料/赵珈绮/20_学生Stage2_B1成果_v1.0.0/WP1/book2_long.csv"
OUT = ROOT / "data/processed/salt_transport/22_legacy_urea_conditions_20261028.csv"
FIELDS = ["scenario", "sample_id", "state", "h_nm", "source_time_min", "scenario_time_min",
          "kappa_uScm", "kappa_initial_uScm", "c_increment_mol_m3", "salt_difference_fraction",
          "urea_receiver_W_mol_m3", "urea_receiver_C_anchor_mol_m3",
          "urea_W_minus_C_anchor_mol_m3", "alpha", "status"]


def concentration(kappa: float) -> float:
    if kappa == 0:
        return 0.0
    if kappa < 10:
        return kappa * concentration(10) / 10
    z = math.log10(kappa)
    return 1000 / 58.44 * 10 ** (0.0165 * z * z + 0.9398 * z - 3.2793)


def main() -> None:
    with BOOK.open(encoding="utf-8-sig", newline="") as stream:
        book = list(csv.DictReader(stream))
    assert len(book) == 41
    baseline = {r["sample_id"]: float(r["kappa_uScm"]) for r in book
                if float(r["time_min"]) == 0}
    assert len(baseline) == 6
    alpha = D_UREA * hindrance_urea() / D_SALT
    output = []
    for conversion in ("C1", "C2"):
        calculated = {}
        for item in book:
            sample = item["sample_id"]
            kappa = float(item["kappa_uScm"])
            initial = baseline[sample]
            increment = (concentration(kappa - initial) if conversion == "C1"
                         else concentration(kappa) - concentration(initial))
            theta = 1 - increment / 300
            if not 0 < theta <= 1:
                raise ArithmeticError("Salt difference fraction outside model interval")
            u_receiver = 5 * (1 - theta ** alpha)  # equal 10 mL pools; initial urea 10 mM
            calculated[(sample, float(item["time_min"]))] = (increment, theta, u_receiver)
        for timing in ("T85", "T80"):
            scenario = f"{conversion}_{timing}"
            for item in book:
                sample = item["sample_id"]
                time = float(item["time_min"])
                h = item["h_nm"]
                increment, theta, u_receiver = calculated[sample, time]
                b_anchor = calculated[f"B{h}", time][2]
                output.append({
                    "scenario": scenario, "sample_id": sample, "state": item["state"], "h_nm": h,
                    "source_time_min": time,
                    "scenario_time_min": 80 if timing == "T80" and time == 85 else time,
                    "kappa_uScm": item["kappa_uScm"], "kappa_initial_uScm": baseline[sample],
                    "c_increment_mol_m3": increment, "salt_difference_fraction": theta,
                    "urea_receiver_W_mol_m3": u_receiver,
                    "urea_receiver_C_anchor_mol_m3": b_anchor,
                    "urea_W_minus_C_anchor_mol_m3": u_receiver - b_anchor,
                    "alpha": alpha,
                    "status": "CONDITIONAL_SHARED_GEOMETRY_EXTERNAL_SCALING_NOT_VERIFIED",
                })
    if len(output) != 164:
        raise ArithmeticError("Unexpected Book2 scenario count")
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(output)
    print(f"rows={len(output)} alpha={alpha:.9g} "
          f"max_W_minus_C={max(abs(r['urea_W_minus_C_anchor_mol_m3']) for r in output):.6g}")


if __name__ == "__main__":
    main()
