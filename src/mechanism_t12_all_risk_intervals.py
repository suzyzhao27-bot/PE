"""Localize the 31/16 concentration-integral difference for all ten unresolved targets."""

from __future__ import annotations

import csv
from pathlib import Path

from mechanism_transport import PoreTable


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/processed/salt_transport"
OUT = DATA / "22_T12_十项未决积分差定位_20261029.csv"
FIELDS = ("a_nm", "sigma_C_m2", "feed_mM", "relative_31_16", "interval_count",
          "largest_lower_mM", "largest_upper_mM", "largest_signed_share",
          "second_lower_mM", "second_upper_mM", "second_signed_share",
          "top2_signed_share", "positive_interval_count", "negative_interval_count",
          "numeric_status")


def read(name: str) -> list[dict[str, str]]:
    with (DATA / name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    main_pb = read("22_pb_main_20261028.csv")
    extras = read("22_pb_target_refinement_20261028.csv")
    targets = [r for r in read("22_pb_target_recheck_20261028.csv")
               if r["status"] == "NUMERIC_UNRESOLVED_OVER_2_PERCENT"]
    if (len(main_pb), len(extras), len(targets)) != (372, 12, 10):
        raise ArithmeticError("Unexpected PB source or target size")
    output = []
    for target in targets:
        a, sigma, feed = (float(target[k]) for k in ("a_nm", "sigma_C_m2", "feed_mM"))
        family = sorted((float(r["c_mol_m3"]), float(r["D_s_local_m2_s"]))
                        for r in main_pb if float(r["a_nm"]) == a
                        and float(r["sigma_C_m2"]) == sigma)
        endpoint = [(float(r["c_mol_m3"]), float(r["D_s_local_m2_s"]))
                    for r in extras if float(r["a_nm"]) == a
                    and float(r["sigma_C_m2"]) == sigma and float(r["c_mol_m3"]) == feed]
        if len(family) != 31 or len(endpoint) != 1:
            raise ArithmeticError("Incomplete PB family or exact endpoint")
        fine, coarse = PoreTable(int(a), sigma), PoreTable(int(a), sigma)
        fine.nodes = sorted([(0.0, 0.0)] + family + endpoint)
        coarse.nodes = sorted([(0.0, 0.0)] + family[::2] + endpoint)
        bounds = [c for c, _ in coarse.nodes if c < feed] + [feed]
        interval = [(low, high, fine.integral(low, high) - coarse.integral(low, high))
                    for low, high in zip(bounds, bounds[1:])]
        total = sum(item[2] for item in interval)
        fine_total = fine.integral(0, feed)
        coarse_total = coarse.integral(0, feed)
        relative = abs(total) / fine_total
        if abs(total - (fine_total - coarse_total)) > 1e-10 * abs(total):
            raise ArithmeticError("Interval sum does not reconcile")
        if abs(relative - float(target["updated_31_16_relative_difference"])) > 1e-10:
            raise ArithmeticError("Target recheck did not reproduce")
        ranked = sorted(interval, key=lambda item: abs(item[2]), reverse=True)
        first, second = ranked[:2]
        output.append({
            "a_nm": a, "sigma_C_m2": sigma, "feed_mM": feed,
            "relative_31_16": relative, "interval_count": len(interval),
            "largest_lower_mM": first[0], "largest_upper_mM": first[1],
            "largest_signed_share": first[2] / total,
            "second_lower_mM": second[0], "second_upper_mM": second[1],
            "second_signed_share": second[2] / total,
            "top2_signed_share": (first[2] + second[2]) / total,
            "positive_interval_count": sum(x[2] > 0 for x in interval),
            "negative_interval_count": sum(x[2] < 0 for x in interval),
            "numeric_status": target["status"],
        })
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(output)
    print(f"targets={len(output)} feed1={sum(r['feed_mM']==1 for r in output)} "
          f"feed10={sum(r['feed_mM']==10 for r in output)} "
          f"top2_share_range={min(r['top2_signed_share'] for r in output):.6g}.."
          f"{max(r['top2_signed_share'] for r in output):.6g}")


if __name__ == "__main__":
    main()
