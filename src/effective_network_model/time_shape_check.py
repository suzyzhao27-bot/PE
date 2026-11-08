"""Check old PE curves against a conditional constant-P finite-bath shape bound.

For nonnegative theta and 0 <= t <= T, the normalized, well-mixed
constant-P solution (1-exp(-theta*t/T))/(1-exp(-theta)) is concave and at
least t/T. This is a shape check, not a dynamic parameter fit or noise test.
"""

import argparse
import csv
import itertools
import math
from pathlib import Path


SOURCE = Path(__file__).resolve().parent / "core" / "inputs" / "Book1__Sheet1.csv"
CURVES = ((65, "col_2_B"), (120, "col_4_D"), (200, "col_5_E"),
          (350, "col_6_F"), (700, "col_7_G"), (2600, "col_8_H"))
COMMON_TIMES = (20, 40, 60, 100)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--common-gate-output", type=Path,
                        help="Optional CSV for the shared time-gate rank test")
    parser.add_argument("--crossing-output", type=Path,
                        help="Optional CSV for scale-invariant normalized-curve crossings")
    args = parser.parse_args()
    with SOURCE.open(encoding="utf-8-sig", newline="") as source:
        rows = list(csv.DictReader(source))
    if len(rows) != 7:
        raise ValueError("Expected seven frozen Book1 rows; do not silently remap time.")
    for index, expected in ((0, 0), (1, 20), (2, 40), (3, 60), (5, 100), (6, 120)):
        if any(float(rows[index][key]) != expected for key in ("col_1_A", "col_3_C")):
            raise ValueError(f"A/C common time labels differ from {expected} min at row {index}.")

    print("h_nm f20 f40 f60 f100 min_minus_t_over_120 at_min late_over_early_60")
    curves = []
    for thickness, column in CURVES:
        start = float(rows[0][column])
        delta = float(rows[6][column]) - start
        if not math.isfinite(delta) or delta <= 0:
            raise ValueError(f"Nonpositive or invalid endpoint increment at {thickness} nm.")
        fractions = [
            (float(rows[index][column]) - start) / delta
            for index in (1, 2, 3, 5)
        ]
        if not all(math.isfinite(value) for value in fractions):
            raise ValueError(f"Invalid normalized value at {thickness} nm.")
        curves.append(dict(h_nm=thickness, delta=delta, fractions=fractions))
        halfway_fraction = fractions[2]
        if not 0 < halfway_fraction < 1:
            raise ValueError(f"Invalid 60-minute fraction at {thickness} nm.")
        # Equal 60-minute windows: ratio of late to early mean increments.
        late_over_early = (1 - halfway_fraction) / halfway_fraction
        differences = [value - time / 120 for value, time in zip(fractions, COMMON_TIMES)]
        worst_index = min(range(len(differences)), key=differences.__getitem__)
        print(
            f"{thickness} "
            + " ".join(f"{value:.5f}" for value in fractions)
            + f" {differences[worst_index]:+.5f} {COMMON_TIMES[worst_index]}"
            + f" {late_over_early:.5f}"
        )

    print("h_nm delta_0_20 delta_20_40 delta_40_60 ratio_20_40_over_0_20 ratio_40_60_over_20_40")
    for thickness, column in CURVES:
        values = [float(rows[index][column]) for index in range(4)]
        increments = [right - left for left, right in zip(values, values[1:])]
        if not all(math.isfinite(value) and value > 0 for value in increments):
            raise ValueError(f"Nonpositive or invalid early 20-minute increment at {thickness} nm.")
        ratios = [later / earlier for earlier, later in zip(increments, increments[1:])]
        print(
            f"{thickness} "
            + " ".join(f"{value:.5f}" for value in increments + ratios)
        )

    if args.crossing_output:
        crossings = []
        for lower, upper in itertools.combinations(curves, 2):
            gaps = [left - right for left, right in
                    zip(lower["fractions"], upper["fractions"])]
            minimum, maximum = min(gaps), max(gaps)
            crossed = minimum < 0 < maximum
            crossings.append(dict(lower_h_nm=lower["h_nm"], upper_h_nm=upper["h_nm"],
                                  **{f"f_gap_at_{time}_min": gap
                                     for time, gap in zip(COMMON_TIMES, gaps)},
                                  min_gap=minimum, max_gap=maximum,
                                  crossing_margin=min(maximum, -minimum) if crossed else 0,
                                  rank_crossing=crossed))
        with args.crossing_output.open("w", encoding="utf-8-sig", newline="") as target:
            writer = csv.DictWriter(target, fieldnames=list(crossings[0]))
            writer.writeheader()
            writer.writerows(crossings)
        print(f"scale_invariant_crossings: "
              f"{sum(row['rank_crossing'] for row in crossings)}/{len(crossings)} pairs; "
              f"margin>=0.02: {sum(row['crossing_margin'] >= .02 for row in crossings)}")

    if args.common_gate_output:
        # Under a common nonnegative gate and common equilibrium amplitude,
        # larger endpoint increments imply larger theta. For 0<u<1,
        # (1-exp(-theta*u))/(1-exp(-theta)) increases with theta.
        # Thus the larger-endpoint curve cannot have a smaller f(t).
        comparisons = []
        for first, second in itertools.combinations(curves, 2):
            larger, smaller = sorted((first, second), key=lambda curve: curve["delta"], reverse=True)
            for index, time in enumerate(COMMON_TIMES):
                gap = larger["fractions"][index] - smaller["fractions"][index]
                comparisons.append(dict(time_min=time,
                                        larger_endpoint_h_nm=larger["h_nm"],
                                        smaller_endpoint_h_nm=smaller["h_nm"],
                                        larger_endpoint_delta=larger["delta"],
                                        smaller_endpoint_delta=smaller["delta"],
                                        larger_endpoint_fraction=larger["fractions"][index],
                                        smaller_endpoint_fraction=smaller["fractions"][index],
                                        fraction_gap=gap,
                                        rank_inversion=gap < 0))
        with args.common_gate_output.open("w", encoding="utf-8-sig", newline="") as target:
            writer = csv.DictWriter(target, fieldnames=list(comparisons[0]))
            writer.writeheader()
            writer.writerows(comparisons)
        for time in COMMON_TIMES:
            relevant = [row for row in comparisons if row["time_min"] == time]
            inversions = sum(row["rank_inversion"] for row in relevant)
            print(f"common_gate {time} min: {inversions}/{len(relevant)} endpoint-rank inversions")


if __name__ == "__main__":
    main()
