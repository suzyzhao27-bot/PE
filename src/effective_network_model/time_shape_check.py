"""Check old PE curves against a conditional constant-P finite-bath shape bound.

For nonnegative theta and 0 <= t <= T, the normalized, well-mixed
constant-P solution (1-exp(-theta*t/T))/(1-exp(-theta)) is concave and at
least t/T. This is a shape check, not a dynamic parameter fit or noise test.
"""

import csv
import math
from pathlib import Path


SOURCE = Path(__file__).resolve().parent / "core" / "inputs" / "Book1__Sheet1.csv"
CURVES = ((65, "col_2_B"), (120, "col_4_D"), (200, "col_5_E"),
          (350, "col_6_F"), (700, "col_7_G"), (2600, "col_8_H"))
COMMON_TIMES = (20, 40, 60, 100)


def main():
    with SOURCE.open(encoding="utf-8-sig", newline="") as source:
        rows = list(csv.DictReader(source))
    if len(rows) != 7:
        raise ValueError("Expected seven frozen Book1 rows; do not silently remap time.")
    for index, expected in ((0, 0), (1, 20), (2, 40), (3, 60), (5, 100), (6, 120)):
        if any(float(rows[index][key]) != expected for key in ("col_1_A", "col_3_C")):
            raise ValueError(f"A/C common time labels differ from {expected} min at row {index}.")

    print("h_nm f20 f40 f60 f100 min_minus_t_over_120 at_min")
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
        differences = [value - time / 120 for value, time in zip(fractions, COMMON_TIMES)]
        worst_index = min(range(len(differences)), key=differences.__getitem__)
        print(
            f"{thickness} "
            + " ".join(f"{value:.5f}" for value in fractions)
            + f" {differences[worst_index]:+.5f} {COMMON_TIMES[worst_index]}"
        )


if __name__ == "__main__":
    main()
