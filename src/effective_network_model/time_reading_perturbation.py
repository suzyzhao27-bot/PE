"""Sharp raw-reading perturbation threshold for the 200 nm time-shape test.

The calculation assumes a one-state well-mixed finite bath with nonnegative,
nondecreasing instantaneous k(t). Epsilon is a hypothetical per-reading bound,
not an instrument uncertainty or a statistical significance estimate.
"""

import argparse
import csv
import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "core/inputs/Book1__Sheet1.csv"
TIMES = (0, 20, 40, 120)
ROW_INDICES = (0, 1, 2, 6)
COLUMN = "col_5_E"


def raw_readings():
    with SOURCE.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 7:
        raise ValueError("Expected seven Book1 observations")
    if any(float(rows[index]["col_3_C"]) != time
           for index, time in zip(ROW_INDICES, TIMES)):
        raise ValueError("200 nm time column C differs from recorded labels")
    return tuple(float(rows[index][COLUMN]) for index in ROW_INDICES)


def best_readings(raw, epsilon):
    """Extrema that maximize the first convexity feasibility condition."""
    y0, y20, y40, y120 = raw
    return (y0 + epsilon, y20 - epsilon, y40 + epsilon, y120 - epsilon)


def increments(readings):
    baseline = readings[0]
    return tuple(value - baseline for value in readings[1:])


def first_condition_margin(readings):
    q20, q40, q120 = increments(readings)
    return q20*q20 - q120*(2*q20 - q40)


def transformed_slopes(readings, equilibrium_increment):
    q20, q40, q120 = increments(readings)
    if not 0 < q20 < q40 < q120 < equilibrium_increment:
        raise ValueError("Need increasing cumulative readings below finite equilibrium")
    x20, x40, x120 = (-math.log1p(-q/equilibrium_increment)
                       for q in (q20, q40, q120))
    return x20/20, (x40-x20)/20, (x120-x40)/80


def normalized_interval(start, middle, end, epsilon):
    """Conservative exact per-time bounds with shared raw baseline/endpoints."""
    if not start + epsilon < middle - epsilon < middle + epsilon < end - epsilon:
        raise ValueError("Perturbation box breaks the observed ordering")
    minimum = ((middle-epsilon)-(start+epsilon))/((end+epsilon)-(start+epsilon))
    maximum = ((middle+epsilon)-(start-epsilon))/((end-epsilon)-(start-epsilon))
    return minimum, maximum


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    raw = raw_readings()
    endpoint_increment = raw[3] - raw[0]
    if tuple(round(x, 2) for x in raw) != (1.53, 2.35, 2.70, 15.45):
        raise AssertionError("Frozen 200 nm source readings changed")
    if abs(endpoint_increment - 13.92) > 1e-12:
        raise AssertionError("Unexpected 200 nm endpoint increment")

    # For all eps <= this search ceiling, 2q20-q40 stays positive for every
    # allowed reading vector. The first nondecreasing-slope condition is then
    # possible only if Qeq <= q20^2/(2q20-q40), while Qeq >= q120.
    ceiling = 0.1175
    if first_condition_margin(best_readings(raw, 0)) >= 0:
        raise AssertionError("Original data unexpectedly feasible")
    if first_condition_margin(best_readings(raw, ceiling)) <= 0:
        raise AssertionError("Perturbation threshold not bracketed")
    low, high = 0.0, ceiling
    for _ in range(70):
        mid = (low + high)/2
        if first_condition_margin(best_readings(raw, mid)) < 0:
            low = mid
        else:
            high = mid
    epsilon_infimum = (low + high)/2
    boundary_readings = best_readings(raw, epsilon_infimum)
    q20, q40, q120 = increments(boundary_readings)
    if not 0 < q20 < q40 < q120 or abs(first_condition_margin(boundary_readings)) > 1e-11:
        raise AssertionError("Analytical lower boundary failed")

    # At the exact infimum Qeq=q120, X(120) diverges. Show a finite feasible
    # witness just above it, with all four raw readings in the same +/-eps box.
    witness_epsilon = epsilon_infimum + 1e-4
    witness_readings = best_readings(raw, witness_epsilon)
    w20, w40, w120 = increments(witness_readings)
    upper_equilibrium = w20*w20/(2*w20-w40)
    if upper_equilibrium <= w120:
        raise AssertionError("No finite equilibrium interval for witness")
    equilibrium = w120 + 0.1*(upper_equilibrium-w120)
    slopes = transformed_slopes(witness_readings, equilibrium)
    if not (0 < slopes[0] <= slopes[1] <= slopes[2]):
        raise AssertionError("Witness does not meet all segment slope conditions")

    with SOURCE.open(encoding="utf-8-sig", newline="") as stream:
        all_rows = list(csv.DictReader(stream))
    crossing_check = []
    for row_index, time in ((1, 20), (2, 40)):
        if any(float(all_rows[row_index][key]) != time for key in ("col_1_A", "col_3_C")):
            raise ValueError("Selected 120/2600 nm time labels are not common")
        a_values = [float(all_rows[index]["col_4_D"]) for index in (0, row_index, 6)]
        b_values = [float(all_rows[index]["col_8_H"]) for index in (0, row_index, 6)]
        a_low, a_high = normalized_interval(*a_values, epsilon_infimum)
        b_low, b_high = normalized_interval(*b_values, epsilon_infimum)
        gap_low, gap_high = a_low-b_high, a_high-b_low
        crossing_check.append(dict(time_min=time, normalized_gap_lower=gap_low,
                                   normalized_gap_upper=gap_high,
                                   guaranteed_positive=gap_low > 0,
                                   guaranteed_negative=gap_high < 0))
    if not crossing_check[0]["guaranteed_positive"] or not crossing_check[1]["guaranteed_negative"]:
        raise AssertionError("120/2600 nm crossing is not robust at this budget")

    result = {
        "source": str(SOURCE.relative_to(HERE.parent.parent)),
        "curve": "Book1 C/E, 200 nm label",
        "times_min": TIMES,
        "raw_conductivity_readings": raw,
        "unit_status": "Book1 E column inherits table conductivity unit; curve-level confirmation pending",
        "original_endpoint_increment": endpoint_increment,
        "assumption": "Fixed-volume well-mixed finite bath; affine conductivity-to-concentration map; no sampling or extra salt; k(t)>=0 and nondecreasing",
        "epsilon_infimum_raw_table_units": epsilon_infimum,
        "epsilon_over_original_endpoint_increment": epsilon_infimum/endpoint_increment,
        "boundary_delta_pattern_0_20_40_120": (1, -1, 1, -1),
        "boundary_perturbed_readings": boundary_readings,
        "boundary_q20_q40_q120": (q20, q40, q120),
        "boundary_first_condition_margin": first_condition_margin(boundary_readings),
        "attainment": "Infimum only: at the exact boundary Qeq=q120 gives divergent X120; any slightly larger epsilon admits a finite-Q witness",
        "finite_witness": {
            "epsilon_raw_table_units": witness_epsilon,
            "readings_0_20_40_120": witness_readings,
            "Qeq_raw_table_increment_units": equilibrium,
            "Qeq_allowed_upper_bound": upper_equilibrium,
            "transformed_segment_slopes_per_min": slopes,
        },
        "conditional_120_2600_crossing_at_same_epsilon": crossing_check,
        "interpretation": "Hypothetical robustness threshold, not measured instrument error, confidence interval, or proof about all dynamic transport models",
    }
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"epsilon_infimum={epsilon_infimum:.9f} raw table units; "
          f"fraction_of_13.92={epsilon_infimum/endpoint_increment:.6%}")
    print(f"finite_witness_epsilon={witness_epsilon:.9f}; "
          f"Qeq={equilibrium:.9f}; segment_slopes={slopes}")


if __name__ == "__main__":
    main()
