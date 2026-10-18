"""Decompose the frozen r=0 model's local thickness response.

No refitting: this uses the v5 shared parameters and the six existing labels.
The output terms add to d ln(P)/d ln(h) at each label.
"""

import argparse
import csv
import json
import math
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "core"))
from theory_model import field, log_resistance_nm, thickness_log_slope  # noqa: E402


def balance(h, ws, ell, a, n=4096):
    # Integrate in u=z/h. The derivative holds u and all fitted parameters fixed.
    weights = []
    derivatives = []
    wmax = field(0, h, ws, ell)
    for i in range(n + 1):
        u = i / n
        first = math.exp(-u * h / ell)
        second = math.exp(-(1 - u) * h / ell)
        w = ws * (first + second)
        dw_dh = -ws / ell * (u * first + (1 - u) * second)
        coeff = 1 if i in (0, n) else (4 if i % 2 else 2)
        weights.append(coeff * math.exp(w - wmax))
        derivatives.append(dw_dh)
    total_weight = sum(weights)
    weighted_derivative = sum(weight * derivative for weight, derivative in zip(weights, derivatives)) / total_weight
    log_membrane = math.log(h * total_weight / (3 * n)) + wmax
    access = math.pi * a / 2
    membrane_fraction = 1 / (1 + math.exp(math.log(access) - log_membrane))
    path_term = -membrane_fraction
    field_relief_term = -membrane_fraction * h * weighted_derivative
    return membrane_fraction, path_term, field_relief_term, path_term + field_relief_term


def local_scale_threshold_rows(ws, ell, a):
    """Conditional 350 nm radius thresholds with 200 nm fixed at a_ref."""
    log_r200 = log_resistance_nm(200, ws, ell, a)
    observed_ratio = 11.89 / 13.92

    def ratio(a350):
        return math.exp(log_r200 - log_resistance_nm(350, ws * a / a350, ell, a350))

    rows = []
    for label, target in (("equal_permeance", 1.0), ("old_endpoint_ratio", observed_ratio)):
        low, high = a * 0.5, a
        if not ratio(low) < target < ratio(high):
            raise AssertionError("The radius threshold is not bracketed")
        for _ in range(60):
            mid = (low + high) / 2
            if ratio(mid) < target:
                low = mid
            else:
                high = mid
        a350 = (low + high) / 2
        predicted = ratio(a350)
        if abs(predicted - target) > 1e-11:
            raise AssertionError("Radius threshold root did not converge")
        rows.append(dict(target=label, h_reference_nm=200, h_test_nm=350,
                         a_reference_nm=a, a_test_nm=a350,
                         relative_a_change=(a350 / a - 1),
                         ws_test_kBT=ws * a / a350, ell_nm=ell,
                         target_P350_over_P200=target,
                         computed_P350_over_P200=predicted,
                         fixed_ws_P350_over_P200=math.exp(
                             log_r200 - log_resistance_nm(350, ws, ell, a350))))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--local-scale-threshold-output", type=Path)
    args = parser.parse_args()
    params = json.loads((ROOT / "core/inputs/model_parameters.json").read_text(encoding="utf-8"))
    ws, ell, a, r = (float(params[key]) for key in ("ws_ref", "ell_nm", "a_ref_nm", "robin_r"))
    if r != 0 or min(ws, ell, a) <= 0:
        raise ValueError("Only the frozen positive-parameter r=0 reference fit is supported")
    reference_log_r = log_resistance_nm(2600, ws, ell, a)
    rows = []
    for h in (65, 120, 200, 350, 700, 2600):
        fraction, path, relief, total = balance(h, ws, ell, a)
        existing = thickness_log_slope(h, ws, ell, a, n=4096)
        step = 1e-4
        finite_difference = -(log_resistance_nm(h * math.exp(step), ws, ell, a) -
                              log_resistance_nm(h * math.exp(-step), ws, ell, a)) / (2 * step)
        if max(abs(total - existing), abs(total - finite_difference)) > 2e-6:
            raise AssertionError(f"Thickness slope inconsistency at {h} nm")
        rows.append(dict(h_nm=h, a_nm=a, ws_ref_kBT=ws, ell_nm=ell, robin_r=r,
                         w_center_kBT=field(h / 2, h, ws, ell),
                         membrane_resistance_fraction=fraction,
                         path_length_term=path, field_relief_term=relief,
                         dlogP_dlogh=total,
                         relative_P_to_2600=math.exp(reference_log_r - log_resistance_nm(h, ws, ell, a))))
    with args.output.open("w", encoding="utf-8-sig", newline="") as target:
        writer = csv.DictWriter(target, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    if args.local_scale_threshold_output:
        threshold_rows = local_scale_threshold_rows(ws, ell, a)
        with args.local_scale_threshold_output.open("w", encoding="utf-8-sig", newline="") as target:
            writer = csv.DictWriter(target, fieldnames=list(threshold_rows[0]))
            writer.writeheader()
            writer.writerows(threshold_rows)
    for row in rows:
        print(f"{row['h_nm']:4d} nm: path={row['path_length_term']:.6f}, "
              f"field={row['field_relief_term']:.6f}, slope={row['dlogP_dlogh']:.6f}")


if __name__ == "__main__":
    main()
