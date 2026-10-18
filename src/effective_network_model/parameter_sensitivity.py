"""Local identifiability of the frozen shared fit after removing one label.

Use centered log-model sensitivities so the unknown common amplitude drops out.
Singular values are local geometry diagnostics, not confidence intervals.
"""

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "core"))
from theory_model import log_resistance_nm  # noqa: E402

LABELS = (65, 120, 200, 350, 700, 2600)


def sensitivities(params, step, n=2048):
    center = np.log([params["ws_ref"], params["ell_nm"]])
    columns = []
    for k in (0, 1):
        offset = np.zeros(2)
        offset[k] = step
        values = []
        for point in (center + offset, center - offset):
            ws, ell = np.exp(point)
            values.append(np.array([
                -log_resistance_nm(h, float(ws), float(ell), params["a_ref_nm"], n=n)
                for h in LABELS
            ]))
        columns.append((values[0] - values[1]) / (2 * step))
    return np.column_stack(columns)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--leverage-output", type=Path,
                        help="Optional held-out contrast along each five-label weak direction")
    args = parser.parse_args()
    params = json.loads((ROOT / "core/inputs/model_parameters.json").read_text(encoding="utf-8"))
    if float(params["robin_r"]) != 0:
        raise ValueError("Only the frozen r=0 fit is evaluated")
    jacobian = sensitivities(params, 1e-4)
    check = sensitivities(params, 5e-5)
    if float(np.max(np.abs(jacobian - check))) > 1e-5:
        raise AssertionError("Finite-difference log sensitivities did not stabilize")
    rows = []
    leverage_rows = []
    for removed in (None, *LABELS):
        indices = [i for i, h in enumerate(LABELS) if h != removed]
        centered = jacobian[indices] - jacobian[indices].mean(axis=0)
        _, singular, right_vectors = np.linalg.svd(centered, full_matrices=False)
        if singular[1] <= 0 or not all(math.isfinite(float(value)) for value in singular):
            raise AssertionError("Invalid centered sensitivity singular values")
        rows.append(dict(removed_h_nm="none" if removed is None else removed,
                         included_label_count=len(indices),
                         strong_singular_value=float(singular[0]),
                         weak_singular_value=float(singular[1]),
                         condition_number=float(singular[0] / singular[1]),
                         parameter_point="frozen_full_six_fit",
                         log_parameters="ln(ws_ref),ln(ell_nm)",
                         common_amplitude="centered_out"))
        if removed is not None:
            weak_direction = right_vectors[-1]
            held_index = LABELS.index(removed)
            held_contrast = float(np.dot(jacobian[held_index] - jacobian[indices].mean(axis=0),
                                         weak_direction))
            train_rms = float(singular[-1] / math.sqrt(len(indices)))
            leverage_rows.append(dict(heldout_h_nm=removed, included_label_count=len(indices),
                                      train_rms_logP_sensitivity=train_rms,
                                      heldout_contrast_logP_sensitivity=abs(held_contrast),
                                      heldout_to_train_ratio=abs(held_contrast) / train_rms,
                                      parameter_point="frozen_full_six_fit",
                                      direction="weak singular vector of five included labels"))
    with args.output.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    for row in rows:
        print(f"remove={row['removed_h_nm']}: weak={row['weak_singular_value']:.6f}, "
              f"condition={row['condition_number']:.6f}")
    if args.leverage_output:
        with args.leverage_output.open("w", encoding="utf-8-sig", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(leverage_rows[0]))
            writer.writeheader()
            writer.writerows(leverage_rows)
        for row in leverage_rows:
            print(f"heldout={row['heldout_h_nm']}: "
                  f"relative_weak_direction_leverage={row['heldout_to_train_ratio']:.6f}")


if __name__ == "__main__":
    main()
