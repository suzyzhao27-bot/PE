"""Map the sign of the frozen v5 model's thickness slope.

This is a conditional model regime map, not a measured PE phase diagram.
All lengths are nm; no old endpoint is used to refit parameters.
"""

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "core"))
from theory_model import log_resistance_nm, thickness_log_slope  # noqa: E402


def slope_scalar(eta, ws, alpha, robin_r, n=4096):
    """Return S, membrane resistance fraction, and field relief term."""
    if not (eta > 0 and ws >= 0 and alpha > 0 and robin_r >= 0):
        raise ValueError("eta, alpha > 0; ws, r >= 0 required")
    if n < 4 or n % 2:
        raise ValueError("Simpson interval count must be even and >=4")
    denominator = 1 + robin_r * math.exp(-eta)
    denominator_derivative = -robin_r * math.exp(-eta)
    total = derivative_total = 0.0
    for index in range(n + 1):
        u = index / n
        left = math.exp(-eta * u)
        right = math.exp(-eta * (1 - u))
        numerator = left + right
        derivative_numerator = -u * left - (1 - u) * right
        w = ws * numerator / denominator
        dw_deta = ws * (derivative_numerator / denominator -
                         numerator * denominator_derivative / denominator**2)
        simpson = 1 if index in (0, n) else (4 if index % 2 else 2)
        weighted = simpson * math.exp(w)
        total += weighted
        derivative_total += weighted * dw_deta
    integral = total / (3 * n)
    average_derivative = derivative_total / total
    membrane_fraction = eta * integral / (eta * integral + alpha)
    field_relief = -membrane_fraction * eta * average_derivative
    return -membrane_fraction + field_relief, membrane_fraction, field_relief


def slope_grid(eta_values, ws_values, alpha, robin_r, n=2048):
    """Vectorize the same Simpson expression over the plotting grid."""
    u = np.linspace(0, 1, n + 1)
    coeff = np.ones(n + 1)
    coeff[1:-1:2] = 4
    coeff[2:-1:2] = 2
    output = np.empty((len(ws_values), len(eta_values)))
    for column, eta in enumerate(eta_values):
        left = np.exp(-eta * u)
        right = np.exp(-eta * (1 - u))
        numerator = left + right
        derivative_numerator = -u * left - (1 - u) * right
        denominator = 1 + robin_r * math.exp(-eta)
        denominator_derivative = -robin_r * math.exp(-eta)
        field = ws_values[:, None] * numerator[None, :] / denominator
        derivative = ws_values[:, None] * (derivative_numerator[None, :] / denominator -
                       numerator[None, :] * denominator_derivative / denominator**2)
        weights = coeff[None, :] * np.exp(field)
        integral = weights.sum(axis=1) / (3 * n)
        mean_derivative = (weights * derivative).sum(axis=1) / weights.sum(axis=1)
        membrane_fraction = eta * integral / (eta * integral + alpha)
        output[:, column] = -membrane_fraction * (1 + eta * mean_derivative)
    return output


def roots(ws, alpha, robin_r, n):
    scan = np.geomspace(0.01, 100, 201)
    values = [slope_scalar(float(eta), ws, alpha, robin_r, n)[0] for eta in scan]
    answer = []
    for index in range(len(scan) - 1):
        if values[index] * values[index + 1] >= 0:
            continue
        low, high = float(scan[index]), float(scan[index + 1])
        sign_low = values[index] < 0
        for _ in range(45):
            mid = (low + high) / 2
            if (slope_scalar(mid, ws, alpha, robin_r, n)[0] < 0) == sign_low:
                low = mid
            else:
                high = mid
        answer.append((low + high) / 2)
    return answer


def write_csv(path, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--figure", type=Path, required=True)
    parser.add_argument("--grid-csv", type=Path, required=True)
    parser.add_argument("--boundaries-csv", type=Path, required=True)
    args = parser.parse_args()
    params = json.loads((HERE / "core/inputs/model_parameters.json").read_text(encoding="utf-8"))
    ws_ref = float(params["ws_ref"])
    ell = float(params["ell_nm"])
    a = float(params["a_ref_nm"])
    alpha = math.pi * a / (2 * ell)
    eta_grid = np.geomspace(0.02, 20, 121)
    ws_grid = np.linspace(0, 15, 121)
    thicknesses = (65, 120, 200, 350, 700, 2600)

    grid_rows = []
    boundary_rows = []
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2), sharey=True, constrained_layout=True)
    for axis, robin_r in zip(axes, (0.0, 0.75)):
        slopes = slope_grid(eta_grid, ws_grid, alpha, robin_r)
        grid_check = max(abs(slopes[i, j] - slope_scalar(float(eta_grid[j]),
                             float(ws_grid[i]), alpha, robin_r, 4096)[0])
                         for i in (0, 60, 120) for j in (0, 60, 120))
        if grid_check > 2e-5:
            raise AssertionError(f"Plot grid quadrature mismatch: {grid_check}")
        for i, ws in enumerate(ws_grid):
            for j, eta in enumerate(eta_grid):
                grid_rows.append(dict(robin_r=robin_r, eta=float(eta), ws_kBT=float(ws),
                                      alpha=alpha, slope=float(slopes[i, j]),
                                      thicker_increases_model_P=bool(slopes[i, j] > 0)))
        axis.contourf(eta_grid, ws_grid, slopes, levels=[-100, 0, 100],
                      colors=["#d6e8f1", "#f3d6ae"], alpha=0.9)
        axis.contour(eta_grid, ws_grid, slopes, levels=[0], colors=["#ad3b2d"], linewidths=2)
        axis.axhline(ws_ref, color="#555555", linestyle="--", linewidth=1)
        for thickness in thicknesses:
            axis.plot(thickness / ell, ws_ref, marker="o", markersize=4, color="black")
            axis.annotate(str(thickness), (thickness / ell, ws_ref), xytext=(2, 4),
                          textcoords="offset points", fontsize=7)
        axis.set_xscale("log")
        axis.set_xlim(0.02, 20)
        axis.set_ylim(0, 15)
        axis.set_xlabel(r"$\eta=h/\ell$  (dimensionless)")
        axis.set_title(f"Fixed Robin r = {robin_r:g}")
        axis.grid(alpha=0.18)
        refined = roots(ws_ref, alpha, robin_r, 4096)
        coarse = roots(ws_ref, alpha, robin_r, 2048)
        if len(refined) != 2 or len(coarse) != 2:
            raise AssertionError(f"Expected two reference roots for r={robin_r}")
        for index, (root, coarse_root) in enumerate(zip(refined, coarse), start=1):
            difference = abs(root - coarse_root)
            if difference > 2e-6 or abs(slope_scalar(root, ws_ref, alpha, robin_r)[0]) > 1e-9:
                raise AssertionError("Boundary failed quadrature refinement")
            boundary_rows.append(dict(scenario="reference_a25", robin_r=robin_r, boundary_index=index,
                                      eta_coarse=coarse_root, eta_refined=root,
                                      n_coarse=2048, n_refined=4096,
                                      absolute_eta_difference=difference,
                                      h_nm_at_frozen_ell=root * ell,
                                      ws_kBT=ws_ref, ell_nm=ell, a_nm=a,
                                      alpha=alpha, not_a_measured_threshold=True))

        # The conditional S/V closure changes both ws and the access term alpha.
        for conditional_a in (20.0, 30.0):
            conditional_ws = ws_ref * a / conditional_a
            conditional_alpha = math.pi * conditional_a / (2 * ell)
            conditional_coarse = roots(conditional_ws, conditional_alpha, robin_r, 4096)
            conditional_refined = roots(conditional_ws, conditional_alpha, robin_r, 8192)
            if len(conditional_coarse) != 2 or len(conditional_refined) != 2:
                raise AssertionError("Conditional radius scan did not find two boundaries")
            for index, (root, coarse_root) in enumerate(zip(conditional_refined, conditional_coarse), start=1):
                difference = abs(root - coarse_root)
                if difference > 2e-6:
                    raise AssertionError("Conditional radius boundary failed refinement")
                boundary_rows.append(dict(scenario="conditional_ws_inverse_a", robin_r=robin_r,
                                          boundary_index=index, eta_coarse=coarse_root, eta_refined=root,
                                          n_coarse=4096, n_refined=8192,
                                          absolute_eta_difference=difference,
                                          h_nm_at_frozen_ell=root * ell, ws_kBT=conditional_ws,
                                          ell_nm=ell, a_nm=conditional_a, alpha=conditional_alpha,
                                          not_a_measured_threshold=True))

    axes[0].set_ylabel(r"$w_s$  ($k_B T$ units)")
    fig.suptitle("Frozen v5 model: where thickness increases relative transfer")
    fig.text(0.5, -0.015, "Orange: d ln P / d ln h > 0; blue: < 0. Dots: old thickness labels at the frozen reference ws, ell.\n"
             "r=0.75 uses the same ws and ell, without refitting. Regions are model scenarios, not measured PE properties.",
             ha="center", va="top", fontsize=8)
    fig.savefig(args.figure, dpi=220, bbox_inches="tight")
    plt.close(fig)
    write_csv(args.grid_csv, grid_rows)
    write_csv(args.boundaries_csv, boundary_rows)

    # Existing r=0 analytical implementation and independent centered log-R differences.
    checks = []
    for robin_r in (0.0, 0.75):
        for h in thicknesses:
            eta = h / ell
            direct = slope_scalar(eta, ws_ref, alpha, robin_r)[0]
            step = 1e-4
            fd = -(log_resistance_nm(h * math.exp(step), ws_ref, ell, a, robin_r) -
                   log_resistance_nm(h * math.exp(-step), ws_ref, ell, a, robin_r)) / (2 * step)
            checks.append(abs(direct - fd))
            if robin_r == 0:
                checks.append(abs(direct - thickness_log_slope(h, ws_ref, ell, a, n=4096)))
    if max(checks) > 2e-6:
        raise AssertionError(f"Thickness slope mismatch: {max(checks)}")
    print(f"grid_rows={len(grid_rows)}; boundaries={len(boundary_rows)}; "
          f"max_slope_check_error={max(checks):.3g}")
    for row in boundary_rows:
        print(f"r={row['robin_r']:g} boundary {row['boundary_index']}: "
              f"eta={row['eta_refined']:.9f}, h={row['h_nm_at_frozen_ell']:.3f} nm, "
              f"refinement={row['absolute_eta_difference']:.3g}")


if __name__ == "__main__":
    main()
