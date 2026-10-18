"""Plot reproduced 24-cell PNP fluxes with an explicit same-sigma baseline.

The transferred-fit point comes from a separate neutral-field calculation.
This figure is a numerical check, not an absolute fit to PE experiments.
"""

import argparse
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


SCENARIOS = ("zero", "rectangle_w3", "two_sided_ws3_ell26p6")
SIGMAS = (0.0, -0.005)


def load_points(path):
    data = json.loads(path.read_text(encoding="utf-8"))
    points = {}
    for entry in data["results"]:
        if entry["cells_per_fiber"] != 24:
            continue
        sigma = round(float(entry["result"]["sigma_c_m2"]), 6)
        key = (sigma, entry["scenario"])
        if key in points:
            raise ValueError(f"Duplicate 24-cell result: {key}")
        points[key] = float(entry["result"]["common_salt_flux_mol_m2_s"])
    expected = {(sigma, scenario) for sigma in SIGMAS for scenario in SCENARIOS}
    if points.keys() != expected or not all(math.isfinite(value) and value > 0 for value in points.values()):
        raise ValueError("Need six positive 24-cell fluxes at sigma=0/-0.005")
    return points


def load_fitted_ratio(path):
    data = json.loads(path.read_text(encoding="utf-8"))
    matches = [row for row in data["rows"] if row["cells_per_fiber"] == 24
               and row["status"] == "converged" and row["sigma_given_c_m2"] == 0
               and row["ws_kbt"] == data["requested_root_coefficients"]["ws_kbt"]]
    if len(matches) != 1:
        raise ValueError("Need exactly one converged neutral 24-cell transferred-fit result")
    ratio = float(matches[0]["flux_ratio_to_same_sigma_zero_energy"])
    if not math.isfinite(ratio) or ratio <= 0:
        raise ValueError("Transferred-fit flux ratio must be finite and positive")
    return ratio


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--representative-json", required=True, type=Path)
    parser.add_argument("--fitted-json", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    points = load_points(args.representative_json)
    fitted = load_fitted_ratio(args.fitted_json)

    fig, ax = plt.subplots(figsize=(8.4, 5.2))
    fig.subplots_adjust(left=0.14, right=0.98, top=0.88, bottom=0.27)
    style = ((0.0, "#17658d", "o", "sigma = 0"),
             (-0.005, "#bd4e37", "s", "sigma = -0.005 C/m2"))
    for sigma, color, marker, label in style:
        baseline = points[(sigma, "zero")]
        values = [points[(sigma, scenario)] / baseline for scenario in SCENARIOS]
        offset = -0.09 if sigma == 0 else 0.09
        ax.scatter([i + offset for i in range(3)], values, s=75, color=color,
                   marker=marker, label=label, zorder=3)
    ax.scatter([3], [fitted], s=95, color="#17658d", marker="D",
               label="Transferred fit, sigma = 0", zorder=3)
    ax.annotate(f"{fitted:.3g}", (3, fitted), xytext=(-12, 11),
                textcoords="offset points", ha="right", fontsize=10)
    ax.set_yscale("log")
    ax.set_ylim(1e-6, 2)
    ax.set_xlim(-0.45, 3.5)
    ax.set_xticks(range(4), ("No entry energy", "Rectangle\nw = 3",
                             "Two-sided\nws = 3", "Transferred\nfitted field"))
    ax.set_ylabel("Salt flux / zero-energy flux at same sigma")
    ax.set_xlabel("Prescribed entry-energy case")
    ax.set_title("Representative 2D PNP transport, 24 cells per fiber")
    ax.grid(axis="y", which="major", color="#d7dee2", linewidth=0.8)
    ax.legend(loc="lower left", frameon=False, fontsize=9)
    fig.text(0.5, 0.055,
             "Five-row idealized network; 30/1 mM; transferred ws = 7.04, ell = 894 nm.",
             ha="center", fontsize=9)
    fig.text(0.5, 0.025, "Numerical comparison only; not measured PE flux.",
             ha="center", fontsize=9)
    fig.savefig(args.output, dpi=200)
    plt.close(fig)
    print(f"Saved {args.output}; transferred-fit ratio={fitted:.9g}")


if __name__ == "__main__":
    main()
