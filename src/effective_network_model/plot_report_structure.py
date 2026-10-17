"""Plot conditional thickness and local-scale scans for report Figure 2."""

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


CORE = Path(__file__).resolve().parent / "core"
sys.path.insert(0, str(CORE))
from theory_model import log_resistance_nm  # noqa: E402


def read_rows(path):
    with path.open(encoding="utf-8-sig", newline="") as source:
        return list(csv.DictReader(source))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--thickness-csv", required=True, type=Path)
    parser.add_argument("--local-csv", required=True, type=Path)
    parser.add_argument("--parameters-json", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    params = json.loads(args.parameters_json.read_text(encoding="utf-8"))
    ws, ell, aref = (float(params[key]) for key in ("ws_ref", "ell_nm", "a_ref_nm"))
    if not all(math.isfinite(value) and value > 0 for value in (ws, ell, aref)):
        raise ValueError("Expected positive frozen model parameters")
    if float(params["robin_r"]) != 0 or aref != 25:
        raise ValueError("This figure is limited to the frozen r=0, a_ref=25 nm fit")

    thickness = read_rows(args.thickness_csv)
    thick_points = [(float(row["h_nm"]), float(row["relative_P_to_2600"]))
                    for row in thickness]
    if len(thick_points) < 2 or any(not math.isfinite(h) or not math.isfinite(p) or h <= 0 or p <= 0
                                    for h, p in thick_points):
        raise ValueError("Thickness scan must contain finite positive points")
    if any(next_h <= h for (h, _), (next_h, _) in zip(thick_points, thick_points[1:])):
        raise ValueError("Thickness scan labels must increase")
    thick_points = [(h, p) for h, p in thick_points if 65 <= h <= 2600]
    if len(thick_points) < 2:
        raise ValueError("Thickness scan does not cover the old label interval")

    local = read_rows(args.local_csv)
    selected = {float(row["a_nm"]): float(row["relative_P_to_reference"])
                for row in local if float(row["h_nm"]) == 65 and float(row["a_nm"]) in (20, 25, 30)}
    if set(selected) != {20, 25, 30} or not all(math.isfinite(value) and value > 0
                                                 for value in selected.values()):
        raise ValueError("Expected positive 65 nm conditional scan at a=20/25/30 nm")
    radius = (20, 25, 30)
    conditional = [selected[a] / selected[25] for a in radius]
    fixed = [math.exp(log_resistance_nm(65, ws, ell, aref) - log_resistance_nm(65, ws, ell, a))
             for a in radius]
    if not math.isclose(conditional[1], fixed[1], rel_tol=1e-12):
        raise ValueError("The two scans must share the a=25 nm reference")

    fig, (left, right) = plt.subplots(1, 2, figsize=(10.7, 4.7))
    fig.subplots_adjust(left=0.09, right=0.98, bottom=0.25, top=0.83, wspace=0.29)
    left.plot([h for h, _ in thick_points], [p for _, p in thick_points], color="#17658d",
              linewidth=2.2)
    left.set_xscale("log")
    left.set_yscale("log")
    left.set_xticks((65, 120, 200, 350, 700, 2600), ("65", "120", "200", "350", "700", "2600"))
    left.set_xlim(60, 2800)
    left.set_xlabel("Thickness label h (nm)")
    left.set_ylabel("Model P(h, 25 nm) / P(2600 nm, 25 nm)")
    left.set_title("A  Thickness at fixed local scale", loc="left", fontweight="bold")
    left.grid(which="major", alpha=0.22)

    right.plot(radius, conditional, color="#c65436", marker="o", linewidth=2,
               label="Conditional ws(a) ∝ 1/a")
    right.plot(radius, fixed, color="#17658d", marker="s", linewidth=2,
               label="Fixed ws; access length only")
    right.axhline(1, color="#777e86", linestyle="--", linewidth=1)
    right.set_yscale("log")
    right.set_xticks(radius)
    right.set_xlabel("Assumed local scale a (nm), h = 65 nm")
    right.set_ylabel("Model P(65 nm, a) / P(65 nm, 25 nm)")
    right.set_title("B  Closure rule controls the response", loc="left", fontweight="bold")
    right.grid(which="major", alpha=0.22)
    right.legend(loc="upper left", fontsize=8, frameon=False)

    fig.text(0.09, 0.92, "Conditional geometry and entry-field sensitivity", fontsize=14,
             fontweight="bold")
    fig.text(0.09, 0.08, "Frozen r=0 shared fit; all values are model ratios."
             "  a and ws(a) are hypotheses, not measured wet-pore properties.",
             fontsize=8, color="#454b52")
    fig.savefig(args.output, dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
