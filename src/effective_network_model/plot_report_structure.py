"""Plot the v5 h/a Figure 2 theme or the earlier closure-rule comparison."""

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
from theory_model import field, log_resistance_nm, relative_permeance  # noqa: E402


def read_rows(path):
    with path.open(encoding="utf-8-sig", newline="") as source:
        return list(csv.DictReader(source))


def plot_mentor_layout(params, output):
    """Reproduce the v5 Figure 2 comparison axes using the frozen core."""
    ws, ell, aref = (float(params[key]) for key in ("ws_ref", "ell_nm", "a_ref_nm"))
    colors = {65: "#21618c", 350: "#148f77", 700: "#c0392b", 2600: "#7d3c98"}
    fig, (left, right) = plt.subplots(1, 2, figsize=(10.8, 4.8))
    fig.subplots_adjust(left=0.09, right=0.98, bottom=0.23, top=0.81, wspace=0.32)

    u = [i / 200 for i in range(201)]
    for h in (65, 350, 700, 2600):
        left.plot(u, [field(position * h, h, ws, ell) for position in u],
                  color=colors[h], linewidth=2, label=f"h = {h} nm")
    left.set_xlabel("Position in membrane z/h")
    left.set_ylabel("Common entry energy w / kBT")
    left.set_xlim(0, 1)
    left.set_ylim(bottom=0)
    left.set_title("A  Double-sided field overlap", loc="left", fontweight="bold")
    left.legend(ncol=2, fontsize=8, frameon=False)
    left.grid(alpha=0.22)

    hs = [65 * (2600 / 65) ** (i / 160) for i in range(161)]
    for a, color in ((20, "#c65436"), (25, "#17658d"), (30, "#148f77")):
        ratios = [relative_permeance(h, a, params, reference_h_nm=2600) for h in hs]
        right.plot(hs, ratios, color=color, linewidth=2, label=f"a = {a} nm")
    right.set_xscale("log")
    right.set_yscale("log")
    right.set_xticks((65, 200, 700, 2600), ("65", "200", "700", "2600"))
    right.set_xlim(65, 2600)
    right.set_xlabel("Thickness label h (nm)")
    right.set_ylabel("Model P(h, a) / P(2600 nm, 25 nm)")
    right.set_title("B  Conditional local-scale coupling", loc="left", fontweight="bold")
    right.legend(fontsize=8, frameon=False)
    right.grid(alpha=0.22)

    fig.text(0.09, 0.91, "Thickness and local scale in the frozen v5 model",
             fontsize=14, fontweight="bold")
    fig.text(0.09, 0.07, "Shared r=0 fit; ws(a) = ws(25 nm) * 25 nm / a."
             "  Curves are conditional model results, not wet-pore measurements.",
             fontsize=8, color="#454b52")
    fig.savefig(output, dpi=200)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--thickness-csv", type=Path)
    parser.add_argument("--local-csv", type=Path)
    parser.add_argument("--parameters-json", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--layout", choices=("closure", "mentor-v5"), default="closure")
    args = parser.parse_args()
    params = json.loads(args.parameters_json.read_text(encoding="utf-8"))
    ws, ell, aref = (float(params[key]) for key in ("ws_ref", "ell_nm", "a_ref_nm"))
    if not all(math.isfinite(value) and value > 0 for value in (ws, ell, aref)):
        raise ValueError("Expected positive frozen model parameters")
    if float(params["robin_r"]) != 0 or aref != 25:
        raise ValueError("This figure is limited to the frozen r=0, a_ref=25 nm fit")

    if args.layout == "mentor-v5":
        plot_mentor_layout(params, args.output)
        return
    if args.thickness_csv is None or args.local_csv is None:
        parser.error("The closure layout requires --thickness-csv and --local-csv")

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
