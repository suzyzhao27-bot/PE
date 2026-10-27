"""Plot frozen full-B training curves and original-unit residuals."""

from __future__ import annotations

import csv
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FIG = ROOT / "data/processed/figures"
assert FIG.is_dir()
os.environ.setdefault("MPLCONFIGDIR", str(FIG))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


MODELS = ("N0", "NS", "ND")
THICKNESSES = (65, 120, 200, 350, 700, 2600)
COLORS = {"N0": "#777777", "NS": "#2767a8", "ND": "#d17925"}
STYLES = {"N0": "--", "NS": "-", "ND": "-."}


def main():
    rows = []
    for model in MODELS:
        path = ROOT / "WP4" / f"v7_{model}_fullB_points.csv"
        with path.open(encoding="utf-8-sig", newline="") as stream:
            selected = list(csv.DictReader(stream))
        if len(selected) != 42 or any(row["model"] != model or
                                       row["series"] != "B" for row in selected):
            raise RuntimeError(f"Unexpected {model} full-B points")
        for row in selected:
            row["residual_pred_minus_obs_uS_cm"] = (
                float(row["predicted_increment_uS_cm"]) -
                float(row["observed_increment_uS_cm"]))
        rows.extend(selected)
    by_key = {}
    for row in rows:
        key = (int(row["h_nm"]), float(row["time_min"]), row["model"])
        if key in by_key:
            raise RuntimeError(f"Duplicate full-B point: {key}")
        by_key[key] = row
    if len(by_key) != 126:
        raise RuntimeError("Full-B comparison requires 126 model-point rows")
    for h in THICKNESSES:
        reference = [r for r in rows if r["model"] == "N0" and
                     int(r["h_nm"]) == h]
        if len(reference) != 7 or sum(r["in_training"] == "True"
                                      for r in reference) != 6:
            raise RuntimeError(f"Full-B training role mismatch: {h} nm")
        for base in reference:
            for model in MODELS[1:]:
                other = by_key[(h, float(base["time_min"]), model)]
                for field in ("curve_id", "conductivity_uS_cm", "baseline_uS_cm",
                              "observed_increment_uS_cm", "in_training"):
                    if base[field] != other[field]:
                        raise RuntimeError(f"Full-B tables disagree on {field}: {h} nm")

    source = FIG / "v7_B_full_comparison_panel.csv"
    with source.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "svg.fonttype": "none"})
    for kind in ("curve", "residual"):
        fig, axes = plt.subplots(2, 3, figsize=(12.6, 6.9), sharex=True)
        for ax, h in zip(axes.flat, THICKNESSES):
            reference = sorted((r for r in rows if r["model"] == "N0" and
                                int(r["h_nm"]) == h),
                               key=lambda r: float(r["time_min"]))
            if kind == "curve":
                ax.scatter([float(r["time_min"]) for r in reference],
                           [float(r["observed_increment_uS_cm"]) for r in reference],
                           s=30, color="black", zorder=5, label="Observed B")
            else:
                ax.axhline(0, color="black", linewidth=0.8)
            for model in MODELS:
                selected = [by_key[(h, float(r["time_min"]), model)]
                            for r in reference]
                field = ("predicted_increment_uS_cm" if kind == "curve" else
                         "residual_pred_minus_obs_uS_cm")
                ax.plot([float(r["time_min"]) for r in selected],
                        [float(r[field]) for r in selected],
                        color=COLORS[model], linewidth=1.7,
                        linestyle=STYLES[model], marker="o" if kind == "residual"
                        else None, markersize=3.3, label=model)
            ax.set_title(f"{h} nm")
            ax.set_xlim(0, 125)
            if kind == "curve":
                ax.set_ylim(bottom=0)
            ax.grid(axis="y", color="#e4e4e4", linewidth=0.6)
        for ax in axes[1]:
            ax.set_xlabel("Time (min)")
        for ax in axes[:, 0]:
            ax.set_ylabel("Conductivity increment (µS/cm)" if kind == "curve"
                          else "Prediction − observation (µS/cm)")
        handles, labels = axes.flat[0].get_legend_handles_labels()
        fig.suptitle("B series: full-data training " +
                     ("curves" if kind == "curve" else "residuals"), y=0.985)
        fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.935),
                   ncol=4, frameon=False)
        fig.subplots_adjust(left=0.075, right=0.985, top=0.84, bottom=0.09,
                            wspace=0.21, hspace=0.25)
        for suffix in ("svg", "png"):
            fig.savefig(FIG / f"v7_B_full_{kind}_panel.{suffix}", dpi=220,
                        bbox_inches="tight")
        plt.close(fig)
    print(f"Saved {source.name} (126 rows) and curve/residual SVG/PNG")


if __name__ == "__main__":
    main()
