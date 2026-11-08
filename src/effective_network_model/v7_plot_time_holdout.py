"""Plot frozen B time-holdout observations and saved N0/NS/ND predictions."""

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


def main():
    rows = []
    for model in MODELS:
        path = ROOT / "WP4" / f"v7_{model}_time_points.csv"
        with path.open(encoding="utf-8-sig", newline="") as stream:
            selected = list(csv.DictReader(stream))
        if len(selected) != 42:
            raise RuntimeError(f"Unexpected {model} time point count: {len(selected)}")
        rows.extend(selected)
    if len(rows) != 126:
        raise RuntimeError("Time comparison requires all three completed models")
    by_key = {}
    for row in rows:
        key = (int(row["h_nm"]), float(row["time_min"]), row["model"])
        if key in by_key:
            raise RuntimeError(f"Duplicate saved prediction: {key}")
        by_key[key] = row
    for h in THICKNESSES:
        reference = [r for r in rows if r["model"] == "N0" and int(r["h_nm"]) == h]
        for base in reference:
            for model in MODELS[1:]:
                other = by_key[(h, float(base["time_min"]), model)]
                for field in ("role", "observed_increment_uS_cm", "baseline_uS_cm"):
                    if base[field] != other[field]:
                        raise RuntimeError(f"Model tables disagree on {field}: h={h}")

    source = FIG / "v7_B_time_holdout_panel.csv"
    with source.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "svg.fonttype": "none"})
    fig, axes = plt.subplots(2, 3, figsize=(12.6, 6.9), sharex=True)
    for ax, h in zip(axes.flat, THICKNESSES):
        reference = [r for r in rows if r["model"] == "N0" and int(r["h_nm"]) == h]
        reference.sort(key=lambda r: float(r["time_min"]))
        for model in MODELS:
            points = [by_key[(h, float(r["time_min"]), model)] for r in reference]
            ax.plot([float(r["time_min"]) for r in points],
                    [float(r["predicted_increment_uS_cm"]) for r in points],
                    color=COLORS[model], linewidth=1.7,
                    linestyle="--" if model == "N0" else "-", label=model)
        for role, marker, face in (("train", "o", "black"),
                                   ("holdout", "o", "white")):
            selected = [r for r in reference if r["role"] == role]
            ax.scatter([float(r["time_min"]) for r in selected],
                       [float(r["observed_increment_uS_cm"]) for r in selected],
                       marker=marker, s=30, facecolors=face, edgecolors="black",
                       linewidths=1, zorder=5, label=f"Observed {role}")
        ax.axvline(60, color="#999999", linestyle=":", linewidth=1)
        ax.set_title(f"{h} nm")
        ax.set_xlim(0, 125)
        ax.set_ylim(bottom=0)
        ax.grid(axis="y", color="#e4e4e4", linewidth=0.6)
    for ax in axes[1]:
        ax.set_xlabel("Time (min)")
    for ax in axes[:, 0]:
        ax.set_ylabel("Conductivity increment (µS/cm)")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.suptitle("B series: train at 20–60 min, predict after 60 min", y=0.985)
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.935),
               ncol=5, frameon=False)
    fig.subplots_adjust(left=0.07, right=0.985, top=0.84, bottom=0.09,
                        wspace=0.19, hspace=0.23)
    for suffix in ("svg", "png"):
        fig.savefig(FIG / f"v7_B_time_holdout_panel.{suffix}", dpi=220,
                    bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {source.name} (126 rows) and SVG/PNG")


if __name__ == "__main__":
    main()
