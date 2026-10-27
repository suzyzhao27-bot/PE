"""Plot frozen L early-fit and late-prediction points without refitting."""

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
THICKNESSES = (200, 350, 700)
COLORS = {"N0": "#777777", "NS": "#2767a8", "ND": "#d17925"}


def main():
    rows = []
    for model in MODELS:
        path = ROOT / "WP4" / f"v7_{model}_L_points.csv"
        with path.open(encoding="utf-8-sig", newline="") as stream:
            selected = list(csv.DictReader(stream))
        if len(selected) != 70 or any(row["model"] != model for row in selected):
            raise RuntimeError(f"Unexpected {model} L point table")
        rows.extend(selected)
    by_key = {}
    for row in rows:
        key = (int(row["h_nm"]), float(row["time_min"]), row["model"])
        if key in by_key:
            raise RuntimeError(f"Duplicate L saved prediction: {key}")
        by_key[key] = row
    if len(by_key) != 210:
        raise RuntimeError("L comparison requires 210 model-point rows")
    for h in THICKNESSES:
        reference = sorted((r for r in rows if r["model"] == "N0" and
                            int(r["h_nm"]) == h),
                           key=lambda r: float(r["time_min"]))
        for base in reference:
            for model in MODELS[1:]:
                other = by_key[(h, float(base["time_min"]), model)]
                for field in ("role", "curve_id", "observed_increment_uS_cm",
                              "baseline_uS_cm"):
                    if base[field] != other[field]:
                        raise RuntimeError(f"L model tables disagree on {field}: {h} nm")
        if sum(r["role"] == "baseline" for r in reference) != 1 or any(
                (r["role"] == "early_train") !=
                (0 < float(r["time_min"]) <= 240) for r in reference):
            raise RuntimeError(f"L training window mismatch: {h} nm")

    source = FIG / "v7_L_conditional_transfer_panel.csv"
    with source.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "svg.fonttype": "none"})
    fig, axes = plt.subplots(3, 2, figsize=(12.6, 9.2))
    for row_index, h in enumerate(THICKNESSES):
        reference = sorted((r for r in rows if r["model"] == "N0" and
                            int(r["h_nm"]) == h),
                           key=lambda r: float(r["time_min"]))
        for column, roles in enumerate(({"baseline", "early_train"},
                                        {"late_prediction"})):
            ax = axes[row_index, column]
            observed = [r for r in reference if r["role"] in roles]
            for model in MODELS:
                selected = [by_key[(h, float(r["time_min"]), model)]
                            for r in observed]
                ax.plot([float(r["time_min"]) for r in selected],
                        [float(r["predicted_increment_uS_cm"]) for r in selected],
                        color=COLORS[model], linewidth=1.7,
                        linestyle={"N0": "--", "NS": "-", "ND": "-."}[model],
                        label=model)
            ax.scatter([float(r["time_min"]) for r in observed],
                       [float(r["observed_increment_uS_cm"]) for r in observed],
                       marker="o", s=29,
                       facecolors="black" if column == 0 else "white",
                       edgecolors="black", linewidths=1, zorder=5,
                       label="Observed fit" if column == 0 else "Observed prediction")
            ax.set_title(f"{h} nm | {'fit ≤240 min' if column == 0 else 'prediction >240 min'}")
            ax.set_ylim(bottom=0)
            ax.grid(axis="y", color="#e4e4e4", linewidth=0.6)
            if column == 0:
                ax.set_xlim(0, 250)
                ax.set_ylabel("Conductivity increment (µS/cm)")
            else:
                times = [float(r["time_min"]) for r in observed]
                ax.set_xlim(min(times) - 20, max(times) + 40)
            if row_index == 2:
                ax.set_xlabel("Time (min)")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    right_handles, right_labels = axes[0, 1].get_legend_handles_labels()
    fig.suptitle("L series: shared early fit and conditional late transfer", y=0.99)
    fig.legend(handles[:3] + right_handles[-1:] + handles[-1:],
               labels[:3] + right_labels[-1:] + labels[-1:],
               loc="upper center", bbox_to_anchor=(0.5, 0.955), ncol=5,
               frameon=False)
    fig.subplots_adjust(left=0.085, right=0.985, top=0.895, bottom=0.075,
                        wspace=0.22, hspace=0.36)
    for suffix in ("svg", "png"):
        fig.savefig(FIG / f"v7_L_conditional_transfer_panel.{suffix}",
                    dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {source.name} (210 rows) and SVG/PNG")


if __name__ == "__main__":
    main()
