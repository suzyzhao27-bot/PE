"""Render E1 equivalent B and conditional L overlays from audited CSV sources."""

from __future__ import annotations

import csv
import os
from collections import defaultdict

from v7_d1_objective import ROOT

os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "WP5"))
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


SOURCE = ROOT / "WP5"
FIGURES = ROOT / "data/processed/figures"
LABELS = ("fullB_best", "equivalent_1", "equivalent_2", "equivalent_3")
COLORS = ("#1b1b1b", "#0077aa", "#dd5e28", "#5a9b44")


def load(name):
    with (SOURCE / name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def draw(name, h_values, source, title, mark_training):
    rows = load(source)
    grouped = defaultdict(list)
    for r in rows:
        grouped[(r["label"], int(r["h_nm"]))].append(r)
    assert len(grouped) == len(LABELS) * len(h_values)
    fig, axes = plt.subplots(1 if mark_training else 2,
                             3, figsize=(12, 3.8 if mark_training else 7.1),
                             squeeze=False, layout="constrained")
    flat = axes.ravel()
    for j, h in enumerate(h_values):
        ax = flat[j]
        observed = sorted(grouped[(LABELS[0], h)],
                          key=lambda r: int(r["point_index"]))
        t = [float(r["time_min"]) for r in observed]
        y_obs = [float(r["observed_increment_uS_cm"]) for r in observed]
        if mark_training:
            ax.axvspan(0, 240, color="#e8e8e8", alpha=0.8, zorder=0)
            ax.axvline(240, color="#999999", linewidth=0.8)
        ax.scatter(t, y_obs, color="black", facecolor="white",
                   s=20, linewidth=0.8, zorder=5, label="observed")
        for label, color in zip(LABELS, COLORS, strict=True):
            group = sorted(grouped[(label, h)],
                           key=lambda r: int(r["point_index"]))
            assert [float(r["time_min"]) for r in group] == t
            y = [float(r["predicted_increment_uS_cm"]) for r in group]
            ax.plot(t, y, color=color, linewidth=1.6,
                    linestyle="--" if label == "fullB_best" else "-",
                    label=label.replace("_", " "))
        ax.set_title(f"{h} nm")
        ax.set_xlabel("Time (min)")
        ax.set_ylabel("Conductivity increment (μS/cm)")
        ax.grid(alpha=0.2, linewidth=0.5)
    handles, labels = flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=5,
               frameon=False)
    fig.suptitle(title)
    for extension in ("svg", "png"):
        fig.savefig(FIGURES / f"{name}.{extension}", dpi=180)
    plt.close(fig)


def main():
    assert FIGURES.is_dir()
    draw("v7_E1_equivalent_B_curves",
         (65, 120, 200, 350, 700, 2600),
         "v7_E1_equivalent_B_curves.csv",
         "NS: full B and three farthest near-equivalent solutions",
         False)
    draw("v7_E1_equivalent_L_late",
         (200, 350, 700),
         "v7_E1_equivalent_L_points.csv",
         "NS: L early calibration (shaded) and late conditional prediction",
         True)
    print("E1 equivalent figures: B and L overlays; SVG/PNG")


if __name__ == "__main__":
    main()
