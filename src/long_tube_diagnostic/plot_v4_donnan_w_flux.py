"""Plot one existing v4 reduced-Donnan scenario without recomputing it."""

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "原始资料/赵珈绮/06_纯计算课题_v4.0.0/results/reduced_partition_scenarios.csv"
OUTPUT = ROOT / "data/processed/figures/v4_donnan_w_flux.png"


def main() -> None:
    with SOURCE.open(newline="", encoding="utf-8-sig") as handle:
        rows = [
            row for row in csv.DictReader(handle)
            if float(row["radius_nm"]) == 25.0
            and float(row["sigma_C_m2"]) == -0.01
            and float(row["donor_M"]) == 0.27
            and float(row["receiver_M"]) == 0.001
        ]
    values = sorted(
        (float(row["common_transfer_energy_kBT"]), float(row["J_over_uncharged_w0"]))
        for row in rows
    )
    assert [w for w, _ in values] == [0.0, 2.0, 4.0, 6.0, 8.0, 12.5]
    expected = {0.0: 1.0020172940320167, 4.0: 0.010286796377280392, 6.0: 0.0002516109941662906}
    assert all(abs(dict(values)[w] - target) < 1e-14 for w, target in expected.items())

    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(figsize=(8.1, 5.3))
    fig.subplots_adjust(left=0.11, right=0.97, bottom=0.24, top=0.80)
    fig.suptitle("Common ion-entry energy suppresses salt transfer", fontsize=14, weight="bold")
    ax.set_title("Reduced Donnan model | a = 25 nm, sigma = -0.01 C/m², donor/receiver = 270/1 mM", fontsize=10)
    w_values, ratios = zip(*values)
    ax.plot(w_values, ratios, color="#1d5d7a", marker="o", linewidth=2, markersize=7)
    ax.axhline(1, color="#888888", linestyle="--", linewidth=1, label="Uncharged, w = 0 reference")
    ax.set_yscale("log")
    ax.set_ylim(1e-10, 4)
    ax.set_xlim(-0.7, 13.2)
    ax.set_xticks(w_values)
    ax.set_xlabel("Common ion-entry energy w (kBT)")
    ax.set_ylabel("Common salt flux J / J0")
    ax.grid(True, which="major", axis="y", alpha=0.25)
    ax.legend(loc="upper right", frameon=True, facecolor="white", edgecolor="none", framealpha=1)
    for w, label, offset in ((0.0, "1.002", (12, -18)), (4.0, "0.01029", (9, 11)), (6.0, "0.0002516", (9, 11))):
        ax.annotate(label, (w, dict(values)[w]), xytext=offset, textcoords="offset points", fontsize=9, weight="bold")
    fig.text(
        0.5, 0.035,
        "Uniform Donnan assumption. Receiver 1 mM is virtual, not the legacy DI experiment.\n"
        "Wide pores and high salt require later spatial validation.",
        ha="center", va="bottom", fontsize=9, color="#444444",
    )
    fig.savefig(OUTPUT, dpi=190)
    plt.close(fig)
    print(f"source={SOURCE}")
    print(f"output={OUTPUT}")
    print(f"rows={len(values)}")
    for w in expected:
        print(f"w={w:g}, J/J0={dict(values)[w]:.15g}")


if __name__ == "__main__":
    main()
