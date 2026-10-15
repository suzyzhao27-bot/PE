"""Plot exact reviewed CSV scenarios; no hardcoded numerical result values."""
import csv
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs"
SOURCE = OUT / "scenario_summary.csv"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    with SOURCE.open(encoding="utf-8-sig", newline="") as stream:
        source_rows = list(csv.DictReader(stream))
    concentrations = [30.0, 270.0]
    sigmas = [-.01, 0., .01]
    review = []
    for c_left in concentrations:
        for sigma in sigmas:
            matched = [row for row in source_rows
                       if float(row["c_left_virtual_mol_m3"]) == c_left
                       and float(row["c_right_virtual_mol_m3"]) == 1.0
                       and float(row["sigma_assumed_c_m2"]) == sigma]
            if len(matched) != 1:
                raise ValueError("Each concentration/charge combination must match exactly one reviewed scenario")
            row = matched[0]
            flux = float(row["pore_area_salt_flux_mol_m2_s"])
            neutral = float(row["neutral_pore_area_salt_flux_mol_m2_s"])
            ratio = float(row["charged_to_neutral_salt_flux_ratio"])
            if abs(ratio - flux / neutral) > 1e-12:
                raise ValueError("Ratio does not reconcile to source flux and neutral denominator")
            review.append(dict(c_left_virtual_mol_m3=c_left, c_right_virtual_mol_m3=1.0,
                               sigma_assumed_c_m2=sigma, pore_area_salt_flux_mol_m2_s=flux,
                               neutral_pore_area_salt_flux_mol_m2_s=neutral,
                               charged_to_neutral_salt_flux_ratio=ratio,
                               radius_assumed_m=float(row["radius_assumed_m"]),
                               length_normalization_m=float(row["length_normalization_m"]),
                               source_status=row["status"]))
    plotdata = OUT / "charge_flux_plotdata.csv"
    with plotdata.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(review[0]))
        writer.writeheader()
        writer.writerows(review)

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "axes.labelsize": 12, "axes.titlesize": 14,
                         "axes.linewidth": .8, "savefig.facecolor": "white"})
    fig, ax = plt.subplots(figsize=(8.6, 5.8))
    fig.subplots_adjust(left=.12, right=.97, bottom=.31, top=.79)
    colors = {-0.01: "#2c6d9e", 0.0: "#8a8a8a", 0.01: "#b48427"}
    width = .21
    for i, sigma in enumerate(sigmas):
        x = [j + (i - 1) * width for j in range(2)]
        values = [next(row["charged_to_neutral_salt_flux_ratio"] for row in review
                       if row["c_left_virtual_mol_m3"] == c and row["sigma_assumed_c_m2"] == sigma)
                  for c in concentrations]
        bars = ax.bar(x, values, width=width*.9, color=colors[sigma], edgecolor="#333333", linewidth=.65, zorder=3)
        for bar, value in zip(bars, values):
            below_reference = value < .98
            ax.text(bar.get_x() + bar.get_width()/2, value - .025 if below_reference else value + .018,
                    f"{value:.5f}", ha="center", va="top" if below_reference else "bottom",
                    color="white" if below_reference else "#222222", fontsize=10,
                    fontweight="bold" if below_reference else "normal")
    ax.axhline(1., color="#222222", lw=1.2, ls="--", zorder=4)
    ax.set_ylim(0, 1.18)  # Absolute ratio bars start at zero; do not exaggerate small changes.
    ax.set_yticks([0, .2, .4, .6, .8, 1.0])
    ax.set_ylabel("Internal salt flux ratio, J/J0")
    ax.set_xticks([0, 1], ["30 to 1 mM", "270 to 1 mM"])
    ax.set_xlabel("Fixed virtual end concentrations (left to right)", labelpad=9)
    ax.set_xlim(-.48, 1.48)
    ax.grid(axis="y", color="#e4e4e4", linewidth=.7, zorder=0)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    handles = [Patch(facecolor=colors[s], edgecolor="#333333", label=f"sigma = {s:+.2f} C/m2" if s else "sigma = 0 C/m2") for s in sigmas]
    handles.append(Line2D([0], [0], color="#222222", ls="--", lw=1.2, label="Neutral reference = 1"))
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.54, .885), ncol=2,
               frameon=False, fontsize=10, handlelength=2, columnspacing=2)
    fig.text(.12, .975, "Charge-sign comparison in a restricted long-tube model", fontsize=14, va="top")
    fig.text(.12, .935, "Hypothetical long-tube PB-NP, not PE measurements", fontsize=12,
             fontweight="bold", va="top")
    caption = [
        "a = 25 nm; L = 350 nm is flux normalization, not long-tube validity certification.",
        "sigma = -0.01, 0, +0.01 C/m2; local fluid velocity u(r,x) = 0; zero net ionic current.",
        "Virtual low end = 1 mM. DI not calculated; no pore entrances, external layers or electrode model.",
        "Deterministic scenarios; no experimental uncertainty bars. Not finite-length PNP certification."
    ]
    for i, line in enumerate(caption):
        fig.text(.12, .19-i*.031, line, fontsize=9, va="top", color="#333333")
    output = OUT / "charge_flux_diagnostic.png"
    fig.savefig(output, dpi=300)
    plt.close(fig)
    report = dict(input_csv=str(SOURCE), input_sha256=sha(SOURCE),
                  script=str(Path(__file__).resolve()), script_sha256=sha(Path(__file__)),
                  diagnostic_code_sha256={name: sha(ROOT / name) for name in
                                          ("pb_v1_readonly_copy.py", "pb_source_identity.json", "long_tube_np.py", "run_diagnostic.py", "test_diagnostic.py")},
                  diagnostic_run_manifest_sha256=sha(OUT / "run_manifest.json"),
                  diagnostic_test_report_sha256=sha(ROOT / "test_report.json"),
                  output_png=str(output), output_sha256=sha(output),
                  plotdata_csv=str(plotdata), plotdata_sha256=sha(plotdata), exact_source_rows=review,
                  plotted_rows=6, duplicates_or_missing=False, ratio_reconciled=True,
                  y_axis_limits=[0,1.18], neutral_reference=1.,
                  uncertainty="No experiment or statistical bounds; numerical convergence recorded separately in test_report.json",
                  caption=caption, visual_qa="pending image inspection", format="PNG only, 300 dpi")
    (OUT / "图数据校验.json").write_text(json.dumps(report, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
    print(json.dumps(dict(output=str(output), rows=6, sha256=report["output_sha256"]), ensure_ascii=False))


if __name__ == "__main__":
    main()
