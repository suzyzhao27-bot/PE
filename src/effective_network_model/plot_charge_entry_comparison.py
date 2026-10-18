"""Plot matched charge/entry-energy PNP cases from the frozen v5 result JSON."""

import argparse
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


SCENARIOS = (
    (0.0, "zero", "Neutral\nno entry energy"),
    (-0.005, "zero", "Charged\nno entry energy"),
    (0.0, "two_sided_ws3_ell26p6", "Neutral\ntwo-sided entry"),
    (-0.005, "two_sided_ws3_ell26p6", "Charged\ntwo-sided entry"),
)


def selected_rows(path):
    report = json.loads(path.read_text(encoding="utf-8"))
    selected = {}
    geometry = None
    for row in report["results"]:
        if row["cells_per_fiber"] != 24:
            continue
        result = row["result"]
        key = (float(result["sigma_c_m2"]), row["scenario"])
        if key not in {(sigma, scenario) for sigma, scenario, _ in SCENARIOS}:
            continue
        if key in selected:
            raise ValueError(f"Duplicate PNP condition: {key}")
        if geometry is None:
            geometry = row["geometry"]
        elif row["geometry"] != geometry:
            raise ValueError("Compared cases must share the same geometry")
        if result["representative_fiber_rows"] != 5 or not result["open_circuit"]:
            raise ValueError("Expected five-row open-circuit calculation")
        if (float(result["c_left_mol_m3"]), float(result["c_right_mol_m3"])) != (30.0, 1.0):
            raise ValueError("Expected common 30/1 mM bath concentrations")
        field = result["common_ion_energy_model"]
        if row["scenario"] == "zero":
            if field is not None:
                raise ValueError("Zero-energy case has an energy field")
        elif (field is None or field["kind"] != "two_sided_exponential"
              or float(field["amplitude_kbt"]) != 3.0 or float(field["decay_nm"]) != 26.6):
            raise ValueError("Expected the prescribed two-sided shared-ion field")
        for name in ("common_salt_flux_mol_m2_s", "right_minus_left_theoretical_aqueous_voltage_mv",
                     "plus_flux_right_mol_m2_s", "minus_flux_right_mol_m2_s"):
            if not math.isfinite(float(result[name])):
                raise ValueError(f"Non-finite {name}: {key}")
        plus = float(result["plus_flux_right_mol_m2_s"])
        minus = float(result["minus_flux_right_mol_m2_s"])
        salt = float(result["common_salt_flux_mol_m2_s"])
        if salt <= 0 or not math.isclose(plus, minus, rel_tol=1e-9) or not math.isclose(
            salt, (plus + minus) / 2, rel_tol=1e-9
        ):
            raise ValueError(f"Open-circuit species fluxes inconsistent: {key}")
        selected[key] = result
    if len(selected) != len(SCENARIOS):
        raise ValueError("Need the four matched 24-cell charge/entry-energy cases")
    return selected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-json", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    selected = selected_rows(args.results_json)
    reference = float(selected[(0.0, "zero")]["common_salt_flux_mol_m2_s"])
    labels = [label for _, _, label in SCENARIOS]
    voltages = [float(selected[(sigma, scenario)]["right_minus_left_theoretical_aqueous_voltage_mv"])
                for sigma, scenario, _ in SCENARIOS]
    fluxes = [float(selected[(sigma, scenario)]["common_salt_flux_mol_m2_s"]) / reference
              for sigma, scenario, _ in SCENARIOS]
    colors = ["#65747f", "#b25c43", "#5391ac", "#17658d"]

    fig, (ax_voltage, ax_flux) = plt.subplots(1, 2, figsize=(10.2, 5.4))
    fig.subplots_adjust(left=0.09, right=0.98, top=0.79, bottom=0.28, wspace=0.27)
    for ax, values in ((ax_voltage, voltages), (ax_flux, fluxes)):
        ax.bar(range(4), values, color=colors, width=0.72)
        ax.set_xticks(range(4), labels, fontsize=8)
        ax.grid(axis="y", alpha=0.18)
        ax.set_axisbelow(True)
    ax_voltage.axhline(0, color="#47515a", linewidth=0.8)
    ax_voltage.set_ylabel("Open-circuit voltage, right minus left (mV)")
    ax_flux.set_ylabel("Common salt flux / neutral zero-energy flux")
    ax_flux.set_ylim(0, 1.25)
    for ax, values in ((ax_voltage, voltages), (ax_flux, fluxes)):
        for i, value in enumerate(values):
            ax.annotate(f"{value:.3g}", (i, value), xytext=(0, 4), textcoords="offset points",
                        ha="center", va="bottom", fontsize=8)
    fig.text(0.09, 0.91, "Supplement | Charge response and common salt transport",
             fontsize=14, fontweight="bold")
    fig.text(0.09, 0.84, "Same five-row virtual network; 30/1 mM NaCl; 24 cells/fiber; σ = −0.005 C/m² when charged.",
             fontsize=9, color="#454b52")
    fig.text(0.09, 0.10, "Two-sided shared-ion entry field: ws = 3 kBT, ell = 26.6 nm. Voltage is a charge-sensitive response,"
             "\nnot a measured ion-selectivity coefficient. Open circuit gives J+ ≈ J−; fluxes are virtual, not PE measurements.",
             fontsize=8, color="#454b52")
    fig.savefig(args.output, dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
