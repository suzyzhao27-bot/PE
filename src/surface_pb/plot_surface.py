from pathlib import Path
import csv
import hashlib
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from fibre_pb import EPS0, E, R


O = Path(__file__).resolve().parent.parent / "outputs"


def build_dielectric_control_figure():
    # Analytical symmetric-valence control, not an atomistic PMF.
    z = np.linspace(0.2, 2, 160)
    epsw, epspe, temperature = 78.3, 2.3, 298.15
    kb = R / 6.02214076e23
    lb = E**2 / (4 * np.pi * EPS0 * epsw * kb * temperature) * 1e9
    eta = (epsw - epspe) / (epsw + epspe)
    dg = lb * eta / (4 * z)

    figure, axis = plt.subplots(figsize=(6.4, 4.35))
    figure.subplots_adjust(left=0.13, right=0.96, bottom=0.23, top=0.86)
    axis.plot(z, dg, color="#237e8a", label="Both monovalent signs (identical q^2)")
    axis.set(
        xlabel="Distance from planar dielectric boundary / nm",
        ylabel="Image self-energy / kBT",
        title="Dielectric contrast alone has no charge-sign preference",
    )
    axis.legend(frameon=False)
    note = figure.text(
        0.1,
        0.025,
        "Assumed eps_water=78.3, eps_PE=2.3; z << fibre radius; no specific affinity.",
        fontsize=8,
    )
    data = {"z": z, "dg": dg, "lb": lb, "eta": eta}
    return figure, axis, note, data


def main():
    rows = list(csv.DictReader((O / "surface_scenarios.csv").open(encoding="utf-8-sig")))
    profiles = list(csv.DictReader((O / "surface_radial_profiles.csv").open(encoding="utf-8-sig")))

    figure, axes = plt.subplots(1, 3, figsize=(13.2, 4.1), constrained_layout=True)
    colors = ["#975330", "#237e8a", "#6e59a0"]
    for pH, color in zip([3, 7, 10], colors):
        profile = [value for value in profiles if float(value["pH"]) == pH]
        axes[0].plot(
            [float(value["r_nm"]) - 13.3 for value in profile],
            [float(value["psi_mV"]) for value in profile],
            label=f"pH {pH}",
            color=color,
        )
    axes[0].axhline(0, color="gray", lw=0.6)
    axes[0].set(
        xlabel="Distance from solid fibre / nm",
        ylabel="Potential relative to reservoir / mV",
        title="A. Charge arises from declared adsorption",
    )
    axes[0].legend(frameon=False)

    profile = [value for value in profiles if float(value["pH"]) == 7]
    axes[1].semilogy(
        [float(value["r_nm"]) - 13.3 for value in profile],
        [float(value["Na_mol_m3"]) for value in profile],
        label="Na+",
        color="#237e8a",
    )
    axes[1].semilogy(
        [float(value["r_nm"]) - 13.3 for value in profile],
        [float(value["Cl_mol_m3"]) for value in profile],
        label="Cl-",
        color="#bb874b",
    )
    axes[1].axhline(1, color="gray", lw=0.6)
    axes[1].set(
        xlabel="Distance from solid fibre / nm",
        ylabel="c(r) / reservoir concentration",
        title="B. Mobile ions around the fibre, pH 7",
    )
    axes[1].legend(frameon=False)

    selected = [
        row
        for row in rows
        if row["model"] == "regulated_H_OH_hypothesis"
        and float(row["pH"]) == 7
        and float(row["assumed_dG_OH_kj"]) == -50
    ]
    salts = [float(row["salt_mM"]) for row in selected]
    axes[2].semilogx(
        salts,
        [float(row["t_cation"]) for row in selected],
        "o-",
        color="#237e8a",
        label="Regulated fibre cell",
    )
    axes[2].semilogx(
        salts,
        [float(row["same_water_uncharged_t_cation"]) for row in selected],
        "--",
        color="gray",
        label="Uncharged, same water",
    )
    axes[2].set(
        xlabel="Background NaCl / mM",
        ylabel="Cation fraction of linear axial current",
        ylim=(0, 1),
        title="C. Screening weakens transport preference",
    )
    axes[2].legend(frameon=False, fontsize=9)
    figure.suptitle(
        "Hypothetical exterior-fibre cell: a=13.3 nm, b=39.5 nm, Gamma=0.05 / nm^2\n"
        "H/OH binding energies -20/-50 kJ/mol are assumptions, not measured PE parameters",
        fontsize=11,
    )
    figure.savefig(O / "surface_regulation_and_selectivity.png", dpi=220)
    plt.close(figure)

    figure, _, _, control = build_dielectric_control_figure()
    z, dg = control["z"], control["dg"]
    with (O / "dielectric_image_control.csv").open(
        "w", encoding="utf-8-sig", newline=""
    ) as output:
        writer = csv.writer(output)
        writer.writerow(["z_nm", "dG_plus_kBT", "dG_minus_kBT", "K_plus", "K_minus"])
        writer.writerows(zip(z, dg, dg, np.exp(-dg), np.exp(-dg)))
    figure.savefig(O / "dielectric_symmetric_control.png", dpi=220)
    plt.close(figure)

    (O / "plot_provenance.json").write_text(
        json.dumps(
            {
                "surface_data_sha256": hashlib.sha256(
                    (O / "surface_scenarios.csv").read_bytes()
                ).hexdigest(),
                "radial_data_sha256": hashlib.sha256(
                    (O / "surface_radial_profiles.csv").read_bytes()
                ).hexdigest(),
                "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                "image_self_energy_scope": (
                    "analytical planar local approximation, hypothetical dielectric constants, "
                    "exactly sign symmetric; not measured PE PMF"
                ),
                "bjerrum_length_nm": control["lb"],
                "dielectric_factor": control["eta"],
            },
            indent=2,
        )
    )
    print("2 PNG, analytical control CSV written")


if __name__ == "__main__":
    main()
