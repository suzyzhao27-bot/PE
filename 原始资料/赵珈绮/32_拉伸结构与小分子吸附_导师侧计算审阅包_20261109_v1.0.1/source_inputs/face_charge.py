"""Flat PE/water face at equilibrium. Run in a local working folder with Python 3."""
import argparse
import csv
import math
from pathlib import Path

R = 8.31446261815324
F = 96485.33212
T = 298.15
EPS = 78.3 * 8.8541878128e-12
E = 1.602176634e-19
VT = R * T / F
SITE_DENSITY = 0.05  # sites per nm2; legacy example, not measured
CAP = E * 1e18 * SITE_DENSITY
PKA = 5.0  # conditional pKa of an already adsorbed monovalent weak acid
SLIP_NM = 0.5  # old comparison convention; not a measured slip plane
H_NM = {"182": 770, "357": 370, "604": 230}


def acid_charge(u, ph):
    alpha = 1 / (1 + 10 ** (PKA - ph) * math.exp(-u))
    return -CAP * alpha


def competitive_charge(u, ph):
    # Effective H(+)/OH(-) association, original v3 example DG=-20/-50 kJ/mol.
    positive = math.exp(20000 / (R*T) - u) * 10 ** (-ph)
    negative = math.exp(50000 / (R*T) + u) * 10 ** (ph-14)
    return CAP * (positive - negative) / (1 + positive + negative)


def plane(charge, ph, ionic_mM, charge_capacity=CAP):
    # ionic_mM is also ionic strength in mol/m3 for the SI formula.
    prefactor = math.sqrt(8 * EPS * R * T * ionic_mM)
    debye_nm = 1e9 * math.sqrt(EPS * R * T / (2 * F**2 * ionic_mM))
    def residual(u):
        return charge(u, ph) - prefactor * math.sinh(u/2)
    lo, hi = -30.0, 30.0
    assert residual(lo) > 0 > residual(hi)
    for _ in range(90):
        mid = (lo + hi) / 2
        if residual(mid) > 0:
            lo = mid
        else:
            hi = mid
    u = (lo + hi) / 2
    sigma = charge(u, ph)
    slip = 4 * math.atanh(math.tanh(u/4) * math.exp(-SLIP_NM / debye_nm))
    assert abs(residual(u)) < 1e-14 and abs(sigma) <= charge_capacity
    assert abs(slip) <= abs(u)
    return 1000*sigma, 1000*VT*u, 1000*VT*slip



def neutral_ph_site_density(row):
    """Infer density from near-neutral measured zeta, keeping chemistry and slip distance fixed."""
    strength = row["conditional_I_mM"]
    debye_nm = 1e9 * math.sqrt(EPS * R * T / (2 * F**2 * strength))
    slip_u = row["measured_zeta_mV"] / (1000 * VT)
    surface_u = 4 * math.atanh(math.tanh(slip_u/4) * math.exp(SLIP_NM / debye_nm))
    required_charge = math.sqrt(8 * EPS * R * T * strength) * math.sinh(surface_u/2)
    return SITE_DENSITY * required_charge / competitive_charge(surface_u, row["actual_pH"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bath_csv", type=Path, help="19 taskbook/inputs/zeta/conditional_bath9.csv")
    parser.add_argument("output_csv", type=Path, help="write to a local output path")
    args = parser.parse_args()
    records = []
    with args.bath_csv.open(encoding="utf-8-sig", newline="") as file:
        for row in csv.DictReader(file):
            ph = float(row["pH"])
            strength = float(row["conditional_I_mM"])
            measured = float(row["zeta_mV"])
            out = dict(draw=row["draw"], h_group_nm=H_NM[row["draw"]], nominal_pH=int(row["pH_group"]),
                       actual_pH=ph, conditional_I_mM=strength, measured_zeta_mV=measured)
            for label, charge in (("acid", acid_charge), ("competitive", competitive_charge)):
                sigma, surface, slip = plane(charge, ph, strength)
                out[label+"_sigma_mC_m2"] = sigma
                out[label+"_surface_mV"] = surface
                out[label+"_at_0p5nm_mV"] = slip
                out[label+"_residual_to_zeta_mV"] = slip - measured
            records.append(out)
    assert len({(r["draw"], r["nominal_pH"]) for r in records}) == len(records) == 9
    assert abs(plane(competitive_charge, 7-30000/(2*R*T*math.log(10)), 1)[0]) < 1e-12
    densities = {r["draw"]: neutral_ph_site_density(r) for r in records if r["nominal_pH"] == 7}
    for r in records:
        density = densities[r["draw"]]
        def local_charge(u, ph):
            return competitive_charge(u, ph) * density / SITE_DENSITY
        sigma, surface, slip = plane(local_charge, r["actual_pH"], r["conditional_I_mM"],
                                     charge_capacity=E * 1e18 * density)
        r["density_from_neutral_pH_per_nm2"] = density
        r["variable_density_role"] = "calibration" if r["nominal_pH"] == 7 else "comparison"
        r["variable_density_sigma_mC_m2"] = sigma
        r["variable_density_surface_mV"] = surface
        r["variable_density_at_0p5nm_mV"] = slip
        r["variable_density_residual_mV"] = slip - r["measured_zeta_mV"]
        if r["nominal_pH"] == 7:
            assert abs(r["variable_density_residual_mV"]) < 1e-9

    records.sort(key=lambda r: (r["nominal_pH"], -r["h_group_nm"]))
    with args.output_csv.open("w", newline="", encoding="utf-8-sig") as file:
        writer = csv.DictWriter(file, fieldnames=records[0].keys())
        writer.writeheader()
        writer.writerows(records)
    for r in records:
        print(r["draw"], r["nominal_pH"], *(round(r[k], 3) for k in (
            "measured_zeta_mV","acid_sigma_mC_m2","acid_surface_mV","acid_at_0p5nm_mV",
            "competitive_sigma_mC_m2","competitive_surface_mV","competitive_at_0p5nm_mV")))


if __name__ == "__main__":
    main()
