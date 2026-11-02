"""Zero-charge, two-reservoir baseline for the frozen Stage 2 B1 scenario.

This is a technical baseline, not the charged PB solver or a material fit.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/processed/salt_transport/22_neutral_baseline_20261028.csv"

AREA_M2 = 1e-4
INITIAL_VOLUME_M3 = 1e-5
SAMPLE_VOLUME_M3 = 2e-7
A_M = 25e-9
TAU = 2.0
D_NA = 1.33e-9
D_CL = 2.03e-9
D_UREA = 1.38e-9
R_UREA_M = 0.18e-9
D_SALT = 2 * D_NA * D_CL / (D_NA + D_CL)
ENTRY_M = math.pi * A_M / (2 * TAU)


@dataclass
class Reservoirs:
    donor_concentration: float
    receiver_concentration: float = 0.0
    donor_volume: float = INITIAL_VOLUME_M3
    receiver_volume: float = INITIAL_VOLUME_M3
    sampled_mass: float = 0.0

    def total_mass(self) -> float:
        return (self.donor_volume * self.donor_concentration
                + self.receiver_volume * self.receiver_concentration
                + self.sampled_mass)

    def advance(self, seconds: float, permeability: float) -> None:
        if seconds < 0 or permeability < 0:
            raise ValueError("Negative time or permeability")
        total_volume = self.donor_volume + self.receiver_volume
        pool_mass = self.total_mass() - self.sampled_mass
        difference = (self.donor_concentration - self.receiver_concentration) * math.exp(
            -AREA_M2 * permeability * (1 / self.donor_volume + 1 / self.receiver_volume) * seconds
        )
        self.donor_concentration = (pool_mass + self.receiver_volume * difference) / total_volume
        self.receiver_concentration = (pool_mass - self.donor_volume * difference) / total_volume

    def sample_both(self) -> None:
        if min(self.donor_volume, self.receiver_volume) <= SAMPLE_VOLUME_M3:
            raise ValueError("Sample exceeds reservoir volume")
        self.sampled_mass += SAMPLE_VOLUME_M3 * (
            self.donor_concentration + self.receiver_concentration
        )
        self.donor_volume -= SAMPLE_VOLUME_M3
        self.receiver_volume -= SAMPLE_VOLUME_M3


def hindrance_urea() -> float:
    ratio = R_UREA_M / A_M
    return (1 - ratio) ** 2 * (1 - 2.104 * ratio + 2.09 * ratio**3 - 0.95 * ratio**5)


def membrane_permeability(h_nm: float, b: float, probe: str) -> float:
    h = h_nm * 1e-9
    if probe == "urea":
        return b * D_UREA / (h / hindrance_urea() + ENTRY_M)
    if probe == "salt_zero_charge":
        return b * D_SALT / (h + ENTRY_M)
    raise ValueError(probe)


def total_permeability(h_nm: float, b: float, probe: str, external: str) -> float:
    membrane = membrane_permeability(h_nm, b, probe)
    if external == "zero":
        resistance = 0.0
    elif external == "reference_2600_neutral":
        diffusivity = D_UREA if probe == "urea" else D_SALT
        resistance = (2600e-9 + ENTRY_M) / (0.01 * diffusivity)
    else:
        raise ValueError(external)
    return 1 / (1 / membrane + resistance)


def run_curve(h_nm: float, b: float, salt_mM: float, probe: str, external: str) -> list[dict[str, object]]:
    initial = 10.0 if probe == "urea" else salt_mM  # 1 mM = 1 mol/m^3
    state = Reservoirs(donor_concentration=initial)
    first_mass = state.total_mass()
    permeability = total_permeability(h_nm, b, probe, external)
    rows = []
    previous_time = 0
    for time_min in (0, 60, 120, 240):
        state.advance((time_min - previous_time) * 60, permeability)
        error = abs(state.total_mass() - first_mass) / first_mass
        if error >= 1e-8:
            raise ArithmeticError(f"Mass conservation error {error}")
        rows.append({
            "h_nm": h_nm, "b": b, "salt_feed_mM": salt_mM, "probe": probe,
            "external": external, "time_min": time_min,
            "donor_mol_m3_before_sample": state.donor_concentration,
            "receiver_mol_m3_before_sample": state.receiver_concentration,
            "permeability_total_m_s": permeability,
            "mass_relative_error": error,
            "status": "MODEL_ONLY_ZERO_CHARGE",
        })
        if time_min in (60, 120, 240):
            old = (state.donor_concentration, state.receiver_concentration)
            state.sample_both()
            assert old == (state.donor_concentration, state.receiver_concentration)
            error = abs(state.total_mass() - first_mass) / first_mass
            if error >= 1e-8:
                raise ArithmeticError(f"Sampling mass conservation error {error}")
        previous_time = time_min
    return rows


def main() -> None:
    if not OUT.parent.is_dir():
        raise FileNotFoundError(OUT.parent)
    rows = []
    for h_nm in (65, 700, 2600):
        for b in (0.0001, 0.01):
            for salt_mM in (1, 10, 600):
                for probe in ("urea", "salt_zero_charge"):
                    for external in ("zero", "reference_2600_neutral"):
                        rows.extend(run_curve(h_nm, b, salt_mM, probe, external))
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"{len(rows)} states; max mass error {max(row['mass_relative_error'] for row in rows):.3g}; {OUT}")


if __name__ == "__main__":
    main()
