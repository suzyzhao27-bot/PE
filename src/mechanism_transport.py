"""Conditional two-reservoir salt/urea transport from the frozen PB table."""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

from mechanism_neutral_baseline import (
    A_M, AREA_M2, D_SALT, D_UREA, ENTRY_M, INITIAL_VOLUME_M3,
    SAMPLE_VOLUME_M3, TAU, membrane_permeability,
)
from mechanism_pb_main_batch import OUT as MAIN_TABLE
from mechanism_pb_target_refinement import OUT as EXTRA_TABLE


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


class PoreTable:
    def __init__(self, a_nm: int = 25, sigma: float = -0.008) -> None:
        rows = _rows(MAIN_TABLE)
        rows += _rows(EXTRA_TABLE)
        points = [(float(r["c_mol_m3"]), float(r["D_s_local_m2_s"])) for r in rows
                  if float(r["a_nm"]) == a_nm and float(r["sigma_C_m2"]) == sigma]
        points.sort()
        if len(points) < 31 or points[0][0] <= 0 or points[-1][0] != 600:
            raise ValueError("Incomplete fixed PB concentration table")
        self.nodes = [(0.0, 0.0)] + points

    def value(self, c: float) -> float:
        if c < 0 or c > self.nodes[-1][0]:
            raise ValueError("Concentration outside PB table")
        for (c0, d0), (c1, d1) in zip(self.nodes, self.nodes[1:]):
            if c <= c1:
                return d0 + (d1 - d0) * (c - c0) / (c1 - c0)
        return self.nodes[-1][1]

    def integral(self, low: float, high: float) -> float:
        if low < 0 or high > self.nodes[-1][0] or high < low:
            raise ValueError("Invalid concentration integral")
        if high == low:
            return 0.0
        knots = [low] + [c for c, _ in self.nodes if low < c < high] + [high]
        total = 0.0
        for c0, c1 in zip(knots, knots[1:]):
            total += (self.value(c0) + self.value(c1)) * (c1 - c0) / 2
        return total


def external_resistance(probe: str, scenario: str) -> float:
    if scenario == "zero":
        return 0.0
    if scenario != "reference_2600_neutral":
        raise ValueError(scenario)
    diffusivity = D_UREA if probe == "urea" else D_SALT
    return (2600e-9 + ENTRY_M) / (0.01 * diffusivity)


def salt_flux(c_donor: float, c_receiver: float, h_nm: float, b: float,
              sigma: float, external: str, table: PoreTable | None = None) -> float:
    if c_donor < c_receiver or c_receiver < 0 or c_donor > 600 or b <= 0:
        raise ValueError("Invalid salt state")
    delta = c_donor - c_receiver
    if delta == 0:
        return 0.0
    if sigma == 0:
        membrane = membrane_permeability(h_nm, b, "salt_zero_charge")
        total = 1 / (1 / membrane + external_resistance("salt", external))
        return total * delta
    if table is None:
        raise ValueError("Charged salt needs fixed PB table")
    epsilon_f = b * TAU
    side = external_resistance("salt", external) / 2 + math.pi * A_M / (4 * epsilon_f * D_SALT)
    upper = delta / (2 * side)
    lower = 0.0
    h = h_nm * 1e-9
    for _ in range(55):
        flux = (lower + upper) / 2
        c_m_d = c_donor - flux * side
        c_m_r = c_receiver + flux * side
        membrane_flux = (b / h) * table.integral(c_m_r, c_m_d)
        if flux > membrane_flux:
            upper = flux
        else:
            lower = flux
    return (lower + upper) / 2


def urea_flux(u_donor: float, u_receiver: float, h_nm: float, b: float,
              external: str) -> float:
    if u_donor < u_receiver or u_receiver < 0:
        raise ValueError("Invalid urea state")
    membrane = membrane_permeability(h_nm, b, "urea")
    total = 1 / (1 / membrane + external_resistance("urea", external))
    return total * (u_donor - u_receiver)


@dataclass
class State:
    donor: float
    receiver: float = 0.0
    volume_donor: float = INITIAL_VOLUME_M3
    volume_receiver: float = INITIAL_VOLUME_M3
    sampled_mass: float = 0.0

    def mass(self) -> float:
        return self.volume_donor * self.donor + self.volume_receiver * self.receiver + self.sampled_mass

    def sample(self) -> None:
        self.sampled_mass += SAMPLE_VOLUME_M3 * (self.donor + self.receiver)
        self.volume_donor -= SAMPLE_VOLUME_M3
        self.volume_receiver -= SAMPLE_VOLUME_M3


def _derivative(donor: float, receiver: float, state: State, probe: str,
                h_nm: float, b: float, sigma: float, external: str,
                table: PoreTable | None) -> tuple[float, float]:
    if donor < receiver and receiver - donor < 1e-12:
        donor = receiver
    if probe == "salt":
        flux = salt_flux(donor, receiver, h_nm, b, sigma, external, table)
    elif probe == "urea":
        flux = urea_flux(donor, receiver, h_nm, b, external)
    else:
        raise ValueError(probe)
    return -AREA_M2 * flux / state.volume_donor, AREA_M2 * flux / state.volume_receiver


def step(state: State, seconds: float, probe: str, h_nm: float, b: float,
         sigma: float, external: str, table: PoreTable | None) -> None:
    d0, r0 = state.donor, state.receiver
    args = (state, probe, h_nm, b, sigma, external, table)
    k1 = _derivative(d0, r0, *args)
    k2 = _derivative(d0 + seconds * k1[0] / 2, r0 + seconds * k1[1] / 2, *args)
    k3 = _derivative(d0 + seconds * k2[0] / 2, r0 + seconds * k2[1] / 2, *args)
    k4 = _derivative(d0 + seconds * k3[0], r0 + seconds * k3[1], *args)
    state.donor += seconds * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6
    state.receiver += seconds * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6
    if min(state.donor, state.receiver) < -1e-12 or state.donor + 1e-12 < state.receiver:
        raise ArithmeticError("Nonphysical reservoir state")


def run_curve(initial: float, probe: str, h_nm: float, b: float,
              sigma: float, external: str, table: PoreTable | None = None,
              step_seconds: float = 60.0) -> list[dict[str, float]]:
    state = State(donor=initial)
    mass_initial = state.mass()
    output = []
    for minute in range(241):
        if minute in (0, 60, 120, 240):
            output.append({"time_min": float(minute), "donor_mol_m3": state.donor,
                           "receiver_mol_m3": state.receiver,
                           "mass_relative_error": abs(state.mass() - mass_initial) / mass_initial})
            if minute in (60, 120, 240):
                state.sample()
                if abs(state.mass() - mass_initial) / mass_initial >= 1e-8:
                    raise ArithmeticError("Sampling mass conservation error")
        if minute < 240:
            if step_seconds != 60:
                parts = int(60 / step_seconds)
                if parts * step_seconds != 60:
                    raise ValueError("Step must divide one minute")
            else:
                parts = 1
            for _ in range(parts):
                step(state, step_seconds, probe, h_nm, b, sigma, external, table)
    if max(row["mass_relative_error"] for row in output) >= 1e-8:
        raise ArithmeticError("Mass conservation error")
    return output
