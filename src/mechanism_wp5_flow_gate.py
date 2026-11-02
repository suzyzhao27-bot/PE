"""Apply the fixed 10%/30% net-volume-flow attribution branches."""

from __future__ import annotations

import math


def flow_gate(delta_volume_m3: float, maximum_concentration_mol_m3: float,
              observed_transfer_mol: float) -> dict[str, float | str]:
    """delta_volume is corrected for withdrawal and evaporation before calling."""
    values = (delta_volume_m3, maximum_concentration_mol_m3, observed_transfer_mol)
    if any(not math.isfinite(x) for x in values) or maximum_concentration_mol_m3 < 0 or observed_transfer_mol <= 0:
        raise ValueError("Need finite concentration and positive observed transfer")
    upper = abs(delta_volume_m3) * maximum_concentration_mol_m3
    fraction = upper / observed_transfer_mol
    below_tenth = fraction < .1 and not math.isclose(fraction, .1, rel_tol=1e-12)
    at_most_three_tenths = fraction < .3 or math.isclose(fraction, .3, rel_tol=1e-12)
    branch = ("DIFFUSION_MAIN_ANALYSIS" if below_tenth else
              "FLOW_ENVELOPE" if at_most_three_tenths else "NO_PRECISE_DIFFUSION_ATTRIBUTION")
    return {"flow_carried_upper_mol": upper, "flow_fraction_upper": fraction,
            "branch": branch, "assumption": "NET_FLOW_REPRESENTS_TRANSMEMBRANE_VOLUME_FLOW"}


def flow_gate_interval(delta_volume_interval_m3: tuple[float, float],
                       maximum_concentration_interval_mol_m3: tuple[float, float],
                       observed_transfer_interval_mol: tuple[float, float]) -> dict[str, float | str]:
    """Use worst-case positive transfer and flow bounds; never substitute LOQ/2."""
    dv_low, dv_high = delta_volume_interval_m3
    c_low, c_high = maximum_concentration_interval_mol_m3
    t_low, t_high = observed_transfer_interval_mol
    if not all(math.isfinite(x) for x in (dv_low, dv_high, c_low, c_high, t_low, t_high)):
        raise ValueError("Finite measured bounds required")
    if dv_low > dv_high or c_low < 0 or c_low > c_high or t_low < 0 or t_low > t_high:
        raise ValueError("Invalid volume, concentration or transfer interval")
    upper = max(abs(dv_low), abs(dv_high)) * c_high
    if t_low == 0:
        return {"flow_carried_upper_mol": upper,
                "branch": "UNDETERMINED_TRANSFER_INTERVAL_INCLUDES_ZERO",
                "assumption": "NET_FLOW_REPRESENTS_TRANSMEMBRANE_VOLUME_FLOW"}
    conservative = flow_gate(max(abs(dv_low), abs(dv_high)), c_high, t_low)
    return {"flow_carried_upper_mol": upper,
            "flow_fraction_upper": conservative["flow_fraction_upper"],
            "branch": conservative["branch"],
            "assumption": "NET_FLOW_REPRESENTS_TRANSMEMBRANE_VOLUME_FLOW"}


def synthetic_branch_check() -> list[str]:
    """Unit and boundary checks; these values are not observations."""
    ratio = [flow_gate(volume, 1.0, 1e-9)["branch"] for volume in (0.05e-9, 0.1e-9, 0.3e-9, 0.4e-9)]
    expected = ["DIFFUSION_MAIN_ANALYSIS", "FLOW_ENVELOPE", "FLOW_ENVELOPE",
                "NO_PRECISE_DIFFUSION_ATTRIBUTION"]
    if ratio != expected:
        raise ArithmeticError("Fixed branch boundary check failed")
    return ratio


def synthetic_interval_check() -> list[str]:
    results = [flow_gate_interval((v, v), (1.0, 1.0), (1e-9, 1.2e-9))["branch"]
               for v in (0.05e-9, 0.2e-9, 0.4e-9)]
    results.append(flow_gate_interval((0, 0.1e-9), (0, 1.0), (0, 1e-9))["branch"])
    expected = ["DIFFUSION_MAIN_ANALYSIS", "FLOW_ENVELOPE",
                "NO_PRECISE_DIFFUSION_ATTRIBUTION", "UNDETERMINED_TRANSFER_INTERVAL_INCLUDES_ZERO"]
    if results != expected:
        raise ArithmeticError("Interval flow branches changed")
    return results


if __name__ == "__main__":
    print("SYNTHETIC_METHOD_CHECK", synthetic_branch_check(), synthetic_interval_check())
