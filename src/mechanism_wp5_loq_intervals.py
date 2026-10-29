"""Preserve [0, LOQ] censoring when reconstructing WP5 inventory."""

from __future__ import annotations

import math

from mechanism_wp5_mass_inventory import EXPECTED_SAMPLE_KEYS, validate_sample_events


Interval = tuple[float, float]


def concentration_interval(measured_mol_m3: float | None, loq_mol_m3: float,
                           below_loq: bool) -> Interval:
    if not math.isfinite(loq_mol_m3) or loq_mol_m3 <= 0:
        raise ValueError("Measured positive LOQ is required")
    if below_loq:
        if measured_mol_m3 is not None:
            raise ValueError("Below-LOQ result must not contain a substituted point value")
        return 0.0, loq_mol_m3
    if measured_mol_m3 is None or not math.isfinite(measured_mol_m3) or measured_mol_m3 < loq_mol_m3:
        raise ValueError("Quantified concentration must be finite and at least LOQ")
    return measured_mol_m3, measured_mol_m3


def inventory_interval(items: list[tuple[float, Interval]]) -> Interval:
    low = high = 0.0
    for volume_m3, (c_low, c_high) in items:
        if not all(math.isfinite(x) for x in (volume_m3, c_low, c_high)):
            raise ValueError("Nonfinite inventory item")
        if volume_m3 < 0 or c_low < 0 or c_high < c_low:
            raise ValueError("Invalid volume or concentration interval")
        low += volume_m3 * c_low
        high += volume_m3 * c_high
    return low, high


def recovery_interval(initial_pools: list[tuple[float, Interval]],
                      final_pools: list[tuple[float, Interval]],
                      removed_samples: list[dict[str, object]]) -> dict[str, float | str]:
    """Use post-240-minute-withdrawal final pools and all six samples exactly once."""
    if len(initial_pools) != 2 or len(final_pools) != 2:
        raise ValueError("Require two pools")
    validate_sample_events(removed_samples)
    initial_low, initial_high = inventory_interval(initial_pools)
    sample_intervals = [(float(s["volume_m3"]), s["concentration_interval_mol_m3"])
                        for s in removed_samples]
    recovered_low, recovered_high = inventory_interval(final_pools + sample_intervals)
    if initial_high <= 0:
        raise ValueError("Initial analyte inventory cannot be zero")
    if initial_low == 0:
        return {"status": "UNDETERMINED_INITIAL_INVENTORY_INCLUDES_ZERO",
                "recovery_lower_percent": 100 * recovered_low / initial_high,
                "recovery_upper_percent": math.inf}
    lower = 100 * recovered_low / initial_high
    upper = 100 * recovered_high / initial_low
    status = ("WHOLE_INTERVAL_WITHIN_90_110" if lower >= 90 and upper <= 110 else
              "WHOLE_INTERVAL_OUTSIDE_90_110" if upper < 90 or lower > 110 else
              "UNDETERMINED_INTERVAL_CROSSES_GATE")
    return {"status": status, "recovery_lower_percent": lower, "recovery_upper_percent": upper}


def synthetic_method_check() -> tuple[str, str]:
    censored = concentration_interval(None, 0.1, True)
    if censored != (0.0, 0.1):
        raise ArithmeticError("LOQ censoring identity failed")
    base = [(1e-5, (1.0, 1.0)), (1e-5, (0.0, 0.0))]
    final = [(9.4e-6, (0.5, 0.5))] * 2
    samples = [{"nominal_time_min": minute, "actual_time_min": minute,
                "side": side, "volume_m3": 0.2e-6,
                "concentration_interval_mol_m3": (0.5, 0.5)}
               for minute, side in sorted(EXPECTED_SAMPLE_KEYS)]
    certain = recovery_interval(base, final, samples)["status"]
    uncertain = recovery_interval(base, [(9.4e-6, (0.0, 1.0))] * 2, samples)["status"]
    if certain != "WHOLE_INTERVAL_WITHIN_90_110" or uncertain != "UNDETERMINED_INTERVAL_CROSSES_GATE":
        raise ArithmeticError("Recovery interval gate failed")
    return certain, uncertain


if __name__ == "__main__":
    print("SYNTHETIC_METHOD_CHECK", synthetic_method_check())
