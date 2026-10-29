"""Reconstruct two-pool solute inventory including every removed sample."""

from __future__ import annotations

import math


def amount(volume_m3: float, concentration_mol_m3: float) -> float:
    if not all(math.isfinite(x) and x >= 0 for x in (volume_m3, concentration_mol_m3)):
        raise ValueError("Volumes and concentrations must be finite and nonnegative")
    return volume_m3 * concentration_mol_m3


def recovery_percent(initial_pools: list[tuple[float, float]],
                     final_pools: list[tuple[float, float]],
                     removed_samples: list[tuple[float, float]]) -> float:
    """Each tuple is (measured volume m3, measured concentration mol/m3).

    final_pools must describe both pools after the 240-minute withdrawal;
    removed_samples includes that withdrawal exactly once.
    """
    if len(initial_pools) != 2 or len(final_pools) != 2 or len(removed_samples) != 6:
        raise ValueError("Require two pools and six measured 60/120/240 min samples")
    initial = sum(amount(*pool) for pool in initial_pools)
    if initial <= 0:
        raise ValueError("Initial analyte inventory must be positive")
    final = sum(amount(*pool) for pool in final_pools)
    samples = sum(amount(*sample) for sample in removed_samples)
    return 100 * (final + samples) / initial


def fixed_protocol_illustration() -> dict[str, float | str]:
    """Algebraic no-loss example, never an experimental observation."""
    initial = [(1e-5, 1.0), (1e-5, 0.0)]
    final = [(9.4e-6, 0.5), (9.4e-6, 0.5)]
    samples = [(0.2e-6, 0.5)] * 6
    full = recovery_percent(initial, final, samples)
    omitted = 100 * sum(amount(*pool) for pool in final) / sum(amount(*pool) for pool in initial)
    if not math.isclose(full, 100, abs_tol=1e-12) or not math.isclose(omitted, 94, abs_tol=1e-12):
        raise ArithmeticError("Fixed protocol inventory identity failed")
    return {"status": "ALGEBRAIC_ILLUSTRATION_ONLY", "with_all_samples_percent": full,
            "if_samples_omitted_percent": omitted,
            "sample_count": 6, "sample_each_side_each_time_mL": 0.2}


if __name__ == "__main__":
    print(fixed_protocol_illustration())
