"""Reconstruct two-pool solute inventory including every removed sample."""

from __future__ import annotations

import math


EXPECTED_SAMPLE_KEYS = {(minute, side) for minute in (60, 120, 240) for side in ("D", "R")}


def validate_sample_events(events: list[dict[str, object]]) -> None:
    """Require each nominal time/side once and retain ordered actual times."""
    if len(events) != 6:
        raise ValueError("Require six fixed sample events")
    keys = {(event.get("nominal_time_min"), event.get("side")) for event in events}
    if keys != EXPECTED_SAMPLE_KEYS:
        raise ValueError("Sample events must cover nominal 60/120/240 min x D/R exactly once")
    for side in ("D", "R"):
        series = sorted((event for event in events if event["side"] == side),
                        key=lambda event: event["nominal_time_min"])
        try:
            actual = [float(event["actual_time_min"]) for event in series]
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError("Measured actual sample times are required") from exc
        if any(not math.isfinite(t) or t <= 0 for t in actual) or not actual[0] < actual[1] < actual[2]:
            raise ValueError("Actual sample times must be finite and increasing on each side")


def amount(volume_m3: float, concentration_mol_m3: float) -> float:
    if not all(math.isfinite(x) and x >= 0 for x in (volume_m3, concentration_mol_m3)):
        raise ValueError("Volumes and concentrations must be finite and nonnegative")
    return volume_m3 * concentration_mol_m3


def recovery_percent(initial_pools: list[tuple[float, float]],
                     final_pools: list[tuple[float, float]],
                     removed_samples: list[dict[str, object]]) -> float:
    """Each tuple is (measured volume m3, measured concentration mol/m3).

    final_pools must describe both pools after the 240-minute withdrawal;
    removed_samples includes that withdrawal exactly once. Each event records
    nominal_time_min, actual_time_min, side, volume_m3 and concentration_mol_m3.
    """
    if len(initial_pools) != 2 or len(final_pools) != 2:
        raise ValueError("Require two pools")
    validate_sample_events(removed_samples)
    initial = sum(amount(*pool) for pool in initial_pools)
    if initial <= 0:
        raise ValueError("Initial analyte inventory must be positive")
    final = sum(amount(*pool) for pool in final_pools)
    samples = sum(amount(float(sample["volume_m3"]), float(sample["concentration_mol_m3"]))
                  for sample in removed_samples)
    return 100 * (final + samples) / initial


def fixed_protocol_illustration() -> dict[str, float | str]:
    """Algebraic no-loss example, never an experimental observation."""
    initial = [(1e-5, 1.0), (1e-5, 0.0)]
    final = [(9.4e-6, 0.5), (9.4e-6, 0.5)]
    samples = [{"nominal_time_min": minute, "actual_time_min": minute, "side": side,
                "volume_m3": 0.2e-6, "concentration_mol_m3": 0.5}
               for minute, side in sorted(EXPECTED_SAMPLE_KEYS)]
    full = recovery_percent(initial, final, samples)
    omitted = 100 * sum(amount(*pool) for pool in final) / sum(amount(*pool) for pool in initial)
    if not math.isclose(full, 100, abs_tol=1e-12) or not math.isclose(omitted, 94, abs_tol=1e-12):
        raise ArithmeticError("Fixed protocol inventory identity failed")
    return {"status": "ALGEBRAIC_ILLUSTRATION_ONLY", "with_all_samples_percent": full,
            "if_samples_omitted_percent": omitted,
            "sample_count": 6, "sample_each_side_each_time_mL": 0.2}


if __name__ == "__main__":
    print(fixed_protocol_illustration())
