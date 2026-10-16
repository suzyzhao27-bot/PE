from __future__ import annotations

from collections.abc import Iterable


def update_volume(
    initial_volume: float,
    sampled_volume: float = 0.0,
    replenished_volume: float = 0.0,
) -> float:
    """Return reservoir volume after sampling and replenishment."""
    values = (initial_volume, sampled_volume, replenished_volume)
    if any(value < 0 for value in values):
        raise ValueError("Volumes must be non-negative")
    final_volume = initial_volume - sampled_volume + replenished_volume
    if final_volume <= 0:
        raise ValueError("Final reservoir volume must be positive")
    return final_volume


def solution_amount(concentration: float, volume: float) -> float:
    """Return amount as concentration times volume in consistent units."""
    if concentration < 0 or volume < 0:
        raise ValueError("Concentration and volume must be non-negative")
    return concentration * volume


def event_amount(events: Iterable[tuple[float, float]]) -> float:
    """Sum concentration-volume pairs for samples or replenishments."""
    return sum(solution_amount(concentration, volume) for concentration, volume in events)


def corrected_system_amount(
    left_concentration: float,
    left_volume: float,
    right_concentration: float,
    right_volume: float,
    removed: Iterable[tuple[float, float]] = (),
    replenished: Iterable[tuple[float, float]] = (),
) -> float:
    """Return current in-cell amount plus removed amount minus added amount."""
    in_cells = solution_amount(left_concentration, left_volume) + solution_amount(
        right_concentration, right_volume
    )
    return in_cells + event_amount(removed) - event_amount(replenished)


def relative_closure_error(current_corrected_amount: float, initial_amount: float) -> float:
    """Return relative salt-balance error; zero means exact closure."""
    if initial_amount <= 0:
        raise ValueError("Initial amount must be positive")
    return (current_corrected_amount - initial_amount) / initial_amount


def interval_flux(
    corrected_amount_change: float,
    effective_area: float,
    elapsed_time: float,
) -> float:
    """Return interval-average flux in units consistent with inputs."""
    if effective_area <= 0 or elapsed_time <= 0:
        raise ValueError("Effective area and elapsed time must be positive")
    return corrected_amount_change / (effective_area * elapsed_time)
