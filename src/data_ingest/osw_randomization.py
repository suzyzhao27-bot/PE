from __future__ import annotations

import random
from collections import Counter


TREATMENTS = ("O", "S", "W")


def assign_treatments(coupon_ids: list[str], seed: int) -> dict[str, str]:
    """Return a reproducible, near-balanced O/S/W assignment for one parent membrane."""
    normalized = [coupon_id.strip() for coupon_id in coupon_ids]
    if len(normalized) < len(TREATMENTS):
        raise ValueError("At least three coupons are required per parent membrane")
    if any(not coupon_id for coupon_id in normalized):
        raise ValueError("coupon_id values must be non-empty")
    if len(normalized) != len(set(normalized)):
        raise ValueError("coupon_id values must be unique")

    assignments = [
        TREATMENTS[index % len(TREATMENTS)] for index in range(len(normalized))
    ]
    generator = random.Random(seed)
    generator.shuffle(assignments)
    return dict(zip(normalized, assignments))


def validate_balance(assignments: dict[str, str]) -> bool:
    """Return true when all groups occur and group counts differ by at most one."""
    if not assignments or any(group not in TREATMENTS for group in assignments.values()):
        return False
    counts = Counter(assignments.values())
    if set(counts) != set(TREATMENTS):
        return False
    return max(counts.values()) - min(counts.values()) <= 1
