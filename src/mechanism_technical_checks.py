"""Read-only numerical audit of frozen 22-series technical checks; no PB solve."""

from __future__ import annotations

import csv
import json
import math
from collections import Counter
from pathlib import Path

from mechanism_neutral_baseline import D_CL, D_NA, D_SALT


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/processed/salt_transport"
OUT = DATA / "22_technical_checks_20261028.json"


def rows(name: str) -> list[dict[str, str]]:
    with (DATA / name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def relative(a: float, b: float) -> float:
    return abs(a - b) / max(abs(a), abs(b), 1e-300)


def main() -> None:
    pb = rows("22_pb_main_20261028.csv")
    neutral = rows("22_neutral_baseline_20261028.csv")
    curves = rows("22_paired_illustrations_20261028.csv")
    effects = rows("22_paired_effects_20261028.csv")
    target = rows("22_pb_target_recheck_20261028.csv")
    ledger = rows("22_pb_call_ledger_20261028.csv")
    assert len(pb) == 372 and len(neutral) == 288 and len(curves) == 864
    assert len(effects) == 54 and len(target) == 12 and len(ledger) == 437

    exchange_error = zero_current_error = naive_zero_current_error = coefficient_error = 0.0
    for row in pb:
        kp, km = float(row["k_plus"]), float(row["k_minus"])
        a, b = D_NA * kp, D_CL * km
        coefficient = 2 * a * b / (a + b)
        coefficient_error = max(coefficient_error, relative(coefficient, float(row["D_s_local_m2_s"])))
        # Opposite charge swaps K+ and K-; exchanging ion diffusivities swaps a,b.
        swapped = 2 * (D_CL * km) * (D_NA * kp) / (D_CL * km + D_NA * kp)
        exchange_error = max(exchange_error, relative(coefficient, swapped))
        # Unit concentration gradient; induced electric term enforces J+ = J-.
        electric = (b - a) / (a + b)
        naive_zero_current_error = max(naive_zero_current_error,
                                       relative(a * (1 + electric), b * (1 - electric)))
        # Avoid cancellation in 1 +/- electric at very asymmetric low-salt points.
        j_plus, j_minus = a * (2 * b / (a + b)), b * (2 * a / (a + b))
        zero_current_error = max(zero_current_error, relative(j_plus, j_minus))

    donnan_error = 0.0
    for psi in (-2.0, -0.5, 0.0, 0.5, 2.0):
        kp, km = math.exp(-psi), math.exp(psi)
        d = D_NA
        numeric = 2 * d * kp * d * km / (d * kp + d * km)
        analytic = d / math.cosh(psi)
        donnan_error = max(donnan_error, relative(numeric, analytic))

    mass_error = max(float(row["mass_relative_error"]) for row in curves)
    neutral_mass_error = max(float(row["mass_relative_error"]) for row in neutral)
    identity_error = max(relative(float(row["E_s_initial"]),
                                  float(row["E_u_initial"]) * float(row["Pi_initial"]))
                         for row in effects)
    known_gain_error = max(
        abs(float(row["E_u_initial"]) - (1 if row["mechanism"] == "C" else 100))
        for row in effects
    )
    statuses = Counter(row["status"] for row in target)
    result = {
        "input_rows": {"pb_main": len(pb), "neutral": len(neutral), "paired_curves": len(curves),
                       "paired_effects": len(effects), "target_recheck": len(target),
                       "pb_calls": len(ledger)},
        "T02_T07_neutral_D_relative_error": relative(D_SALT, 2 * D_NA * D_CL / (D_NA + D_CL)),
        "T05_exchange_max_relative_error_372": exchange_error,
        "T06_zero_current_max_relative_error_372": zero_current_error,
        "T06_naive_electric_reconstruction_max_relative_error_372": naive_zero_current_error,
        "T06_stored_coefficient_max_relative_error_372": coefficient_error,
        "T08_uniform_donnan_max_relative_error_5": donnan_error,
        "T09_paired_max_mass_relative_error_864": mass_error,
        "T09_neutral_max_mass_relative_error_288": neutral_mass_error,
        "T11_Es_Eu_Pi_max_relative_error_54": identity_error,
        "T11_known_Eu_max_absolute_error_54": known_gain_error,
        "T12_target_status_counts": dict(statuses),
        "T12_target_max_updated_31_16_relative_difference": max(
            float(row["updated_31_16_relative_difference"]) for row in target),
        "note": "Existing artifacts only; no PB calls. Algebraic checks are module checks, not material validation.",
    }
    print(json.dumps(result, ensure_ascii=False))
    assert exchange_error < 1e-12 and zero_current_error < 1e-6
    assert donnan_error < 1e-12 and mass_error < 1e-8 and identity_error < 1e-12
    assert known_gain_error < 1e-10
    assert statuses["NUMERIC_UNRESOLVED_OVER_2_PERCENT"] == 10
    with OUT.open("w", encoding="utf-8") as stream:
        json.dump(result, stream, ensure_ascii=False, indent=2)


if __name__ == "__main__":
    main()
