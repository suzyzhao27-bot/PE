"""Frozen one-at-a-time instantaneous flux sensitivity at saved pool states."""

from __future__ import annotations

import bisect
import csv
import json
import math
from pathlib import Path

from mechanism_neutral_baseline import D_CL, D_NA, D_SALT, D_UREA
from mechanism_transport import external_resistance, salt_flux, urea_flux, PoreTable


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/processed/salt_transport"
MATRIX = ROOT / "原始资料/赵珈绮/22_机制预测与实验验证方案_v1.0.0/matrices/paired_predictions.json"
OUT = DATA / "22_algebraic_sensitivity_20261028.csv"
CASES = ("base", "a_2", "a_10", "a_50", "tau_1", "tau_4", "ru_015", "ru_025",
         "Dplus_080", "Dplus_120", "Dminus_080", "Dminus_120", "Dhigh_060",
         "Dhigh_140", "Rext_0_100", "Rext_100_0", "ionr_015", "ionr_035")


def read(name: str) -> list[dict[str, str]]:
    with (DATA / name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def renkin(radius_nm: float, pore_nm: float) -> float:
    x = radius_nm / pore_nm
    return (1 - x) ** 2 * (1 - 2.104 * x + 2.09 * x**3 - 0.95 * x**5)


class CoefficientTable:
    def __init__(self, rows: list[dict[str, str]], pore_nm: float,
                 d_plus: float, d_minus: float) -> None:
        source = sorted((float(r["c_mol_m3"]), float(r["k_plus"]), float(r["k_minus"]))
                        for r in rows if float(r["a_nm"]) == pore_nm
                        and float(r["sigma_C_m2"]) == -0.008)
        if len(source) < 31 or source[-1][0] != 600:
            raise ValueError(f"Incomplete PB table for a={pore_nm}")
        self.x = [0.0] + [v[0] for v in source]
        self.y = [0.0] + [2 * d_plus * d_minus * v[1] * v[2] /
                          (d_plus * v[1] + d_minus * v[2]) for v in source]
    def value(self, c: float) -> float:
        if not 0 <= c <= 600:
            raise ValueError(c)
        if c == 0:
            return 0.0
        i = min(bisect.bisect_left(self.x, c), len(self.x) - 1)
        dx = c - self.x[i-1]
        slope = (self.y[i] - self.y[i-1]) / (self.x[i] - self.x[i-1])
        return self.y[i-1] + slope * dx

    def integral(self, low: float, high: float) -> float:
        if low == high:
            return 0.0
        knots = [low] + [c for c in self.x if low < c < high] + [high]
        return sum((self.value(a) + self.value(b)) * (b - a) / 2
                   for a, b in zip(knots, knots[1:]))


def parameters(case: str, pore_nm: float, tau: float, ru_nm: float,
               dp: float, dm: float, feed: float) -> tuple[float, float, float, float, float, float]:
    if case.startswith("a_"):
        pore_nm = float(case.split("_")[1])
    elif case.startswith("tau_"):
        tau = float(case.split("_")[1])
    elif case == "ru_015":
        ru_nm = 0.15
    elif case == "ru_025":
        ru_nm = 0.25
    elif case == "Dplus_080":
        dp *= .8
    elif case == "Dplus_120":
        dp *= 1.2
    elif case == "Dminus_080":
        dm *= .8
    elif case == "Dminus_120":
        dm *= 1.2
    elif case == "Dhigh_060" and feed == 600:
        dp *= .6
        dm *= .6
    elif case == "Dhigh_140" and feed == 600:
        dp *= 1.4
        dm *= 1.4
    elif case == "ionr_015" and pore_nm == 2:
        dp *= renkin(.15, 2)
        dm *= renkin(.15, 2)
    elif case == "ionr_035" and pore_nm == 2:
        dp *= renkin(.35, 2)
        dm *= renkin(.35, 2)
    return pore_nm, tau, ru_nm, dp, dm, 2 * dp * dm / (dp + dm)


def salt_instant(donor: float, receiver: float, h_nm: float, b: float,
                 charged: bool, external: str, pore_nm: float, tau: float,
                 dp: float, dm: float, table: CoefficientTable | None,
                 left_fraction: float) -> float:
    ds = 2 * dp * dm / (dp + dm)
    ext = external_resistance("salt", external)  # Keep independently specified total external resistance fixed.
    entry = math.pi * pore_nm * 1e-9 / (4 * b * tau * ds)
    if not charged:
        return (donor - receiver) / (h_nm * 1e-9 / (b * ds) + 2 * entry + ext)
    if table is None:
        raise ValueError("Missing charged PB table")
    side_l, side_r = entry + left_fraction * ext, entry + (1 - left_fraction) * ext
    lo, hi = 0.0, (donor - receiver) / (side_l + side_r)
    for _ in range(55):
        flux = (lo + hi) / 2
        cm_d, cm_r = donor - flux * side_l, receiver + flux * side_r
        internal = b / (h_nm * 1e-9) * table.integral(cm_r, cm_d)
        if flux > internal:
            hi = flux
        else:
            lo = flux
    return (lo + hi) / 2


def urea_instant(donor: float, receiver: float, h_nm: float, b: float,
                 external: str, pore_nm: float, tau: float, ru_nm: float) -> float:
    h_u = renkin(ru_nm, pore_nm)
    p_mem = b * D_UREA / (h_nm * 1e-9 / h_u + math.pi * pore_nm * 1e-9 / (2 * tau))
    return (donor - receiver) / (1 / p_mem + external_resistance("urea", external))


def main() -> None:
    pairs = {p["pair_id"]: p for p in json.loads(MATRIX.read_text(encoding="utf-8-sig"))}
    curves = read("22_paired_illustrations_20261028.csv")
    if len(pairs) != 54 or len(curves) != 864 or len(read("22_pb_call_ledger_20261028.csv")) != 437:
        raise ArithmeticError("Frozen input or PB budget mismatch")
    pb = read("22_pb_main_20261028.csv") + read("22_pb_target_refinement_20261028.csv")
    integral_status = {(float(r["a_nm"]), float(r["sigma_C_m2"]), float(r["feed_mM"])):
                       r["status"] for r in read("22_pb_integral_check_20261028.csv")}
    integral_status.update({(float(r["a_nm"]), float(r["sigma_C_m2"]), float(r["feed_mM"])):
                            r["status"] for r in read("22_pb_target_recheck_20261028.csv")})
    cache: dict[tuple[float, float, float], CoefficientTable] = {}
    output = []
    for row in curves:
        if float(row["time_min"]) == 0:
            continue
        pair = pairs[row["pair_id"]]
        probe, state = row["probe"], row["state"]
        feed = float(pair["c_feed_mM"])
        h, b = float(pair["h_nm"]), float(pair[f"b_{state}"])
        charged = float(pair[f"sigma_{state}_C_m2"]) != 0
        donor, receiver = float(row["donor_mol_m3"]), float(row["receiver_mol_m3"])
        values: dict[str, float | str] = {}
        unresolved_cases = []
        for case in CASES:
            if case.startswith("Dhigh_") and feed != 600:
                continue
            if case.startswith("ionr_") and probe != "salt":
                continue
            if case.startswith("ionr_") and probe == "salt":
                # Ion-size endpoints accompany a=2 nm as one declared conditional subcase.
                pore, tau, ru, dp, dm, _ = parameters("a_2", 25, 2, .18, D_NA, D_CL, feed)
            else:
                pore, tau, ru, dp, dm, _ = parameters(case, 25, 2, .18, D_NA, D_CL, feed)
            if case.startswith("ionr_"):
                radius = .15 if case == "ionr_015" else .35
                dp *= renkin(radius, 2)
                dm *= renkin(radius, 2)
            if probe == "urea":
                if case.startswith(("Dplus_", "Dminus_", "Dhigh_", "ionr_", "Rext_")):
                    continue
                flux = urea_instant(donor, receiver, h, b, pair["external"], pore, tau, ru)
            else:
                if case.startswith("ru_"):
                    continue
                if charged and integral_status[(pore, -0.008, feed)] == "NUMERIC_UNRESOLVED_OVER_2_PERCENT":
                    unresolved_cases.append(case)
                key = (pore, dp, dm)
                table = None
                if charged:
                    if key not in cache:
                        cache[key] = CoefficientTable(pb, pore, dp, dm)
                    table = cache[key]
                left = 0.0 if case == "Rext_0_100" else 1.0 if case == "Rext_100_0" else .5
                flux = salt_instant(donor, receiver, h, b, charged, pair["external"],
                                    pore, tau, dp, dm, table, left)
            if flux < 0 or not math.isfinite(flux):
                raise ArithmeticError(f"Invalid flux: {row['pair_id']} {case}")
            values[case] = flux
        if probe == "salt":
            base = salt_flux(donor, receiver, h, b,
                             float(pair[f"sigma_{state}_C_m2"]), pair["external"],
                             PoreTable() if charged else None)
        else:
            base = urea_flux(donor, receiver, h, b, pair["external"])
        if abs(values["base"] - base) > 1e-8 * max(abs(base), 1e-30):
            raise ArithmeticError(f"Baseline mismatch: {row['pair_id']} {state} {probe} "
                                  f"new={values['base']:.12g} old={base:.12g}")
        candidates = values
        lower = min(candidates, key=candidates.get)
        upper = max(candidates, key=candidates.get)
        output.append({"pair_id": row["pair_id"], "h_nm": h, "feed_mM": feed,
                       "mechanism": pair["mechanism"], "external": pair["external"],
                       "state": state, "probe": probe, "time_min": row["time_min"],
                       "base_numeric_status": row["numeric_status"],
                       "variant_numeric_status": ("ALGEBRAIC_VARIANTS_T12_NOT_RECHECKED" if charged and probe == "salt"
                                                  else "ANALYTIC_NEUTRAL"),
                       "known_unresolved_cases": "|".join(unresolved_cases) or "NONE",
                       **{f"J_{k}_mol_m2_s": values.get(k, "NA") for k in CASES},
                       "min_case": lower, "min_J_mol_m2_s": candidates[lower],
                       "max_case": upper, "max_J_mol_m2_s": candidates[upper]})
    keys = {(r["pair_id"], r["state"], r["probe"], r["time_min"]) for r in output}
    if len(output) != 648 or len(keys) != 648:
        raise ArithmeticError(f"Expected 648 states, got {len(output)}")
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output[0]))
        writer.writeheader()
        writer.writerows(output)
    print(f"rows={len(output)} pb_calls_added=0 min_J={min(r['min_J_mol_m2_s'] for r in output):.6g}")


if __name__ == "__main__":
    main()
