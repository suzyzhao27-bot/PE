"""Compare frozen sampled curves with a no-withdrawal counterfactual."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from mechanism_transport import PoreTable, State, step


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/processed/salt_transport"
MATRIX = ROOT / "原始资料/赵珈绮/22_机制预测与实验验证方案_v1.0.0/matrices/paired_predictions.json"
FROZEN = DATA / "22_paired_illustrations_20261028.csv"
OUT = {
    "neutral": DATA / "22_中性曲线取样省略影响_20261029.csv",
    "charged": DATA / "22_带电盐曲线取样省略影响_20261029.csv",
}
FIELDS = ("pair_id", "h_nm", "feed_mM", "mechanism", "external", "state", "probe",
          "time_min", "frozen_sampled_receiver_uM", "no_withdrawal_receiver_uM",
          "sampled_minus_no_withdrawal_uM", "relative_to_no_withdrawal",
          "no_withdrawal_mass_relative_error", "interpretation")


def no_withdrawal(initial: float, probe: str, h: float, b: float, sigma: float,
                  external: str, table: PoreTable | None) -> list[tuple[int, float, float]]:
    state = State(donor=initial)
    mass0 = state.mass()
    output = []
    for minute in range(241):
        if minute in (60, 120, 240):
            output.append((minute, state.receiver, abs(state.mass() - mass0) / mass0))
        if minute < 240:
            step(state, 60, probe, h, b, sigma, external, table)
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=OUT, required=True)
    mode = parser.parse_args().mode
    pairs = json.loads(MATRIX.read_text(encoding="utf-8-sig"))
    with FROZEN.open(encoding="utf-8-sig", newline="") as stream:
        frozen_rows = list(csv.DictReader(stream))
    if len(pairs) != 54 or len(frozen_rows) != 864:
        raise ArithmeticError("Fixed source size changed")
    frozen = {(r["pair_id"], r["state"], r["probe"], int(float(r["time_min"]))): r
              for r in frozen_rows}
    if len(frozen) != 864:
        raise ArithmeticError("Duplicate frozen state")
    table = PoreTable() if mode == "charged" else None
    rows = []
    curves = 0
    for pair in pairs:
        for state in ("B", "E"):
            sigma = float(pair[f"sigma_{state}_C_m2"])
            for probe in ("salt", "urea"):
                charged = probe == "salt" and sigma != 0
                if charged != (mode == "charged"):
                    continue
                curves += 1
                initial = float(pair["c_feed_mM"]) if probe == "salt" else 10.0
                for minute, receiver, mass_error in no_withdrawal(
                    initial, probe, float(pair["h_nm"]), float(pair[f"b_{state}"]),
                    sigma, pair["external"], table
                ):
                    sampled = float(frozen[pair["pair_id"], state, probe, minute]["receiver_mol_m3"])
                    if minute == 60 and abs(sampled - receiver) > 1e-12:
                        raise ArithmeticError("Sampling changed pre-withdrawal 60-minute state")
                    rows.append({"pair_id": pair["pair_id"], "h_nm": pair["h_nm"],
                                 "feed_mM": pair["c_feed_mM"], "mechanism": pair["mechanism"],
                                 "external": pair["external"], "state": state, "probe": probe,
                                 "time_min": minute, "frozen_sampled_receiver_uM": sampled * 1000,
                                 "no_withdrawal_receiver_uM": receiver * 1000,
                                 "sampled_minus_no_withdrawal_uM": (sampled - receiver) * 1000,
                                 "relative_to_no_withdrawal": (sampled - receiver) / receiver if receiver else 0,
                                 "no_withdrawal_mass_relative_error": mass_error,
                                 "interpretation": "NO_WITHDRAWAL_COUNTERFACTUAL_NOT_PROTOCOL"})
    expected = 36 if mode == "charged" else 180
    if curves != expected or len(rows) != expected * 3:
        raise ArithmeticError("Curve coverage changed")
    with OUT[mode].open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    later = [r for r in rows if r["time_min"] > 60]
    print(f"mode={mode} curves={curves} rows={len(rows)} "
          f"max_abs_uM={max(abs(r['sampled_minus_no_withdrawal_uM']) for r in later):.8g} "
          f"max_abs_relative={max(abs(r['relative_to_no_withdrawal']) for r in later):.8g} "
          f"max_mass_error={max(r['no_withdrawal_mass_relative_error'] for r in rows):.8g}")


if __name__ == "__main__":
    main()
