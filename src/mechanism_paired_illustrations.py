"""Run 54 frozen illustrative pairs; preserve conditional numeric-risk labels."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from mechanism_transport import PoreTable, run_curve


ROOT = Path(__file__).resolve().parents[1]
MATRIX = ROOT / "原始资料/赵珈绮/22_机制预测与实验验证方案_v1.0.0/matrices/paired_predictions.json"
OUT = ROOT / "data/processed/salt_transport/22_paired_illustrations_20261028.csv"
FIELDS = ["pair_id", "h_nm", "c_feed_mM", "mechanism", "external", "state",
          "probe", "time_min", "donor_mol_m3", "receiver_mol_m3",
          "mass_relative_error", "numeric_status"]


def main() -> None:
    with MATRIX.open(encoding="utf-8-sig") as stream:
        pairs = json.load(stream)
    assert len(pairs) == 54 and len({p["pair_id"] for p in pairs}) == 54
    table = PoreTable(a_nm=25, sigma=-0.008)
    rows = []
    for pair in pairs:
        for state in ("B", "E"):
            b = float(pair[f"b_{state}"])
            sigma = float(pair[f"sigma_{state}_C_m2"])
            for probe in ("salt", "urea"):
                initial = float(pair["c_feed_mM"]) if probe == "salt" else 10.0
                curve = run_curve(initial, probe, float(pair["h_nm"]), b, sigma,
                                  pair["external"], table if sigma else None)
                status = ("NUMERIC_UNRESOLVED_T12" if probe == "salt" and sigma
                          and float(pair["c_feed_mM"]) == 10
                          else "MODEL_ONLY_TARGET_CHECK")
                for item in curve:
                    rows.append({
                        "pair_id": pair["pair_id"], "h_nm": pair["h_nm"],
                        "c_feed_mM": pair["c_feed_mM"], "mechanism": pair["mechanism"],
                        "external": pair["external"], "state": state, "probe": probe,
                        "time_min": item["time_min"], "donor_mol_m3": item["donor_mol_m3"],
                        "receiver_mol_m3": item["receiver_mol_m3"],
                        "mass_relative_error": item["mass_relative_error"],
                        "numeric_status": status,
                    })
    if not OUT.parent.is_dir():
        raise FileNotFoundError(OUT.parent)
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"pairs={len(pairs)} single_ODE={len(pairs)*2} rows={len(rows)} "
          f"max_mass_error={max(r['mass_relative_error'] for r in rows):.3g} "
          f"numeric_unresolved_rows={sum(r['numeric_status']=='NUMERIC_UNRESOLVED_T12' for r in rows)}")


if __name__ == "__main__":
    main()
