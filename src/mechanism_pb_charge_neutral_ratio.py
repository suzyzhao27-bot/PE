"""Compare fixed charged PB coefficients and integrals with the neutral salt limit."""

from __future__ import annotations

import csv
from pathlib import Path

from mechanism_neutral_baseline import D_SALT
from mechanism_transport import PoreTable


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/processed/salt_transport"
RECHECK = DATA / "22_pb_target_recheck_20261028.csv"
INITIAL = DATA / "22_pb_integral_check_20261028.csv"
OUT = DATA / "22_PB荷电与中性盐系数比_20261029.csv"
FIELDS = ("a_nm", "sigma_mC_m2", "feed_mM", "charged_local_D_m2_s",
          "neutral_local_D_m2_s", "charged_to_neutral_local_ratio",
          "charged_integral_mol_m_s", "neutral_integral_mol_m_s",
          "charged_to_neutral_integral_ratio", "T12_status", "interpretation")


def read(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    initial = read(INITIAL)
    recheck = read(RECHECK)
    if len(initial) != 36 or len(recheck) != 12:
        raise ArithmeticError("Fixed PB audit tables changed")
    risk = {(int(float(r["a_nm"])), round(float(r["sigma_C_m2"]) * 1000),
             int(float(r["feed_mM"]))): r["status"] for r in recheck}
    rows = []
    for a in (2, 10, 25, 50):
        for sigma_mc in (-2, -8, -20):
            table = PoreTable(a_nm=a, sigma=sigma_mc / 1000)
            for feed in (1, 10, 600):
                local = table.value(feed)
                integral = table.integral(0, feed)
                source = next(r for r in initial if int(float(r["a_nm"])) == a
                              and round(float(r["sigma_C_m2"]) * 1000) == sigma_mc
                              and int(float(r["feed_mM"])) == feed)
                if source["status"] == "PASS":
                    status = "GRID_CHECK_PASS"
                else:
                    revised = risk[a, sigma_mc, feed]
                    status = ("T12_UNRESOLVED" if revised == "NUMERIC_UNRESOLVED_OVER_2_PERCENT"
                              else "GRID_CHECK_PASS_AFTER_ENDPOINT")
                rows.append({"a_nm": a, "sigma_mC_m2": sigma_mc, "feed_mM": feed,
                             "charged_local_D_m2_s": local,
                             "neutral_local_D_m2_s": D_SALT,
                             "charged_to_neutral_local_ratio": local / D_SALT,
                             "charged_integral_mol_m_s": integral,
                             "neutral_integral_mol_m_s": D_SALT * feed,
                             "charged_to_neutral_integral_ratio": integral / (D_SALT * feed),
                             "T12_status": status,
                             "interpretation": "FIXED_PB_MODEL_ONLY"})
    if len(rows) != 36:
        raise ArithmeticError("Expected 12 PB families x three feeds")
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    high_main = [r for r in rows if r["feed_mM"] == 600 and r["a_nm"] >= 10]
    print(f"rows={len(rows)} high_main_local_ratio={min(r['charged_to_neutral_local_ratio'] for r in high_main):.8g}.."
          f"{max(r['charged_to_neutral_local_ratio'] for r in high_main):.8g} "
          f"high_main_integral_ratio={min(r['charged_to_neutral_integral_ratio'] for r in high_main):.8g}.."
          f"{max(r['charged_to_neutral_integral_ratio'] for r in high_main):.8g}")


if __name__ == "__main__":
    main()
