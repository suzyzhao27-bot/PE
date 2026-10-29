"""Cross-check fixed WP6 candidates against short-pore and T12 flags."""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/processed/salt_transport"
SHORT = DATA / "22_WP6候选短孔适用范围_20261029.csv"
T12 = DATA / "22_WP6候选T12低盐影响范围_20261029.csv"
OUT = DATA / "22_WP6短孔与T12风险交集_20261029.csv"
FIELDS = ("candidate_id", "candidate_scope", "low_salt_choice_mM", "a_nm", "tau",
          "sigma_B_mC_m2", "sigma_E_mC_m2", "thin_nominal_tau_h_over_a",
          "charged_short_pore_flag", "low_salt_T12_flag", "joint_flag",
          "other_model_limit", "actual_experiment_status")


def read(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    geometry = read(SHORT)
    t12 = read(T12)
    if len(geometry) != 576 or len(t12) != 576:
        raise ArithmeticError("Fixed audit inputs changed")
    thin = {r["candidate_id"]: r for r in geometry if r["thickness_level"] == "thin"}
    if len(thin) != 192:
        raise ArithmeticError("Unexpected thin-candidate count")
    output = []
    for status in t12:
        feed = int(status["salt_feed_mM"])
        if feed == 600:
            continue
        candidate = thin[status["candidate_id"]]
        if any(candidate[k] != status[k] for k in ("candidate_scope", "a_nm", "tau",
                                                    "sigma_B_mC_m2", "sigma_E_mC_m2")):
            raise ArithmeticError("Candidate axis mismatch")
        charge_active = int(candidate["sigma_B_mC_m2"]) != 0 or int(candidate["sigma_E_mC_m2"]) != 0
        short = candidate["geometry_status"] == "SHORT_PORE_CONDITIONAL" and charge_active
        numeric = status["candidate_T12_status"] == "T12_UNRESOLVED"
        label = "BOTH_FLAGS" if short and numeric else ("SHORT_PORE_ONLY" if short else
                "T12_ONLY" if numeric else "NEITHER_OF_THESE_TWO_FLAGS")
        output.append({"candidate_id": status["candidate_id"],
                       "candidate_scope": status["candidate_scope"],
                       "low_salt_choice_mM": feed, "a_nm": status["a_nm"],
                       "tau": status["tau"], "sigma_B_mC_m2": status["sigma_B_mC_m2"],
                       "sigma_E_mC_m2": status["sigma_E_mC_m2"],
                       "thin_nominal_tau_h_over_a": candidate["tau_h_over_a"],
                       "charged_short_pore_flag": int(short), "low_salt_T12_flag": int(numeric),
                       "joint_flag": label,
                       "other_model_limit": "UNCONFIRMED_2NM_CONTINUUM" if int(status["a_nm"]) == 2 else "NOT_ASSESSED_HERE",
                       "actual_experiment_status": "NOT_STARTED"})
    if len(output) != 384 or len({(r["candidate_id"], r["low_salt_choice_mM"]) for r in output}) != 384:
        raise ArithmeticError("Expected 192 candidates x two low-salt choices")
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(output)
    counts = Counter((r["candidate_scope"], r["low_salt_choice_mM"], r["joint_flag"]) for r in output)
    for feed in (1, 10):
        main = {flag: counts["MAIN_144", feed, flag] for flag in
                ("BOTH_FLAGS", "SHORT_PORE_ONLY", "T12_ONLY", "NEITHER_OF_THESE_TWO_FLAGS")}
        if sum(main.values()) != 144:
            raise ArithmeticError("Main candidate partition failed")
        print(f"low={feed}mM main={main}")


if __name__ == "__main__":
    main()
