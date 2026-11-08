"""Map existing PB concentration-grid checks onto fixed WP6 B/E charge candidates."""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/processed/salt_transport"
SHORT = DATA / "22_WP6候选短孔适用范围_20261029.csv"
INITIAL = DATA / "22_pb_integral_check_20261028.csv"
RECHECK = DATA / "22_pb_target_recheck_20261028.csv"
OUT = DATA / "22_WP6候选T12低盐影响范围_20261029.csv"
FIELDS = ("candidate_id", "candidate_scope", "a_nm", "tau", "sigma_B_mC_m2",
          "sigma_E_mC_m2", "salt_feed_mM", "B_grid_status", "E_grid_status",
          "candidate_T12_status", "postfreeze_25nm_check_note", "measured_status")


def read(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def key(a: float, sigma_c: float, feed: float) -> tuple[int, int, int]:
    return int(a), round(sigma_c * 1000), int(feed)


def main() -> None:
    short = read(SHORT)
    initial = read(INITIAL)
    recheck = read(RECHECK)
    if len(short) != 576 or len(initial) != 36 or len(recheck) != 12:
        raise ArithmeticError("Fixed inputs changed")
    base = {key(float(r["a_nm"]), float(r["sigma_C_m2"]), float(r["feed_mM"])): r["status"]
            for r in initial}
    updated = {key(float(r["a_nm"]), float(r["sigma_C_m2"]), float(r["feed_mM"])): r["status"]
               for r in recheck}
    if len(base) != 36 or len(updated) != 12:
        raise ArithmeticError("Duplicate PB target")
    candidates = [r for r in short if r["thickness_level"] == "thin"]
    if len(candidates) != 192 or len({r["candidate_id"] for r in candidates}) != 192:
        raise ArithmeticError("Candidate list changed")

    def status(a: int, sigma: int, feed: int) -> str:
        if sigma == 0:
            return "ZERO_CHARGE_ANALYTIC"
        k = a, sigma, feed
        old = base[k]
        if old == "PASS":
            return "GRID_CHECK_PASS"
        newer = updated[k]
        if newer == "PASS_AFTER_ENDPOINT":
            return "GRID_CHECK_PASS_AFTER_ENDPOINT"
        if newer == "NUMERIC_UNRESOLVED_OVER_2_PERCENT":
            return "T12_UNRESOLVED"
        raise ArithmeticError(f"Unknown PB status {old}/{newer}")

    rows = []
    for candidate in candidates:
        a = int(candidate["a_nm"])
        sb, se = int(candidate["sigma_B_mC_m2"]), int(candidate["sigma_E_mC_m2"])
        for feed in (1, 10, 600):
            b_status, e_status = status(a, sb, feed), status(a, se, feed)
            unresolved = "T12_UNRESOLVED" in (b_status, e_status)
            rows.append({"candidate_id": candidate["candidate_id"],
                         "candidate_scope": candidate["candidate_scope"], "a_nm": a,
                         "tau": candidate["tau"], "sigma_B_mC_m2": sb,
                         "sigma_E_mC_m2": se, "salt_feed_mM": feed,
                         "B_grid_status": b_status, "E_grid_status": e_status,
                         "candidate_T12_status": "T12_UNRESOLVED" if unresolved else "NO_KNOWN_T12_FLAG",
                         "postfreeze_25nm_check_note": "LOCAL_EXTRA_POINTS_DO_NOT_CLEAR_T12" if a == 25 and feed == 10 and -8 in (sb, se) else "NA",
                         "measured_status": "NOT_STARTED"})
    if len(rows) != 576:
        raise ArithmeticError("Expected 192 candidates x 3 feeds")
    counts = Counter((r["candidate_scope"], int(r["salt_feed_mM"]))
                     for r in rows if r["candidate_T12_status"] == "T12_UNRESOLVED")
    if any(r["candidate_T12_status"] == "T12_UNRESOLVED" for r in rows if int(r["salt_feed_mM"]) == 600):
        raise ArithmeticError("Unexpected high-salt T12 flag")
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"candidate_feed_rows={len(rows)} main_1mM_unresolved={counts['MAIN_144', 1]} "
          f"main_10mM_unresolved={counts['MAIN_144', 10]} "
          f"exploratory_1mM_unresolved={counts['UNCONFIRMED_2NM_EXPLORATORY', 1]} "
          f"exploratory_10mM_unresolved={counts['UNCONFIRMED_2NM_EXPLORATORY', 10]}")


if __name__ == "__main__":
    main()
