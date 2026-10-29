"""Derive the fixed three-parent paired-log 95% t-interval precision boundary."""

from __future__ import annotations

import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/processed/salt_transport/22_三母膜配对区间精度边界_20261029.json"


def main() -> None:
    n = 3
    # For df=2, CDF(t)=1/2+t/(2*sqrt(t*t+2)).
    q = 0.95  # 2*CDF(t_0.975)-1
    t95_df2 = math.sqrt(2 * q * q / (1 - q * q))
    halfwidth_per_sample_sd = t95_df2 / math.sqrt(n)
    max_log_pi_sd_when_mean_one = math.log(2) / halfwidth_per_sample_sd
    result = {
        "source_rule": "22_04_section_6_parent_film_n3_log_effect_95pct_t_df2",
        "n_independent_parent_films": n,
        "df": 2,
        "two_sided_95pct_t_critical": t95_df2,
        "log_interval_halfwidth_per_sample_sd": halfwidth_per_sample_sd,
        "positive_log_gain_lower_bound_condition": "mean_log_gain > halfwidth_per_sample_sd * sample_sd_log_gain",
        "pi_band": [0.5, 2.0],
        "log_pi_band_halfwidth": math.log(2),
        "max_sample_sd_log_pi_for_95pct_interval_inside_band_when_mean_pi_is_1": max_log_pi_sd_when_mean_one,
        "general_pi_band_condition": "abs(mean_log_pi) + halfwidth_per_sample_sd * sample_sd_log_pi < log(2)",
        "note": "Analytic design boundary only. No experimental SD, gain, or confidence interval estimated.",
    }
    if not math.isclose(0.5 + t95_df2 / (2 * math.sqrt(t95_df2**2 + 2)),
                        0.975, abs_tol=1e-14):
        raise ArithmeticError("Incorrect df=2 t critical")
    OUT.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"n={n} tcrit={t95_df2:.10f} halfwidth_per_sd={halfwidth_per_sample_sd:.10f} "
          f"pi_log_sd_limit_at_mean_one={max_log_pi_sd_when_mean_one:.10f}")


if __name__ == "__main__":
    main()
