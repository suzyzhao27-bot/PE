"""Prepare the twelve frozen E3 raw-reading perturbation inputs."""

from __future__ import annotations

import csv
import json

import numpy as np

from v7_d1_objective import H_VALUES, ROOT, curves


OUT = ROOT / "WP5" / "v7_E3_perturbed_B_readings.csv"
FROZEN = (ROOT / "原始资料/赵珈绮/"
          "14_论文补强学生任务书_v7.0.0/configs/frozen_spec.json")


def main() -> None:
    spec = json.loads(FROZEN.read_text(encoding="utf-8-sig"))["matrices"]["E3"]
    assert (spec["eta"], spec["sigma_floor_raw"], spec["seeds"],
            spec["max_jobs"], spec["perturb_raw_then_subtract_shared_baseline"]) == (
                [0.01, 0.05], 0.01, [8101, 8102, 8103, 8104, 8105, 8106],
                12, True)
    if OUT.exists():
        raise FileExistsError(f"E3 input already exists; no silent overwrite: {OUT}")

    observed = curves()
    rows = []
    for eta in spec["eta"]:
        for seed in spec["seeds"]:
            # Reuse the same standard-normal realization across eta at each seed.
            rng = np.random.Generator(np.random.PCG64(seed))
            for h in H_VALUES:
                points = observed[h]
                raw = np.array([reading for _, reading, _ in points], dtype=float)
                sigma = np.maximum(spec["sigma_floor_raw"], eta * np.abs(raw))
                noisy = raw + rng.normal(0.0, sigma)
                baseline = float(noisy[0])
                for (time_min, raw_input, curve_id), raw_noisy, scale in zip(
                        points, noisy, sigma, strict=True):
                    rows.append({
                        "eta": eta, "seed": seed, "h_nm": h,
                        "curve_id": curve_id, "time_min": time_min,
                        "raw_input_uS_cm": raw_input,
                        "sigma_raw_uS_cm": scale,
                        "raw_perturbed_uS_cm": float(raw_noisy),
                        "baseline_perturbed_uS_cm": baseline,
                        "increment_perturbed_uS_cm": float(raw_noisy - baseline),
                    })
    assert len(rows) == 12 * 6 * 7
    with OUT.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    negatives = sum(row["increment_perturbed_uS_cm"] < 0 for row in rows)
    print(f"E3 prepared: jobs=12; rows={len(rows)}; negative increments={negatives}")


if __name__ == "__main__":
    main()
