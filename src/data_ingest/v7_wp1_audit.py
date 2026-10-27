"""Rebuild the v7 WP1 data, geometry, and descriptive time audits.

The v7 taskbook under 原始资料 is read only. No conductivity is converted to flux.
"""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
TASK = ROOT / "原始资料/赵珈绮/14_论文补强学生任务书_v7.0.0"
OUT = ROOT / "WP1"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, fields: list[str], rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def checksum(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    spec_path = TASK / "configs/frozen_spec.json"
    spec = json.loads(spec_path.read_text(encoding="utf-8-sig"))
    source = TASK / "source_manifest.csv"
    manifest = read_csv(source)
    for item in manifest:
        path = TASK / item["package_relative"]
        assert path.is_file(), path
        assert path.stat().st_size == int(item["bytes"]), path
        assert checksum(path) == item["sha256"], path

    salt_path = TASK / spec["observations"]["file"]
    assert checksum(salt_path) == spec["observations"]["sha256"]
    salt = read_csv(salt_path)
    roles = read_csv(TASK / spec["observations"]["roles_file"])
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for point in salt:
        grouped[point["curve_id"]].append(point)
    assert len(salt) == 167 and len(grouped) == 17
    assert set(grouped) == {r["curve_id"] for r in roles}
    assert Counter({r["role"]: sum(len(grouped[x["curve_id"]]) for x in roles if x["role"] == r["role"]) for r in roles}) == {"B": 42, "L": 70, "C": 14, "U": 20, "D": 21}

    role_rows = []
    for role in roles:
        points = sorted(grouped[role["curve_id"]], key=lambda x: int(x["point_index"]))
        assert len(points) == int(role["n_points"])
        assert all(x["eligible_for_mechanistic_fit"] == role["original_eligibility"] for x in points)
        assert all(x["time_unit"] == "min" for x in points)
        t = [float(x["time"]) for x in points]
        assert t[0] == 0 and all(b > a for a, b in zip(t, t[1:]))
        assert max(t) == float(role["t_max_min"])
        role_rows.append({**role, "source_opju": points[0]["source_opju"], "condition": points[0]["condition"],
                          "time_column_evidence": points[0]["time_column_evidence"],
                          "unit_evidence": points[0]["unit_evidence"],
                          "conductivity_unit": points[0]["conductivity_unit"],
                          "identity_level": "A_source_curve" if role["role"] in {"B", "D"} else "U_pairing" if role["role"] == "U" else "B_same_material_family" if role["role"] == "L" else "A_control_curve",
                          "same_film_geometry_pairing": "unconfirmed",
                          "main_fit": "yes_36_nonzero_points" if role["role"] == "B" else "no",
                          "late_transfer": "conditional_after_240_min" if role["role"] == "L" else "no"})
    write_csv(OUT / "data_role_table.csv", list(role_rows[0]), role_rows)

    duplicate_rows = []
    for thickness in (65, 700, 2600):
        book1 = next(grouped[r["curve_id"]] for r in roles if r["role"] == "B" and r["sample_id"] == f"{thickness} nm")
        book2 = next(grouped[r["curve_id"]] for r in roles if r["role"] == "D" and r["sample_id"] == f"{thickness} nm")
        for p1, p2 in zip(book1, book2):
            assert float(p1["conductivity"]) == float(p2["conductivity"])
            duplicate_rows.append({"sample_id": f"{thickness} nm", "point_index": p1["point_index"],
                                   "book1_time_min": p1["time"], "book2_time_min": p2["time"],
                                   "book1_kappa": p1["conductivity"], "book2_kappa": p2["conductivity"],
                                   "conductivity_equal": "true", "time_equal": str(p1["time"] == p2["time"]).lower()})
    assert sum(x["time_equal"] == "false" for x in duplicate_rows) == 2
    write_csv(OUT / "book2_duplicate_check.csv", list(duplicate_rows[0]), duplicate_rows)

    geometry_rows = []
    geom_dir = TASK / "inputs/geometry/geometry_raw"
    for anchor in spec["geometry"]["anchors"]:
        draw = anchor["draw"]
        means = {}
        counts = {}
        for version in (["cal2", "base_historical"] if draw == 604 else ["main"]):
            paths = sorted(geom_dir.glob(f"Polymer_SEM_{draw}_view*{'_' + version if version != 'main' else ''}.csv"))
            assert len(paths) == 2, (draw, version, paths)
            values = []
            for path in paths:
                for item in read_csv(path):
                    if item[" "].strip().isdigit():
                        values.append(float(item["Length"]))
            assert len(values) == 100, (draw, version, len(values))
            means[version] = float(np.mean(values))
            counts[version] = len(values)
        primary = "cal2" if draw == 604 else "main"
        assert abs(means[primary] - anchor["df_nm"]) < 1e-5, (draw, means[primary], anchor["df_nm"])
        geometry_rows.append({"draw": draw, "h_nm": anchor["h_nm"], "SEM_count": counts[primary],
                              "SEM_mean_nm": f"{means[primary]:.5f}", "SEM_version": primary,
                              "SEM_base_historical_nm": f"{means['base_historical']:.5f}" if draw == 604 else "",
                              "dp_nm": anchor["dp_nm"], "dp_low_nm": anchor["dp_interval_nm"][0],
                              "dp_high_nm": anchor["dp_interval_nm"][1],
                              "phi_A": anchor["phi_A"], "phi_low": anchor["phi_interval"][0],
                              "phi_high": anchor["phi_interval"][1],
                              "source_type": anchor["source_type"],
                              "digitization_file": "inputs/geometry/Polymer2023_p7.png Fig.9c" if draw != 604 else "Polymer2023 main text and SEM cal2",
                              "changed_before_fit": "false", "reason": "source-based midpoint retained; not fitted to salt"})
    write_csv(OUT / "geometry_anchor_audit.csv", list(geometry_rows[0]), geometry_rows)

    interval_rows = []
    fit_rows = []
    for role in roles:
        points = sorted(grouped[role["curve_id"]], key=lambda x: int(x["point_index"]))
        t = np.array([float(x["time"]) for x in points])
        k = np.array([float(x["conductivity"]) for x in points])
        y = k - k[0]
        avg_rate = (y[-1] - y[0]) / (t[-1] - t[0])
        for i, point in enumerate(points):
            rate = (y[i] - y[i-1]) / (t[i] - t[i-1]) if i else float("nan")
            interval_rows.append({"curve_id": role["curve_id"], "role": role["role"], "point_index": point["point_index"],
                                  "time_min": t[i], "kappa_raw": k[i], "baseline": k[0], "y_increment": y[i],
                                  "interval_start_min": t[i-1] if i else "", "interval_rate_raw_per_min": rate if i else "",
                                  "whole_span_mean_rate": avg_rate, "rate_over_mean": rate / avg_rate if i and avg_rate != 0 else "",
                                  "conductivity_unit": point["conductivity_unit"], "time_unit": "min"})
        linear = np.column_stack([np.ones(len(t)), t])
        coeff = np.linalg.lstsq(linear, y, rcond=None)[0]
        pred = linear @ coeff
        sse = float(np.sum((y-pred)**2))
        best = None
        for j in range(2, len(t)-3):
            knot = t[j]
            design = np.column_stack([np.ones(len(t)), t, np.maximum(0, t-knot)])
            c = np.linalg.lstsq(design, y, rcond=None)[0]
            trial = design @ c
            trial_sse = float(np.sum((y-trial)**2))
            if best is None or trial_sse < best[0]:
                best = (trial_sse, knot, c, trial)
        # At least three observed points on each side of the interior knot.
        assert best is not None
        positive_rates = (np.diff(y) / np.diff(t))
        positive_rates = positive_rates[positive_rates > 0]
        shape_ratio = float(np.max(positive_rates) / np.min(positive_rates)) if len(positive_rates) >= 2 else ""
        late = np.where(t > 240)[0]
        early = np.where(t <= 240)[0]
        fit_rows.append({"curve_id": role["curve_id"], "role": role["role"], "n_points": len(t),
                         "linear_sse_raw": sse, "piecewise_sse_raw": best[0],
                         "sse_reduction": (sse-best[0])/sse if sse else "",
                         "knot_min": best[1], "linear_rmse_asinh": float(np.sqrt(np.mean((np.arcsinh(pred)-np.arcsinh(y))**2))),
                         "piecewise_rmse_asinh": float(np.sqrt(np.mean((np.arcsinh(best[3])-np.arcsinh(y))**2))),
                         "positive_interval_rate_ratio": shape_ratio,
                         "n_early_le_240": len(early) if role["role"] == "L" else "",
                         "n_late_gt_240": len(late) if role["role"] == "L" else "",
                         "early_mean_rate": (y[early[-1]]-y[early[0]])/(t[early[-1]]-t[early[0]]) if role["role"] == "L" and len(early) > 1 else "",
                         "late_mean_rate": (y[late[-1]]-y[late[0]])/(t[late[-1]]-t[late[0]]) if role["role"] == "L" and len(late) > 1 else "",
                         "interpretation": "descriptive_only_not_wetting_time"})
    write_csv(OUT / "time_shape_diagnostics.csv", list(interval_rows[0]), interval_rows)
    write_csv(OUT / "time_shape_fit_summary.csv", list(fit_rows[0]), fit_rows)

    actual = {"source_config": str(spec_path.relative_to(ROOT)).replace('\\', '/'),
              "source_config_sha256": checksum(spec_path), "difference_from_frozen": {},
              "frozen_spec": spec, "execution": {"WP1_data_audit": "completed", "N0_NS_ND": "not_started"}}
    (OUT / "actual_config.json").write_text(json.dumps(actual, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"{len(manifest)} source copies verified; {len(grouped)} curves, {len(salt)} points; {len(interval_rows)} time rows; {len(fit_rows)} shape fits")
    print(f"Book2 conductivity equal 21/21; differing time bindings {sum(x['time_equal']=='false' for x in duplicate_rows)}")
    print("SEM primary means:", ", ".join(f"{x['draw']}={x['SEM_mean_nm']} nm" for x in geometry_rows))


if __name__ == "__main__":
    main()
