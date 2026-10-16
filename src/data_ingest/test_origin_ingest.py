import csv
import json
import unittest
from pathlib import Path

from analyze_salt_transport import structural_qc, write_exact_duplicates
from normalize_salt_transport import CURVES, EXPORT_ROOT, OUTPUT_ROOT, normalize_curve
from validate_experiment_records import load_tracker, validate_tracker


class OriginIngestTests(unittest.TestCase):
    def test_source_hashes_match(self):
        manifest = json.loads((EXPORT_ROOT / "manifest.json").read_text(encoding="utf-8-sig"))
        self.assertEqual(len(manifest), 5)
        for project in manifest:
            self.assertEqual(project["source_sha256_before"], project["source_sha256_after"])

    def test_curve_mapping_and_point_count(self):
        self.assertEqual(len(CURVES), 17)
        rows = [row for curve in CURVES for row in normalize_curve(curve)]
        self.assertEqual(len(rows), 167)
        self.assertTrue(all(row["eligible_for_mechanistic_fit"] is False for row in rows))

    def test_long_duration_file_is_not_six_hours(self):
        rows = [
            row
            for curve in CURVES
            if curve.condition == "dry_6h"
            for row in normalize_curve(curve)
        ]
        self.assertEqual(max(row["time"] for row in rows), 2831.0)
        self.assertGreater(max(row["time"] for row in rows), 6 * 60)

    def test_exact_copies_are_not_replicates(self):
        with (OUTPUT_ROOT / "salt_transport_tidy.csv").open(
            "r", encoding="utf-8-sig", newline=""
        ) as handle:
            rows = list(csv.DictReader(handle))
        duplicates = write_exact_duplicates(rows)
        self.assertEqual(len(duplicates), 3)

    def test_record_tracker_covers_all_curves_without_blank_fields(self):
        report = validate_tracker(load_tracker())
        self.assertTrue(report["valid"], report["errors"])
        self.assertEqual(report["curve_count"], 17)
        self.assertEqual(report["status_counts"]["pending"], 17)
        self.assertEqual(report["unresolved_curve_count"], 17)
        self.assertFalse(report["mechanistic_fit_authorized"])

    def test_complete_record_cannot_retain_unknown_fields(self):
        rows = load_tracker()
        rows[0] = dict(rows[0], status="complete")
        report = validate_tracker(rows)
        self.assertFalse(report["valid"])
        self.assertTrue(
            any("complete row still has unknown" in error for error in report["errors"])
        )

    def test_structural_qc_preserves_the_single_observed_decrease(self):
        with (OUTPUT_ROOT / "salt_transport_tidy.csv").open(
            "r", encoding="utf-8-sig", newline=""
        ) as handle:
            qc_rows = structural_qc(list(csv.DictReader(handle)))
        self.assertEqual(len(qc_rows), 17)
        self.assertTrue(all(row["time_strictly_increasing"] for row in qc_rows))
        self.assertTrue(all(row["nonpositive_conductivity_count"] == 0 for row in qc_rows))
        decreases = [row for row in qc_rows if row["conductivity_decrease_count"]]
        self.assertEqual(len(decreases), 1)
        self.assertTrue(decreases[0]["curve_id"].endswith(":700 nm"))
        self.assertEqual(decreases[0]["largest_conductivity_decrease"], -20.0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
