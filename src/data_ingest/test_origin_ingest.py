import csv
import json
import unittest
from pathlib import Path

from analyze_salt_transport import (structural_qc, write_exact_duplicates,
                                    write_matching_conductivity_sequences)
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

    def test_book1_times_follow_saved_origin_graph_bindings(self):
        evidence = json.loads(
            (OUTPUT_ROOT / "polymer_salt_2h_plot_bindings.json").read_text(encoding="utf-8")
        )
        graph1 = {plot["dataset"]: plot["x_range"] for plot in evidence["plots"]
                  if plot["graph"] == "Graph1"}
        expected = {"65 nm": ("B", "A"), "120 nm": ("D", "C"),
                    "200 nm": ("E", "C"), "350 nm": ("F", "C"),
                    "700 nm": ("G", "C"), "2600 nm": ("H", "C"),
                    "PC filter, 75 μm thick, 200 nm pore": ("I", "C")}
        for sample, (y_column, x_column) in expected.items():
            curve = next(curve for curve in CURVES
                         if curve.project == "Polymer_salt_2h"
                         and curve.workbook == "Book1" and curve.sample_id == sample)
            self.assertEqual(curve.time_column, f"col_{ord(x_column) - 64}_{x_column}")
            self.assertIn(f"[Book1]Sheet1!{x_column}", graph1[f"Book1_{y_column}"])
            self.assertEqual(normalize_curve(curve)[4]["time"],
                             85.0 if sample == "65 nm" else 80.0)

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
        self.assertEqual(len(duplicates), 1)
        copies = write_matching_conductivity_sequences(rows)
        self.assertEqual(len(copies), 3)
        self.assertEqual(sum(not copy["same_time_grid"] for copy in copies), 2)
        self.assertTrue(all(copy["mismatched_point_indices"] == "5"
                            for copy in copies if not copy["same_time_grid"]))

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
