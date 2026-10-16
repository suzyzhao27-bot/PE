import csv
import json
import unittest
from pathlib import Path

from analyze_salt_transport import write_exact_duplicates
from normalize_salt_transport import CURVES, EXPORT_ROOT, OUTPUT_ROOT, normalize_curve


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


if __name__ == "__main__":
    unittest.main(verbosity=2)
