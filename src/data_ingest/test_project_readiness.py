import unittest

from assess_project_readiness import assess, compatibility_links_ready
from test_liquid_compatibility import row as compatibility_row
from test_osw_protocol import complete_row


class ProjectReadinessTests(unittest.TestCase):
    def test_current_repository_state_is_reported_without_authorization(self):
        report = assess()
        self.assertEqual(report["salt_transport"]["curve_count"], 17)
        self.assertEqual(report["salt_transport"]["unresolved_curve_count"], 17)
        self.assertEqual(report["liquid_compatibility"]["row_count"], 0)
        self.assertEqual(report["liquid_compatibility_thresholds"]["row_count"], 0)
        self.assertEqual(report["osw_protocol"]["row_count"], 0)
        self.assertEqual(report["conductivity_calibration"]["row_count"], 0)
        self.assertEqual(report["salt_transport_runs"]["row_count"], 0)
        self.assertEqual(report["blanks"]["row_count"], 0)
        self.assertEqual(report["quality_thresholds"]["row_count"], 0)
        self.assertEqual(report["independent_readouts"]["row_count"], 0)
        self.assertFalse(report["mechanistic_fit_authorized"])
        self.assertFalse(report["experiment_authorized"])

    def test_current_repository_state_lists_nine_next_actions(self):
        report = assess()
        self.assertEqual(len(report["next_actions"]), 9)

    def test_osw_requires_matching_ready_compatible_liquid(self):
        compatibility = compatibility_row()
        compatibility.update(
            {
                "protocol_liquid_identity": "water",
                "protocol_liquid_batch": "batch-1",
            }
        )
        osw_rows = [
            complete_row(parent, treatment, f"{parent}-{treatment}")
            for parent in ("M1", "M2")
            for treatment in ("O", "S", "W")
        ]
        for item in osw_rows:
            item.update(
                {
                    "compatibility_record_id": "compat-1",
                    "protocol_liquid_identity": "water",
                    "protocol_liquid_batch": "batch-1",
                }
            )
        self.assertTrue(
            compatibility_links_ready([compatibility], ["compat-1"], osw_rows)
        )

        osw_rows[0]["protocol_liquid_batch"] = "wrong-batch"
        self.assertFalse(
            compatibility_links_ready([compatibility], ["compat-1"], osw_rows)
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
