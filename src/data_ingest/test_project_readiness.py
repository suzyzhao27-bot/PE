import unittest

from assess_project_readiness import assess


class ProjectReadinessTests(unittest.TestCase):
    def test_current_repository_state_is_reported_without_authorization(self):
        report = assess()
        self.assertEqual(report["salt_transport"]["curve_count"], 17)
        self.assertEqual(report["salt_transport"]["unresolved_curve_count"], 17)
        self.assertEqual(report["osw_protocol"]["row_count"], 0)
        self.assertEqual(report["conductivity_calibration"]["row_count"], 0)
        self.assertEqual(report["salt_transport_runs"]["row_count"], 0)
        self.assertEqual(report["blanks"]["row_count"], 0)
        self.assertEqual(report["quality_thresholds"]["row_count"], 0)
        self.assertEqual(report["independent_readouts"]["row_count"], 0)
        self.assertFalse(report["mechanistic_fit_authorized"])
        self.assertFalse(report["experiment_authorized"])

    def test_current_repository_state_lists_seven_next_actions(self):
        report = assess()
        self.assertEqual(len(report["next_actions"]), 7)


if __name__ == "__main__":
    unittest.main(verbosity=2)
