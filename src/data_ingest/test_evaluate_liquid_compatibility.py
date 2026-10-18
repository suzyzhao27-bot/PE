import unittest

from evaluate_liquid_compatibility import absolute_percent_change, evaluate_records
from test_liquid_compatibility import row as compatibility_row
from test_liquid_compatibility_thresholds import threshold_row


def matching_rows() -> tuple[dict[str, str], dict[str, str]]:
    record = compatibility_row()
    threshold = threshold_row()
    shared = {
        "protocol_liquid_identity": "water",
        "protocol_liquid_batch": "batch-1",
        "exposure_duration": "60",
        "duration_unit": "min",
        "temperature": "25",
        "temperature_unit": "degC",
        "agitation": "static",
        "measurement_protocol_id": "method-1",
    }
    record.update(shared)
    threshold.update(shared)
    record["compatibility_threshold_id"] = threshold["compatibility_threshold_id"]
    return record, threshold


class EvaluateLiquidCompatibilityTests(unittest.TestCase):
    def test_absolute_percent_change_uses_pre_exposure_denominator(self):
        self.assertAlmostEqual(absolute_percent_change("10", "9.5"), 5.0)

    def test_matching_record_within_limits_is_compatible(self):
        record, threshold = matching_rows()
        report = evaluate_records([record], [threshold])
        self.assertTrue(report["valid"], report["errors"])
        self.assertEqual(report["evaluated_ready_compatible_ids"], ["compat-1"])

    def test_failed_limit_uses_predeclared_failure_decision(self):
        record, threshold = matching_rows()
        record.update({"post_mass": "1.2", "decision": "review"})
        report = evaluate_records([record], [threshold])
        self.assertTrue(report["valid"], report["errors"])
        self.assertEqual(report["results"][0]["derived_decision"], "review")
        self.assertFalse(report["compatibility_ready"])

    def test_condition_mismatch_is_rejected(self):
        record, threshold = matching_rows()
        record["protocol_liquid_batch"] = "wrong-batch"
        report = evaluate_records([record], [threshold])
        self.assertFalse(report["valid"])

    def test_recorded_decision_cannot_override_thresholds(self):
        record, threshold = matching_rows()
        record.update({"post_width": "12", "decision": "compatible"})
        report = evaluate_records([record], [threshold])
        self.assertFalse(report["valid"])
        self.assertTrue(any("does not match" in item for item in report["errors"]))


if __name__ == "__main__":
    unittest.main(verbosity=2)
