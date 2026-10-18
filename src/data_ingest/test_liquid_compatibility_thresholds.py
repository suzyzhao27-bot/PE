import unittest

from validate_liquid_compatibility_thresholds import (
    FIELDS,
    TEMPLATE_PATH,
    load_registry,
    validate_registry,
)


def threshold_row() -> dict[str, str]:
    row = {field: "recorded" for field in FIELDS}
    row.update(
        {
            "compatibility_threshold_id": "compat-threshold-1",
            "exposure_duration": "60",
            "temperature": "25",
            "max_abs_mass_change_pct": "5",
            "max_abs_length_change_pct": "2",
            "max_abs_width_change_pct": "2",
            "visual_damage_allowed": "no",
            "failure_decision": "review",
            "input_level": "provisional",
            "status": "ready",
        }
    )
    return row


class LiquidCompatibilityThresholdTests(unittest.TestCase):
    def test_template_has_frozen_columns_and_no_fake_rows(self):
        self.assertEqual(load_registry(TEMPLATE_PATH), [])

    def test_ready_threshold_is_reported(self):
        report = validate_registry([threshold_row()])
        self.assertTrue(report["valid"], report["errors"])
        self.assertTrue(report["thresholds_ready"])

    def test_ready_threshold_rejects_unknown_limit(self):
        item = threshold_row()
        item["max_abs_mass_change_pct"] = "unknown"
        report = validate_registry([item])
        self.assertFalse(report["valid"])

    def test_negative_change_limit_is_rejected(self):
        item = threshold_row()
        item["max_abs_width_change_pct"] = "-1"
        report = validate_registry([item])
        self.assertFalse(report["valid"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
