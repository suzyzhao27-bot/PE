import unittest

from validate_liquid_compatibility import (
    FIELDS,
    TEMPLATE_PATH,
    load_registry,
    validate_registry,
)


def row() -> dict[str, str]:
    result = {field: "recorded" for field in FIELDS}
    result.update(
        {
            "compatibility_record_id": "compat-1",
            "exposure_duration": "60",
            "temperature": "25",
            "pre_mass": "1",
            "post_mass": "1",
            "pre_length": "10",
            "post_length": "10",
            "pre_width": "10",
            "post_width": "10",
            "visual_damage": "no",
            "input_level": "provisional",
            "status": "ready",
            "decision": "compatible",
            "exclusion_reason": "not_applicable",
        }
    )
    return result


class LiquidCompatibilityTests(unittest.TestCase):
    def test_template_has_frozen_columns_and_no_fake_rows(self):
        self.assertEqual(load_registry(TEMPLATE_PATH), [])

    def test_ready_compatible_record_is_reported(self):
        report = validate_registry([row()])
        self.assertTrue(report["valid"], report["errors"])
        self.assertEqual(report["ready_compatible_ids"], ["compat-1"])
        self.assertFalse(report["experiment_authorized"])

    def test_ready_record_requires_decision(self):
        item = row()
        item["decision"] = "unknown"
        report = validate_registry([item])
        self.assertFalse(report["valid"])

    def test_planned_record_may_keep_unknown_measurements(self):
        item = row()
        for field in ("pre_mass", "post_mass", "pre_length", "post_length"):
            item[field] = "unknown"
        item.update({"status": "planned", "decision": "unknown"})
        report = validate_registry([item])
        self.assertTrue(report["valid"], report["errors"])
        self.assertFalse(report["compatibility_ready"])

    def test_duplicate_record_id_is_rejected(self):
        report = validate_registry([row(), row()])
        self.assertFalse(report["valid"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
