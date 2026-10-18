import unittest

from validate_blank_registry import FIELDS, TEMPLATE_PATH, load_registry, validate_registry


def row(blank_id: str, blank_type: str, applies_to: str) -> dict[str, str]:
    result = {field: "recorded" for field in FIELDS}
    result.update(
        {
            "blank_id": blank_id,
            "blank_type": blank_type,
            "applies_to": applies_to,
            "left_volume": "10",
            "right_volume": "10",
            "temperature": "25",
            "duration": "120",
            "input_level": "provisional",
            "status": "ready",
            "exclusion_reason": "not_applicable",
        }
    )
    return result


class BlankRegistryTests(unittest.TestCase):
    def test_template_has_frozen_columns_and_no_fake_rows(self):
        self.assertEqual(load_registry(TEMPLATE_PATH), [])

    def test_minimum_salt_blank_set_is_ready(self):
        rows = [
            row("blank-open", "no_membrane", "salt_transport"),
            row("blank-sealed", "impermeable", "salt_transport"),
        ]
        report = validate_registry(rows)
        self.assertTrue(report["valid"], report["errors"])
        self.assertTrue(report["salt_blank_set_ready"])
        self.assertFalse(report["experiment_authorized"])

    def test_missing_impermeable_blank_is_reported(self):
        report = validate_registry(
            [row("blank-open", "no_membrane", "salt_transport")]
        )
        self.assertFalse(report["salt_blank_set_ready"])
        self.assertTrue(any("impermeable" in item for item in report["warnings"]))

    def test_ready_blank_cannot_retain_unknown_fields(self):
        item = row("blank-open", "no_membrane", "salt_transport")
        item["fixture_id"] = "unknown"
        report = validate_registry([item])
        self.assertFalse(report["valid"])
        self.assertTrue(any("ready row still has unknown" in item for item in report["errors"]))

    def test_planned_blank_may_keep_unknown_numeric_fields(self):
        item = row("blank-open", "no_membrane", "salt_transport")
        item.update(
            {
                "left_volume": "unknown",
                "right_volume": "unknown",
                "temperature": "unknown",
                "duration": "unknown",
                "status": "planned",
            }
        )
        report = validate_registry([item])
        self.assertTrue(report["valid"], report["errors"])
        self.assertFalse(report["salt_blank_set_ready"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
