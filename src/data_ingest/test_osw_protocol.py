import unittest

from validate_osw_protocol import FIELDS, TEMPLATE_PATH, load_registry, validate_registry


def complete_row(parent: str, treatment: str, coupon: str) -> dict[str, str]:
    row = {field: "recorded" for field in FIELDS}
    row.update(
        {
            "parent_membrane_id": parent,
            "coupon_id": coupon,
            "treatment_group": treatment,
            "input_level": "provisional",
            "status": "ready",
            "exclusion_reason": "not_applicable",
        }
    )
    return row


class OswProtocolTests(unittest.TestCase):
    def test_template_has_frozen_columns_and_no_fake_rows(self):
        self.assertEqual(TEMPLATE_PATH.name, "osw_sample_registry_template.csv")
        rows = load_registry(TEMPLATE_PATH)
        self.assertEqual(rows, [])

    def test_two_parent_balanced_design_is_structurally_ready(self):
        rows = [
            complete_row(parent, treatment, f"{parent}-{treatment}")
            for parent in ("M1", "M2")
            for treatment in ("O", "S", "W")
        ]
        report = validate_registry(rows)
        self.assertTrue(report["valid"], report["errors"])
        self.assertTrue(report["protocol_ready"], report["warnings"])
        self.assertFalse(report["experiment_authorized"])

    def test_missing_treatment_is_reported(self):
        rows = [
            complete_row("M1", treatment, f"M1-{treatment}")
            for treatment in ("O", "S")
        ]
        report = validate_registry(rows)
        self.assertFalse(report["protocol_ready"])
        self.assertTrue(any("missing active treatment" in item for item in report["warnings"]))

    def test_ready_row_cannot_retain_unknown(self):
        row = complete_row("M1", "O", "M1-O")
        row["temperature"] = "unknown"
        report = validate_registry([row])
        self.assertFalse(report["valid"])
        self.assertTrue(any("ready row still has unknown" in item for item in report["errors"]))


if __name__ == "__main__":
    unittest.main(verbosity=2)
