import unittest

from validate_independent_readouts import (
    FIELDS,
    TEMPLATE_PATH,
    load_registry,
    validate_registry,
)


def complete_row(role: str) -> dict[str, str]:
    row = {field: "recorded" for field in FIELDS}
    row.update(
        {
            "method_id": f"method-{role}",
            "readout_role": role,
            "input_level": "provisional",
            "status": "ready",
            "exclusion_reason": "not_applicable",
        }
    )
    return row


class IndependentReadoutTests(unittest.TestCase):
    def test_template_has_frozen_columns_and_no_fake_rows(self):
        self.assertEqual(load_registry(TEMPLATE_PATH), [])

    def test_minimum_readout_set_is_structurally_ready(self):
        rows = [
            complete_row("independent_na"),
            complete_row("independent_cl"),
            complete_row("neutral_tracer"),
        ]
        report = validate_registry(rows)
        self.assertTrue(report["valid"], report["errors"])
        self.assertTrue(report["selection_ready"], report["warnings"])
        self.assertFalse(report["experiment_authorized"])

    def test_missing_tracer_is_reported(self):
        rows = [complete_row("independent_na"), complete_row("independent_cl")]
        report = validate_registry(rows)
        self.assertFalse(report["selection_ready"])
        self.assertTrue(any("neutral_tracer" in item for item in report["warnings"]))

    def test_ready_tracer_requires_charge_and_recovery_checks(self):
        row = complete_row("neutral_tracer")
        row["charge_state_basis"] = "unknown"
        report = validate_registry([row])
        self.assertFalse(report["valid"])
        self.assertTrue(any("ready neutral tracer" in item for item in report["errors"]))


if __name__ == "__main__":
    unittest.main(verbosity=2)
