import unittest

from validate_salt_transport_runs import (
    FIELDS,
    TEMPLATE_PATH,
    load_registry,
    validate_registry,
)


def row(time: str) -> dict[str, str]:
    result = {field: "recorded" for field in FIELDS}
    result.update(
        {
            "run_id": "run-1",
            "coupon_id": "M1-O-1",
            "parent_membrane_id": "M1",
            "treatment_group": "O",
            "receiver_side": "right",
            "time": time,
            "time_unit": "min",
            "record_phase": "pre_event",
            "left_conductivity": "100",
            "right_conductivity": "2",
            "conductivity_unit": "uS/cm",
            "calibration_id": "cal-1",
            "left_volume": "10",
            "right_volume": "10",
            "volume_unit": "mL",
            "left_sampled_volume": "0",
            "right_sampled_volume": "0",
            "left_replenished_volume": "0",
            "right_replenished_volume": "0",
            "left_replenishment_concentration": "0",
            "right_replenishment_concentration": "0",
            "concentration_unit": "mmol/L",
            "temperature": "25",
            "temperature_unit": "degC",
            "effective_area": "1",
            "area_unit": "cm2",
            "blank_id": "blank-1",
            "fixture_id": "fixture-1",
            "leakage_flag": "no",
            "input_level": "provisional",
        }
    )
    return result


class SaltTransportRunTests(unittest.TestCase):
    def test_template_has_frozen_columns_and_no_fake_rows(self):
        self.assertEqual(load_registry(TEMPLATE_PATH), [])

    def test_two_point_run_is_structurally_ready(self):
        report = validate_registry([row("0"), row("10")])
        self.assertTrue(report["run_log_ready"], report["errors"])
        self.assertFalse(report["mechanistic_fit_authorized"])

    def test_duplicate_time_is_rejected(self):
        report = validate_registry([row("0"), row("0")])
        self.assertFalse(report["valid"])
        self.assertTrue(any("duplicate time" in item for item in report["errors"]))

    def test_mixed_area_is_rejected(self):
        first = row("0")
        second = row("10")
        second["effective_area"] = "2"
        report = validate_registry([first, second])
        self.assertFalse(report["valid"])
        self.assertTrue(any("mixed effective_area" in item for item in report["errors"]))


if __name__ == "__main__":
    unittest.main(verbosity=2)
