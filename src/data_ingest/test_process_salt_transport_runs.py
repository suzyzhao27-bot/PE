import unittest

from conductivity_calibration import FIELDS as CALIBRATION_FIELDS
from process_salt_transport_runs import process_runs
from validate_salt_transport_runs import FIELDS as RUN_FIELDS


def calibration_row(standard: str, concentration: str, conductivity: str) -> dict[str, str]:
    result = {field: "recorded" for field in CALIBRATION_FIELDS}
    result.update(
        {
            "calibration_id": "cal-1",
            "standard_id": standard,
            "concentration": concentration,
            "concentration_unit": "mmol/L",
            "conductivity": conductivity,
            "conductivity_unit": "uS/cm",
            "temperature": "25",
            "temperature_unit": "degC",
            "temperature_compensation": "none",
            "input_level": "provisional",
        }
    )
    return result


def run_row(time: str, left: str, right: str) -> dict[str, str]:
    result = {field: "recorded" for field in RUN_FIELDS}
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
            "left_conductivity": left,
            "right_conductivity": right,
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


class ProcessSaltTransportRunTests(unittest.TestCase):
    def setUp(self):
        self.calibration = [
            calibration_row("s0", "0", "0"),
            calibration_row("s1", "10", "100"),
        ]

    def test_conversion_closes_synthetic_two_pool_balance(self):
        output = process_runs(
            [run_row("0", "100", "0"), run_row("10", "90", "10")],
            self.calibration,
        )
        self.assertEqual(len(output), 2)
        self.assertAlmostEqual(output[1]["relative_closure_error"], 0.0)
        self.assertAlmostEqual(output[1]["corrected_receiver_transfer_mol"], 1e-5)
        self.assertAlmostEqual(
            output[1]["cumulative_average_flux_mol_m2_s"], 1.0 / 6000.0
        )
        self.assertFalse(output[1]["mechanistic_fit_authorized"])

    def test_temperature_must_match_calibration(self):
        rows = [run_row("0", "100", "0"), run_row("10", "90", "10")]
        rows[1]["temperature"] = "24"
        with self.assertRaisesRegex(ValueError, "temperature differs"):
            process_runs(rows, self.calibration)

    def test_sampling_is_added_back_at_next_timepoint(self):
        rows = [
            run_row("0", "100", "0"),
            run_row("10", "90", "10"),
            run_row("20", "90", "10"),
        ]
        rows[1]["right_sampled_volume"] = "1"
        rows[2]["right_volume"] = "9"
        output = process_runs(rows, self.calibration)
        self.assertAlmostEqual(
            output[2]["corrected_receiver_transfer_mol"],
            output[1]["corrected_receiver_transfer_mol"],
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
