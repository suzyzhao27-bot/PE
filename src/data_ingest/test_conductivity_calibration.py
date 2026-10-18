import unittest

from conductivity_calibration import (
    FIELDS,
    TEMPLATE_PATH,
    CalibrationPoint,
    interpolate_concentration,
    load_registry,
    validate_registry,
)


def row(standard_id: str, concentration: str, conductivity: str) -> dict[str, str]:
    result = {field: "recorded" for field in FIELDS}
    result.update(
        {
            "calibration_id": "cal-1",
            "standard_id": standard_id,
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


class ConductivityCalibrationTests(unittest.TestCase):
    def test_template_has_frozen_columns_and_no_fake_rows(self):
        self.assertEqual(load_registry(TEMPLATE_PATH), [])

    def test_monotonic_table_is_ready_and_interpolates(self):
        rows = [row("s0", "0", "1"), row("s1", "10", "3")]
        report = validate_registry(rows)
        self.assertTrue(report["conversion_ready"], report["errors"])
        points = [CalibrationPoint(0.0, 1.0), CalibrationPoint(10.0, 3.0)]
        self.assertAlmostEqual(interpolate_concentration(2.0, points), 5.0)

    def test_nonmonotonic_table_is_rejected(self):
        rows = [row("s0", "10", "1"), row("s1", "5", "3")]
        report = validate_registry(rows)
        self.assertFalse(report["valid"])

    def test_interpolation_rejects_extrapolation(self):
        points = [CalibrationPoint(0.0, 1.0), CalibrationPoint(10.0, 3.0)]
        with self.assertRaisesRegex(ValueError, "outside the calibration range"):
            interpolate_concentration(4.0, points)


if __name__ == "__main__":
    unittest.main(verbosity=2)
