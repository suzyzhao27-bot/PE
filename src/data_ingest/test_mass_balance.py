import unittest

from mass_balance import (
    corrected_system_amount,
    interval_flux,
    relative_closure_error,
    update_volume,
)


class MassBalanceTests(unittest.TestCase):
    def test_volume_tracks_sampling_and_replenishment(self):
        self.assertAlmostEqual(update_volume(10.0, 1.5, 1.0), 9.5)

    def test_corrected_amount_restores_removed_sample(self):
        current = corrected_system_amount(
            left_concentration=2.0,
            left_volume=4.0,
            right_concentration=1.0,
            right_volume=5.0,
            removed=[(1.0, 1.0)],
        )
        self.assertAlmostEqual(current, 14.0)
        self.assertAlmostEqual(relative_closure_error(current, 14.0), 0.0)

    def test_replenishment_amount_is_not_counted_as_membrane_transfer(self):
        current = corrected_system_amount(
            left_concentration=2.0,
            left_volume=4.0,
            right_concentration=1.0,
            right_volume=6.0,
            replenished=[(1.0, 1.0)],
        )
        self.assertAlmostEqual(current, 13.0)

    def test_flux_requires_positive_area_and_time(self):
        self.assertAlmostEqual(interval_flux(6.0, 2.0, 3.0), 1.0)
        with self.assertRaises(ValueError):
            interval_flux(6.0, 0.0, 3.0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
