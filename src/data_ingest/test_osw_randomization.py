import unittest

from osw_randomization import assign_treatments, validate_balance


class OswRandomizationTests(unittest.TestCase):
    def test_assignment_is_reproducible_and_balanced(self):
        coupons = [f"M1-{index}" for index in range(1, 7)]
        first = assign_treatments(coupons, seed=20261016)
        second = assign_treatments(coupons, seed=20261016)
        self.assertEqual(first, second)
        self.assertTrue(validate_balance(first))

    def test_assignment_does_not_depend_on_input_results(self):
        coupons = ["M1-left", "M1-center", "M1-right"]
        assignment = assign_treatments(coupons, seed=7)
        self.assertEqual(set(assignment), set(coupons))
        self.assertEqual(set(assignment.values()), {"O", "S", "W"})

    def test_assignment_rejects_duplicate_or_too_few_coupons(self):
        with self.assertRaises(ValueError):
            assign_treatments(["a", "a", "b"], seed=1)
        with self.assertRaises(ValueError):
            assign_treatments(["a", "b"], seed=1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
