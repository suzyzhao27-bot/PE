import unittest

from validate_quality_thresholds import (
    FIELDS,
    REQUIRED_ROLES,
    TEMPLATE_PATH,
    evaluate_metrics,
    evaluate_threshold,
    load_registry,
    validate_registry,
)


def row(role: str) -> dict[str, str]:
    result = {field: "recorded" for field in FIELDS}
    result.update(
        {
            "threshold_id": f"threshold-{role}",
            "metric_role": role,
            "metric_name": role,
            "comparison": "abs_le",
            "lower_bound": "not_applicable",
            "upper_bound": "1",
            "unit": "declared_unit",
            "action": "review",
            "input_level": "provisional",
            "status": "ready",
        }
    )
    return result


class QualityThresholdTests(unittest.TestCase):
    def test_template_has_frozen_columns_and_no_fake_rows(self):
        self.assertEqual(load_registry(TEMPLATE_PATH), [])

    def test_minimum_threshold_set_is_ready(self):
        report = validate_registry([row(role) for role in sorted(REQUIRED_ROLES)])
        self.assertTrue(report["valid"], report["errors"])
        self.assertTrue(report["threshold_set_ready"], report["warnings"])
        self.assertFalse(report["experiment_authorized"])

    def test_missing_role_is_reported(self):
        report = validate_registry([row("mass_closure_error")])
        self.assertFalse(report["threshold_set_ready"])
        self.assertTrue(report["warnings"])

    def test_planned_threshold_may_keep_unknown_bounds(self):
        item = row("mass_closure_error")
        item.update({"lower_bound": "unknown", "upper_bound": "unknown", "status": "planned"})
        report = validate_registry([item])
        self.assertTrue(report["valid"], report["errors"])

    def test_between_rejects_reversed_bounds(self):
        item = row("temperature_deviation")
        item.update({"comparison": "between", "lower_bound": "2", "upper_bound": "1"})
        report = validate_registry([item])
        self.assertFalse(report["valid"])
        self.assertTrue(any("exceeds" in item for item in report["errors"]))

    def test_absolute_limit_is_applied_inclusively(self):
        self.assertTrue(evaluate_threshold(-1.0, "abs_le", None, 1.0))
        self.assertFalse(evaluate_threshold(1.1, "abs_le", None, 1.0))

    def test_between_limit_is_applied_inclusively(self):
        self.assertTrue(evaluate_threshold(2.0, "between", 1.0, 2.0))
        self.assertFalse(evaluate_threshold(2.1, "between", 1.0, 2.0))

    def test_threshold_application_rejects_missing_bound(self):
        with self.assertRaisesRegex(ValueError, "requires upper_bound"):
            evaluate_threshold(1.0, "le", None, None)

    def test_ready_threshold_set_evaluates_all_metrics(self):
        rows = [row(role) for role in sorted(REQUIRED_ROLES)]
        metrics = {role: 0.5 for role in REQUIRED_ROLES}
        results = evaluate_metrics(metrics, rows)
        self.assertEqual(len(results), len(REQUIRED_ROLES))
        self.assertTrue(all(result["decision"] == "pass" for result in results))

    def test_failed_metric_uses_predeclared_action(self):
        rows = [row(role) for role in sorted(REQUIRED_ROLES)]
        target = next(item for item in rows if item["metric_role"] == "mass_closure_error")
        target["action"] = "exclude"
        metrics = {role: 0.5 for role in REQUIRED_ROLES}
        metrics["mass_closure_error"] = 2.0
        results = evaluate_metrics(metrics, rows)
        decision = next(
            item["decision"]
            for item in results
            if item["metric_role"] == "mass_closure_error"
        )
        self.assertEqual(decision, "exclude")

    def test_metric_evaluation_rejects_missing_role_value(self):
        rows = [row(role) for role in sorted(REQUIRED_ROLES)]
        with self.assertRaisesRegex(ValueError, "Missing metric value"):
            evaluate_metrics({}, rows)


if __name__ == "__main__":
    unittest.main(verbosity=2)
