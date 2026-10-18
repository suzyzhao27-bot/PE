"""Test whether one r=0 field can satisfy selected old thickness ratios.

Solve only the 350/200 and 2600/700 endpoint ratios for shared ws and ell.
Other ratios are diagnostics, not fitted targets. No point-specific parameters.
"""

import argparse
import csv
import json
import math
import sys
from pathlib import Path

from scipy.optimize import least_squares


ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "core"))
from theory_model import log_resistance_nm  # noqa: E402


def observations():
    source = ROOT / "core/inputs/Book1__Sheet1.csv"
    with source.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 7:
        raise ValueError("Expected the frozen seven-row Book1 export")
    mapping = ((65, "col_2_B"), (120, "col_4_D"), (200, "col_5_E"),
               (350, "col_6_F"), (700, "col_7_G"), (2600, "col_8_H"))
    result = {h: float(rows[-1][column]) - float(rows[0][column]) for h, column in mapping}
    if any(not math.isfinite(value) or value <= 0 for value in result.values()):
        raise ValueError("All six 0–120 min increments must be positive and finite")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tradeoff-output", type=Path,
                        help="Optional comparison of three shared-parameter choices")
    args = parser.parse_args()
    params = json.loads((ROOT / "core/inputs/model_parameters.json").read_text(encoding="utf-8"))
    if float(params["robin_r"]) != 0:
        raise ValueError("This check preserves the frozen r=0 field family")
    a = float(params["a_ref_nm"])
    observed = observations()
    targets = ((200, 350), (700, 2600))

    def residual(log_parameters, n):
        ws, ell = map(math.exp, log_parameters)
        log_r = {h: log_resistance_nm(h, ws, ell, a, n=n)
                 for pair in targets for h in pair}
        return [log_r[left] - log_r[right] - math.log(observed[right] / observed[left])
                for left, right in targets]

    fit = least_squares(lambda x: residual(x, 2048),
                        [math.log(params["ws_ref"]), math.log(params["ell_nm"])],
                        bounds=([math.log(.01), math.log(.1)],
                                [math.log(100), math.log(10000)]),
                        xtol=1e-12, ftol=1e-12, gtol=1e-12)
    if not fit.success or max(abs(value) for value in residual(fit.x, 4096)) > 1e-8:
        raise RuntimeError("The targeted two-ratio solution did not pass fine integration")
    ws, ell = map(math.exp, fit.x)
    log_r = {h: log_resistance_nm(h, ws, ell, a, n=4096) for h in observed}
    reference_log_r = log_r[200]
    rows = []
    for h, increment in observed.items():
        predicted_ratio = math.exp(reference_log_r - log_r[h])
        observed_ratio = increment / observed[200]
        rows.append(dict(h_nm=h, observed_delta_kappa=increment,
                         observed_ratio_to_200=observed_ratio,
                         model_ratio_to_200=predicted_ratio,
                         model_ratio_over_observed_ratio=predicted_ratio / observed_ratio,
                         role="ratio_constraint" if h in (200, 350, 700, 2600) else "unused_endpoint",
                         ws_kBT=ws, ell_nm=ell, robin_r=0, a_nm=a))
    with args.output.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    middle_model = math.exp(log_r[350] - log_r[700])
    middle_observed = observed[700] / observed[350]
    print(f"ws={ws:.9f}, ell_nm={ell:.9f}")
    for left, right in targets:
        print(f"{right}/{left}: model={math.exp(log_r[left]-log_r[right]):.9f}, "
              f"old={observed[right]/observed[left]:.9f}")
    print(f"700/350 unused: model={middle_model:.9f}, old={middle_observed:.9f}, "
          f"model/old={middle_model/middle_observed:.9f}")

    if args.tradeoff_output:
        pairs = ((65, 120), (200, 350), (350, 700), (700, 2600))
        fitted_pairs = pairs[1:]

        def three_ratio_residual(log_parameters, n):
            fit_ws, fit_ell = map(math.exp, log_parameters)
            resistances = {h: log_resistance_nm(h, fit_ws, fit_ell, a, n=n)
                           for pair in fitted_pairs for h in pair}
            return [resistances[left] - resistances[right] -
                    math.log(observed[right] / observed[left])
                    for left, right in fitted_pairs]

        starts = ([math.log(params["ws_ref"]), math.log(params["ell_nm"])],
                  list(fit.x), [math.log(3), math.log(100)])
        candidates = [least_squares(lambda x: three_ratio_residual(x, 2048), start,
                                    bounds=([math.log(.01), math.log(.1)],
                                            [math.log(100), math.log(10000)]),
                                    xtol=1e-12, ftol=1e-12, gtol=1e-12)
                      for start in starts]
        successful = [candidate for candidate in candidates if candidate.success]
        if not successful:
            raise RuntimeError("The three-ratio fit did not converge")
        compromise = min(successful, key=lambda candidate: candidate.cost)
        three_ws, three_ell = map(math.exp, compromise.x)
        choices = (("reference_six_endpoint", float(params["ws_ref"]), float(params["ell_nm"])),
                   ("two_ratio_constraints", ws, ell),
                   ("three_ratio_compromise", three_ws, three_ell))
        diagnostics = []
        for name, choice_ws, choice_ell in choices:
            resistances = {h: log_resistance_nm(h, choice_ws, choice_ell, a, n=4096)
                           for pair in pairs for h in pair}
            for left, right in pairs:
                old_ratio = observed[right] / observed[left]
                model_ratio = math.exp(resistances[left] - resistances[right])
                diagnostics.append(dict(model_choice=name, left_h_nm=left, right_h_nm=right,
                                        ratio_role="unused" if (left, right) == (65, 120) else "three_adjacent",
                                        old_ratio=old_ratio, model_ratio=model_ratio,
                                        model_over_old=model_ratio / old_ratio,
                                        log_ratio_residual=math.log(model_ratio / old_ratio),
                                        ws_kBT=choice_ws, ell_nm=choice_ell, robin_r=0, a_nm=a))
        with args.tradeoff_output.open("w", encoding="utf-8-sig", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(diagnostics[0]))
            writer.writeheader()
            writer.writerows(diagnostics)
        for name, _, _ in choices:
            relevant = [row for row in diagnostics if row["model_choice"] == name and
                        row["ratio_role"] == "three_adjacent"]
            norm = math.sqrt(sum(row["log_ratio_residual"] ** 2 for row in relevant))
            ratios = ",".join(f"{row['model_over_old']:.6f}" for row in relevant)
            print(f"{name}: three-adjacent log-ratio norm={norm:.9f}; "
                  f"model/old={ratios}")


if __name__ == "__main__":
    main()
