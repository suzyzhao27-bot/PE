#!/usr/bin/env python3
"""Run one reviewable scenario from JSON. No experimental fit or claim upgrade."""
import argparse
import hashlib
import json
import platform
import time
from pathlib import Path
from models import (Geometry, ambipolar_diffusivity, cylindrical_pb,
                    neutral_curve, neutral_resistances, reduced_curve)
from run_demo import write_csv


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    config_path = args.config.resolve()
    config = json.loads(config_path.read_text(encoding="utf-8"))
    if config["evidence_status"] not in ("synthetic", "hypothesis_calculation"):
        raise ValueError("evidence_status must be synthetic or hypothesis_calculation; calculation cannot claim experiment validation")
    out = args.out.resolve()
    # Require a fresh directory so a previous reviewable scenario is not overwritten.
    out.mkdir(parents=True, exist_ok=False)
    started = time.perf_counter()
    g = Geometry(**config["geometry"]).validate()
    dp, dm = config["d_plus_m2_s"], config["d_minus_m2_s"]
    d = ambipolar_diffusivity(dp, dm)
    cl, cr = config["c_left0_mol_m3"], config["c_right0_mol_m3"]
    times = config["times_s"]
    sigma = config["sigma_c_m2"]
    pb = cylindrical_pb(config["pb_reference_bulk_mol_m3"], sigma, g.radius_m,
                        config["pb_cells"], config["temperature_k"], config["epsilon_r"])
    neutral = neutral_curve(g, cl, cr, times, d)
    reduced = reduced_curve(g, cl, cr, times, sigma, dp, dm,
                            max_step_s=config["rk4_max_step_s"])
    write_csv(out / "neutral_exact_curve.csv", neutral)
    write_csv(out / "donnan_reduced_curve.csv", reduced)
    write_csv(out / "radial_pb_profile.csv", pb.pop("profile"))
    (out / "input_snapshot.json").write_bytes(config_path.read_bytes())
    here = Path(__file__).resolve().parent
    result = dict(scenario_id=config["scenario_id"], evidence_status=config["evidence_status"],
                  scope="model scenario only; NOT experiment validation/finite-length PNP",
                  python=platform.python_version(), runtime_seconds=time.perf_counter() - started,
                  config_sha256=sha256(config_path),
                  code_sha256={name: sha256(here / name) for name in ("models.py", "run_scenario.py", "run_demo.py")},
                  neutral_resistances=neutral_resistances(g, d), radial_pb_summary=pb)
    (out / "scenario_manifest.json").write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(dict(output=str(out), evidence_status=config["evidence_status"],
                          config_sha256=result["config_sha256"], runtime_seconds=result["runtime_seconds"]), ensure_ascii=False))


if __name__ == "__main__":
    main()
