import csv
import hashlib
import json
import platform
import time
from pathlib import Path
from long_tube_np import integrate_tube

ROOT = Path(__file__).resolve().parent


def write_csv(path, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    started = time.perf_counter()
    out = ROOT / "outputs"
    out.mkdir(exist_ok=True)
    summaries, locals_ = [], []
    for c_left in (30.0, 270.0):
        for sigma in (0.0, -.01, .01):
            result = integrate_tube(c_left, 1.0, sigma, cells=512, intervals=256)
            profiles = result.pop("local_coefficients")
            summaries.append(result)
            for row in profiles:
                locals_.append(dict(c_left_virtual_mol_m3=c_left, sigma_c_m2=sigma, **row))
    write_csv(out / "scenario_summary.csv", summaries)
    write_csv(out / "local_transport_coefficients.csv", locals_)
    hashes = {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
              for name in ("pb_v1_readonly_copy.py", "long_tube_np.py", "run_diagnostic.py")}
    source_identity = json.loads((ROOT / "pb_source_identity.json").read_text(encoding="utf-8"))
    registered_hash = source_identity["registered_source_sha256"]
    if hashes["pb_v1_readonly_copy.py"] != registered_hash:
        raise ValueError("Local PB copy differs from its fixed registered v1 source SHA256")
    manifest = dict(status="restricted physical diagnostic; no real fit; no finite-length certification",
                    python=platform.python_version(), runtime_seconds=time.perf_counter() - started,
                    copied_v1_pb_source_provenance_only=source_identity["original_path_provenance_only"],
                    registered_pb_source_sha256=registered_hash,
                    copy_hash_matches_registered=registered_hash == hashes["pb_v1_readonly_copy.py"],
                    original_absolute_file_is_not_a_runtime_dependency=True,
                    source_identity_sha256=hashlib.sha256((ROOT / "pb_source_identity.json").read_bytes()).hexdigest(),
                    code_sha256=hashes, assumptions=dict(radius_nm=25, length_for_flux_normalization_nm=350,
                    length_is_not_validity_certification=True, no_flow=True, zero_net_current=True,
                    no_pore_ends=True, no_external_layers=True, no_membrane_faces=True,
                    endpoint_concentrations_are_virtual_not_measured_reservoir_after_access=True,
                    sigma_is_assumed_not_inferred_from_zeta_KPFM=True), scenarios=summaries)
    (out / "run_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(dict(runtime_seconds=manifest["runtime_seconds"], scenarios=summaries), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
