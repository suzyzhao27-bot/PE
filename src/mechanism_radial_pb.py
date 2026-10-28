"""First fixed radial PB grid point, with a persistent actual-call ledger.

Finite-volume cylindrical PB, nonuniform radial grid, damped Newton steps.
No dependency on the archived planar/fiber geometry solver.
"""

from __future__ import annotations

import csv
import json
import math
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MATRIX = ROOT / "原始资料/赵珈绮/22_机制预测与实验验证方案_v1.0.0/matrices/pb_main.json"
OUT = ROOT / "data/processed/salt_transport/22_pb_first_point_20261028.json"
LEDGER = ROOT / "data/processed/salt_transport/22_pb_call_ledger_20261028.csv"
F = 96485.33212
R = 8.314462618
T = 298.15
EPSILON_W = 78.3 * 8.8541878128e-12


def tridiagonal(lower: list[float], diagonal: list[float], upper: list[float], rhs: list[float]) -> list[float]:
    n = len(rhs)
    d = diagonal.copy()
    y = rhs.copy()
    for i in range(1, n):
        factor = lower[i] / d[i - 1]
        d[i] -= factor * upper[i - 1]
        y[i] -= factor * y[i - 1]
    result = [0.0] * n
    result[-1] = y[-1] / d[-1]
    for i in range(n - 2, -1, -1):
        result[i] = (y[i] - upper[i] * result[i + 1]) / d[i]
    return result


def geometry(cells: int) -> tuple[list[float], list[float], list[float]]:
    faces = [1 - (1 - i / cells) ** 3 for i in range(cells + 1)]
    centers = [math.sqrt((faces[i] ** 2 + faces[i + 1] ** 2) / 2) for i in range(cells)]
    volumes = [(faces[i + 1] ** 2 - faces[i] ** 2) / 2 for i in range(cells)]
    return faces, centers, volumes


def system(psi: list[float], faces: list[float], centers: list[float], volumes: list[float],
           kappa_a_squared: float, wall_slope: float) -> tuple[list[float], list[float], list[float], list[float]]:
    n = len(psi)
    residual = [0.0] * n
    lower = [0.0] * n
    diagonal = [0.0] * n
    upper = [0.0] * n
    for i in range(n):
        source = kappa_a_squared * math.sinh(psi[i]) * volumes[i]
        residual[i] = -source
        diagonal[i] = -kappa_a_squared * math.cosh(psi[i]) * volumes[i]
        if i:
            conductance = faces[i] / (centers[i] - centers[i - 1])
            residual[i] -= conductance * (psi[i] - psi[i - 1])
            lower[i] = conductance
            diagonal[i] -= conductance
        if i < n - 1:
            conductance = faces[i + 1] / (centers[i + 1] - centers[i])
            residual[i] += conductance * (psi[i + 1] - psi[i])
            upper[i] = conductance
            diagonal[i] -= conductance
        else:
            residual[i] += wall_slope
    return residual, lower, diagonal, upper


def solve(point: dict[str, object]) -> dict[str, object]:
    a = float(point["a_nm"]) * 1e-9
    sigma = float(point["sigma_C_m2"])
    c = float(point["c_mol_m3"])
    cells = int(point["cells"])
    if not (a > 0 and c > 0 and sigma != 0 and cells >= 4):
        raise ValueError("Invalid fixed grid point")
    kappa_a_squared = 2 * F * F * c * a * a / (EPSILON_W * R * T)
    wall_slope = F * a * sigma / (EPSILON_W * R * T)
    faces, centers, volumes = geometry(cells)
    initial = math.asinh(2 * wall_slope / kappa_a_squared)
    psi = [initial] * cells
    iterations = 0
    for iterations in range(1, 61):
        residual, lower, diagonal, upper = system(
            psi, faces, centers, volumes, kappa_a_squared, wall_slope
        )
        size = max(abs(v) for v in residual)
        if size < 1e-8 * max(1.0, abs(wall_slope)):
            break
        step = tridiagonal(lower, diagonal, upper, [-v for v in residual])
        accepted = False
        for reduction in range(25):
            scale = 2.0 ** -reduction
            proposal = [x + scale * dx for x, dx in zip(psi, step)]
            try:
                trial = system(proposal, faces, centers, volumes, kappa_a_squared, wall_slope)[0]
            except OverflowError:
                continue
            if max(abs(v) for v in trial) < size:
                psi = proposal
                accepted = True
                break
        if not accepted:
            raise ArithmeticError(f"Newton line search stalled at residual {size:.3g}")
    else:
        raise ArithmeticError("PB did not converge in 60 Newton steps")
    residual = system(psi, faces, centers, volumes, kappa_a_squared, wall_slope)[0]
    max_residual = max(abs(v) for v in residual)
    k_plus = 2 * sum(math.exp(-p) * v for p, v in zip(psi, volumes))
    k_minus = 2 * sum(math.exp(p) * v for p, v in zip(psi, volumes))
    gauss = F * c * (k_plus - k_minus) + 2 * sigma / a
    gauss_relative = abs(gauss) / abs(2 * sigma / a)
    d_s_local = 2 * 1.33e-9 * 2.03e-9 * k_plus * k_minus / (
        1.33e-9 * k_plus + 2.03e-9 * k_minus
    )
    if not (max_residual < 1e-8 * max(1.0, abs(wall_slope))
            and gauss_relative < 1e-4 and max(psi) <= 1e-12
            and psi[-1] < 0
            and k_plus > 1 and k_minus < 1):
        raise ArithmeticError(
            f"PB check failed: residual={max_residual:.3g}, gauss={gauss_relative:.3g}, "
            f"psi_min={min(psi):.3g}, psi_max={max(psi):.3g}, "
            f"k_plus={k_plus:.9g}, k_minus={k_minus:.9g}"
        )
    return {
        "id": point["id"], "a_nm": point["a_nm"], "sigma_C_m2": sigma,
        "c_mM": point["c_mM"], "c_mol_m3": c, "cells": cells,
        "psi_axis": psi[0], "psi_wall_cell": psi[-1],
        "k_plus": k_plus, "k_minus": k_minus, "D_s_local_m2_s": d_s_local,
        "max_fv_residual": max_residual, "gauss_relative": gauss_relative,
        "newton_iterations": iterations, "status": "MODEL_ONLY_NUMERIC_POINT",
        "note": ("Authorized feed endpoint supplement" if "_feed" in str(point["id"])
                 else "Fixed PB grid point; no material validation"),
    }


def main() -> None:
    with MATRIX.open(encoding="utf-8-sig") as stream:
        point = json.load(stream)[0]
    if not OUT.parent.is_dir():
        raise FileNotFoundError(OUT.parent)
    prior = []
    if LEDGER.exists():
        with LEDGER.open(encoding="utf-8-sig", newline="") as stream:
            prior = list(csv.DictReader(stream))
    if any(row["point_id"] == point["id"] and row["status"] == "PASS" for row in prior):
        print(f"{point['id']} already passed; no new PB solve")
        return
    if len(prior) >= 444:
        raise RuntimeError("PB solve cap reached")
    entry = {"call_number": len(prior) + 1, "utc_time": datetime.now(timezone.utc).isoformat(),
             "point_id": point["id"], "cells": point["cells"], "status": "STARTED",
             "command": "python src/mechanism_radial_pb.py", "exit_code": "", "reason": ""}
    prior.append(entry)
    with LEDGER.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(entry))
        writer.writeheader()
        writer.writerows(prior)
    try:
        result = solve(point)
    except Exception as error:
        prior[-1]["status"] = "FAILED"
        prior[-1]["exit_code"] = "1"
        prior[-1]["reason"] = str(error)
        raise
    else:
        prior[-1]["status"] = "PASS"
        prior[-1]["exit_code"] = "0"
        with OUT.open("w", encoding="utf-8") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
        print(json.dumps(result, ensure_ascii=False))
    finally:
        with LEDGER.open("w", encoding="utf-8-sig", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(entry))
            writer.writeheader()
            writer.writerows(prior)


if __name__ == "__main__":
    main()
