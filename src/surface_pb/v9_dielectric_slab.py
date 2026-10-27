"""Symmetric and asymmetric planar water/PE/water dielectric reference.

The PE interior solves Laplace exactly: its potential is linear, so the
displacement coupling between the two charge-regulating faces is explicit.
This is the frozen dense-slab check, not a 2D fiber dielectric solver.
"""

from dataclasses import dataclass

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

try:
    from .v9_analytic import state
    from .v9_pb_1d import Chemical, EPS0, F, R, T, Z, reservoir
except ImportError:
    from v9_analytic import state
    from v9_pb_1d import Chemical, EPS0, F, R, T, Z, reservoir


@dataclass
class SlabResult:
    h_nm: float
    sigma_left_C_m2: float
    sigma_right_C_m2: float
    psi_left_mV: float
    psi_right_mV: float
    E_PE_V_m: float
    displacement_PE_C_m2: float
    total_Gauss_absolute_C_m2: float
    newton_iterations: int
    scaled_residual: float


def solve_slab(h_nm, left: Chemical, right: Chemical,
               salt_mM=1.0, pH=7.0, L_nm=120.0, n=512,
               epsilon_water_r=78.3, epsilon_PE_r=2.3,
               tol=1e-11, max_iter=80):
    if min(h_nm, L_nm, epsilon_water_r, epsilon_PE_r) <= 0 or n < 8:
        raise ValueError("positive physical slab and water domains required")
    eps_w = EPS0 * epsilon_water_r
    eps_pe = EPS0 * epsilon_PE_r
    vt = R * T / F
    h = h_nm * 1e-9
    dx = L_nm * 1e-9 / n
    bath = reservoir(salt_mM, pH)
    poisson = F * dx * dx / (eps_w * vt)
    surface_factor = dx / (eps_w * vt)
    beta = surface_factor / 2
    solid_coupling = eps_pe * vt / h
    size = 2 * n + 2
    s_left, s_right = 2 * n, 2 * n + 1

    def evaluate(y, jacobian=False):
        us = (float(y[s_left]), float(y[s_right]))
        sig = []
        ds = []
        for surface_u, chem in zip(us, (left, right)):
            charge = state(surface_u, pH, chem.gamma_nm2, chem.pI, chem.b)
            sig.append(charge[3])
            ds.append(charge[4])
        displacement = solid_coupling * (us[1] - us[0])
        effective = (sig[0] + displacement, sig[1] - displacement)
        eff_deriv = ((ds[0] - solid_coupling, solid_coupling),
                     (solid_coupling, ds[1] - solid_coupling))
        residual = np.zeros(size)
        rows, cols, values = [], [], []
        concentrations = []
        for side in range(2):
            offset = side * n
            u = y[offset:offset + n]
            exponent = -Z[:, None] * u[None, :]
            if not np.all(np.isfinite(exponent)) or np.max(exponent) > 700:
                raise FloatingPointError("dielectric PB exponent outside range")
            conc = bath[:, None] * np.exp(exponent)
            concentrations.append(conc)
            rho = Z @ conc
            screen = (Z * Z) @ conc
            residual[offset] = u[0] - u[1] - poisson * rho[0] - surface_factor * effective[side]
            residual[offset + 1:offset + n - 1] = (2 * u[1:-1] - u[:-2] - u[2:]
                                                    - poisson * rho[1:-1])
            residual[offset + n - 1] = 3 * u[-1] - u[-2] - poisson * rho[-1]
            if jacobian:
                for i in range(n):
                    diagonal = (1 if i == 0 else 3 if i == n - 1 else 2) + poisson * screen[i]
                    rows.append(offset + i); cols.append(offset + i); values.append(diagonal)
                    if i:
                        rows.append(offset + i); cols.append(offset + i - 1); values.append(-1)
                    if i < n - 1:
                        rows.append(offset + i); cols.append(offset + i + 1); values.append(-1)
                rows.extend((offset, offset))
                cols.extend((s_left, s_right))
                values.extend((-surface_factor * eff_deriv[side][0],
                               -surface_factor * eff_deriv[side][1]))
        residual[s_left] = us[0] - y[0] - beta * effective[0]
        residual[s_right] = us[1] - y[n] - beta * effective[1]
        if jacobian:
            for side, row in enumerate((s_left, s_right)):
                rows.extend((row, row, row))
                cols.extend((side * n, s_left, s_right))
                values.extend((-1, (1 if side == 0 else 0) - beta * eff_deriv[side][0],
                               (1 if side == 1 else 0) - beta * eff_deriv[side][1]))
        if not jacobian:
            return residual
        matrix = coo_matrix((values, (rows, cols)), shape=(size, size)).tocsr()
        return residual, matrix, sig, displacement, concentrations

    y = np.zeros(size)
    for iteration in range(max_iter + 1):
        residual, jac, sigma, displacement, concentrations = evaluate(y, True)
        norm = float(np.max(np.abs(residual)))
        if norm <= tol:
            break
        if iteration == max_iter:
            raise RuntimeError(f"dielectric slab Newton failed {norm:g}")
        step = spsolve(jac, -residual)
        if not np.all(np.isfinite(step)):
            raise FloatingPointError("nonfinite dielectric Newton step")
        maximum = float(np.max(np.abs(step)))
        if maximum > 2:
            step *= 2 / maximum
        alpha = 1.0
        for _ in range(30):
            trial = y + alpha * step
            if float(np.max(np.abs(evaluate(trial)))) < norm:
                y = trial
                break
            alpha *= 0.5
        else:
            raise RuntimeError("dielectric slab line search stalled")
    residual, jac, sigma, displacement, concentrations = evaluate(y, True)
    # Sum both aqueous charge balances; the PE displacement cancels exactly.
    mobile = F * dx * sum(float(np.sum(Z @ c)) for c in concentrations)
    outer = eps_w * vt * 2 / dx * (y[n - 1] + y[2 * n - 1])
    gauss = abs(mobile + sigma[0] + sigma[1] - outer)
    return SlabResult(h_nm, float(sigma[0]), float(sigma[1]),
                      float(y[s_left] * vt * 1000), float(y[s_right] * vt * 1000),
                      float(-vt * (y[s_right] - y[s_left]) / h),
                      float(displacement), gauss, iteration,
                      float(np.max(np.abs(residual))))
