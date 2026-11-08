"""Linear pressure flow and short-circuit mobile-current observation for v9."""

from dataclasses import dataclass

import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
from scipy.linalg import solve_banded

try:
    from .v9_pb_1d import EPS0, F, Z
except ImportError:
    from v9_pb_1d import EPS0, F, Z


@dataclass
class Observation:
    H_um: float
    I_over_Q_C_m3: float
    zeta_app_mV: float
    I_half_A: float
    Q_half_m3_s: float
    Q_near_m3_s: float
    Q_outer_m3_s: float
    velocity_m_s: np.ndarray
    v_top_mean_m_s: float
    v_top_relative_spread: float
    flow_residual_scaled: float


def solve_velocity(geometry, field, H_um, G_Pa_m=1.0, mu_Pa_s=0.001):
    if geometry.kind != "backed":
        raise ValueError("the declared two-wall observation uses backed geometry")
    if H_um <= 0 or G_Pa_m <= 0 or mu_Pa_s <= 0:
        raise ValueError("positive gap, pressure gradient and viscosity required")
    mesh = field.mesh
    H_m = H_um * 1e-6
    L_m = mesh.z_edges_nm[-1] * 1e-9
    if L_m <= 0 or L_m >= H_m / 2:
        raise ValueError("near-field outer boundary must lie inside half-gap")
    nz, nx = mesh.water.shape
    dx = np.diff(mesh.x_edges_nm) * 1e-9
    dz = np.diff(mesh.z_edges_nm) * 1e-9
    count = len(field.u)
    diagonal = np.zeros(count)
    np.add.at(diagonal, mesh.outer_cell, -mesh.outer_conductance)
    eps_vt = EPS0 * 78.3 * (8.314462618 * 298.15 / F)
    wall_cond = mesh.wall_length_m / (mesh.wall_beta_m2_C * eps_vt)
    np.add.at(diagonal, mesh.wall_cell, wall_cond)
    # The support floor is no-slip in the water gaps between the bottom fibers.
    for xi in range(nx):
        cell = int(mesh.ids[0, xi])
        if cell >= 0:
            diagonal[cell] += dx[xi] / (dz[0] / 2)
    matrix = mesh.laplacian + diags(diagonal)
    source = (G_Pa_m / mu_Pa_s) * mesh.volumes_m2.copy()
    top_gradient = G_Pa_m * (H_m / 2 - L_m) / mu_Pa_s
    for xi in range(nx):
        cell = int(mesh.ids[-1, xi])
        if cell < 0:
            raise AssertionError("top clean-water boundary intersects solid")
        source[cell] += top_gradient * dx[xi]
    velocity = spsolve(matrix, source)
    if not np.all(np.isfinite(velocity)) or np.min(velocity) < -1e-16:
        raise RuntimeError("nonphysical pressure velocity")
    residual = matrix @ velocity - source
    v_at_L = np.array([velocity[int(mesh.ids[-1, xi])] + top_gradient * dz[-1] / 2
                       for xi in range(nx)])
    v_mean = float(np.average(v_at_L, weights=dx))
    spread = float(np.ptp(v_at_L) / max(abs(v_mean), 1e-30))
    Q_near = float(velocity @ mesh.volumes_m2)
    a = H_m / 2 - L_m
    width = geometry.period_nm * 1e-9
    Q_outer = width * (v_mean * a + G_Pa_m / mu_Pa_s * a**3 / 3)
    Q_half = Q_near + Q_outer
    rho_mobile = F * (Z @ field.concentrations_mol_m3)
    I_half = float((rho_mobile * velocity) @ mesh.volumes_m2)
    ratio = I_half / Q_half
    zeta = -(H_m**2 / (12 * EPS0 * 78.3)) * ratio * 1000
    return Observation(H_um, ratio, zeta, I_half, Q_half, Q_near, Q_outer,
                       velocity, v_mean, spread,
                       float(np.max(np.abs(residual)) / max(float(np.max(np.abs(source))), 1e-30)))


def flat_reference_zeta(pb_result, H_um, G_Pa_m=1.0, mu_Pa_s=0.001):
    """Independent planar PB + exact Poiseuille reference for HS checking."""
    if pb_result.geometry != "planar":
        raise ValueError("planar PB required")
    H_m = H_um * 1e-6
    centers = pb_result.center_nm * 1e-9
    if H_m / 2 <= centers[-1]:
        raise ValueError("half-gap must exceed PB reference domain")
    dx = centers[1] - centers[0]
    v = G_Pa_m / (2 * mu_Pa_s) * centers * (H_m - centers)
    rho = F * (Z @ pb_result.concentrations_mol_m3)
    I = float(np.sum(rho * v) * dx)
    Q = G_Pa_m * H_m**3 / (24 * mu_Pa_s)
    return -H_m**2 / (12 * EPS0 * 78.3) * I / Q * 1000


def empty_slit_half_reference(H_um=100.0, L_nm=380.0, n=256,
                              G_Pa_m=1.0, mu_Pa_s=0.001):
    """No-fiber no-slip half-slit FVM versus Poiseuille at cell centers."""
    H = H_um * 1e-6
    L = L_nm * 1e-9
    if not 0 < L < H / 2:
        raise ValueError("L inside half-gap required")
    dz = L / n
    diagonal = np.full(n, 2 / dz)
    diagonal[0] += 1 / dz  # wall distance is half a cell
    diagonal[-1] -= 1 / dz  # top has prescribed derivative, not Dirichlet
    off = np.full(n - 1, -1 / dz)
    rhs = np.full(n, G_Pa_m / mu_Pa_s * dz)
    rhs[-1] += G_Pa_m / mu_Pa_s * (H / 2 - L)
    ab = np.zeros((3, n))
    ab[0, 1:] = off
    ab[1] = diagonal
    ab[2, :-1] = off
    v = solve_banded((1, 1), ab, rhs)
    z = (np.arange(n) + 0.5) * dz
    exact = G_Pa_m / (2 * mu_Pa_s) * z * (H - z)
    relative = float(np.max(np.abs(v - exact)) / np.max(exact))
    return relative, float(np.min(v))
