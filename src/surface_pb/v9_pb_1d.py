"""One-dimensional four-species PB reference for v9 technical checks.

Planar and radial cells share conservative faces. This is a local reference,
not the finite-layer membrane model or a calibration target.
"""

import math
from dataclasses import dataclass

import numpy as np
from scipy.linalg import solve_banded
from scipy.optimize import brentq

try:
    from .v9_analytic import state
except ImportError:
    from v9_analytic import state

F = 96485.33212
R = 8.314462618
T = 298.15
EPS0 = 8.8541878128e-12
E = 1.602176634e-19
Z = np.array([1.0, -1.0, 1.0, -1.0])


def reservoir(salt_mM, pH):
    """Na, Cl, H, OH in mol/m3 for the declared NaCl+acid/base scenario."""
    if salt_mM <= 0:
        raise ValueError("positive background NaCl required")
    h = 1000 * 10 ** (-pH)
    oh = 1000 * 10 ** (pH - 14)
    return np.array([salt_mM + max(oh - h, 0),
                     salt_mM + max(h - oh, 0), h, oh])


@dataclass(frozen=True)
class Chemical:
    gamma_nm2: float
    pI: float
    b: float


@dataclass
class Result:
    geometry: str
    sigma_C_m2: float
    surface_psi_mV: float
    surface_u: float
    max_species_M: float
    min_species_M: float
    gauss_absolute: float
    gauss_relative: float | None
    newton_iterations: int
    scaled_residual: float
    center_nm: np.ndarray
    u: np.ndarray
    concentrations_mol_m3: np.ndarray
    outer_electric_flux: float
    mobile_charge: float
    bound_charge: float


def solve(kind="planar", inner_nm=0.0, outer_nm=120.0, n=512,
          salt_mM=1.0, pH=7.0, chemical=None, fixed_sigma=None,
          outer_condition="bulk", eps_r=78.3, tolerance=1e-11,
          max_iter=100):
    if kind not in ("planar", "radial") or outer_condition not in ("bulk", "cell"):
        raise ValueError("unsupported one-dimensional reference geometry or outer boundary")
    if outer_nm <= inner_nm or (kind == "radial" and inner_nm <= 0) or n < 8:
        raise ValueError("invalid domain")
    if (chemical is None) == (fixed_sigma is None):
        raise ValueError("select chemical regulation or fixed sigma, not both")
    if chemical is not None and chemical.gamma_nm2 < 0:
        raise ValueError("nonnegative Gamma required")
    edges_nm = np.linspace(inner_nm, outer_nm, n + 1)
    centers_nm = (edges_nm[:-1] + edges_nm[1:]) / 2
    edges = edges_nm * 1e-9
    centers = centers_nm * 1e-9
    dx = edges[1] - edges[0]
    if kind == "planar":
        volume = np.diff(edges)
        interior_face = np.ones(n - 1)
        inner_face = outer_face = 1.0
    else:
        volume = np.diff(edges * edges) / 2
        interior_face = edges[1:-1]
        inner_face, outer_face = edges[0], edges[-1]
    conduct = interior_face / dx
    diagonal = np.r_[conduct, 0.0] + np.r_[0.0, conduct]
    boundary_conduct = 2 * outer_face / dx if outer_condition == "bulk" else 0.0
    diagonal[-1] += boundary_conduct
    epsilon = EPS0 * eps_r
    vt = R * T / F
    coeff = F / (epsilon * vt)
    bath = reservoir(salt_mM, pH)
    if abs(float(Z @ bath)) > 1e-11:
        raise AssertionError("reservoir is not electroneutral")
    capacity = E * chemical.gamma_nm2 * 1e18 if chemical is not None else None
    wall_beta = (dx / 2) / (epsilon * vt)

    def wall(u0):
        if chemical is None:
            return fixed_sigma, 0.0, u0 + wall_beta * fixed_sigma
        if capacity == 0:
            return 0.0, 0.0, u0
        def residual(sig):
            return sig - state(u0 + wall_beta * sig, pH, chemical.gamma_nm2,
                               chemical.pI, chemical.b)[3]
        sigma = brentq(residual, -capacity, capacity, xtol=1e-15)
        surface_u = u0 + wall_beta * sigma
        derivative = state(surface_u, pH, chemical.gamma_nm2,
                           chemical.pI, chemical.b)[4]
        return sigma, derivative / (1 - wall_beta * derivative), surface_u

    def evaluate(u, jacobian=False):
        exponent = -Z[:, None] * u[None, :]
        if not np.all(np.isfinite(exponent)) or np.max(exponent) > 700:
            raise FloatingPointError("PB exponent outside representable range")
        concentrations = bath[:, None] * np.exp(exponent)
        charge = Z @ concentrations
        screen = (Z * Z) @ concentrations
        sigma, derivative, surface_u = wall(float(u[0]))
        residual = diagonal * u
        residual[1:] -= conduct * u[:-1]
        residual[:-1] -= conduct * u[1:]
        residual -= coeff * volume * charge
        residual[0] -= inner_face / (epsilon * vt) * sigma
        if not jacobian:
            return residual
        jac_diagonal = diagonal + coeff * volume * screen
        jac_diagonal[0] -= inner_face / (epsilon * vt) * derivative
        scale = max(1.0, abs(inner_face * sigma / (epsilon * vt)),
                    float(np.max(coeff * volume * np.sum(concentrations, axis=0))))
        return residual, jac_diagonal, concentrations, sigma, surface_u, scale

    u = np.zeros(n)
    for iteration in range(max_iter + 1):
        residual, jd, conc, sigma, surface_u, scale = evaluate(u, True)
        norm = float(np.max(np.abs(residual)) / scale)
        if norm <= tolerance:
            break
        if iteration == max_iter:
            raise RuntimeError(f"PB Newton did not converge: {norm:g}")
        ab = np.zeros((3, n))
        ab[0, 1:] = -conduct
        ab[1] = jd
        ab[2, :-1] = -conduct
        step = solve_banded((1, 1), ab, -residual)
        if not np.all(np.isfinite(step)):
            raise FloatingPointError("nonfinite Newton step")
        step_max = float(np.max(np.abs(step)))
        if step_max > 2:
            step *= 2 / step_max
        old = float(np.max(np.abs(residual)))
        alpha = 1.0
        for _ in range(32):
            trial = u + alpha * step
            if float(np.max(np.abs(evaluate(trial)))) < old:
                u = trial
                break
            alpha *= 0.5
        else:
            raise RuntimeError(f"PB Newton line search stalled: {norm:g}")
    residual, jd, conc, sigma, surface_u, scale = evaluate(u, True)
    mobile = F * float((Z @ conc) @ volume)
    bound = inner_face * sigma
    outer_flux = epsilon * vt * boundary_conduct * float(u[-1])
    # The outward displacement is -epsilon*grad(psi); with bulk u=0,
    # its signed value is epsilon*Vt*conduct*u_last.
    imbalance = mobile + bound - outer_flux
    denominator = max(abs(bound),
                      (capacity * inner_face * 1e-6) if capacity is not None else 0.0)
    relative = abs(imbalance) / denominator if denominator else None
    return Result(kind, sigma, surface_u * vt * 1000, surface_u,
                  float(np.max(conc)) / 1000, float(np.min(conc)) / 1000,
                  abs(imbalance), relative, iteration,
                  float(np.max(np.abs(residual)) / scale), centers_nm, u, conc,
                  outer_flux, mobile, bound)
