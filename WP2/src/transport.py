"""v6 M0/M1/M2 conditional salt model; concentrations are scaled by C0.

Lengths and diffusivities are SI. Flux is per connected water area and is
reported in units of C0 m/s. This is a model calculation, not calibrated flux.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from math import exp, log, pi

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.sparse import csc_matrix, lil_matrix


D_NA = 1.33e-9
D_CL = 2.03e-9
D_B_REF = 2 * D_NA * D_CL / (D_NA + D_CL)


@dataclass(frozen=True)
class Config:
    h_m: float
    a_m: float = 25e-9
    ell_m: float = 893.953885e-9
    ws_ref: float = 7.039790538
    a_ref_m: float = 25e-9
    robin_r: float = 0.0
    dm_m2_s: float = D_B_REF
    db_m2_s: float = D_B_REF
    gamma_per_m: float = 1e4
    volume_ratio: float = 1.0
    chi: float = 0.0
    tau_s: float = 1.0
    q_conductivity: float = 1.0
    access: bool = True
    boundary_layer_m: float = 0.0
    fixed_pools: bool = False

    def __post_init__(self) -> None:
        positive = (self.h_m, self.a_m, self.ell_m, self.a_ref_m,
                    self.dm_m2_s, self.db_m2_s, self.gamma_per_m,
                    self.volume_ratio, self.tau_s, self.q_conductivity)
        if any(x <= 0 for x in positive) or self.ws_ref < 0 or self.chi < 0 or self.boundary_layer_m < 0:
            raise ValueError("positive scales, nonnegative field and capacity required")
        if 1 + self.robin_r * exp(-self.h_m / self.ell_m) <= 0:
            raise ValueError("invalid Robin denominator")

    def with_nm(self, h_nm: float) -> Config:
        return replace(self, h_m=h_nm * 1e-9)


def field(x_m: float | np.ndarray, cfg: Config) -> float | np.ndarray:
    x = np.asarray(x_m)
    ws = cfg.ws_ref * cfg.a_ref_m / cfg.a_m
    numerator = np.exp(-x / cfg.ell_m) + np.exp(-(cfg.h_m - x) / cfg.ell_m)
    result = ws * numerator / (1 + cfg.robin_r * exp(-cfg.h_m / cfg.ell_m))
    return float(result) if result.ndim == 0 else result


_GAUSS_X, _GAUSS_W = np.polynomial.legendre.leggauss(8)


def field_resistance(x0: float, x1: float, cfg: Config) -> float:
    """Integral exp(w)/D_m dx on a short interval, in s/m."""
    center = (x0 + x1) / 2
    half = (x1 - x0) / 2
    return float(half * np.dot(_GAUSS_W, np.exp(field(center + half * _GAUSS_X, cfg))) / cfg.dm_m2_s)


def steady_resistance(cfg: Config) -> float:
    access = (pi * cfg.a_m / (2 * cfg.db_m2_s) if cfg.access else 0.0)
    access += 2 * cfg.boundary_layer_m / cfg.db_m2_s
    # Adaptive integration resolves narrow surface layers when ell << h.
    points = sorted({x for x in (cfg.ell_m, cfg.h_m - cfg.ell_m)
                     if 0 < x < cfg.h_m})
    integral = quad(lambda x: exp(field(x, cfg)) / cfg.dm_m2_s,
                    0, cfg.h_m, points=points, epsrel=1e-10, limit=300)[0]
    return integral + access


def steady_permeance(cfg: Config) -> float:
    """Connected-area permeance p_h, m/s; whole-area value is epsilon*p_h."""
    return 1 / steady_resistance(cfg)


def m0_right_fraction(times_s: np.ndarray, cfg: Config) -> np.ndarray:
    t = np.asarray(times_s, dtype=float)
    if np.any(t < 0):
        raise ValueError("nonnegative times required")
    rate = cfg.gamma_per_m * steady_permeance(cfg) * (1 + 1 / cfg.volume_ratio)
    return cfg.volume_ratio / (1 + cfg.volume_ratio) * (-np.expm1(-rate * t))


@dataclass
class Operator:
    matrix: csc_matrix
    centers_m: np.ndarray
    field_centers: np.ndarray
    dx_m: float
    left_flux: tuple[float, float]
    right_flux: tuple[float, float]
    n_cells: int
    has_storage: bool
    config: Config

    def initial(self, left: float = 1.0, right: float = 0.0,
                mobile: np.ndarray | None = None,
                storage: np.ndarray | None = None) -> np.ndarray:
        n = self.n_cells
        y = np.zeros(2 + n + (n if self.has_storage else 0))
        y[0], y[n + 1] = left, right
        if mobile is not None:
            y[1:n + 1] = mobile
        if self.has_storage and storage is not None:
            y[n + 2:] = storage
        return y

    def fluxes(self, y: np.ndarray) -> tuple[float, float]:
        al, bl = self.left_flux
        ar, br = self.right_flux
        n = self.n_cells
        return al * y[0] - bl * y[1], ar * y[n] - br * y[n + 1]

    def mass(self, y: np.ndarray) -> float:
        cfg = self.config
        n = self.n_cells
        total = cfg.volume_ratio * y[0] + y[n + 1]
        total += cfg.gamma_per_m * self.dx_m * np.sum(y[1:n + 1])
        if self.has_storage:
            total += cfg.gamma_per_m * self.dx_m * np.sum(y[n + 2:])
        return float(total)

    def free_energy(self, y: np.ndarray) -> float:
        """Closed-system F/(RT V_R C0), with cref=C0 and x log x -> 0."""
        if self.config.fixed_pools:
            raise ValueError("open fixed-pool system needs external chemical work")
        if np.min(y) < 0:
            raise ValueError("free energy requires nonnegative state")

        def ideal(z: np.ndarray) -> np.ndarray:
            return np.where(z > 0, z * np.log(np.maximum(z, np.finfo(float).tiny)) - z, 0)

        n = self.n_cells
        cfg = self.config
        mobile = y[1:n + 1]
        f = cfg.volume_ratio * ideal(np.asarray(y[0])) + ideal(np.asarray(y[n + 1]))
        f += cfg.gamma_per_m * self.dx_m * np.sum(ideal(mobile) + mobile * self.field_centers)
        if self.has_storage:
            stored = y[n + 2:]
            f += cfg.gamma_per_m * self.dx_m * np.sum(
                ideal(stored) - stored * log(cfg.chi) + stored * self.field_centers
            )
        return float(f)


def build_operator(cfg: Config, n_cells: int = 100) -> Operator:
    if n_cells < 3:
        raise ValueError("at least three cells required")
    n = n_cells
    dx = cfg.h_m / n
    centers = (np.arange(n) + 0.5) * dx
    wc = np.asarray(field(centers, cfg))
    has_storage = cfg.chi > 0
    size = n + 2 + (n if has_storage else 0)
    mat = lil_matrix((size, size), dtype=float)

    def put_flux(left_idx: int, right_idx: int, a: float, b: float,
                 left_factor: float, right_factor: float) -> None:
        # j = a*y[left] - b*y[right]; divergence sends +j rightward.
        for row, factor in ((left_idx, left_factor), (right_idx, right_factor)):
            mat[row, left_idx] += factor * a
            mat[row, right_idx] -= factor * b

    entry_half = field_resistance(0, centers[0], cfg)
    exit_half = field_resistance(centers[-1], cfg.h_m, cfg)
    access_half = (pi * cfg.a_m / (4 * cfg.db_m2_s) if cfg.access else 0.0)
    access_half += cfg.boundary_layer_m / cfg.db_m2_s
    left_r = entry_half + access_half
    right_r = exit_half + access_half
    al, bl = 1 / left_r, exp(wc[0]) / left_r
    ar, br = exp(wc[-1]) / right_r, 1 / right_r
    put_flux(0, 1, al, bl,
             0 if cfg.fixed_pools else -cfg.gamma_per_m / cfg.volume_ratio,
             1 / dx)
    for i in range(n - 1):
        r = field_resistance(centers[i], centers[i + 1], cfg)
        a, b = exp(wc[i]) / r, exp(wc[i + 1]) / r
        put_flux(i + 1, i + 2, a, b, -1 / dx, 1 / dx)
    put_flux(n, n + 1, ar, br, -1 / dx,
             0 if cfg.fixed_pools else cfg.gamma_per_m)
    if has_storage:
        for i in range(n):
            c_idx, s_idx = 1 + i, n + 2 + i
            mat[c_idx, c_idx] -= cfg.chi / cfg.tau_s
            mat[c_idx, s_idx] += 1 / cfg.tau_s
            mat[s_idx, c_idx] += cfg.chi / cfg.tau_s
            mat[s_idx, s_idx] -= 1 / cfg.tau_s
    return Operator(mat.tocsc(), centers, wc, dx, (al, bl), (ar, br), n, has_storage, cfg)


@dataclass
class Trajectory:
    times_s: np.ndarray
    states: np.ndarray
    operator: Operator
    solver_nfev: int
    solver_status: str

    @property
    def right_fraction(self) -> np.ndarray:
        return self.states[self.operator.n_cells + 1]

    @property
    def conductivity_increment(self) -> np.ndarray:
        return self.operator.config.q_conductivity * self.right_fraction

    @property
    def mass_error(self) -> float:
        if self.operator.config.fixed_pools:
            raise ValueError("fixed reservoirs exchange salt with their environment")
        masses = np.array([self.operator.mass(self.states[:, i]) for i in range(self.states.shape[1])])
        return float(np.max(np.abs(masses - masses[0])) / max(abs(masses[0]), 1e-300))


def solve(cfg: Config, times_s: np.ndarray, n_cells: int = 100,
          *, left: float = 1.0, right: float = 0.0,
          mobile: np.ndarray | None = None, storage: np.ndarray | None = None,
          rtol: float = 1e-9, atol: float = 1e-12) -> Trajectory:
    times = np.asarray(times_s, dtype=float)
    if times.ndim != 1 or len(times) < 2 or times[0] != 0 or np.any(np.diff(times) <= 0):
        raise ValueError("times must start at zero and increase strictly")
    op = build_operator(cfg, n_cells)
    y0 = op.initial(left, right, mobile, storage)
    result = solve_ivp(lambda _t, y: op.matrix @ y, (0, times[-1]), y0,
                       method="BDF", jac=op.matrix, t_eval=times,
                       rtol=rtol, atol=atol)
    if not result.success:
        raise RuntimeError(result.message)
    return Trajectory(times, result.y, op, result.nfev, result.message)


def minutes_to_seconds(time_min: np.ndarray) -> np.ndarray:
    return np.asarray(time_min, dtype=float) * 60
