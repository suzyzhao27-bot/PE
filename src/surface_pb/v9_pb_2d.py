"""Conservative four-species PB on a finite periodic square-fiber cell.

The mesh aligns each fiber face exactly. This equilibrium solver has no
hydrodynamic observation or calibration logic.
"""

import math
from dataclasses import dataclass

import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import spsolve

try:
    from .v9_geometry import Geometry
    from .v9_pb_1d import Chemical, E, EPS0, F, R, T, Z, reservoir
except ImportError:
    from v9_geometry import Geometry
    from v9_pb_1d import Chemical, E, EPS0, F, R, T, Z, reservoir


@dataclass
class Mesh:
    x_edges_nm: np.ndarray
    z_edges_nm: np.ndarray
    water: np.ndarray
    ids: np.ndarray
    layer_grid: np.ndarray
    volumes_m2: np.ndarray
    laplacian: object
    wall_cell: np.ndarray
    wall_layer: np.ndarray
    wall_length_m: np.ndarray
    wall_beta_m2_C: np.ndarray
    outer_cell: np.ndarray
    outer_conductance: np.ndarray
    target_nm: float
    L_nm: float


@dataclass
class Result:
    sigma_mean_C_m2: float
    sigma_outer_C_m2: float
    sigma_core_C_m2: float | None
    Q_proj_C_m2: float
    layer_sigma_C_m2: dict
    Gauss_absolute_C_m: float
    Gauss_relative: float | None
    max_species_M: float
    min_species_M: float
    physical_PB_warning: bool
    scaled_residual: float
    Newton_iterations: int
    mesh: Mesh
    u: np.ndarray
    concentrations_mol_m3: np.ndarray
    wall_sigma_C_m2: np.ndarray
    wall_surface_u: np.ndarray
    mobile_charge_C_m: float
    bound_charge_C_m: float
    outer_displacement_C_m: float


def _split_segments(boundaries, step_selector):
    edges = [boundaries[0]]
    for a, b in zip(boundaries[:-1], boundaries[1:]):
        count = max(1, math.ceil((b - a) / step_selector((a + b) / 2)))
        edges.extend(float(a + (b - a) * j / count) for j in range(1, count + 1))
    return np.asarray(edges)


def _unique(values):
    return sorted({round(float(v), 11) for v in values})


def build_mesh(geometry: Geometry, bath, eps_r=78.3, step_nm=None, outer_factor=1.0):
    ionic = float((Z * Z) @ bath)
    epsilon = EPS0 * eps_r
    lambda_nm = 1e9 * math.sqrt(epsilon * R * T / (F * F * ionic))
    width = geometry.period_nm
    if outer_factor < 1:
        raise ValueError("outer domain factor must be at least one")
    L_nm = outer_factor * max(5 * width, 12 * lambda_nm)
    target = min(geometry.d_nm / 12, geometry.g_nm / 20, lambda_nm / 8)
    if step_nm is not None:
        if step_nm <= 0:
            raise ValueError("positive requested grid step required")
        target = step_nm
    x_boundaries = [0.0, width]
    for layer in geometry.layers:
        x_boundaries.extend((layer.x_start_periodic_nm,
                             (layer.x_start_periodic_nm + geometry.d_nm) % width))
    x_edges_nm = _split_segments(_unique(x_boundaries), lambda _: target)
    film_bottom = min(layer.z_bottom_nm for layer in geometry.layers)
    film_top = max(layer.z_top_nm for layer in geometry.layers)
    z_low = film_bottom - L_nm if geometry.kind == "free" else film_bottom
    z_high = film_top + L_nm
    z_boundaries = [z_low, z_high, film_bottom, film_top]
    for layer in geometry.layers:
        z_boundaries.extend((layer.z_bottom_nm, layer.z_top_nm))
    if geometry.kind == "free":
        z_boundaries.append(film_bottom - min(12 * lambda_nm, L_nm))
    z_boundaries.append(film_top + min(12 * lambda_nm, L_nm))
    def z_step(mid):
        far = mid < film_bottom - 12 * lambda_nm or mid > film_top + 12 * lambda_nm
        return 4 * target if far else target
    z_edges_nm = _split_segments(_unique(z_boundaries), z_step)
    xc = (x_edges_nm[:-1] + x_edges_nm[1:]) / 2
    zc = (z_edges_nm[:-1] + z_edges_nm[1:]) / 2
    nz, nx = len(zc), len(xc)
    layer_grid = np.zeros((nz, nx), dtype=np.int32)
    for layer in geometry.layers:
        in_z = (zc > layer.z_bottom_nm) & (zc < layer.z_top_nm)
        in_x = ((xc - layer.x_start_periodic_nm) % width) < geometry.d_nm
        layer_grid[np.ix_(in_z, in_x)] = layer.layer_id
    water = layer_grid == 0
    ids = np.full((nz, nx), -1, dtype=np.int32)
    ids[water] = np.arange(int(water.sum()), dtype=np.int32)
    count = int(water.sum())
    dx = np.diff(x_edges_nm) * 1e-9
    dz = np.diff(z_edges_nm) * 1e-9
    volume = (dz[:, None] * dx[None, :])[water]
    rows, columns, data = [], [], []
    diagonal = np.zeros(count)
    wall_cell, wall_layer, wall_length, wall_beta = [], [], [], []
    outer_cell, outer_conduct = [], []
    vt = R * T / F
    for zi, xi in zip(*np.nonzero(water)):
        cell = int(ids[zi, xi])
        for direction, xn, zn in (("right", (xi + 1) % nx, zi),
                                  ("up", xi, zi + 1)):
            if zn >= nz:
                g = dx[xi] / (dz[zi] / 2)
                diagonal[cell] += g
                outer_cell.append(cell)
                outer_conduct.append(g)
                continue
            length = dz[zi] if direction == "right" else dx[xi]
            if water[zn, xn]:
                distance = ((dx[xi] + dx[xn]) / 2 if direction == "right"
                            else (dz[zi] + dz[zn]) / 2)
                conduct = length / distance
                other = int(ids[zn, xn])
                diagonal[cell] += conduct
                diagonal[other] += conduct
                rows.extend((cell, other))
                columns.extend((other, cell))
                data.extend((-conduct, -conduct))
            else:
                wall_cell.append(cell)
                wall_layer.append(int(layer_grid[zn, xn]))
                wall_length.append(length)
                wall_beta.append((dx[xi] if direction == "right" else dz[zi]) / (2 * epsilon * vt))
        for direction, xn, zn in (("left", (xi - 1) % nx, zi),
                                  ("down", xi, zi - 1)):
            if zn < 0:
                if geometry.kind == "free":
                    g = dx[xi] / (dz[zi] / 2)
                    diagonal[cell] += g
                    outer_cell.append(cell)
                    outer_conduct.append(g)
                continue  # backed bottom is neutral support
            if not water[zn, xn]:
                length = dz[zi] if direction == "left" else dx[xi]
                wall_cell.append(cell)
                wall_layer.append(int(layer_grid[zn, xn]))
                wall_length.append(length)
                wall_beta.append((dx[xi] if direction == "left" else dz[zi]) / (2 * epsilon * vt))
    rows.extend(range(count))
    columns.extend(range(count))
    data.extend(diagonal.tolist())
    laplacian = coo_matrix((data, (rows, columns)), shape=(count, count)).tocsr()
    mesh = Mesh(x_edges_nm, z_edges_nm, water, ids, layer_grid, volume,
                laplacian, np.asarray(wall_cell, dtype=np.int32),
                np.asarray(wall_layer, dtype=np.int32), np.asarray(wall_length),
                np.asarray(wall_beta), np.asarray(outer_cell, dtype=np.int32),
                np.asarray(outer_conduct), target, L_nm)
    # Each exact face must be counted once, including periodically wrapped faces.
    expected = geometry.wet_area_m2_per_m
    actual = float(np.sum(mesh.wall_length_m))
    if not math.isclose(actual, expected, rel_tol=1e-10, abs_tol=1e-17):
        raise AssertionError(f"wet perimeter mismatch {actual} versus {expected}")
    return mesh


def _wall_regulation(u_cell, beta, pH, chemical):
    cap = E * chemical.gamma_nm2 * 1e18
    if cap == 0:
        zeros = np.zeros_like(u_cell)
        return zeros, zeros, u_cell
    sigma = np.zeros_like(u_cell)
    l10 = math.log(10)
    for _ in range(35):
        us = u_cell + beta * sigma
        a = l10 * (chemical.b + chemical.pI - pH) - us
        b = l10 * (chemical.b - chemical.pI + pH) + us
        anchor = np.maximum(0, np.maximum(a, b))
        w0, wp, wm = np.exp(-anchor), np.exp(a - anchor), np.exp(b - anchor)
        den = w0 + wp + wm
        mean = (wp - wm) / den
        var = (wp + wm) / den - mean * mean
        derivative = -cap * var
        correction = (sigma - cap * mean) / (1 - beta * derivative)
        sigma -= correction
        if float(np.max(np.abs(correction))) < 1e-14:
            break
    else:
        raise RuntimeError("local surface-regulation root did not converge")
    us = u_cell + beta * sigma
    return sigma, derivative / (1 - beta * derivative), us


def solve_equilibrium(geometry: Geometry, salt_mM, pH, chemical=None,
                      fixed_sigma=None, mesh_step_nm=None, outer_factor=1.0, eps_r=78.3,
                      tolerance=1e-9, max_iter=60):
    if (chemical is None) == (fixed_sigma is None):
        raise ValueError("select chemical regulation or fixed surface charge")
    bath = reservoir(salt_mM, pH)
    mesh = build_mesh(geometry, bath, eps_r=eps_r, step_nm=mesh_step_nm,
                      outer_factor=outer_factor)
    epsilon = EPS0 * eps_r
    vt = R * T / F
    coeff = F / (epsilon * vt)
    face_coefficient = mesh.wall_length_m / (epsilon * vt)
    count = len(mesh.volumes_m2)

    def evaluate(u, jacobian=False):
        exponent = -Z[:, None] * u[None, :]
        if not np.all(np.isfinite(exponent)) or np.max(exponent) > 700:
            raise FloatingPointError("PB exponent outside representable range")
        conc = bath[:, None] * np.exp(exponent)
        rho = Z @ conc
        screen = (Z * Z) @ conc
        if chemical is None:
            sigma = np.full(len(mesh.wall_cell), fixed_sigma)
            derivative = np.zeros_like(sigma)
            surface_u = u[mesh.wall_cell] + mesh.wall_beta_m2_C * sigma
        else:
            sigma, derivative, surface_u = _wall_regulation(
                u[mesh.wall_cell], mesh.wall_beta_m2_C, pH, chemical)
        wall = np.bincount(mesh.wall_cell, weights=face_coefficient * sigma, minlength=count)
        residual = mesh.laplacian @ u - coeff * mesh.volumes_m2 * rho - wall
        if not jacobian:
            return residual
        wall_derivative = np.bincount(mesh.wall_cell,
                                      weights=face_coefficient * derivative, minlength=count)
        jac = mesh.laplacian + diags(coeff * mesh.volumes_m2 * screen - wall_derivative)
        scale = max(1.0, float(np.max(np.abs(wall))),
                    float(np.max(coeff * mesh.volumes_m2 * np.sum(conc, axis=0))))
        return residual, jac, conc, sigma, surface_u, scale

    u = np.zeros(count)
    for iteration in range(max_iter + 1):
        residual, jac, conc, sigma, surface_u, scale = evaluate(u, True)
        norm = float(np.max(np.abs(residual)) / scale)
        if norm <= tolerance:
            break
        if iteration == max_iter:
            raise RuntimeError(f"finite-fiber PB Newton failed: {norm:g}")
        step = spsolve(jac, -residual)
        if not np.all(np.isfinite(step)):
            raise FloatingPointError("nonfinite PB step")
        step_max = float(np.max(np.abs(step)))
        if step_max > 2:
            step *= 2 / step_max
        old = float(np.max(np.abs(residual)))
        alpha = 1.0
        for _ in range(30):
            trial = u + alpha * step
            if float(np.max(np.abs(evaluate(trial)))) < old:
                u = trial
                break
            alpha *= 0.5
        else:
            raise RuntimeError(f"finite-fiber PB line search stalled: {norm:g}")
    residual, jac, conc, sigma, surface_u, scale = evaluate(u, True)
    if not np.all(np.isfinite(sigma)) or not np.all(np.isfinite(conc)):
        raise FloatingPointError("nonfinite finite-fiber PB state")
    if chemical is not None:
        capacity = E * chemical.gamma_nm2 * 1e18
        if np.max(np.abs(sigma)) > capacity * (1 + 1e-12):
            raise AssertionError("surface site capacity exceeded")
    mobile = F * float((Z @ conc) @ mesh.volumes_m2)
    bound = float(np.sum(sigma * mesh.wall_length_m))
    outer = epsilon * vt * float(np.sum(mesh.outer_conductance * u[mesh.outer_cell]))
    imbalance = mobile + bound - outer
    cap = E * chemical.gamma_nm2 * 1e18 if chemical is not None else abs(fixed_sigma)
    denominator = max(float(np.sum(np.abs(sigma) * mesh.wall_length_m)),
                      cap * geometry.wet_area_m2_per_m * 1e-6)
    by_layer = {}
    for layer in geometry.layers:
        mask = mesh.wall_layer == layer.layer_id
        area = float(np.sum(mesh.wall_length_m[mask]))
        assert math.isclose(area, layer.wet_perimeter_nm * 1e-9, rel_tol=1e-10)
        by_layer[layer.layer_id] = float(np.sum(sigma[mask] * mesh.wall_length_m[mask]) / area)
    outer_mask = np.isin(mesh.wall_layer, [layer.layer_id for layer in geometry.layers if layer.is_outer])
    outer_area = float(np.sum(mesh.wall_length_m[outer_mask]))
    sigma_outer = float(np.sum(sigma[outer_mask] * mesh.wall_length_m[outer_mask]) / outer_area)
    core_ids = [layer.layer_id for layer in geometry.layers if layer.is_core]
    if core_ids:
        core_mask = np.isin(mesh.wall_layer, core_ids)
        core_area = float(np.sum(mesh.wall_length_m[core_mask]))
        sigma_core = float(np.sum(sigma[core_mask] * mesh.wall_length_m[core_mask]) / core_area)
    else:
        sigma_core = None
    max_m = float(np.max(conc)) / 1000
    return Result(bound / geometry.wet_area_m2_per_m, sigma_outer, sigma_core,
                  bound / geometry.projected_area_m2_per_m, by_layer,
                  abs(imbalance), abs(imbalance) / denominator if denominator else None,
                  max_m, float(np.min(conc)) / 1000, max_m > 0.5,
                  float(np.max(np.abs(residual)) / scale), iteration,
                  mesh, u, conc, sigma, surface_u, mobile, bound, outer)
