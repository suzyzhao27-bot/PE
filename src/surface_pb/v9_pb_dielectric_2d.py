"""Finite-water/PE dielectric PB on the same exact square-fiber mesh.

Every water/solid face carries one regulated surface charge. The face
potential and PE displacement are eliminated locally; the global system
retains the solid Laplace degrees of freedom and water PB degrees of freedom.
"""

import math
from dataclasses import dataclass

import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import spsolve

try:
    from .v9_pb_1d import Chemical, E, EPS0, F, R, T, Z, reservoir
    from .v9_pb_2d import build_mesh, _wall_regulation
except ImportError:
    from v9_pb_1d import Chemical, E, EPS0, F, R, T, Z, reservoir
    from v9_pb_2d import build_mesh, _wall_regulation


@dataclass
class DielectricResult:
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
    mesh: object
    u: np.ndarray
    concentrations_mol_m3: np.ndarray
    wall_sigma_C_m2: np.ndarray
    wall_surface_u: np.ndarray
    mobile_charge_C_m: float
    bound_charge_C_m: float
    outer_displacement_C_m: float
    solid_u: np.ndarray
    eps_PE_r: float


def solve_dielectric(geometry, salt_mM, pH, chemical: Chemical,
                     eps_PE_r=2.3, eps_water_r=78.3,
                     mesh_step_nm=None, outer_factor=1.0,
                     tolerance=1e-9, max_iter=60):
    if eps_PE_r <= 0 or eps_water_r <= 0 or chemical.gamma_nm2 < 0:
        raise ValueError("positive dielectric constants and site capacity required")
    bath = reservoir(salt_mM, pH)
    mesh = build_mesh(geometry, bath, eps_r=eps_water_r,
                      step_nm=mesh_step_nm, outer_factor=outer_factor)
    water = mesh.water
    nz, nx = water.shape
    count = nx * nz
    water_ids = np.flatnonzero(water.ravel())
    solid_ids = np.flatnonzero(~water.ravel())
    dx = np.diff(mesh.x_edges_nm) * 1e-9
    dz = np.diff(mesh.z_edges_nm) * 1e-9
    relative_PE = eps_PE_r / eps_water_r
    eps_vt = EPS0 * eps_water_r * R * T / F
    diagonal = np.zeros(count)
    rows, cols, values = [], [], []
    face_water, face_solid, face_layer = [], [], []
    face_length, face_alpha, face_beta, face_regulation_beta = [], [], [], []
    outer_ids, outer_conductance = [], []

    def connect(first, second, conductance):
        diagonal[first] += conductance
        diagonal[second] += conductance
        rows.extend((first, second))
        cols.extend((second, first))
        values.extend((-conductance, -conductance))

    def edge(za, xa, zb, xb, length, da, db):
        first, second = za * nx + xa, zb * nx + xb
        first_water, second_water = bool(water[za, xa]), bool(water[zb, xb])
        if first_water == second_water:
            relative = 1.0 if first_water else relative_PE
            connect(first, second, relative * length / ((da + db) / 2))
            return
        w, s = (first, second) if first_water else (second, first)
        dw, ds = (da, db) if first_water else (db, da)
        A, B = 2 / dw, 2 * relative_PE / ds
        factor = A + B
        connect(w, s, A * B / factor * length)
        face_water.append(w)
        face_solid.append(s)
        face_layer.append(int(mesh.layer_grid.ravel()[s]))
        face_length.append(length)
        face_alpha.append(A / factor)
        face_beta.append(B / factor)
        face_regulation_beta.append(1 / (eps_vt * factor))

    for z in range(nz):
        for x in range(nx):
            edge(z, x, z, (x + 1) % nx, dz[z], dx[x], dx[(x + 1) % nx])
            if z + 1 < nz:
                edge(z, x, z + 1, x, dx[x], dz[z], dz[z + 1])
            elif water[z, x]:
                cell = z * nx + x
                g = dx[x] / (dz[z] / 2)
                diagonal[cell] += g
                outer_ids.append(cell)
                outer_conductance.append(g)
            if z == 0 and geometry.kind == "free" and water[z, x]:
                cell = x
                g = dx[x] / (dz[z] / 2)
                diagonal[cell] += g
                outer_ids.append(cell)
                outer_conductance.append(g)
    rows.extend(range(count))
    cols.extend(range(count))
    values.extend(diagonal.tolist())
    matrix = coo_matrix((values, (rows, cols)), shape=(count, count)).tocsr()
    face_water = np.asarray(face_water, dtype=np.int32)
    face_solid = np.asarray(face_solid, dtype=np.int32)
    face_layer = np.asarray(face_layer, dtype=np.int32)
    face_length = np.asarray(face_length)
    face_alpha = np.asarray(face_alpha)
    face_beta = np.asarray(face_beta)
    face_regulation_beta = np.asarray(face_regulation_beta)
    outer_ids = np.asarray(outer_ids, dtype=np.int32)
    outer_conductance = np.asarray(outer_conductance)
    if not math.isclose(float(sum(face_length)), geometry.wet_area_m2_per_m,
                        rel_tol=1e-10, abs_tol=1e-17):
        raise AssertionError("dielectric interface perimeter mismatch")
    volume = np.zeros(count)
    volume[water_ids] = mesh.volumes_m2
    mobile_coefficient = F * volume[water_ids] / eps_vt
    face_coefficient = face_length / eps_vt

    def evaluate(u, jacobian=False):
        exponent = -Z[:, None] * u[water_ids][None, :]
        if not np.all(np.isfinite(exponent)) or float(np.max(exponent)) > 700:
            raise FloatingPointError("dielectric PB exponent outside representable range")
        concentrations = bath[:, None] * np.exp(exponent)
        charge = Z @ concentrations
        screen = (Z * Z) @ concentrations
        mix = face_alpha * u[face_water] + face_beta * u[face_solid]
        sigma, derivative, surface_u = _wall_regulation(
            mix, face_regulation_beta, pH, chemical)
        residual = matrix @ u
        residual[water_ids] -= mobile_coefficient * charge
        np.add.at(residual, face_water, -face_coefficient * face_alpha * sigma)
        np.add.at(residual, face_solid, -face_coefficient * face_beta * sigma)
        if not jacobian:
            return residual
        w_w = -face_coefficient * face_alpha * derivative * face_alpha
        w_s = -face_coefficient * face_alpha * derivative * face_beta
        s_w = -face_coefficient * face_beta * derivative * face_alpha
        s_s = -face_coefficient * face_beta * derivative * face_beta
        row = np.concatenate((face_water, face_water, face_solid, face_solid))
        col = np.concatenate((face_water, face_solid, face_water, face_solid))
        val = np.concatenate((w_w, w_s, s_w, s_s))
        jac = matrix + coo_matrix((val, (row, col)), shape=(count, count)).tocsr()
        jac += coo_matrix((mobile_coefficient * screen, (water_ids, water_ids)),
                          shape=(count, count)).tocsr()
        scale = max(1.0, float(np.max(np.abs(face_coefficient * sigma))),
                    float(np.max(mobile_coefficient * np.sum(concentrations, axis=0))))
        return residual, jac, concentrations, sigma, surface_u, scale

    u = np.zeros(count)
    for iteration in range(max_iter + 1):
        residual, jac, concentrations, sigma, surface_u, scale = evaluate(u, True)
        norm = float(np.max(np.abs(residual)) / scale)
        if norm <= tolerance:
            break
        if iteration == max_iter:
            raise RuntimeError(f"dielectric fiber Newton failed: {norm:g}")
        step = spsolve(jac, -residual)
        if not np.all(np.isfinite(step)):
            raise FloatingPointError("nonfinite dielectric PB step")
        maximum = float(np.max(np.abs(step)))
        if maximum > 2:
            step *= 2 / maximum
        old_norm = float(np.max(np.abs(residual)))
        alpha = 1.0
        for _ in range(30):
            trial = u + alpha * step
            if float(np.max(np.abs(evaluate(trial)))) < old_norm:
                u = trial
                break
            alpha *= 0.5
        else:
            raise RuntimeError(f"dielectric fiber line search stalled: {norm:g}")
    residual, jac, concentrations, sigma, surface_u, scale = evaluate(u, True)
    capacity = E * chemical.gamma_nm2 * 1e18
    if not np.all(np.isfinite(sigma)) or float(np.max(np.abs(sigma))) > capacity * (1 + 1e-12):
        raise AssertionError("dielectric surface capacity exceeded")
    mobile = F * float((Z @ concentrations) @ mesh.volumes_m2)
    bound = float(sigma @ face_length)
    outer = eps_vt * float(np.sum(outer_conductance * u[outer_ids]))
    imbalance = mobile + bound - outer
    denominator = max(float(np.sum(np.abs(sigma) * face_length)),
                      capacity * geometry.wet_area_m2_per_m * 1e-6)
    by_layer = {}
    for layer in geometry.layers:
        mask = face_layer == layer.layer_id
        area = float(np.sum(face_length[mask]))
        if not math.isclose(area, layer.wet_perimeter_nm * 1e-9, rel_tol=1e-10):
            raise AssertionError("dielectric layer area mismatch")
        by_layer[layer.layer_id] = float((sigma[mask] @ face_length[mask]) / area)
    outer_mask = np.isin(face_layer, [l.layer_id for l in geometry.layers if l.is_outer])
    outer_area = float(np.sum(face_length[outer_mask]))
    sigma_outer = float((sigma[outer_mask] @ face_length[outer_mask]) / outer_area)
    core_mask = np.isin(face_layer, [l.layer_id for l in geometry.layers if l.is_core])
    core_area = float(np.sum(face_length[core_mask]))
    sigma_core = float((sigma[core_mask] @ face_length[core_mask]) / core_area) if core_area else None
    max_species_M = float(np.max(concentrations)) / 1000
    return DielectricResult(bound / geometry.wet_area_m2_per_m, sigma_outer, sigma_core,
                            bound / geometry.projected_area_m2_per_m, by_layer,
                            abs(imbalance), abs(imbalance) / denominator if denominator else None,
                            max_species_M, float(np.min(concentrations)) / 1000,
                            max_species_M > 0.5, float(np.max(np.abs(residual)) / scale),
                            iteration, mesh, u[water_ids], concentrations, sigma, surface_u,
                            mobile, bound, outer, u[solid_ids], eps_PE_r)
