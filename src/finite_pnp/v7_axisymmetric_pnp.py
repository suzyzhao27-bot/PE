"""Axisymmetric aqueous cylinder PNP for the frozen v7 local reference.

The domain is 0<=r<=a, 0<=z<=L. No external bath or access resistance.
Unknown concentrations are logarithmic, so Newton iterates remain positive.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import spsolve

from solver import bernoulli


F = 96485.33212
R = 8.31446261815324
EPS0 = 8.8541878128e-12
D_NA = 1.33e-9
D_CL = 2.03e-9
D_REF = 1e-9


@dataclass(frozen=True)
class Cylinder:
    a_m: float
    length_m: float
    nr: int = 24
    nz: int = 240

    def build(self):
        if self.a_m <= 0 or self.length_m <= 0 or self.nr < 2 or self.nz < 2:
            raise ValueError("positive cylinder dimensions and at least two cells per axis")
        nr, nz = self.nr, self.nz
        dr, dz = self.a_m / nr, self.length_m / nz
        radial_face = np.arange(nr + 1) * dr
        ring_area = (radial_face[1:]**2 - radial_face[:-1]**2) / 2  # per 2*pi
        volume = np.tile(ring_area * dz, nz)
        edge_i, edge_j, weights = [], [], []
        for z in range(nz):
            for r in range(nr):
                i = z * nr + r
                if r + 1 < nr:
                    edge_i.append(i); edge_j.append(i + 1)
                    weights.append(radial_face[r + 1] * dz / dr)
                if z + 1 < nz:
                    edge_i.append(i); edge_j.append(i + nr)
                    weights.append(ring_area[r] / dz)
        left = np.arange(nr)
        right = (nz - 1) * nr + np.arange(nr)
        boundary_nodes = np.r_[left, right]
        boundary_weights = np.r_[2 * ring_area / dz, 2 * ring_area / dz]
        ei, ej = np.array(edge_i, dtype=int), np.array(edge_j, dtype=int)
        weights = np.array(weights) / self.a_m
        boundary_weights /= self.a_m
        diagonal = np.bincount(np.r_[ei, ej], weights=np.r_[weights, weights], minlength=nr * nz)
        np.add.at(diagonal, boundary_nodes, boundary_weights)
        laplacian = coo_matrix((np.r_[diagonal, -weights, -weights],
                                (np.r_[np.arange(nr * nz), ei, ej],
                                 np.r_[np.arange(nr * nz), ej, ei])),
                               shape=(nr * nz, nr * nz)).tocsr()
        return Mesh(self, volume, ei, ej, weights, boundary_nodes,
                    boundary_weights, laplacian, diagonal)


@dataclass
class Mesh:
    cylinder: Cylinder
    volume: np.ndarray
    ei: np.ndarray
    ej: np.ndarray
    weights: np.ndarray
    boundary_nodes: np.ndarray
    boundary_weights: np.ndarray
    laplacian: object
    diagonal: np.ndarray

    @property
    def n(self):
        return self.cylinder.nr * self.cylinder.nz


def initial_guess(mesh: Mesh, c_left: float, c_right: float, right_u: float, open_circuit: bool):
    z = np.repeat((np.arange(mesh.cylinder.nz) + 0.5) / mesh.cylinder.nz,
                  mesh.cylinder.nr)
    c = c_left * np.exp(z * math.log(c_right / c_left))
    y = np.r_[right_u * z, np.log(c), np.log(c)]
    return np.r_[y, right_u] if open_circuit else y


def assemble(mesh: Mesh, y: np.ndarray, c_left: float, c_right: float,
             sigma: float, temperature: float, eps_r: float,
             open_circuit: bool, fixed_right_u: float, jacobian: bool = True):
    n, nr = mesh.n, mesh.cylinder.nr
    u = y[:n]
    cp, cm = np.exp(y[n:2*n]), np.exp(y[2*n:3*n])
    right_u = float(y[-1]) if open_circuit else fixed_right_u
    eps = EPS0 * eps_r
    charge = F*F / (eps*R*temperature) * mesh.volume / mesh.cylinder.a_m
    wall = np.zeros(n)
    wall[np.arange(nr - 1, n, nr)] = (
        F * sigma / (eps*R*temperature) * mesh.cylinder.a_m
        * (mesh.cylinder.length_m / mesh.cylinder.nz) / mesh.cylinder.a_m)
    boundary_u = np.r_[np.zeros(nr), np.full(nr, right_u)]
    boundary_c = np.r_[np.full(nr, c_left), np.full(nr, c_right)]
    rhs_u = np.bincount(mesh.boundary_nodes,
                        weights=mesh.boundary_weights * boundary_u, minlength=n)
    residual = np.zeros(3*n + int(open_circuit))
    residual[:n] = mesh.laplacian @ u - charge*(cp-cm) - wall - rhs_u
    jr, jc, jv = [], [], []

    def add(rows, cols, vals):
        jr.extend(np.asarray(rows, dtype=int).ravel().tolist())
        jc.extend(np.asarray(cols, dtype=int).ravel().tolist())
        jv.extend(np.asarray(vals, dtype=float).ravel().tolist())

    if jacobian:
        l = mesh.laplacian.tocoo()
        add(l.row, l.col, l.data)
        add(np.arange(n), n + np.arange(n), -charge*cp)
        add(np.arange(n), 2*n + np.arange(n), charge*cm)
        if open_circuit:
            add(mesh.boundary_nodes[nr:], np.full(nr, 3*n), -mesh.boundary_weights[nr:])

    boundary_fluxes = []
    for sign, c, offset, d_ratio in ((1, cp, n, D_NA/D_REF),
                                      (-1, cm, 2*n, D_CL/D_REF)):
        ei, ej, w = mesh.ei, mesh.ej, mesh.weights
        v = sign * (u[ej] - u[ei])
        bv, dbv = bernoulli(v)
        bm, dbm = bernoulli(-v)
        flux = w * (bv*c[ei] - bm*c[ej])
        q = w * (dbv*c[ei] + dbm*c[ej])
        np.add.at(residual, offset + ei, flux)
        np.add.at(residual, offset + ej, -flux)
        if jacobian:
            add(offset + ei, offset + ei, w*bv*c[ei])
            add(offset + ei, offset + ej, -w*bm*c[ej])
            add(offset + ej, offset + ei, -w*bv*c[ei])
            add(offset + ej, offset + ej, w*bm*c[ej])
            add(offset + ei, ei, -sign*q)
            add(offset + ei, ej, sign*q)
            add(offset + ej, ei, sign*q)
            add(offset + ej, ej, -sign*q)
        bn, bw = mesh.boundary_nodes, mesh.boundary_weights
        vb = sign * (boundary_u - u[bn])
        bb, dbb = bernoulli(vb)
        bm, dbm = bernoulli(-vb)
        bflux = bw * (bb*c[bn] - bm*boundary_c)
        bq = bw * (dbb*c[bn] + dbm*boundary_c)
        np.add.at(residual, offset + bn, bflux)
        boundary_fluxes.append(bflux)
        if jacobian:
            add(offset + bn, offset + bn, bw*bb*c[bn])
            add(offset + bn, bn, -sign*bq)
            if open_circuit:
                right_nodes = bn[nr:]
                add(offset + right_nodes, np.full(nr, 3*n), sign*bq[nr:])
                current_sign = sign
                add(np.full(nr, 3*n), offset + right_nodes,
                    current_sign*d_ratio*bw[nr:]*bb[nr:]*c[right_nodes])
                add(np.full(nr, 3*n), right_nodes,
                    -current_sign*d_ratio*sign*bq[nr:])
                add(np.full(nr, 3*n), np.full(nr, 3*n),
                    current_sign*d_ratio*sign*bq[nr:])
    if open_circuit:
        residual[-1] = (D_NA/D_REF)*sum(boundary_fluxes[0][nr:]) - (
            D_CL/D_REF)*sum(boundary_fluxes[1][nr:])
    scale = np.r_[1 + mesh.diagonal + charge*(cp+cm) + abs(wall),
                  1 + mesh.diagonal*cp, 1 + mesh.diagonal*cm]
    if open_circuit:
        scale = np.r_[scale, 1 + max(c_left,c_right)*sum(mesh.boundary_weights[nr:])
                      * max(D_NA,D_CL)/D_REF]
    jac = coo_matrix((jv, (jr,jc)), shape=(len(y), len(y))).tocsr() if jacobian else None
    return residual, jac, scale, (cp, cm, u, right_u, boundary_fluxes, wall, charge)


def solve(c_left: float, c_right: float, sigma: float, cylinder: Cylinder,
          temperature: float = 298.15, eps_r: float = 78.5,
          open_circuit: bool = True, right_voltage_v: float = 0.0,
          initial: np.ndarray | None = None, tol: float = 2e-10, max_iter: int = 60):
    if min(c_left, c_right, temperature, eps_r) <= 0:
        raise ValueError("positive aqueous concentrations and conditions required")
    mesh = cylinder.build()
    vt = R*temperature/F
    right_u0 = (D_CL-D_NA)/(D_CL+D_NA)*math.log(c_right/c_left) if open_circuit else right_voltage_v/vt
    y = initial.copy() if initial is not None else initial_guess(mesh,c_left,c_right,right_u0,open_circuit)
    if len(y) != 3*mesh.n + int(open_circuit):
        raise ValueError("initial state belongs to another mesh or boundary mode")
    history = []
    for iteration in range(max_iter+1):
        residual, jac, scale, extra = assemble(mesh,y,c_left,c_right,sigma,temperature,eps_r,
                                                open_circuit,right_u0)
        norm = float(np.max(np.abs(residual)/scale))
        history.append(norm)
        if norm < tol:
            break
        if iteration == max_iter:
            raise RuntimeError(f"axisymmetric Newton max_iter: residual={norm:g}")
        step = spsolve(diags(1/scale)@jac, -residual/scale)
        if not np.all(np.isfinite(step)):
            raise RuntimeError("nonfinite Newton step")
        factor = min(1.0, 1.5/max(1.5,float(np.max(abs(step[mesh.n:3*mesh.n])))),
                     2.0/max(2.0,float(np.max(abs(step[:mesh.n])))))
        for _ in range(28):
            candidate = y + factor*step
            if np.max(abs(candidate[mesh.n:3*mesh.n])) < 100:
                test,_,_,_ = assemble(mesh,candidate,c_left,c_right,sigma,temperature,eps_r,
                                      open_circuit,right_u0,jacobian=False)
                if float(np.max(np.abs(test)/scale)) < norm:
                    y = candidate
                    break
            factor *= 0.5
        else:
            raise RuntimeError(f"axisymmetric Newton line search stalled: residual={norm:g}")
    cp, cm, u, right_u, boundary_fluxes, wall, charge = extra
    plus, minus = boundary_fluxes
    rate_plus_r = 2*math.pi*cylinder.a_m*D_NA*float(sum(plus[cylinder.nr:]))
    rate_minus_r = 2*math.pi*cylinder.a_m*D_CL*float(sum(minus[cylinder.nr:]))
    rate_plus_l = -2*math.pi*cylinder.a_m*D_NA*float(sum(plus[:cylinder.nr]))
    rate_minus_l = -2*math.pi*cylinder.a_m*D_CL*float(sum(minus[:cylinder.nr]))
    denom = max(abs(rate_plus_r)+abs(rate_minus_r),1e-300)
    charge_rhs = float(sum(charge*(cp-cm)) + sum(wall))
    boundary_field = float(sum(mesh.boundary_weights*(u[mesh.boundary_nodes]-
                                 np.r_[np.zeros(cylinder.nr),np.full(cylinder.nr,right_u)])))
    volume = math.pi*cylinder.a_m**2*cylinder.length_m
    inventory = 2*math.pi*float(np.dot(mesh.volume,(cp+cm)/2))
    result = {"a_m":cylinder.a_m,"L_m":cylinder.length_m,"nr":cylinder.nr,"nz":cylinder.nz,
              "c_L_mol_m3":c_left,"c_R_mol_m3":c_right,"sigma_C_m2":sigma,
              "delta_phi_V":right_u*vt,"open_circuit":open_circuit,
              "F_plus_mol_s":rate_plus_r,"F_minus_mol_s":rate_minus_r,
              "F_salt_mol_s":(rate_plus_r+rate_minus_r)/2,
              "plus_conservation_relative":abs(rate_plus_l-rate_plus_r)/denom,
              "minus_conservation_relative":abs(rate_minus_l-rate_minus_r)/denom,
              "current_relative":abs(rate_plus_r-rate_minus_r)/denom,
              "poisson_gauss_relative":abs(boundary_field-charge_rhs)/
              max(1.0,abs(boundary_field),abs(charge_rhs)),
              "salt_inventory_mol":inventory,"mean_salt_concentration_mol_m3":inventory/volume,
              "min_c_mol_m3":float(min(cp.min(),cm.min())),
              "residual_scaled_max":norm,"newton_iterations":iteration}
    return result, y, history
