"""Restricted PB + local-radial-equilibrium, no-flow long-tube NP diagnostic.

Not a finite-length/entrance PNP solver. Inputs use SI. Axial flux is per
WATER-FILLED PORE cross-sectional area, not nominal membrane area.
"""
import math
from functools import lru_cache
from pb_v1_readonly_copy import cylindrical_pb, donnan_partition, ambipolar_diffusivity, R, F


def _coefficients(k_plus, k_minus, d_plus, d_minus):
    gp = d_plus * k_plus
    gm = d_minus * k_minus
    return dict(k_plus=k_plus, k_minus=k_minus,
                effective_salt_D_m2_s=2 * gp * gm / (gp + gm),
                d_reduced_virtual_potential_d_log_c=(gm - gp) / (gp + gm))


@lru_cache(maxsize=10000)
def local_pb_transport(c, sigma, radius, cells=256, temperature=298.15,
                       eps_r=78.3, d_plus=1.33e-9, d_minus=2.03e-9):
    if c <= 0 or radius <= 0 or min(d_plus, d_minus) <= 0:
        raise ValueError("positive virtual concentration/radius/diffusivities required")
    pb = cylindrical_pb(c, sigma, radius, cells, temperature, eps_r)
    values = _coefficients(pb["c_plus_partition"], pb["c_minus_partition"], d_plus, d_minus)
    values.update(pb_residual_scaled=pb["residual_scaled"], pb_gauss_error_scaled=pb["gauss_error_scaled"],
                  kappa_a=pb["kappa_a"], a_over_gc=pb["a_over_gouy_chapman"])
    return values


def local_donnan_transport(c, sigma, radius, d_plus=1.33e-9, d_minus=2.03e-9):
    partition = donnan_partition(c, sigma, radius)
    return _coefficients(partition["c_plus_mol_m3"] / c,
                         partition["c_minus_mol_m3"] / c, d_plus, d_minus)


def integrate_tube(c_left, c_right, sigma, radius=25e-9, length=350e-9,
                   cells=256, intervals=128, temperature=298.15, eps_r=78.3,
                   d_plus=1.33e-9, d_minus=2.03e-9, approximation="pb"):
    """Composite Simpson in y=ln(c), including the dc=c dy Jacobian.

    x goes left->right. Positive common salt flux goes left->right.
    Output delta virtual V = Phi_v(right)-Phi_v(left), not electrode voltage.
    """
    if min(c_left, c_right, length, radius) <= 0:
        raise ValueError("positive endpoint virtual concentrations, length and radius required; DI not silently regularized")
    if intervals < 2 or intervals % 2:
        raise ValueError("positive even Simpson intervals required")
    if approximation not in ("pb", "donnan"):
        raise ValueError("approximation must be pb or donnan")
    yl, yr = math.log(c_left), math.log(c_right)
    ylo, yhi = min(yl, yr), max(yl, yr)
    dy = (yhi - ylo) / intervals
    sums = [0.0, 0.0]
    rows = []
    for i in range(intervals + 1):
        y = ylo + i * dy
        c = math.exp(y)
        vals = local_pb_transport(c, sigma, radius, cells, temperature, eps_r, d_plus, d_minus) if approximation == "pb" else local_donnan_transport(c, sigma, radius, d_plus, d_minus)
        w = 1 if i in (0, intervals) else 4 if i % 2 else 2
        sums[0] += w * vals["effective_salt_D_m2_s"] * c
        sums[1] += w * vals["d_reduced_virtual_potential_d_log_c"]
        rows.append(dict(c_virtual_mol_m3=c, **vals))
    int_d_dc = sums[0] * dy / 3
    int_h_dlogc = sums[1] * dy / 3
    direction = 1 if c_left >= c_right else -1
    salt_flux = direction * int_d_dc / length
    delta_v = -direction * int_h_dlogc * R * temperature / F
    neutral_flux = ambipolar_diffusivity(d_plus, d_minus) * (c_left - c_right) / length
    return dict(c_left_virtual_mol_m3=c_left, c_right_virtual_mol_m3=c_right,
                sigma_assumed_c_m2=sigma, radius_assumed_m=radius, length_normalization_m=length,
                pb_cells=cells, log_c_Simpson_intervals=intervals, approximation=approximation,
                pore_area_salt_flux_mol_m2_s=salt_flux,
                neutral_pore_area_salt_flux_mol_m2_s=neutral_flux,
                charged_to_neutral_salt_flux_ratio=salt_flux / neutral_flux if neutral_flux else None,
                virtual_delta_phi_right_minus_left_mv=delta_v * 1000,
                max_pb_residual_scaled=max((r.get("pb_residual_scaled", 0) for r in rows)),
                max_pb_gauss_error_scaled=max((r.get("pb_gauss_error_scaled", 0) for r in rows)),
                status="hypothetical long-tube no-flow local-radial-equilibrium NP; NOT finite-length/entrance PNP or PE validation",
                local_coefficients=rows)
