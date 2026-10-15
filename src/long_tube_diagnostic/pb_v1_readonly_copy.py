"""Teaching models for PE salt diffusion. Python standard library only.

All numerical inputs use SI: m, s, K, mol/m^3, C/m^2. These are
synthetic/ideal models, not an experimental fit or a finite-length PNP solver.
"""
import math
from dataclasses import dataclass

F = 96485.33212                 # C/mol
R = 8.31446261815324            # J/(mol K)
EPS0 = 8.8541878128e-12         # F/m
NACL_MOLAR_MASS = 58.44277       # g/mol


def molar_to_si(c_mol_l):
    return c_mol_l * 1000.0


def si_to_mg_l(c_mol_m3, molar_mass=NACL_MOLAR_MASS):
    # 1 mg/L = 1 g/m^3, hence no extra 1000 factor.
    return c_mol_m3 * molar_mass


def mg_l_to_si(c_mg_l, molar_mass=NACL_MOLAR_MASS):
    return c_mg_l / molar_mass


def debye_length(c_mol_m3, temperature=298.15, eps_r=78.3):
    """Ideal monovalent symmetric electrolyte; c is each ion's bulk molarity."""
    if c_mol_m3 <= 0 or temperature <= 0 or eps_r <= 0:
        raise ValueError("Debye/PB requires positive salt, T and permittivity")
    return math.sqrt(EPS0 * eps_r * R * temperature / (2 * F * F * c_mol_m3))


def ambipolar_diffusivity(d_plus=1.33e-9, d_minus=2.03e-9):
    if d_plus <= 0 or d_minus <= 0:
        raise ValueError("diffusivities must be positive")
    return 2 * d_plus * d_minus / (d_plus + d_minus)


@dataclass(frozen=True)
class Geometry:
    radius_m: float = 25e-9
    length_m: float = 350e-9
    porosity: float = 0.40
    wetted_fraction: float = 1.0
    area_m2: float = 1e-4
    volume_left_m3: float = 10e-6
    volume_right_m3: float = 10e-6
    layer_left_m: float = 100e-6
    layer_right_m: float = 100e-6

    def validate(self):
        for name in ("radius_m", "length_m", "area_m2", "volume_left_m3", "volume_right_m3"):
            if getattr(self, name) <= 0:
                raise ValueError(name + " must be positive")
        if not 0 < self.porosity <= 1:
            raise ValueError("porosity must be in (0, 1]")
        if not 0 <= self.wetted_fraction <= 1:
            raise ValueError("wetted_fraction must be in [0, 1]")
        if self.layer_left_m < 0 or self.layer_right_m < 0:
            raise ValueError("external layer lengths cannot be negative")
        return self


def neutral_resistances(g, diffusivity):
    """Projected-area resistances [s/m]; circular independent-pore access.

    Per pore: internal L/(D*pi*a^2), two-side access 1/(2*D*a).
    Effective pore density fw*phi/(pi*a^2). Outer planar layers use the
    nominal projected area and are NOT silently rescaled by fw.
    """
    g.validate()
    if diffusivity <= 0:
        raise ValueError("diffusivity must be positive")
    if g.wetted_fraction == 0:
        return dict(internal_s_m=math.inf, access_s_m=math.inf,
                    external_s_m=(g.layer_left_m + g.layer_right_m) / diffusivity,
                    total_s_m=math.inf, permeability_m_s=0.0)
    fraction = g.porosity * g.wetted_fraction
    internal = g.length_m / (diffusivity * fraction)
    access = math.pi * g.radius_m / (2 * diffusivity * fraction)
    external = (g.layer_left_m + g.layer_right_m) / diffusivity
    total = internal + access + external
    return dict(internal_s_m=internal, access_s_m=access,
                external_s_m=external, total_s_m=total,
                permeability_m_s=1 / total)


def neutral_curve(g, c_left0, c_right0, times_s, diffusivity):
    """Exact well-mixed finite-reservoir result; salt-formula-unit balance."""
    if min(c_left0, c_right0) < 0 or any(t < 0 for t in times_s):
        raise ValueError("concentrations/times must be nonnegative")
    p = neutral_resistances(g, diffusivity)["permeability_m_s"]
    vl, vr = g.volume_left_m3, g.volume_right_m3
    eq = (vl * c_left0 + vr * c_right0) / (vl + vr)
    rate = g.area_m2 * p * (1 / vl + 1 / vr)
    rows = []
    for t in times_s:
        delta = (c_left0 - c_right0) * math.exp(-rate * t)
        left = eq + vr / (vl + vr) * delta
        right = eq - vl / (vl + vr) * delta
        rows.append(dict(time_s=t, c_left_mol_m3=left, c_right_mol_m3=right,
                         c_right_mg_l=si_to_mg_l(right),
                         salt_flux_mol_m2_s=p * delta,
                         total_salt_mol=vl * left + vr * right))
    return rows


def donnan_partition(c_bulk, sigma_c_m2, radius_m, temperature=298.15):
    """Uniform-potential/electroneutral overlap limit. Not general pore PB."""
    if c_bulk < 0 or radius_m <= 0:
        raise ValueError("bad Donnan concentration/radius")
    x = 2 * sigma_c_m2 / (F * radius_m)  # signed fixed charge, mol/m^3
    a = abs(x)
    # Stable small coion concentration, including c=0 at finite charge.
    co = c_bulk if a == 0 else c_bulk * c_bulk / (math.hypot(c_bulk, a / 2) + a / 2)
    counter = co + a
    cp, cm = (counter, co) if x < 0 else (co, counter)
    if x == 0:
        cp = cm = c_bulk
    u = math.asinh(x / (2 * c_bulk)) if c_bulk > 0 else (
        0.0 if x == 0 else math.copysign(math.inf, x))
    return dict(fixed_charge_signed_mol_m3=x, c_plus_mol_m3=cp,
                c_minus_mol_m3=cm, c_coion_mol_m3=co,
                coion_partition=co / c_bulk if c_bulk else None,
                reduced_potential=u, potential_v=u * R * temperature / F)


def _tridiagonal(lower, diagonal, upper, rhs):
    """Thomas solver; PB Newton matrices are symmetric positive definite."""
    n = len(diagonal)
    d, b = diagonal[:], rhs[:]
    for i in range(1, n):
        factor = lower[i] / d[i - 1]
        d[i] -= factor * upper[i - 1]
        b[i] -= factor * b[i - 1]
    out = [0.0] * n
    out[-1] = b[-1] / d[-1]
    for i in range(n - 2, -1, -1):
        out[i] = (b[i] - upper[i] * out[i + 1]) / d[i]
    return out


def cylindrical_pb(c_bulk, sigma_c_m2, radius_m=25e-9, cells=256,
                   temperature=298.15, eps_r=78.3, tol=1e-11, max_iter=80):
    """Nonlinear radial PB, annular finite volumes + damped Newton.

    u=F*psi/(R*T), x=r/a: (x*u')'= (kappa*a)^2*x*sinh(u).
    At x=0, x*u'=0; at x=1, u'=sigma*a/(eps*V_T).
    The cross-section is in radial equilibrium with a virtual bulk salt.
    It has no axial flux, no pore ends, and no charge-regulation chemistry.
    """
    if radius_m <= 0 or cells < 8 or int(cells) != cells:
        raise ValueError("radius positive; integer cells >=8")
    lam = debye_length(c_bulk, temperature, eps_r)
    vt = R * temperature / F
    eps = EPS0 * eps_r
    k2 = (radius_m / lam) ** 2
    wall_flux = sigma_c_m2 * radius_m / (eps * vt)
    dx = 1.0 / cells
    weights = [((i + 1) ** 2 - i ** 2) * dx * dx / 2 for i in range(cells)]
    left = [float(i) for i in range(cells)]
    right = [float(i + 1) if i + 1 < cells else 0.0 for i in range(cells)]
    u = [math.asinh(2 * wall_flux / k2)] * cells

    def residual(profile):
        ans = []
        for i in range(cells):
            value = k2 * weights[i] * math.sinh(profile[i])
            if i > 0:
                value += left[i] * (profile[i] - profile[i - 1])
            if i + 1 < cells:
                value += right[i] * (profile[i] - profile[i + 1])
            else:
                value -= wall_flux
            ans.append(value)
        return ans

    scale = 1 + abs(wall_flux)
    iterations = 0
    for iterations in range(max_iter + 1):
        residuals = residual(u)
        norm = max(abs(v) for v in residuals) / scale
        if norm <= tol:
            break
        if iterations == max_iter:
            raise RuntimeError("PB Newton did not converge")
        diag = [left[i] + right[i] + k2 * weights[i] * math.cosh(u[i])
                for i in range(cells)]
        delta = _tridiagonal([-v for v in left], diag, [-v for v in right],
                             [-v for v in residuals])
        step = 1.0
        while step >= 2 ** -24:
            candidate = [u[i] + step * delta[i] for i in range(cells)]
            try:
                new_norm = max(abs(v) for v in residual(candidate)) / scale
            except OverflowError:
                new_norm = math.inf
            if new_norm < norm:
                u = candidate
                break
            step /= 2
        else:
            raise RuntimeError("PB Newton line search stalled")

    kp = 2 * sum(w * math.exp(-v) for w, v in zip(weights, u))
    km = 2 * sum(w * math.exp(v) for w, v in zip(weights, u))
    integrated_rhs = k2 * sum(w * math.sinh(v) for w, v in zip(weights, u))
    gauss_error = abs(integrated_rhs - wall_flux) / scale
    # Boundary-face extrapolation; report its numerical approximation explicitly.
    u_wall = u[-1] + wall_flux * dx / 2
    u_center = u[0] - k2 * math.sinh(u[0]) * dx * dx / 16
    rows = [dict(r_over_a=(i + 0.5) * dx, r_nm=(i + 0.5) * dx * radius_m * 1e9,
                 potential_mv=u[i] * vt * 1000, c_plus_over_bulk=math.exp(-u[i]),
                 c_minus_over_bulk=math.exp(u[i])) for i in range(cells)]
    return dict(profile=rows, cells=cells, iterations=iterations,
                residual_scaled=norm, gauss_error_scaled=gauss_error,
                debye_nm=lam * 1e9, kappa_a=radius_m / lam,
                a_over_gouy_chapman=abs(sigma_c_m2) * radius_m / (2 * eps * vt),
                center_potential_mv=u_center * vt * 1000,
                wall_potential_extrapolated_mv=u_wall * vt * 1000,
                c_plus_partition=kp, c_minus_partition=km,
                coion_partition=km if sigma_c_m2 < 0 else kp if sigma_c_m2 > 0 else 1.0,
                mobile_charge_mean_c_m3=F * c_bulk * (kp - km),
                wall_charge_per_pore_volume_c_m3=2 * sigma_c_m2 / radius_m)


def donnan_transport_integral(c_left, c_right, sigma, radius, d_plus, d_minus):
    """Integral of ideal electroneutral zero-current 1D Donnan diffusivity.

    J+=J-; fixed charge constant; no convection. Returns mol/(m s).
    This approximation does not import the PB profile as an axial PNP result.
    """
    if min(c_left, c_right) < 0 or radius <= 0:
        raise ValueError("bad transport concentration/radius")
    d_salt = ambipolar_diffusivity(d_plus, d_minus)
    charge = 2 * abs(sigma) / (F * radius)
    if charge == 0:
        return d_salt * (c_left - c_right)
    q_left = donnan_partition(c_left, sigma, radius)["c_coion_mol_m3"]
    q_right = donnan_partition(c_right, sigma, radius)["c_coion_mol_m3"]
    d_counter = d_plus if sigma < 0 else d_minus
    ds = d_plus + d_minus
    b = d_counter * charge
    log_ratio = math.log1p(ds * (q_left - q_right) / (ds * q_right + b))
    return d_plus * d_minus * (2 * (q_left - q_right) / ds +
                              (charge - 2 * b / ds) / ds * log_ratio)


def reduced_donnan_flux(g, c_left, c_right, sigma, d_plus=1.33e-9, d_minus=2.03e-9):
    """Donnan membrane in series with NEUTRAL access/external layers.

    Scalar bisection enforces the same salt flux through each region.
    Boundary concentration matching is a teaching approximation, not a
    charged-access solution; neither electric current nor capacitance is evolved.
    """
    g.validate()
    if min(c_left, c_right) < 0:
        raise ValueError("negative reservoir concentration")
    if g.wetted_fraction == 0 or c_left == c_right:
        return 0.0
    d_salt = ambipolar_diffusivity(d_plus, d_minus)
    fwphi = g.wetted_fraction * g.porosity
    half_access = math.pi * g.radius_m / (4 * d_salt * fwphi)
    rl = half_access + g.layer_left_m / d_salt
    rr = half_access + g.layer_right_m / d_salt
    prefactor = fwphi / g.length_m

    def membrane(j):
        bl = max(0.0, c_left - j * rl)
        br = c_right + j * rr
        return prefactor * donnan_transport_integral(bl, br, sigma, g.radius_m, d_plus, d_minus)

    bound = (c_left - c_right) / (rl + rr)
    low, high = (0.0, bound) if bound > 0 else (bound, 0.0)
    for _ in range(70):
        mid = (low + high) / 2
        if mid < membrane(mid):
            low = mid
        else:
            high = mid
    return (low + high) / 2


def reduced_curve(g, c_left0, c_right0, times_s, sigma,
                  d_plus=1.33e-9, d_minus=2.03e-9, max_step_s=30.0):
    """RK4 in reservoir concentration difference; salt balance exact by reconstruction.

    Step halving protects positivity, but accuracy must be checked with another
    max_step_s. This is not an adaptive error-certified PNP transient.
    """
    g.validate()
    if min(c_left0, c_right0) < 0 or max_step_s <= 0:
        raise ValueError("bad initial concentration/time step")
    times = list(times_s)
    if any(t < 0 for t in times) or times != sorted(times):
        raise ValueError("times must be nonnegative and sorted")
    vl, vr = g.volume_left_m3, g.volume_right_m3
    eq = (vl * c_left0 + vr * c_right0) / (vl + vr)
    factor = g.area_m2 * (1 / vl + 1 / vr)
    delta = c_left0 - c_right0

    def concentrations(diff):
        return eq + vr / (vl + vr) * diff, eq - vl / (vl + vr) * diff

    def rhs(diff):
        l, r = concentrations(diff)
        if l < 0 or r < 0:
            raise ValueError("RK stage left physical concentration interval")
        return -factor * reduced_donnan_flux(g, l, r, sigma, d_plus, d_minus)

    current, rows = 0.0, []
    for target in times:
        while current < target - 1e-12:
            dt = min(max_step_s, target - current)
            while True:
                try:
                    k1 = rhs(delta)
                    k2 = rhs(delta + dt * k1 / 2)
                    k3 = rhs(delta + dt * k2 / 2)
                    k4 = rhs(delta + dt * k3)
                    candidate = delta + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
                    ok = candidate * delta >= 0 and abs(candidate) <= abs(delta) + 1e-12
                except ValueError:
                    ok = False
                if ok:
                    break
                dt /= 2
                if dt < 1e-10:
                    raise RuntimeError("reservoir integration stalled")
            delta = candidate
            current += dt
        left, right = concentrations(delta)
        rows.append(dict(time_s=target, c_left_mol_m3=left, c_right_mol_m3=right,
                         c_right_mg_l=si_to_mg_l(right),
                         salt_flux_mol_m2_s=reduced_donnan_flux(g, left, right, sigma, d_plus, d_minus),
                         total_salt_mol=vl * left + vr * right))
    return rows
