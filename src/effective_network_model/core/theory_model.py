"""Shared effective network field and spatial partition resistance.

The field is a frozen structural hypothesis. It is not a microscopic dielectric
self-energy calculation. Core calculations require only Python's standard library.
Lengths are nm here. Absolute transport requires independent prefactors.
"""
import math

D_NA = 1.33e-9
D_CL = 2.03e-9
D_AMB = 2*D_NA*D_CL/(D_NA+D_CL)
FARADAY = 96485.33212


def field(z_nm, h_nm, ws, ell_nm, robin_r=0.0):
    if h_nm <= 0 or ell_nm <= 0 or not 0 <= z_nm <= h_nm:
        raise ValueError('Positive lengths and 0 <= z <= h are required')
    denominator = 1 + robin_r*math.exp(-h_nm/ell_nm)
    if denominator <= 0:
        raise ValueError('Invalid Robin denominator')
    return ws*(math.exp(-z_nm/ell_nm)+math.exp(-(h_nm-z_nm)/ell_nm))/denominator


def log_integral_exp_field(h_nm, ws, ell_nm, robin_r=0.0, n=4096):
    """log integral exp(w) dz, Simpson quadrature with exponential scaling."""
    if n < 4 or n % 2:
        raise ValueError('Simpson n must be even and >= 4')
    wmax = max(field(0,h_nm,ws,ell_nm,robin_r),
               field(h_nm/2,h_nm,ws,ell_nm,robin_r))
    weighted = 0.0
    for i in range(n+1):
        weight = 1 if i in (0,n) else (4 if i % 2 else 2)
        weighted += weight*math.exp(field(h_nm*i/n,h_nm,ws,ell_nm,robin_r)-wmax)
    return wmax+math.log(h_nm*weighted/(3*n))


def logadd(a, b):
    m = max(a,b)
    return m+math.log(math.exp(a-m)+math.exp(b-m))


def log_resistance_nm(h_nm, ws, ell_nm, a_nm=25.0, robin_r=0.0,
                      access=True, n=4096):
    if a_nm <= 0:
        raise ValueError('Local radius must be positive')
    logr = log_integral_exp_field(h_nm,ws,ell_nm,robin_r,n)
    # Independent dilute-aperture geometric reference; not a measured outer layer.
    return logadd(logr,math.log(math.pi*a_nm/2)) if access else logr


def relative_permeance(h_nm, a_nm, params, reference_h_nm=2600.0):
    """Conditional b/a coupling inspired by interfacial area per volume.

    This is an explicit effective S/V closure, not a measured radius dependence.
    The fit identifies ws at a_ref only, not b and radius separately.
    """
    aref = params.get('a_ref_nm',25.0)
    ws_ref = params['ws_ref']
    ell = params['ell_nm']
    r = params.get('robin_r',0.0)
    lr = log_resistance_nm(h_nm,ws_ref*aref/a_nm,ell,a_nm,r)
    lref = log_resistance_nm(reference_h_nm,ws_ref,ell,aref,r)
    return math.exp(lref-lr)


def thickness_log_slope(h_nm,ws,ell_nm,a_nm=25.,n=1024):
    """d ln(P)/d ln(h) from the exact derivative of h*int_0^1 exp(w)du.

    Main r=0 boundary choice. A positive value means permeability increases
    when the film gets thicker at fixed local radius and shared field parameters.
    """
    wmax=field(0,h_nm,ws,ell_nm)
    weight_sum,derivative_sum=0.,0.
    for i in range(n+1):
        u=i/n
        first=math.exp(-u*h_nm/ell_nm)
        second=math.exp(-(1-u)*h_nm/ell_nm)
        w=ws*(first+second)
        dh=-ws/ell_nm*(u*first+(1-u)*second)
        coeff=1 if i in (0,n) else (4 if i%2 else 2)
        e=coeff*math.exp(w-wmax)
        weight_sum+=e;derivative_sum+=e*dh
    log_internal=math.log(h_nm*weight_sum/(3*n))+wmax
    frac_internal=1/(1+math.exp(math.log(math.pi*a_nm/2)-log_internal))
    return -frac_internal*(1+h_nm*derivative_sum/weight_sum)


def neutral_profile(h_nm, ws, ell_nm, c_left, c_right, n=400):
    """Exact equal-D, uncharged 1D profile with chemical-potential matching.

    Values at the membrane sides are c_bath*exp(-w_side), not c_bath.
    Integrals are cumulative trapezoids for the displayed profile only.
    """
    dz = h_nm/n
    zz = [i*dz for i in range(n+1)]
    ww = [field(z,h_nm,ws,ell_nm) for z in zz]
    shift = max(ww)
    ee = [math.exp(w-shift) for w in ww]
    ccum = [0.0]
    for i in range(n):
        ccum.append(ccum[-1]+dz*(ee[i]+ee[i+1])/2)
    return [dict(z_nm=z,w_kBT=w,
                 concentration=math.exp(-w)*(c_left-(c_left-c_right)*cum/ccum[-1]))
            for z,w,cum in zip(zz,ww,ccum)]


def normalized_finite_bath(time_ratio, theta):
    """C_R(t)/C_R(T), theta=k*T. Constant P, equal well-mixed reservoirs."""
    if theta < 0 or time_ratio < 0:
        raise ValueError('Nonnegative time and theta required')
    if theta < 1e-10:
        return time_ratio
    return -math.expm1(-theta*time_ratio)/(-math.expm1(-theta))


def finite_bath(t_s,p_m_s,area_m2,v_left_m3,v_right_m3,c_left0,c_right0):
    """Mass-conserving constant-P finite reservoirs, concentrations mol/m3."""
    if min(p_m_s,area_m2,v_left_m3,v_right_m3) <= 0 or t_s < 0:
        raise ValueError('Positive model inputs required')
    vtot=v_left_m3+v_right_m3
    ceq=(v_left_m3*c_left0+v_right_m3*c_right0)/vtot
    delta=(c_left0-c_right0)*math.exp(-area_m2*p_m_s*(1/v_left_m3+1/v_right_m3)*t_s)
    return ceq+v_right_m3/vtot*delta,ceq-v_left_m3/vtot*delta
