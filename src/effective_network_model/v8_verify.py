"""Eight frozen technical checks for the v8.0.1 implementation."""

import csv
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.linalg import expm

from v7_network_geometry import generate, morphology
from v8_prescreen import DB, RS, conductance
from v8_transport import baths, steady, weighted_graph

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / "原始资料/赵珈绮/16_物理诊断与接续任务_v8.0.1"
if not PACKAGE.exists():
    PACKAGE = ROOT / "taskbook"
OUT = ROOT / "WP3"


def add(rows, number, measured, limit, detail):
    status = "PASS_TECHNICAL" if measured <= limit else "NUMERICALLY_UNRESOLVED"
    rows.append({"check_id": number, "status": status, "measured": measured,
                 "limit": limit, "detail": detail})
    if status != "PASS_TECHNICAL":
        raise AssertionError(rows[-1])


def main():
    from sys import path
    path.insert(0, str(ROOT))
    from src.data_ingest.v8_step1_units import concentration, conductivity
    rows = []

    # 1. Conversion at all seven frozen points, including the branch join.
    samples = [0, 0.08, 1, 10, 100, 1363.25, 6872.54]
    errors = [abs(conductivity(concentration(x)) - x) / x if x else
              abs(conductivity(concentration(x))) for x in samples]
    add(rows, 1, max(errors), 1e-8, "7 frozen uS/cm points; zero absolute <=1e-12; 1 uS/cm=1e-4 S/m")
    assert errors[0] <= 1e-12
    assert math.isclose(DB, 1.6070833333333333e-9, rel_tol=1e-12)

    # 2. Conditional local kernel and series port resistance.
    lambdas = [0, 0.001, 0.1, 0.5, 0.9, 0.999999]
    kernels = [(1 - z) ** 2 * (1 - 2.104*z + 2.09*z**3 - 0.95*z**5) for z in lambdas]
    assert all(0 <= k <= 1 for k in kernels)
    assert conductance(RS, 1e-8, False) == 0
    assert conductance(RS / 2, 1e-8, True) == 0
    a = 2e-9
    length = 5e-9
    k = kernels[2]  # separate explicit value below, not a fitted quantity
    z = RS/a
    k = (1-z)**2 * (1-2.104*z+2.09*z**3-0.95*z**5)
    expected_port = 1 / (length/(DB*math.pi*a*a*k)+1/(4*DB*(a-RS)))
    add(rows, 2, abs(conductance(a,length,True)/expected_port-1), 1e-12,
        "K_R in [0,1]; a<=r_s closes edge; port has half-throat plus one access")

    # 3. One-layer, one-node graph has exactly two half-layer ports.
    h = 10e-9
    simple = {"h_nm": 10, "nodes": [{"r_m": 2e-9}],
              "edges": [{"i":"L","j":0,"L_m":h/2,"port_flag":True},
                        {"i":0,"j":"R","L_m":h/2,"port_flag":True}]}
    geom = {"s_nm":100,"nz":1,"epsilon":0.3,"nx":1,"ny":1}
    values = []
    for model, beta in (("Tcal", 0), ("Scal", 0.5)):
        weighted, area = weighted_graph(simple, geom, 1.0, model, beta)
        result = steady(1, weighted, area)
        a = 2e-9
        if model == "Tcal":
            z = RS/a
            kval = (1-z)**2*(1-2.104*z+2.09*z**3-0.95*z**5)
            a_access = a-RS
        else:
            kval = math.exp(-beta*25e-9/a)
            a_access = a
        analytical = 1/(h/(DB*math.pi*a*a*kval)+2/(4*DB*a_access))
        values.append(abs(result["G_m3s"]/analytical-1))
    add(rows, 3, max(values), 1e-8, "nz=1: both h/2 throats and both one-side access terms")

    # 4. Two parallel edges and two serial edges.
    parallel = steady(0, [(0,1,2.0),(0,1,3.0)], 1)
    serial = steady(1, [(1,0,2.0),(0,2,3.0)], 1)
    add(rows, 4, max(abs(parallel["G_m3s"]/5-1), abs(serial["G_m3s"]/1.2-1)),
        1e-8, "parallel 2+3=5; series (1/2+1/3)^-1=1.2")

    # 5. Ion-disconnected geometry and detached branch cannot pollute the solve.
    disconnected = steady(2, [(2,0,2.0),(1,3,3.0)], 1)
    branch = steady(2, [(2,0,2.0),(0,3,3.0),(0,1,7.0)], 1)
    add(rows, 5, max(disconnected["G_m3s"], abs(branch["G_m3s"]/1.2-1)),
        1e-8, "ion-disconnected zero; dangling branch no singular pollution")

    # 6. Three-bath conservation, nonnegative state and tau=0 analytic limit.
    times = np.array([0, 20, 60, 120]) * 60.0
    cases = [baths(1e-8, 10, tau, times) for tau in (0, 300, 1200)]
    analytic = 600 * (-np.expm1(-2*10*1e-8*times)) / 2
    error = max([x["mass_relative_error"] for x in cases] +
                [float(np.max(np.abs(cases[0]["receiver_c_molm3"] - analytic))/600)])
    assert min(x["minimum_u"] for x in cases) >= -1e-12
    add(rows, 6, error, 1e-8, "tau=0/300/1200 s; conserved D+phi near+(1-phi) bulk")

    # 7. Nonuniform two-node storage identity, independent quadrature of flux deficit.
    g1, g2, g3 = 2.0, 1.3, 3.0
    lap = np.array([[g1+g2,-g2],[-g2,g2+g3]])
    m = np.diag([0.7,1.1])
    v = np.linalg.solve(lap, np.array([g1,0.0]))
    G = g3*v[1]
    predicted_lag = (0.7*v[0]*(1-v[0])+1.1*v[1]*(1-v[1]))/G
    decay = -np.linalg.solve(m,lap)
    integrated_lag = quad(lambda t:g3*(expm(decay*t)@v)[1]/G, 0, 30,
                          epsabs=1e-11, epsrel=1e-11)[0]
    add(rows, 7, abs(integrated_lag/predicted_lag-1), 0.01,
        "independent flux-deficit integral vs sum M_i v_i(1-v_i)/G")

    # 8. One main candidate graph: refined linear tolerance and time integration.
    spec = json.loads((PACKAGE / "inputs/legacy_v7_spec.json").read_text(encoding="utf-8"))
    spec["geometry"] = dict(spec["geometry"])
    spec["geometry"].update(nx=12, ny=12, p_z=0.65, sigma_ln=0.35, q=0.5)
    hs = [65,120,200,350,700,2600]
    max_nz = max(math.ceil(h / sum(morphology(h,"S",spec["geometry"]["anchors"])[:2])) for h in hs)
    graph, meta = generate(350,"S",7001,spec,max_nz)
    ginfo = {"s_nm":meta["s_nm"],"nz":meta["nz"],"epsilon":meta["epsilon"],"nx":12,"ny":12}
    weighted, area = weighted_graph(graph,ginfo,0.5,"Tcal")
    direct = steady(len(graph["nodes"]),weighted,area)
    refined = steady(len(graph["nodes"]),weighted,area,rtol=1e-11,solver="cg")
    loose = baths(direct["P_m_s"],100,300,times)
    precise = baths(refined["P_m_s"],100,300,times)
    delta_direct = np.array([conductivity(c) for c in loose["receiver_c_molm3"]])
    delta_refined = np.array([conductivity(c) for c in precise["receiver_c_molm3"]])
    # Independent high-accuracy solve_ivp on the same fixed 3-bath equations.
    rate = 100*direct["P_m_s"]
    phi, tau = 0.05, 300.0
    matrix = np.array([[-rate,rate,0],[rate/phi,-rate/phi-(1-phi)/tau,(1-phi)/tau],
                       [0,phi/tau,-phi/tau]])
    ivp = solve_ivp(lambda _t,u:matrix@u,(0,float(times[-1])),[1,0,0],
                    t_eval=times,method="BDF",jac=matrix,rtol=1e-10,atol=1e-13)
    assert ivp.success
    delta_ivp = np.array([conductivity(max(0,600*u)) for u in ivp.y[2]])
    difference = max(float(np.max(np.abs(delta_direct-delta_refined))),
                     float(np.max(np.abs(delta_direct-delta_ivp))))
    tolerance = max(0.1, 0.01*float(np.max(delta_direct)))
    add(rows, 8, difference, tolerance,
        f"350 nm seed7001 q=0.5; direct vs CG 1e-11 and matrix-exponential vs BDF 1e-10; max delta={float(np.max(delta_direct)):.9g} uS/cm")

    probe = {"model":"Tcal","h_nm":350,"seed":7001,"q":0.5,"r_s_nm":0.35,
             "alpha_m_inv":100,"tau_s":300,"C0_molm3":600,"times_s":times.tolist(),
             "G_m3s":direct["G_m3s"],"P_m_s":direct["P_m_s"],
             "spanning_nodes":direct["spanning_nodes"],
             "linear_residual_relative":direct["residual_relative"],
             "conservation_relative":direct["conservation_relative"],
             "mass_relative_error":loose["mass_relative_error"],
             "minimum_activity":loose["minimum_u"],
             "delta_uScm":delta_direct.tolist(),
             "refinement_max_abs_uScm":difference,
             "observation_status":"synthetic model probe; not fitted to B350"}
    (OUT / "v8_技术探针.json").write_text(json.dumps(probe,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")

    with (OUT / "v8_技术检查.csv").open("w",encoding="utf-8-sig",newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"v8 technical checks: {len(rows)}/8 PASS; max probe delta={float(np.max(delta_direct)):.6g} uS/cm")


if __name__ == "__main__":
    main()
