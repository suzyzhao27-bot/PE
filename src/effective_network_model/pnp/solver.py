"""Finite 2D aqueous Poisson-Nernst-Planck, periodic transverse fiber cell.

Water domain includes both entrances and finite external aqueous regions.
Square impermeable fiber cross-section is an explicitly representative local
layer, not the full 230nm or real 3D PE network. Local fluid velocity = 0.
NumPy/SciPy are existing dependencies; nothing is installed by this module.
"""
import math
from dataclasses import dataclass, asdict
import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import spsolve

F = 96485.33212
RG = 8.31446261815324
EPS0 = 8.8541878128e-12
E_CHARGE = 1.602176634e-19
SCALE_M = 1e-9
CSCALE = 1.0  # 1 mol/m3 = 1 mM


@dataclass(frozen=True)
class CommonIonEnergy:
    """Prescribed shared Na/Cl partition energy w=U/(k_B T).

    x is the axial coordinate (the array's first index, called z in the v3
    implementation). The default membrane boundaries are the outer faces of
    the first/last fiber rows. Energy is assigned only to aqueous cell centers
    in that span; both finite external baths have w=0. This is a phenomenological
    structure/interface entry energy, not a computed dielectric self-energy.
    """
    kind: str = "rectangle"
    amplitude_kbt: float = 0.
    decay_nm: float = 10.
    x_left_nm: float | None = None
    x_right_nm: float | None = None

    def build(self, mesh):
        g = mesh.geometry
        left = g.external_fiber_sides_each*g.fiber_side_nm if self.x_left_nm is None else float(self.x_left_nm)
        right = left+(2*g.fiber_rows-1)*g.fiber_side_nm if self.x_right_nm is None else float(self.x_right_nm)
        if self.kind not in ("rectangle", "two_sided_exponential"):
            raise ValueError("energy kind must be rectangle or two_sided_exponential")
        if not np.isfinite([left,right,self.amplitude_kbt,self.decay_nm]).all() or not 0 < left < right < mesh.height_nm:
            raise ValueError("finite energy parameters and strictly internal membrane boundaries required")
        if self.kind == "two_sided_exponential" and self.decay_nm <= 0:
            raise ValueError("positive decay_nm required")
        x = (np.nonzero(mesh.mask)[0]+.5)*mesh.h
        inside = (x >= left) & (x < right)
        w = np.zeros(mesh.n)
        if self.kind == "rectangle":
            w[inside] = self.amplitude_kbt
        else:
            w[inside] = self.amplitude_kbt*(np.exp(-(x[inside]-left)/self.decay_nm)
                                          + np.exp(-(right-x[inside])/self.decay_nm))
        return w, dict(**asdict(self), resolved_x_left_nm=left,
                       resolved_x_right_nm=right, membrane_span_nm=right-left,
                       interpretation="prescribed effective shared structure/interface entry energy U/kBT; not microscopic dielectric self-energy",
                       discontinuity="outside bath energy is zero; same SG electrochemical-potential difference is used across every face")


def bernoulli(v):
    v = np.asarray(v, dtype=float)
    b, db = np.empty_like(v), np.empty_like(v)
    small, high, low = np.abs(v) < 1e-3, v > 50, v < -50
    middle = ~(small | high | low)
    q = v[small]
    b[small] = 1-q/2+q*q/12-q**4/720+q**6/30240
    db[small] = -.5+q/6-q**3/180+q**5/5040
    q = v[high]
    b[high] = q*np.exp(-q)
    db[high] = (1-q)*np.exp(-q)
    b[low], db[low] = -v[low], -1.
    q = v[middle]
    em = np.expm1(q)
    b[middle] = q/em
    db[middle] = (em-q*(em+1))/(em*em)
    return b, db


@dataclass(frozen=True)
class Geometry:
    fiber_side_nm: float = 26.6
    cells_per_fiber: int = 8
    period_fiber_sides: int = 3
    external_fiber_sides_each: int = 2
    fiber_present: bool = True
    fiber_rows: int = 1
    staggered_rows: bool = True

    def build(self):
        m = self.cells_per_fiber
        if m < 2 or m % 2 or self.external_fiber_sides_each < 1 or self.period_fiber_sides < 2 or self.fiber_rows<1:
            raise ValueError("even cells_per_fiber>=2; positive external water; period>=2 sides")
        nx = self.period_fiber_sides*m
        ny = (2*self.fiber_rows-1+2*self.external_fiber_sides_each)*m
        h = self.fiber_side_nm/m
        mask = np.ones((ny,nx), dtype=bool)
        if self.fiber_present:
            x0=(nx-m)//2
            for layer in range(self.fiber_rows):
                z0=(self.external_fiber_sides_each+2*layer)*m
                shift=nx//2 if self.staggered_rows and layer%2 else 0
                columns=(np.arange(x0,x0+m)+shift)%nx
                mask[z0:z0+m,columns] = False
        ids = np.full((ny,nx), -1, dtype=int)
        ids[mask] = np.arange(mask.sum())
        edge_i, edge_j = [], []
        wall_faces = np.zeros(mask.sum())
        for z,x in zip(*np.nonzero(mask)):
            i = ids[z,x]
            # Once per periodic undirected x-face and once per axial z-face.
            xr = (x+1)%nx
            if mask[z,xr]:
                edge_i.append(i);edge_j.append(ids[z,xr])
            if z+1 < ny and mask[z+1,x]:
                edge_i.append(i);edge_j.append(ids[z+1,x])
            for zz,xx in ((z,(x-1)%nx),(z,xr),(z-1,x),(z+1,x)):
                if 0 <= zz < ny and not mask[zz,xx]:
                    wall_faces[i] += 1
        return Mesh(self, nx,ny,h,mask,ids,np.array(edge_i),np.array(edge_j),wall_faces)


class Mesh:
    def __init__(self,g,nx,ny,h,mask,ids,ei,ej,wall_faces):
        self.geometry,self.nx,self.ny,self.h = g,nx,ny,h
        self.mask,self.ids,self.ei,self.ej,self.wall_faces = mask,ids,ei,ej,wall_faces
        self.n = int(mask.sum())
        self.left = ids[0,:]
        self.right = ids[-1,:]
        self.bnodes = np.r_[self.left,self.right]
        self.bsides = np.r_[np.zeros(nx,dtype=int),np.ones(nx,dtype=int)]
        n=self.n
        diagonal=np.bincount(np.r_[ei,ej],minlength=n).astype(float)
        np.add.at(diagonal,self.bnodes,2.)
        self.L = coo_matrix((np.r_[diagonal,-np.ones(len(ei)),-np.ones(len(ei))],
                             (np.r_[np.arange(n),ei,ej],np.r_[np.arange(n),ej,ei])),shape=(n,n)).tocsr()
        self.width_m = nx*h*SCALE_M
        self.height_nm = ny*h


def initial_guess(mesh,c_left,c_right,right_u,open_circuit):
    z,x = np.nonzero(mesh.mask)
    fraction=(z+.5)/mesh.ny
    c=c_left+(c_right-c_left)*fraction
    u=right_u*fraction
    y=np.r_[u,np.log(c),np.log(c)]
    return np.r_[y,right_u] if open_circuit else y


def adsorption_wall(mesh,cp,cm,params,T,eps_r):
    """Hypothetical single-site Na/Cl adsorption, chemical Gauss/Robin wall.

    Gamma [sites/nm2], K [1/M]. No H/OH chemistry or electronic charge source.
    Wall concentration uses zero normal electrochemical flux over half a cell.
    Its implicit sigma correction is solved locally and differentiated exactly.
    """
    gamma=float(params["site_density_nm2"])
    kp=float(params.get("K_Na_per_M",0))/1000
    km=float(params.get("K_Cl_per_M",0))/1000
    if gamma<0 or min(kp,km)<0:
        raise ValueError("nonnegative site density/adsorption constants")
    cap=E_CHARGE*gamma*1e18
    distance=mesh.h*SCALE_M/2
    q=distance/(EPS0*eps_r*RG*T/F)
    sigma=cap*(kp*cp-km*cm)/(1+kp*cp+km*cm)
    for _ in range(20):
        a=kp*cp*np.exp(-sigma*q)
        b=km*cm*np.exp(sigma*q)
        den=1+a+b
        residual=sigma-cap*(a-b)/den
        slope=1+cap*q*(a+b+4*a*b)/(den*den)
        delta=residual/slope
        sigma-=delta
        if float(np.max(abs(delta)))<1e-14:
            break
    else:
        raise RuntimeError("local adsorption charge root failed")
    a=kp*cp*np.exp(-sigma*q);b=km*cm*np.exp(sigma*q);den=1+a+b
    slope=1+cap*q*(a+b+4*a*b)/(den*den)
    dsplus=cap*a*(1+2*b)/(den*den)/slope
    dsminus=-cap*b*(1+2*a)/(den*den)/slope
    return sigma,dsplus,dsminus


def assemble(mesh,y,c_left,c_right,sigma,dp,dm,T,eps_r,open_circuit,fixed_right_u,jacobian=True,surface_model=None,common_energy=None):
    n=mesh.n
    u=y[:n];cp=np.exp(y[n:2*n]);cm=np.exp(y[2*n:3*n])
    w = np.zeros(n) if common_energy is None else np.asarray(common_energy,dtype=float)
    if w.shape != (n,) or not np.isfinite(w).all():
        raise ValueError("common_energy must contain one finite dimensionless value per water cell")
    right_u=float(y[-1]) if open_circuit else fixed_right_u
    beta=F*F*CSCALE*SCALE_M*SCALE_M/(EPS0*eps_r*RG*T)
    wall_coefficient=SCALE_M*F/(EPS0*eps_r*RG*T)*mesh.h*mesh.wall_faces
    if surface_model is None:
        sigma_local=np.full(n,sigma);dsplus=dsminus=np.zeros(n)
    else:
        sigma_local,dsplus,dsminus=adsorption_wall(mesh,cp,cm,surface_model,T,eps_r)
    wall=sigma_local*wall_coefficient
    boundary_u=right_u*mesh.bsides
    dir_rhs=np.zeros(n)
    np.add.at(dir_rhs,mesh.bnodes,2.*boundary_u)
    residual=np.zeros(3*n+int(open_circuit))
    residual[:n]=mesh.L@u-beta*mesh.h**2*(cp-cm)-wall-dir_rhs
    jr,jc,jv=[],[],[]
    def add(rows,cols,vals):
        jr.extend(np.asarray(rows,dtype=int).ravel().tolist())
        jc.extend(np.asarray(cols,dtype=int).ravel().tolist())
        jv.extend(np.asarray(vals,dtype=float).ravel().tolist())
    if jacobian:
        l=mesh.L.tocoo()
        add(l.row,l.col,l.data)
        add(np.arange(n),n+np.arange(n),-beta*mesh.h**2*cp-wall_coefficient*dsplus)
        add(np.arange(n),2*n+np.arange(n),beta*mesh.h**2*cm-wall_coefficient*dsminus)
        if open_circuit:
            add(mesh.right,np.full(mesh.nx,3*n),-2*np.ones(mesh.nx))
    ei,ej=mesh.ei,mesh.ej
    boundary_c=np.where(mesh.bsides==0,c_left,c_right)
    boundary_fluxes=[]
    for z,c,offset,dr in ((1,cp,n,dp/1e-9),(-1,cm,2*n,dm/1e-9)):
        # SG is fitted to mu_i/kBT=log(c_i)+z_i*u+w_i. Using the same
        # potential difference at the membrane/bath jump preserves Boltzmann
        # partitioning and zero equilibrium flux without an extra jump source.
        v=z*(u[ej]-u[ei])+w[ej]-w[ei]
        bv,dbv=bernoulli(v);bm,dbm=bernoulli(-v)
        flux=bv*c[ei]-bm*c[ej]
        q=dbv*c[ei]+dbm*c[ej]
        np.add.at(residual,offset+ei,flux)
        np.add.at(residual,offset+ej,-flux)
        if jacobian:
            add(offset+ei,offset+ei,bv*c[ei]);add(offset+ei,offset+ej,-bm*c[ej])
            add(offset+ej,offset+ei,-bv*c[ei]);add(offset+ej,offset+ej,bm*c[ej])
            add(offset+ei,ei,-z*q);add(offset+ei,ej,z*q)
            add(offset+ej,ei,z*q);add(offset+ej,ej,-z*q)
        bn=mesh.bnodes
        # Reservoir boundary partition energy is zero.
        v=z*(boundary_u-u[bn])-w[bn]
        bv,dbv=bernoulli(v);bm,dbm=bernoulli(-v)
        bflux=2*(bv*c[bn]-bm*boundary_c)
        bq=2*(dbv*c[bn]+dbm*boundary_c)
        np.add.at(residual,offset+bn,bflux)
        boundary_fluxes.append(bflux)
        if jacobian:
            add(offset+bn,offset+bn,2*bv*c[bn]);add(offset+bn,bn,-z*bq)
            if open_circuit:
                right_mask=mesh.bsides==1
                right_nodes=bn[right_mask]
                add(offset+right_nodes,np.full(mesh.nx,3*n),z*bq[right_mask])
                sign=1 if z==1 else -1
                add(np.full(mesh.nx,3*n),offset+right_nodes,sign*dr*2*bv[right_mask]*c[right_nodes])
                add(np.full(mesh.nx,3*n),right_nodes,-sign*dr*z*bq[right_mask])
                add(np.full(mesh.nx,3*n),np.full(mesh.nx,3*n),sign*dr*z*bq[right_mask])
    if open_circuit:
        residual[-1]=(dp/1e-9)*sum(boundary_fluxes[0][mesh.nx:])-(dm/1e-9)*sum(boundary_fluxes[1][mesh.nx:])
    # Fixed per-current-iterate scaling for Newton/line search; dimensions are documented.
    scale=np.r_[1+np.abs(wall)+beta*mesh.h**2*(cp+cm),1+cp,1+cm]
    if open_circuit:
        scale=np.r_[scale,1+max(c_left,c_right)*mesh.nx/mesh.ny*max(dp,dm)/1e-9]
    matrix=coo_matrix((jv,(jr,jc)),shape=(len(y),len(y))).tocsr() if jacobian else None
    extras=dict(cp=cp,cm=cm,u=u,right_u=right_u,wall_rhs=wall,beta=beta,boundary_fluxes=boundary_fluxes,sigma_local=sigma_local,common_energy_kbt=w)
    return residual,matrix,scale,extras


def solve(c_left=30.,c_right=1.,sigma=0.,geometry=Geometry(),dp=1.33e-9,dm=2.03e-9,
          temperature=298.15,eps_r=78.3,open_circuit=True,right_voltage_v=0.,
          initial=None,tol=2e-10,max_iter=80,surface_model=None,ion_energy_model=None):
    if min(c_left,c_right,dp,dm,temperature,eps_r)<=0:
        raise ValueError("Positive known aqueous concentrations/parameters; no arbitrary DI regularization")
    mesh=geometry.build()
    if ion_energy_model is None:
        energy = np.zeros(mesh.n)
        energy_metadata = None
    else:
        if isinstance(ion_energy_model,dict):
            ion_energy_model = CommonIonEnergy(**ion_energy_model)
        energy,energy_metadata = ion_energy_model.build(mesh)
    if surface_model is not None and np.any(energy):
        raise ValueError("shared entry-energy extension is validated with fixed sigma only; do not combine with the unvalidated v3 adsorption chemistry")
    vt=RG*temperature/F
    h=(dm-dp)/(dm+dp)
    right_u=h*math.log(c_right/c_left) if open_circuit else right_voltage_v/vt
    y=initial.copy() if initial is not None else initial_guess(mesh,c_left,c_right,right_u,open_circuit)
    if initial is None:
        y[mesh.n:2*mesh.n] -= energy
        y[2*mesh.n:3*mesh.n] -= energy
    if len(y)!=3*mesh.n+int(open_circuit):
        raise ValueError("Initial vector belongs to a different mesh or voltage mode")
    history=[]
    for iteration in range(max_iter+1):
        res,jac,scale,extra=assemble(mesh,y,c_left,c_right,sigma,dp,dm,temperature,eps_r,open_circuit,right_u,surface_model=surface_model,common_energy=energy)
        norm=float(np.max(np.abs(res)/scale))
        history.append(norm)
        if norm<tol:
            break
        if iteration==max_iter:
            raise RuntimeError(f"Newton failed: residual={norm:g}, sigma={sigma:g}")
        step=spsolve(diags(1/scale)@jac,-res/scale)
        if not np.all(np.isfinite(step)):
            raise RuntimeError("Nonfinite Newton step")
        alpha=min(1.,1.5/max(1.5,float(np.max(np.abs(step[mesh.n:3*mesh.n])))),
                  2./max(2.,float(np.max(np.abs(step[:mesh.n])))))
        accepted=False
        for _ in range(28):
            candidate=y+alpha*step
            if np.max(np.abs(candidate[mesh.n:3*mesh.n]))<100:
                rr,_,_,_=assemble(mesh,candidate,c_left,c_right,sigma,dp,dm,temperature,eps_r,open_circuit,right_u,False,surface_model=surface_model,common_energy=energy)
                nn=float(np.max(np.abs(rr)/scale))
                if nn<norm:
                    y=candidate;accepted=True;break
            alpha*=.5
        if not accepted:
            raise RuntimeError(f"Newton line search stalled at {norm:g}")
    plus,minus=extra["boundary_fluxes"]
    flux_plus_right=dp*CSCALE*float(sum(plus[mesh.nx:]))/mesh.width_m
    flux_minus_right=dm*CSCALE*float(sum(minus[mesh.nx:]))/mesh.width_m
    flux_plus_left=-dp*CSCALE*float(sum(plus[:mesh.nx]))/mesh.width_m
    flux_minus_left=-dm*CSCALE*float(sum(minus[:mesh.nx]))/mesh.width_m
    net_i=F*(flux_plus_right-flux_minus_right)
    boundary_field=float(np.sum(2*(extra["u"][mesh.bnodes]-extra["right_u"]*mesh.bsides)))
    charge_rhs=float(extra["beta"]*mesh.h**2*np.sum(extra["cp"]-extra["cm"])+np.sum(extra["wall_rhs"]))
    denom=max(abs(flux_plus_right),abs(flux_minus_right),1e-30)
    result=dict(c_left_mol_m3=c_left,c_right_mol_m3=c_right,sigma_c_m2=sigma,
                surface_boundary="fixed given sigma" if surface_model is None else "hypothetical Na/Cl single-site chemical Gauss/Robin",
                surface_model_parameters=surface_model,
                common_ion_energy_model=energy_metadata,
                common_ion_energy_min_kbt=float(np.min(energy)),
                common_ion_energy_max_kbt=float(np.max(energy)),
                wall_sigma_mean_c_m2=float(np.sum(extra["sigma_local"]*mesh.wall_faces)/np.sum(mesh.wall_faces)) if np.sum(mesh.wall_faces) else 0.,
                wall_sigma_min_c_m2=float(np.min(extra["sigma_local"][mesh.wall_faces>0])) if np.sum(mesh.wall_faces) else 0.,
                wall_sigma_max_c_m2=float(np.max(extra["sigma_local"][mesh.wall_faces>0])) if np.sum(mesh.wall_faces) else 0.,
                open_circuit=open_circuit,local_fluid_velocity_zero=True,
                right_minus_left_theoretical_aqueous_voltage_mv=extra["right_u"]*vt*1000,
                plus_flux_right_mol_m2_s=flux_plus_right,minus_flux_right_mol_m2_s=flux_minus_right,
                plus_flux_left_mol_m2_s=flux_plus_left,minus_flux_left_mol_m2_s=flux_minus_left,
                common_salt_flux_mol_m2_s=(flux_plus_right+flux_minus_right)/2,
                right_ionic_current_a_m2=net_i,
                relative_net_current_error=abs(flux_plus_right-flux_minus_right)/denom if open_circuit else None,
                plus_conservation_relative=abs(flux_plus_left-flux_plus_right)/denom,
                minus_conservation_relative=abs(flux_minus_left-flux_minus_right)/denom,
                poisson_gauss_balance_scaled=abs(boundary_field-charge_rhs)/(1+abs(boundary_field)+abs(charge_rhs)),
                max_scaled_equation_residual=norm,newton_iterations=iteration,
                water_cells=mesh.n,grid_nm=mesh.h,width_nm=mesh.width_m/SCALE_M,
                total_aqueous_domain_height_nm=mesh.height_nm,
                represented_local_fiber_layer_thickness_nm=geometry.fiber_side_nm,
                representative_fiber_rows=geometry.fiber_rows,
                representative_network_span_nm=(2*geometry.fiber_rows-1)*geometry.fiber_side_nm,
                nominal_open_gap_nm=(geometry.period_fiber_sides-1)*geometry.fiber_side_nm,
                external_water_each_nm=geometry.external_fiber_sides_each*geometry.fiber_side_nm,
                full_230nm_membrane_not_reconstructed=True,
                status="finite 2D representative local fiber-layer PNP; assumed charge; not PE measurements or 3D network validation")
    return result,y,extra,mesh,history


def linear_conductance(c=1.,sigma=0.,geometry=Geometry(),dp=1.33e-9,dm=2.03e-9,bias_v=1e-4,surface_model=None,ion_energy_model=None):
    equilibrium=solve(c,c,sigma,geometry,dp,dm,open_circuit=False,right_voltage_v=0.,surface_model=surface_model,ion_energy_model=ion_energy_model)
    pos=solve(c,c,sigma,geometry,dp,dm,open_circuit=False,right_voltage_v=bias_v,initial=equilibrium[1],surface_model=surface_model,ion_energy_model=ion_energy_model)
    neg=solve(c,c,sigma,geometry,dp,dm,open_circuit=False,right_voltage_v=-bias_v,initial=equilibrium[1],surface_model=surface_model,ion_energy_model=ion_energy_model)
    rp,rn=pos[0],neg[0]
    derivative_i=(rp["right_ionic_current_a_m2"]-rn["right_ionic_current_a_m2"])/(2*bias_v)
    derivative_jp=(rp["plus_flux_right_mol_m2_s"]-rn["plus_flux_right_mol_m2_s"])/(2*bias_v)
    derivative_jm=(rp["minus_flux_right_mol_m2_s"]-rn["minus_flux_right_mol_m2_s"])/(2*bias_v)
    tplus=derivative_jp/(derivative_jp-derivative_jm)
    return dict(c_bulk_mol_m3=c,sigma_c_m2=sigma,bias_v=bias_v,
                wall_sigma_mean_equilibrium_c_m2=equilibrium[0]["wall_sigma_mean_c_m2"],
                surface_boundary=equilibrium[0]["surface_boundary"],surface_model_parameters=surface_model,
                projected_area_ionic_conductance_s_m2=-derivative_i,
                linear_current_cation_fraction_t_plus=tplus,
                equilibrium_mobile_charge_mol_per_m_depth=CSCALE*SCALE_M**2*geometry.build().h**2*float(np.sum(equilibrium[2]["cp"]-equilibrium[2]["cm"])),
                static_charge_is_given_boundary_not_origin_prediction=surface_model is None,
                static_sigma_calculated_from_assumed_chemical_model=surface_model is not None,
                chemical_identity_and_constants_are_unvalidated=surface_model is not None,
                max_scaled_equation_residual=max(x[0]["max_scaled_equation_residual"] for x in (equilibrium,pos,neg)),
                scope="local u=0, equal known salt, finite representative water/fiber cell; conductance per nominal width and unit out-of-plane depth"),equilibrium
