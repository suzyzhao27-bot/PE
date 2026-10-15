"""Equilibrium outside a solid fibre: nonlinear PB + optional competitive adsorption.

Water is a<r<b, never the inside of a pore. Units SI except mesh in nm.
Cell boundary du/dr(b)=0; u is relative to chemical reservoir, not set zero at b.
Bound adsorption charge and all diffuse mobile ions balance. No electron-transfer model.
"""
from dataclasses import dataclass
import math
import numpy as np
from scipy.linalg import solve_banded
from scipy.optimize import brentq

F=96485.33212; R=8.31446261815324; EPS0=8.8541878128e-12; E=1.602176634e-19
D={'Na':1.33e-9,'Cl':2.03e-9,'H':9.31e-9,'OH':5.27e-9}
Z={'Na':1.,'Cl':-1.,'H':1.,'OH':-1.}


def reservoir(salt_mM, pH=7.):
    # HCl or NaOH establishes pH. Assumed ideal Kw=1e-14 at 298.15K.
    h=10**(-pH);oh=10**(pH-14);c=salt_mM/1000
    return {'Na':1000*(c+max(oh-h,0)), 'Cl':1000*(c+max(h-oh,0)),
            'H':1000*h,'OH':1000*oh}

@dataclass
class Adsorption:
    density_nm2:float=0.05
    # Standard association free energies relative to 1 M; these are hypotheses.
    delta_g_kj:dict=None
    def state(self,u_surface, bulk,T=298.15):
        if self.delta_g_kj is None:return 0.,0.,{}
        names=list(self.delta_g_kj)
        logs=np.array([-self.delta_g_kj[n]*1000/(R*T)+math.log(bulk[n]/1000)-Z[n]*u_surface for n in names])
        # Include neutral empty-site state with log-weight0.
        maxlog=max(0.,float(logs.max()));weights=np.exp(logs-maxlog);empty=math.exp(-maxlog)
        probs=weights/(empty+weights.sum());zs=np.array([Z[n] for n in names])
        mean=float(probs@zs);second=float(probs@(zs*zs));scale=E*self.density_nm2*1e18
        return scale*mean,scale*(mean*mean-second),dict(zip(names,map(float,probs)))


@dataclass
class AcidSites:
    density_nm2:float=0.1
    pKa:float=7.
    def state(self,u_surface,bulk,T=298.15):
        pH=-math.log10(bulk['H']/1000)
        logq=math.log(10)*(self.pKa-pH)-u_surface
        alpha=1/(1+math.exp(min(700,max(-700,logq))))
        scale=E*self.density_nm2*1e18
        return -scale*alpha,-scale*alpha*(1-alpha),{'deprotonated':alpha}


def solve(salt_mM=1., pH=7.,radius_nm=13.3,cell_nm=39.5,n=512,
          sigma=None,adsorption=None,T=298.15,eps_r=78.3,tol=2e-11,outer='cell'):
    if not 0<radius_nm<cell_nm:raise ValueError('0<a<b')
    if (sigma is None)==(adsorption is None):raise ValueError('Choose fixed sigma OR adsorption')
    if outer not in ('cell','bulk'):raise ValueError('cell or bulk')
    bulk=reservoir(salt_mM,pH);names=list(bulk);cs=np.array([bulk[x] for x in names]);zs=np.array([Z[x] for x in names])
    edges=np.linspace(radius_nm,cell_nm,n+1);rs=(edges[:-1]+edges[1:])/2;dr=edges[1]-edges[0]
    vol=(edges[1:]**2-edges[:-1]**2)/2;conduct=edges[1:-1]/dr
    eps=EPS0*eps_r;chi=F*F/(eps*R*T)*1e-18;wallfac=radius_nm*F/(eps*R*T)*1e-9
    beta=0.5*dr*F/(eps*R*T)*1e-9
    diag=np.r_[conduct,0.]+np.r_[0.,conduct];lower=-conduct;upper=-conduct
    outercoef=cell_nm/(dr/2) if outer=='bulk' else 0.
    diag[-1]+=outercoef
    def boundary(u0):
        if adsorption is None:return sigma,0.,u0+beta*sigma,{}
        maxs=E*adsorption.density_nm2*1e18
        def eq(s):return s-adsorption.state(u0+beta*s,bulk,T)[0]
        s=brentq(eq,-maxs,maxs,xtol=1e-15);us=u0+beta*s
        _,ds,theta=adsorption.state(us,bulk,T)
        return s,ds/(1-beta*ds),us,theta
    def evaluate(u,jac=False):
        conc=cs[:,None]*np.exp(np.clip(-zs[:,None]*u[None,:],-80,80))
        rho=zs@conc;screen=(zs*zs)@conc
        sig,ds,us,theta=boundary(float(u[0]))
        res=diag*u;res[1:]+=lower*u[:-1];res[:-1]+=upper*u[1:]
        res-=chi*vol*rho;res[0]-=wallfac*sig
        if not jac:return res
        jd=diag+chi*vol*screen;jd[0]-=wallfac*ds
        return res,jd,conc,sig,us,theta
    u=np.zeros(n);history=[]
    for it in range(150):
        res,jd,conc,sig,us,theta=evaluate(u,True)
        scale=max(1.,abs(wallfac*sig),float(np.max(chi*vol*conc.sum(axis=0))))
        norm=float(np.max(abs(res))/scale);history.append(norm)
        if norm<tol:break
        ab=np.zeros((3,n));ab[0,1:]=upper;ab[1]=jd;ab[2,:-1]=lower
        step=solve_banded((1,1),ab,-res)
        if np.max(abs(step))>2:step*=2/np.max(abs(step))
        alpha=1.;old=float(np.max(abs(res)))
        while alpha>1e-8:
            trial=u+alpha*step
            if np.max(abs(evaluate(trial)))<old:break
            alpha*=.5
        if alpha<=1e-8:raise RuntimeError(f'Newton stalled {norm}')
        u=trial
    else:raise RuntimeError('No convergence')
    res,jd,conc,sig,us,theta=evaluate(u,True)
    mean=conc@vol/vol.sum();charge_water=F*float((zs@conc)@vol)*1e-18
    charge_surface=sig*radius_nm*1e-9
    # Per 2*pi*one metre fibre length. Outer electric-flux charge closes truncated bulk cells.
    outer_charge=eps*R*T/F*outercoef*u[-1] if outer=='bulk' else 0.
    gauss=charge_water+charge_surface-outer_charge
    area_water_nm2=2*vol.sum();area_total_nm2=cell_nm**2
    g=np.array([D[k] for k in names])*mean
    cation=float(g[zs>0].sum()/g.sum())
    sodium_proxy=float(g[0]/(g[0]+g[1]))
    result={'salt_mM':salt_mM,'pH':pH,'radius_nm_assumed':radius_nm,'cell_radius_nm_assumed':cell_nm,
            'water_fraction_cell':1-radius_nm**2/cell_nm**2,'outer_boundary':outer,'cells':n,
            'surface_sigma_C_m2':sig,'surface_u':us,'surface_psi_mV':us*R*T/F*1000,
            'cell_edge_psi_mV':u[-1]*R*T/F*1000,'slip_0p5nm_psi_mV':float(np.interp(radius_nm+.5,rs,u))*R*T/F*1000,
            'mobile_mean_mol_m3':dict(zip(names,map(float,mean))),
            'mobile_partition':dict(zip(names,map(float,mean/cs))),
            'linear_axial_current_fraction_cations':cation,'Na_fraction_NaCl_only_proxy':sodium_proxy,
            'relative_Na_Cl_partition_ratio':float((mean[0]/cs[0])/(mean[1]/cs[1])),
            'same_reservoir_uncharged_cation_fraction':float(sum(D[k]*bulk[k] for k in names if Z[k]>0)/sum(D[k]*bulk[k] for k in names)),
            'surface_adsorption_occupancies':theta,'surface_bound_charge_per_2pi_length_C_m':charge_surface,
            'mobile_charge_per_2pi_length_C_m':charge_water,'outer_electric_flux_C_m':outer_charge,
            'relative_Gauss_error':abs(gauss)/max(abs(charge_surface),1e-28),
            'Newton_iterations':it,'residual_scaled':history[-1],
            'scope':'hypothetical charged/charge-regulating exterior-fibre cell; axial linear conduction only, no finite-film/percolation/hydrodynamics/charge-origin proof',
            'radius_nm':rs,'u_profile':u,'concentration_profiles_mol_m3':conc,'species':names}
    return result

if __name__=='__main__':
    v=solve(adsorption=Adsorption(.05,{'H':-20,'OH':-50}))
    print({k:v[k] for k in ('surface_sigma_C_m2','surface_psi_mV','cell_edge_psi_mV','linear_axial_current_fraction_cations','relative_Gauss_error')})
