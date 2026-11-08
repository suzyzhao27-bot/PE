"""Fit a shared phenomenological overlap field to six PE endpoint increments.

No source writes or installation. Finite-reservoir/dynamic effects are not fitted.
The common amplitude absorbs calibration, A/V, drive, time and D.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import sys

import numpy as np
import scipy
from scipy.optimize import minimize, minimize_scalar
from numpy.polynomial.legendre import leggauss

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--input',type=Path,default=HERE/'inputs'/'Book1__Sheet1.csv')
parser.add_argument('--output',type=Path,default=HERE/'results',
    help='Write only to an independently copied local working directory, never the formal release.')
args=parser.parse_args()
SOURCE=args.input.resolve()
OUT=args.output.resolve()
if any(part.lower() in ('onedrive','dropbox') for part in OUT.parts):
    raise ValueError('Copy the package to a non-sync local working directory before running this optional refit.')
SOURCE_HASH_BEFORE=hashlib.sha256(SOURCE.read_bytes()).hexdigest()
OUT.mkdir(parents=True,exist_ok=True)
A_NM = 25.0
ACCESS_NM = math.pi*A_NM/2
THICK = np.array([65.,120.,200.,350.,700.,2600.])
with SOURCE.open(encoding='utf-8-sig',newline='') as f:
    raw=list(csv.DictReader(f))
if len(raw)<2 or any(float(raw[0][x])!=0 or float(raw[-1][x])!=120 for x in ['col_1_A','col_3_C']):
    raise ValueError('This endpoint analysis requires observed 0 and 120 min in both explicit X columns.')
Y=np.array([float(raw[-1][c])-float(raw[0][c]) for c in ['col_2_B','col_4_D','col_5_E','col_6_F','col_7_G','col_8_H']])
if not np.all(np.isfinite(Y)) or not np.all(Y>0):
    raise ValueError('Log residual fitting requires six finite positive endpoint increments.')
LY=np.log(Y)
WS_BOUNDS=(0.01,100.)
ELL_BOUNDS=(0.1,10000.)
BOUNDS=[tuple(np.log(WS_BOUNDS)),tuple(np.log(ELL_BOUNDS))]
GX,GW=leggauss(48)


def log_resistance(h,ws,ell,access=ACCESS_NM,nodes=None,robin_r=0):
    """R/D in nm: h + 2ell integral_0^(h/2ell) expm1(w(x)) dx + access.

    Piecewise quadrature resolves the boundary width ell/ws. Tail x>32 is
    negligible for ws<=100 and h/ell>=64, and is checked independently below.
    """
    gx,gw=(GX,GW) if nodes is None else nodes
    h=np.atleast_1d(h)
    res=[]
    for hh in h:
        H=hh/ell
        xmax=min(H/2,32.)
        cuts=sorted(set([0.,xmax]+[v for v in [.01,.05,.2,.5,1.,2.,4.,8.,16.] if v<xmax]))
        val=0.
        for lo,hi in zip(cuts,cuts[1:]):
            x=lo+(gx+1)*(hi-lo)/2
            w=ws*(np.exp(-x)+np.exp(-(H-x)))/(1+robin_r*np.exp(-H))
            val+=(hi-lo)/2*np.dot(gw,np.expm1(w))
        res.append(math.log(hh+2*ell*val+access))
    return np.array(res)


def evaluate_shape(logr,train):
    logamp=float(np.mean(LY[train]+logr[train]))
    residual=logamp-logr-LY
    return logamp,residual


def fit_overlap(train,access=ACCESS_NM,global_scan=True,starts=None,robin_r=0):
    train=np.array(train,dtype=int)
    def objective(x):
        logr=log_resistance(THICK,*np.exp(x),access,robin_r=robin_r)
        _,r=evaluate_shape(logr,train)
        return float(np.dot(r[train],r[train]))
    grid=[]
    if global_scan:
        for ws in np.geomspace(*WS_BOUNDS,41):
            for ell in np.geomspace(*ELL_BOUNDS,61):
                x=np.log([ws,ell]);v=objective(x)
                grid.append((v,ws,ell))
        starts=[np.log([r[1],r[2]]) for r in sorted(grid)[:12]]
        starts += [np.log([1.,100.]),np.log([10.,1000.]),np.log([80.,8000.])]
    best=None
    for start in starts:
        opt=minimize(objective,start,method='L-BFGS-B',bounds=BOUNDS,
            options=dict(ftol=1e-13,gtol=1e-9,maxiter=400))
        if best is None or opt.fun<best.fun:best=opt
    ws,ell=np.exp(best.x)
    logr=log_resistance(THICK,ws,ell,access,robin_r=robin_r)
    logamp,r=evaluate_shape(logr,train)
    return dict(ws=float(ws),ell_nm=float(ell),amplitude=float(np.exp(logamp)),log_amplitude=logamp,robin_r=robin_r,
        train_indices=train.tolist(),train_SSE_log=float(np.dot(r[train],r[train])),
        log_resistance=logr.tolist(),predicted=np.exp(logamp-logr).tolist(),residual_log=r.tolist(),
        optimizer_success=bool(best.success),optimizer_message=str(best.message),
        _grid=grid)


def metrics(r):
    r=np.asarray(r)
    return dict(log_RMSE=float(np.sqrt(np.mean(r*r))),log_MAE=float(np.mean(abs(r))),
        geometric_mean_absolute_fold_error=float(np.exp(np.mean(abs(r)))),
        max_fold_error=float(np.exp(np.max(abs(r)))),SSE_log=float(np.dot(r,r)))


def fit_baseline(train,kind):
    train=np.asarray(train,dtype=int)
    w=None
    if kind=='D_over_h':lr=np.log(THICK)
    elif kind=='D_over_h_plus_fixed_access':lr=np.log(THICK+ACCESS_NM)
    elif kind=='shared_constant_positive_w':
        def obj(ww):
            lr=np.log(THICK*np.exp(ww)+ACCESS_NM)
            _,r=evaluate_shape(lr,train)
            return float(np.dot(r[train],r[train]))
        op=minimize_scalar(obj,bounds=(0,100),method='bounded',options=dict(xatol=1e-10))
        w=float(op.x);lr=np.log(THICK*np.exp(w)+ACCESS_NM)
    else:raise ValueError(kind)
    la,r=evaluate_shape(lr,train)
    return dict(model=kind,constant_w=w,amplitude=float(np.exp(la)),predicted=np.exp(la-lr).tolist(),
        residual_log=r.tolist(),metrics=metrics(r))


def csvsave(name,rows):
    with (OUT/name).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


full=fit_overlap(range(6))
grid=full.pop('_grid')
full['metrics']=metrics(full['residual_log'])
print('REPRESENTATIVE_FIT',json.dumps({k:full[k] for k in ['ws','ell_nm','amplitude','metrics']}),flush=True)
baselines=[fit_baseline(range(6),k) for k in ['D_over_h','D_over_h_plus_fixed_access','shared_constant_positive_w']]

loo=[]
for hold in range(6):
    train=[i for i in range(6) if i!=hold]
    ff=fit_overlap(train)
    ff.pop('_grid')
    row=dict(held_h_nm=float(THICK[hold]),observed=float(Y[hold]),predicted=ff['predicted'][hold],
        prediction_over_observed=ff['predicted'][hold]/Y[hold],held_log_residual=ff['residual_log'][hold],
        ws=ff['ws'],ell_nm=ff['ell_nm'],train_SSE_log=ff['train_SSE_log'])
    for k in ['D_over_h','D_over_h_plus_fixed_access','shared_constant_positive_w']:
        bb=fit_baseline(train,k)
        row[k+'_held_log_residual']=bb['residual_log'][hold]
    loo.append(row)
    print('LOO',json.dumps(row),flush=True)

split=fit_overlap([0,2,4,5]);split.pop('_grid')
split['held_indices']=[1,3]
split['held_metrics']=metrics(np.array(split['residual_log'])[[1,3]])

noaccess=fit_overlap(range(6),access=0.);noaccess.pop('_grid');noaccess['metrics']=metrics(noaccess['residual_log'])
rplus=fit_overlap(range(6),robin_r=1);rplus.pop('_grid');rplus['metrics']=metrics(rplus['residual_log'])
print('RPLUS1',json.dumps(rplus),flush=True)

profile=[]
profile_ells=sorted(set(np.geomspace(*ELL_BOUNDS,151).tolist()+
    np.geomspace(full['ell_nm']*.65,full['ell_nm']*1.5,121).tolist()+[full['ell_nm'],1.,25.,100.,300.,500.]))
for ell in profile_ells:
    def obj(lws):
        lr=log_resistance(THICK,float(np.exp(lws)),float(ell))
        _,r=evaluate_shape(lr,range(6))
        return float(np.dot(r,r))
    opt=minimize_scalar(obj,bounds=BOUNDS[0],method='bounded',options=dict(xatol=1e-8))
    profile.append(dict(ell_nm=float(ell),ws=float(np.exp(opt.x)),SSE_log=float(opt.fun),
        log_RMSE=math.sqrt(opt.fun/6)))
profile_ranges={}
for multiplier in [1.01,1.1,1.25]:
    accepted=[r for r in profile if r['SSE_log']<=multiplier*full['train_SSE_log']]
    profile_ranges[str(multiplier)]=dict(criterion='profile SSE <= multiplier * best SSE; descriptive, not confidence interval',
        ell_nm=[min(r['ell_nm'] for r in accepted),max(r['ell_nm'] for r in accepted)] if accepted else [],
        ws=[min(r['ws'] for r in accepted),max(r['ws'] for r in accepted)] if accepted else [])

ws=full['ws'];ell=full['ell_nm']
field=[]
for i,h in enumerate(THICK):
    field.append(dict(h_nm=float(h),observed_delta_kappa=float(Y[i]),predicted_delta_kappa=full['predicted'][i],
        prediction_over_observed=full['predicted'][i]/Y[i],log_residual=full['residual_log'][i],
        w_surface=ws*(1+math.exp(-h/ell)),w_center=2*ws*math.exp(-h/(2*ell)),
        w_average=2*ws*ell/h*(-math.expm1(-h/ell)),
        log_effective_resistance_nm=full['log_resistance'][i],
        resistance_nm=math.exp(full['log_resistance'][i])))

curve=[]
for h in np.geomspace(20,10000,1001):
    lr=log_resistance([h],ws,ell)[0]
    rplr=log_resistance([h],rplus['ws'],rplus['ell_nm'],robin_r=1)[0]
    curve.append(dict(h_nm=float(h),predicted_delta_kappa=float(np.exp(full['log_amplitude']-lr)),
        r_plus1_predicted_delta_kappa=float(np.exp(rplus['log_amplitude']-rplr)),
        w_surface=ws*(1+math.exp(-h/ell)),w_center=2*ws*math.exp(-h/(2*ell))))
mincurve=min(curve,key=lambda r:r['predicted_delta_kappa'])

lr48=log_resistance(THICK,ws,ell)
lr96=log_resistance(THICK,ws,ell,nodes=leggauss(96))
checks=dict(quadrature_48_vs96_max_logR_difference=float(np.max(abs(lr48-lr96))),
    neutral_ws0_max_logR_error=float(np.max(abs(log_resistance(THICK,0,ell)-np.log(THICK+ACCESS_NM)))),
    source_sha256_after=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
    source_unchanged=SOURCE_HASH_BEFORE==hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
    physical_ell_labels='network/interfacial structural correlation depth; not Debye length or molecular self-energy range')
result=dict(model='w(z,h)=ws*(exp(-z/ell)+exp(-(h-z)/ell)); R=integral exp(w)dz+pi*a/2',
    a_nm=A_NM,access_nm=ACCESS_NM,ws_bounds=WS_BOUNDS,ell_bounds_nm=ELL_BOUNDS,
    loss='natural-log residual; amplitude analytically profiled in every training set',
    source=str(SOURCE),source_sha256=SOURCE_HASH_BEFORE,
    environment=dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__),
    thickness_nm=THICK.tolist(),observed_delta_kappa=Y.tolist(),full_fit=full,baselines=baselines,
    leave_one_out=loo,LOO_metrics=metrics([r['held_log_residual'] for r in loo]),
    LOO_baseline_metrics={k:metrics([r[k+'_held_log_residual'] for r in loo]) for k in ['D_over_h','D_over_h_plus_fixed_access','shared_constant_positive_w']},
    heldout_120_350=split,no_access_sensitivity=noaccess,profile_ranges=profile_ranges,
    fixed_surface_r_plus1=rplus,
    fixed_surface_formula='w=ws*(exp(-z/ell)+exp(-(h-z)/ell))/(1+exp(-h/ell)); w(0)=w(h)=ws',
    amplitude_units='delta_kappa times nm; not an absolute D or permeance',
    short_ell_best_profile={str(cap):min((r for r in profile if r['ell_nm']<=cap),key=lambda r:r['SSE_log']) for cap in [25.,100.,300.,500.]},
    representative_fields=field,curve_minimum_grid=mincurve,checks=checks,
    limitations=[
      'Six labels are not six verified independent membrane replicates; no inferential confidence interval or sigma identification.',
      'Endpoint fit is a shared phenomenological shape comparison, not independent prediction or a dynamic finite-reservoir fit.',
      'Amplitude absorbs D, conductivity calibration, area/volume, driving force and endpoint time; those are not separately identified.',
      'A common affine conductivity relation and measurement geometry are conditional assumptions.',
      'No measured per-point errors; log residuals weight multiplicative errors equally.',
      'w is an effective neutral entry field; fit alone does not establish dielectric, molecular or electrostatic origin.',
      'The 80/85 min ambiguity does not change 0/120 min endpoint targets.'
    ])
csvsave('endpoint_fit.csv',field)
csvsave('leave_one_out.csv',loo)
csvsave('ell_profile.csv',profile)
csvsave('parameter_grid.csv',[dict(SSE_log=r[0],ws=r[1],ell_nm=r[2]) for r in grid])
csvsave('thickness_curve.csv',curve)
(OUT/'fit_results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
print('FINAL',json.dumps(dict(full=full,baselines=baselines,LOO=result['LOO_metrics'],
    LOO_baselines=result['LOO_baseline_metrics'],heldout=split,noaccess=noaccess,
    profile_ranges=profile_ranges,r_plus1=rplus,fields=field,minimum=mincurve,checks=checks),ensure_ascii=False,indent=2),flush=True)
