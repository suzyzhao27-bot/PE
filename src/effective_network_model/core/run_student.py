"""Reproduce the frozen shared model. Only the Python standard library is used."""
import argparse
import csv
import hashlib
import json
import math
import platform
import sys
from pathlib import Path
from theory_model import (field,log_resistance_nm,relative_permeance,thickness_log_slope,
                          neutral_profile,normalized_finite_bath)

ROOT=Path(__file__).resolve().parent


def write_rows(path,rows):
    with path.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)


def run(output):
    output.mkdir(parents=True,exist_ok=True)
    params=json.loads((ROOT/'inputs'/'model_parameters.json').read_text(encoding='utf-8'))
    source=ROOT/'inputs'/'Book1__Sheet1.csv'
    with source.open(encoding='utf-8-sig',newline='') as f:
        raw=list(csv.reader(f))[1:]
    mapping=[(65,1),(120,3),(200,4),(350,5),(700,6),(2600,7)]
    rows=[]
    for h,col in mapping:
        y=float(raw[-1][col])-float(raw[0][col])
        logr=log_resistance_nm(h,params['ws_ref'],params['ell_nm'],params['a_ref_nm'])
        predicted=math.exp(params['log_amplitude']-logr)
        rows.append(dict(h_label_nm=h,observed_delta_kappa=y,predicted_delta_kappa=predicted,
                         predicted_over_observed=predicted/y,log_residual=math.log(predicted/y),
                         w_surface_kBT=field(0,h,params['ws_ref'],params['ell_nm']),
                         w_center_kBT=field(h/2,h,params['ws_ref'],params['ell_nm'])))
    write_rows(output/'shared_model_comparison.csv',rows)
    thick=[]
    for i in range(241):
        h=20*(10000/20)**(i/240)
        logr=log_resistance_nm(h,params['ws_ref'],params['ell_nm'],params['a_ref_nm'])
        thick.append(dict(h_nm=h,predicted_delta_kappa=math.exp(params['log_amplitude']-logr),
                          relative_P_to_2600=relative_permeance(h,25,params),
                          w_center_kBT=field(h/2,h,params['ws_ref'],params['ell_nm']),
                          dlogP_dlogh=thickness_log_slope(h,params['ws_ref'],params['ell_nm'])))
    write_rows(output/'thickness_scan.csv',thick)
    local=[]
    for a in (5,10,15,20,25,30,40,60,100):
        for h in (65,120,200,350,700,1200,2600):
            local.append(dict(a_nm=a,h_nm=h,ws_kBT=params['ws_ref']*25/a,
                              relative_P_to_reference=relative_permeance(h,a,params),
                              evidence='conditional S/V closure; no radius fit'))
    write_rows(output/'local_confinement_scan.csv',local)
    profiles=[]
    for h in (65,350,700,2600):
        for row in neutral_profile(h,params['ws_ref'],params['ell_nm'],270.,1.,200):
            profiles.append(dict(h_nm=h,**row))
    write_rows(output/'neutral_profiles.csv',profiles)
    shapes=[]
    for h,col in mapping:
        timecol=2 if h==120 else 0
        base=float(raw[0][col]);delta=float(raw[-1][col])-base
        for row in raw:
            t=float(row[timecol])
            shapes.append(dict(h_nm=h,time_min=t,observed_normalized=(float(row[col])-base)/delta,
                               constant_P_sink=t/120,
                               constant_P_theta_0p223=normalized_finite_bath(t/120,-math.log(.8)),
                               note='time mapping inherited; endpoints unambiguous'))
    write_rows(output/'time_shapes.csv',shapes)
    summary=dict(version='5.0.0',python=sys.version,platform=platform.platform(),
                 source_csv_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                 shared_parameters=params,
                 log_rmse=math.sqrt(sum(r['log_residual']**2 for r in rows)/len(rows)),
                 max_fit_error_factor=max(max(r['predicted_over_observed'],1/r['predicted_over_observed']) for r in rows),
                 endpoint_count=len(rows),
                 interpretation='Shared effective-field fit, not independent PE microscopic mechanism validation',
                 scope='neutral spatial resistance plus fixed S/V closure; finite PNP is separate optional module',
                 no_new_experiments=True)
    (output/'run_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(summary,ensure_ascii=False,indent=2))


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--output',type=Path,default=ROOT/'student_results')
    run(ap.parse_args().output)
