"""Descriptive neutral P fit/freeze/holdout. Standard library, SI only.

This is not a mechanism selector. Shared P requires the same material,
geometry/pre-treatment/chemistry domain and independently valid concentration
calibration. Synthetic data generated with the same model test software only.
"""
from pathlib import Path
import argparse, csv, hashlib, json, math, statistics
import models as m

ROOT=Path(__file__).resolve().parent
FIELDS=['source_type','source_ref','mother_id','run_id','time_s','area_m2',
        'volume_left_m3','volume_right_m3','c_left0_mol_m3','c_right0_mol_m3',
        'c_left_mol_m3','c_right_mol_m3']
NUMERIC=FIELDS[4:]
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def model_identity():
    return {p.name:sha(p) for p in [ROOT/'models.py',Path(__file__)]}
def write_json(p,obj):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding='utf8')
def read_data(path):
    with open(path,encoding='utf-8-sig',newline='') as f:
        reader=csv.DictReader(f)
        if not reader.fieldnames or any(k not in reader.fieldnames for k in FIELDS):
            raise ValueError('Required schema missing: '+','.join(FIELDS))
        rows=[]
        for row in reader:
            if row['source_type'] not in ['synthetic','real']:
                raise ValueError('source_type must be synthetic or real')
            if any(not row[k].strip() for k in FIELDS[:4]):
                raise ValueError('source reference/mother/run identity required')
            for k in NUMERIC:
                row[k]=float(row[k])
                if not math.isfinite(row[k]) or row[k]<0:
                    raise ValueError('Invalid nonnegative SI value: '+k)
            if min(row[k] for k in ['area_m2','volume_left_m3','volume_right_m3'])<=0:
                raise ValueError('Area and volumes positive')
            if row['c_left0_mol_m3']<=row['c_right0_mol_m3']:
                raise ValueError('This teaching pipeline requires positive left-to-right initial gradient')
            rows.append(row)
    if not rows: raise ValueError('Empty data')
    types=set(r['source_type'] for r in rows)
    if len(types)!=1: raise ValueError('Do not combine real and synthetic data')
    groups={}
    for r in rows: groups.setdefault(r['run_id'],[]).append(r)
    for name,run in groups.items():
        if len(run)<3: raise ValueError('At least three time points per whole run')
        times=[r['time_s'] for r in run]
        if len(times)!=len(set(times)): raise ValueError('Duplicate time in run '+name)
        for k in ['source_type','source_ref','mother_id','area_m2','volume_left_m3',
                  'volume_right_m3','c_left0_mol_m3','c_right0_mol_m3']:
            if len(set(r[k] for r in run))!=1:
                raise ValueError('Constant-volume teaching contract violated: '+k)
    return rows
def predict(row,P):
    vl,vr=row['volume_left_m3'],row['volume_right_m3']
    cl,cr=row['c_left0_mol_m3'],row['c_right0_mol_m3']
    eq=(vl*cl+vr*cr)/(vl+vr)
    delta=(cl-cr)*math.exp(-row['area_m2']*P*(1/vl+1/vr)*row['time_s'])
    return eq+vr*delta/(vl+vr),eq-vl*delta/(vl+vr)
def objective(rows,P):
    groups={}
    for r in rows:
        pred=predict(r,P)[1]
        error=(pred-r['c_right_mol_m3'])/(r['c_left0_mol_m3']-r['c_right0_mol_m3'])
        groups.setdefault((r['mother_id'],r['run_id']),[]).append(error*error)
    mothers={}
    for (mother,run),errors in groups.items():
        mothers.setdefault(mother,[]).append(statistics.mean(errors))
    return statistics.mean(statistics.mean(e) for e in mothers.values())
def fit(data,params):
    rows=read_data(data)
    a,b=-14.,-1.;rho=(math.sqrt(5)-1)/2
    c,d=b-rho*(b-a),a+rho*(b-a)
    fc,fd=objective(rows,10**c),objective(rows,10**d)
    for _ in range(100):
        if fc<fd:
            b,d,fd=d,c,fc;c=b-rho*(b-a);fc=objective(rows,10**c)
        else:
            a,c,fc=c,d,fd;d=a+rho*(b-a);fd=objective(rows,10**d)
    logP=(a+b)/2;P=10**logP
    if not -13.9<logP<-.1:
        raise ValueError('Fit at search boundary; P not resolved under this contract')
    result={'model':'neutral_shared_effective_P','parameter_scope':'descriptive total resistance',
            'permeability_m_s':P,'source_type':rows[0]['source_type'],
            'status':'synthetic_software_training' if rows[0]['source_type']=='synthetic' else 'descriptive_real_fit_not_mechanism_validation',
            'training_data_sha256':sha(data),'training_mother_ids':sorted(set(r['mother_id'] for r in rows)),
            'training_run_ids':sorted(set(r['run_id'] for r in rows)),
            'objective_equal_mother_and_run_weight':objective(rows,P),
            'model_source_sha256':model_identity(),
            'contract':'constant reservoirs/no sampling/common geometric and pretreatment domain; not a charge/wetting attribution'}
    write_json(params,result);return result
def evaluate(params,data,output):
    p=json.loads(Path(params).read_text());rows=read_data(data)
    if p['model_source_sha256']!=model_identity():
        raise ValueError('Model source changed after freezing; save new version explicitly')
    if rows[0]['source_type']!=p['source_type']:
        raise ValueError('Cannot treat synthetic fit as real validation')
    overlap=set(r['mother_id'] for r in rows)&set(p['training_mother_ids'])
    if overlap: raise ValueError('Held-out validation requires new mother IDs: '+str(overlap))
    if set(r['run_id'] for r in rows)&set(p['training_run_ids']):
        raise ValueError('Validation contains training run')
    preds=[];byrun={};mass=[]
    for r in rows:
        cl,cr=predict(r,p['permeability_m_s'])
        preds.append(dict(r,predicted_left_mol_m3=cl,predicted_right_mol_m3=cr,
                          right_error_mol_m3=cr-r['c_right_mol_m3']))
        byrun.setdefault(r['run_id'],[]).append((cr-r['c_right_mol_m3'])**2)
        before=r['volume_left_m3']*r['c_left0_mol_m3']+r['volume_right_m3']*r['c_right0_mol_m3']
        after=r['volume_left_m3']*r['c_left_mol_m3']+r['volume_right_m3']*r['c_right_mol_m3']
        mass.append(abs(after-before)/before)
    out=Path(output);out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('w',encoding='utf8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(preds[0]));writer.writeheader();writer.writerows(preds)
    summary={'frozen_parameter_sha256':sha(params),'test_data_sha256':sha(data),
             'source_type':rows[0]['source_type'],'new_mother_n':len(set(r['mother_id'] for r in rows)),
             'whole_run_n':len(byrun),'run_RMSE_mol_m3':{k:math.sqrt(statistics.mean(v)) for k,v in byrun.items()},
             'max_observed_relative_mass_error':max(mass),
             'experimental_mechanism_validated':False,
             'note':'Observed mass error is a diagnostic, not automatically dismissed or corrected; detection/calibration uncertainty requires separate records.'}
    write_json(out.with_suffix('.summary.json'),summary);return summary
def demo(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    times=[0,300,600,1200,2400,3600,7200]
    def generate(mothers,concs,target):
        rows=[]
        for i,mother in enumerate(mothers):
            for j,c in enumerate(concs):
                for t in times:
                    r={'source_type':'synthetic','source_ref':'generated_by_neutral_same_model_not_PE',
                       'mother_id':mother,'run_id':mother+'_'+str(j),'time_s':t,
                       'area_m2':1e-4,'volume_left_m3':10e-6,'volume_right_m3':10e-6,
                       'c_left0_mol_m3':c,'c_right0_mol_m3':0.}
                    cl,cr=predict(r,3e-7*(1+.04*(i-1)))
                    rows.append(dict(r,c_left_mol_m3=cl,c_right_mol_m3=cr))
        with target.open('w',encoding='utf8',newline='') as f:
            w=csv.DictWriter(f,fieldnames=FIELDS);w.writeheader();w.writerows(rows)
    generate(['SYN_T1','SYN_T2','SYN_T3'],[1.,100.],out/'train.csv')
    generate(['SYN_V1','SYN_V2'],[10.,30.],out/'validation.csv')
    frozen=fit(out/'train.csv',out/'frozen.json')
    score=evaluate(out/'frozen.json',out/'validation.csv',out/'predictions.csv')
    try: evaluate(out/'frozen.json',out/'train.csv',out/'forbidden.csv')
    except ValueError as e:
        if 'mother IDs' not in str(e): raise
        refusal=str(e)
    else: raise AssertionError('Training reuse not rejected')
    result={'fit':frozen,'validation':score,'training_reuse_correctly_rejected':refusal,
            'synthetic_does_not_prove_model_choice':True}
    write_json(out/'pipeline_receipt.json',result)
    return result
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='cmd',required=True)
    d=sub.add_parser('demo');d.add_argument('--out',default='outputs/pipeline')
    f=sub.add_parser('fit');f.add_argument('data');f.add_argument('params')
    e=sub.add_parser('evaluate');e.add_argument('params');e.add_argument('data');e.add_argument('output')
    args=parser.parse_args()
    result=demo(args.out) if args.cmd=='demo' else fit(args.data,args.params) if args.cmd=='fit' else evaluate(args.params,args.data,args.output)
    print(json.dumps(result,ensure_ascii=False,indent=2))
