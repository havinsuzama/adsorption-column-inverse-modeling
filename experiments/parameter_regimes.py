"""Four predeclared synthetic regimes; resumable references and 64 inverse fits."""
import argparse,hashlib,json,platform,sys,time
from dataclasses import replace
from pathlib import Path
import numpy as np
import scipy
from scipy.optimize import least_squares
from physics_model import Column,solve
from fast_column_sensitivity import sensitive

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'runs_parameter_regimes'/'pilot_v1'
NAMES=['k_min','density_capacity_product','KL_l_mg','porosity']
CASES={'slow_kinetics':[.0005,48,.04,.4],'fast_kinetics':[.002,48,.04,.4],'low_affinity':[.001,48,.01,.4],'high_affinity':[.001,48,.16,.4]}
DESIGNS={'one_condition':[(100,40)],'two_concentrations':[(100,20),(200,20)]}
STARTS=[[.0007,38.4,.024,.32],[.0015,57.6,.072,.5]]

def save(path,value):
    tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(value,indent=2),encoding='utf-8');tmp.replace(path)

def column(v,feed):
    return replace(Column(),inlet_mg_l=feed,k_min=v[0],qmax_mg_g=v[1]/.8,KL_l_mg=v[2],porosity=v[3])

def setup():
    OUT.mkdir(parents=True,exist_ok=True)
    provenance={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in ['parameter_regimes.py','physics_model.py','fast_column_sensitivity.py','FARKLI_PARAMETRELER_PROTOKOLU.md']}
    if (OUT/'provenance.json').exists():assert json.loads((OUT/'provenance.json').read_text())==provenance,'Source changed: preserve this run and use a new directory.'
    else:save(OUT/'provenance.json',provenance)
    save(OUT/'protocol.json',dict(cases=CASES,designs=DESIGNS,starts=STARTS,concentration_seeds=[20261121,20261122],porosity_seeds=[20261131,20261132],sigma_c=.02,sigma_log_porosity=.05,density_observed=.8254263138594846,max_nfev=30,reference_cells=1024,check_cells=2048,fit_cells=256,reference_gate=dict(rmse=.001,max_error=.002,relative_mass_error=1e-6)))
    save(OUT/'environment.json',dict(python=sys.version,platform=platform.platform(),numpy=np.__version__,scipy=scipy.__version__))

def references():
    times=np.unique(np.r_[0,np.linspace(.005,1,20),np.linspace(.005,1,40),np.linspace(0,1,101)])
    for case,v in CASES.items():
        folder=OUT/case;folder.mkdir(exist_ok=True)
        arrays={};checks={}
        for feed in [100,150,200]:
            curves=[]
            for n in [1024,2048]:
                cache=folder/f'reference_{feed}_{n}.npz'
                if cache.exists():
                    cached=np.load(cache);curves.append(cached['outlet']);checks[f'{feed}_{n}']=json.loads(str(cached['check']));continue
                tic=time.perf_counter();p=column(v,feed);r=solve(p,n,times=times,rtol=2e-9,atol=1e-11)
                flux=r['integrated_net_flux'];storage=np.mean(r['c']+p.beta*r['s'],axis=1)
                check=dict(relative_mass_error=float(max(abs(storage-flux))/max(max(abs(flux)),1e-12)),c_min=float(r['c'].min()),c_max=float(r['c'].max()),seconds=time.perf_counter()-tic,nfev=r['nfev'])
                curves.append(r['c'][:,-1]);checks[f'{feed}_{n}']=check
                np.savez_compressed(cache,times=times,outlet=curves[-1],check=json.dumps(check));print('REFERENCE',case,feed,n,check,flush=True)
            delta=curves[0]-curves[1];grid=dict(rmse=float(np.sqrt(np.mean(delta**2))),max_error=float(max(abs(delta))))
            grid['passed']=bool(grid['rmse']<.001 and grid['max_error']<.002 and all(checks[f'{feed}_{n}']['relative_mass_error']<1e-6 for n in [1024,2048]));checks[f'{feed}_grid']=grid;arrays[str(feed)]=curves[0]
            save(folder/'reference_checks.json',checks)
            if not grid['passed']:raise RuntimeError(f'Reference gate failed: {case} {feed}: {grid}; retain case and revise resolution before fitting.')
        np.savez_compressed(folder/'reference.npz',times=times,**arrays)
    save(OUT/'reference_status.json',dict(all_four_cases_passed=True))

def fit_all():
    assert json.loads((OUT/'reference_status.json').read_text())['all_four_cases_passed']
    for case,truth in CASES.items():
        folder=OUT/case;ref=np.load(folder/'reference.npz');ts=ref['times']
        if not (folder/'observations.npz').exists():
            arrays={}
            for i,seed in enumerate([20261121,20261122]):
                z=np.random.default_rng(seed).normal(size=40)
                for design,conditions in DESIGNS.items():
                    clean=np.concatenate([np.interp(np.linspace(.005,1,count),ts,ref[str(feed)]) for feed,count in conditions]);arrays[f'{design}_{i}']=clean+.02*z;arrays[f'{design}_{i}_clean']=clean
            arrays['porosity_observed']=np.array([.4*np.exp(.05*np.random.default_rng(seed).normal()) for seed in [20261131,20261132]])
            np.savez_compressed(folder/'observations.npz',**arrays)
        obs=np.load(folder/'observations.npz');rows=json.loads((folder/'fits.json').read_text()) if (folder/'fits.json').exists() else []
        for design,conditions in DESIGNS.items():
            for i in range(2):
                for j in range(2):
                    for si,start in enumerate(STARTS):
                        key=f'{design}_{i}_{j}_{si}'
                        if any(r['key']==key for r in rows):continue
                        lastx=None;cached=None;tic=time.perf_counter()
                        def evaluate(x):
                            nonlocal lastx,cached
                            if lastx is None or not np.array_equal(lastx,x):
                                pairs=[sensitive(x,feed,np.linspace(.005,1,count),256) for feed,count in conditions];pred=np.concatenate([p[0] for p in pairs]);J=np.vstack([p[1] for p in pairs])/.02
                                cached=(np.r_[(pred-obs[f'{design}_{i}'])/.02,(x[3]-np.log(obs['porosity_observed'][j]))/.05],np.vstack([J,[0,0,0,20]]));lastx=x.copy()
                            return cached
                        print('FIT',case,key,flush=True)
                        try:
                            fit=least_squares(lambda x:evaluate(x)[0],np.log(start),jac=lambda x:evaluate(x)[1],bounds=np.log([[1e-4,5,.005,.15],[.01,200,.4,.75]]),max_nfev=30,ftol=1e-8,xtol=1e-8,gtol=1e-6);v=np.exp(fit.x)
                            row=dict(key=key,design=design,concentration_index=i,porosity_index=j,start_index=si,estimate=dict(zip(NAMES,v.tolist())),relative_errors=(v/np.array(truth)-1).tolist(),cost=float(fit.cost),nfev=fit.nfev,success=bool(fit.success),status=int(fit.status),message=fit.message,optimality=float(fit.optimality),active_bounds=fit.active_mask.tolist())
                        except Exception as exc:
                            row=dict(key=key,design=design,concentration_index=i,porosity_index=j,start_index=si,success=False,error=repr(exc))
                        row['seconds']=time.perf_counter()-tic;rows.append(row);save(folder/'fits.json',rows);print('SAVED',case,key,row.get('cost',row.get('error')),flush=True)
        selected=[]
        for design in DESIGNS:
            for i in range(2):
                for j in range(2):
                    candidates=[r for r in rows if r['design']==design and r['concentration_index']==i and r['porosity_index']==j and 'cost' in r]
                    if not candidates:
                        selected.append(dict(design=design,concentration_index=i,porosity_index=j,error='No finite fit'));continue
                    best=min(candidates,key=lambda r:r['cost']);v=[best['estimate'][n] for n in NAMES];times=np.linspace(0,1,101);curve=solve(column(v,150),1024,times=times,rtol=2e-9,atol=1e-11)['c'][:,-1];target=np.interp(times,ts,ref['150']);selected.append(dict(best,common_150_rmse=float(np.sqrt(np.mean((curve-target)**2))),start_cost_gap=float(max(r['cost'] for r in candidates)-best['cost'])));save(folder/'selected.json',selected)
        save(folder/'result.json',dict(case=case,truth=dict(zip(NAMES,truth)),fits=rows,selected=selected));print('CASE DONE',case,flush=True)
    save(OUT/'completion.json',dict(fits_complete=True,cases=list(CASES)))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--references-only',action='store_true');args=parser.parse_args();setup();references()
    if not args.references_only:fit_all()
