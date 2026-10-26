"""Resumable crossed concentration/porosity noise pilot; synthetic data only."""
import hashlib,json,sys,time
from pathlib import Path
from dataclasses import replace
import numpy as np
from scipy.optimize import least_squares
from physics_model import Column,solve
from fast_column_sensitivity import sensitive

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'runs_crossed_replication'/'pilot_v1'
TRUTH=np.array([.001,48.,.04,.4])
NAMES=['k_min','density_capacity_product','KL_l_mg','porosity']
DESIGNS={'one_condition':[(100,40)],'two_concentrations':[(100,20),(200,20)]}

def save(name,value):
    tmp=OUT/(name+'.tmp');tmp.write_text(json.dumps(value,indent=2),encoding='utf-8');tmp.replace(OUT/name)

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    provenance={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in ['crossed_replication.py','physics_model.py','fast_column_sensitivity.py','CAPRAZ_TEKRAR_PROTOKOLU.md']}
    if (OUT/'provenance.json').exists():
        assert json.loads((OUT/'provenance.json').read_text())==provenance,'Sources changed; use a new run directory.'
    else:save('provenance.json',provenance)
    protocol=dict(concentration_seeds=[20261101,20261102,20261103],porosity_seeds=[20261111,20261112,20261113],sigma_c=.02,sigma_log_porosity=.05,density_observed=.8254263138594846,fit_cells=256,reference_cells=1024,check_cells=2048,max_nfev=30,starts=[[.0007,38.4,.024,.32],[.0015,57.6,.072,.5]],reference_gate=dict(outlet_rmse=.001,outlet_max=.002,relative_mass_error=1e-6),scope='Three concentration noise vectors crossed with three auxiliary observations; not nine independent datasets. Same FV family; synthetic only.')
    save('protocol.json',protocol)
    times=np.unique(np.r_[0,np.linspace(.005,1,20),np.linspace(.005,1,40),np.linspace(0,1,101)])
    if not (OUT/'reference.npz').exists():
        arrays={};checks={}
        for feed in [100,150,200]:
            curves=[]
            for n in [1024,2048]:
                tic=time.perf_counter();p=replace(Column(),inlet_mg_l=feed);r=solve(p,n,times=times,rtol=2e-9,atol=1e-11)
                storage=np.mean(r['c']+p.beta*r['s'],axis=1);flux=r['integrated_net_flux'];mass=float(np.max(abs(storage-flux))/max(np.max(abs(flux)),1e-12))
                curves.append(r['c'][:,-1]);checks[f'{feed}_{n}']=dict(relative_mass_error=mass,seconds=time.perf_counter()-tic,nfev=r['nfev'],c_min=float(r['c'].min()),c_max=float(r['c'].max()))
                print('REFERENCE',feed,n,checks[f'{feed}_{n}'],flush=True)
                save('reference_checks.json',checks)
            delta=curves[0]-curves[1];rmse=float(np.sqrt(np.mean(delta**2)));maximum=float(max(abs(delta)));checks[f'{feed}_grid']=dict(rmse=rmse,max_error=maximum)
            assert rmse<.001 and maximum<.002 and all(checks[f'{feed}_{n}']['relative_mass_error']<1e-6 for n in [1024,2048]),checks
            arrays[str(feed)]=curves[0]
        save('reference_checks.json',checks);np.savez_compressed(OUT/'reference.npz',times=times,**arrays)
    ref=np.load(OUT/'reference.npz');ts=ref['times']
    if not (OUT/'observations.npz').exists():
        obs={}
        for i,seed in enumerate(protocol['concentration_seeds']):
            z=np.random.default_rng(seed).normal(size=40)
            for design,conditions in DESIGNS.items():
                clean=np.concatenate([np.interp(np.linspace(.005,1,count),ts,ref[str(feed)]) for feed,count in conditions]);obs[f'{design}_{i}']=clean+.02*z;obs[f'{design}_{i}_clean']=clean
        obs['porosity_observed']=np.array([.4*np.exp(.05*np.random.default_rng(seed).normal()) for seed in protocol['porosity_seeds']]);np.savez_compressed(OUT/'observations.npz',**obs)
    obs=np.load(OUT/'observations.npz');rows=json.loads((OUT/'fits.json').read_text()) if (OUT/'fits.json').exists() else []
    for design,conditions in DESIGNS.items():
        for i in range(3):
            for j in range(3):
                for start_index,start in enumerate(protocol['starts']):
                    key=f'{design}_{i}_{j}_{start_index}'
                    if any(r['key']==key for r in rows):continue
                    cachedx=None;cached=None;tic=time.perf_counter()
                    def evaluate(x):
                        nonlocal cachedx,cached
                        if cachedx is None or not np.array_equal(x,cachedx):
                            pairs=[sensitive(x,feed,np.linspace(.005,1,count),256) for feed,count in conditions];pred=np.concatenate([p[0] for p in pairs]);J=np.vstack([p[1] for p in pairs])/.02
                            cached=(np.r_[(pred-obs[f'{design}_{i}'])/.02,(x[3]-np.log(obs['porosity_observed'][j]))/.05],np.vstack([J,[0,0,0,20]]));cachedx=x.copy()
                        return cached
                    print('FIT',key,flush=True)
                    fit=least_squares(lambda x:evaluate(x)[0],np.log(start),jac=lambda x:evaluate(x)[1],bounds=np.log([[1e-4,5,.005,.15],[.01,200,.4,.75]]),max_nfev=30,ftol=1e-8,xtol=1e-8,gtol=1e-6)
                    v=np.exp(fit.x);row=dict(key=key,design=design,concentration_index=i,porosity_index=j,start_index=start_index,estimate=dict(zip(NAMES,v.tolist())),relative_errors=(v/TRUTH-1).tolist(),cost=float(fit.cost),nfev=fit.nfev,success=bool(fit.success),status=int(fit.status),message=fit.message,optimality=float(fit.optimality),active_bounds=fit.active_mask.tolist(),seconds=time.perf_counter()-tic)
                    rows.append(row);save('fits.json',rows);print('SAVED',key,'cost',fit.cost,'seconds',row['seconds'],flush=True)
    selected=[]
    for design in DESIGNS:
        for i in range(3):
            for j in range(3):
                candidates=[r for r in rows if r['design']==design and r['concentration_index']==i and r['porosity_index']==j];best=min(candidates,key=lambda r:r['cost']);v=np.array([best['estimate'][n] for n in NAMES]);p=replace(Column(),inlet_mg_l=150,k_min=v[0],qmax_mg_g=v[1]/.8,KL_l_mg=v[2],porosity=v[3]);curve=solve(p,1024,times=np.linspace(0,1,101),rtol=2e-9,atol=1e-11)['c'][:,-1];target=np.interp(np.linspace(0,1,101),ts,ref['150']);selected.append(dict(best,common_150_rmse=float(np.sqrt(np.mean((curve-target)**2))),start_cost_gap=float(max(r['cost'] for r in candidates)-best['cost'])));save('selected.json',selected)
    save('result.json',dict(protocol=protocol,truth=dict(zip(NAMES,TRUTH.tolist())),fits=rows,selected=selected));print('DONE',OUT,flush=True)

if __name__=='__main__':main()
