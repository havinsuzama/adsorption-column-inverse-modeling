"""Sparse nonlinear KL profile diagnostic; not calibrated confidence intervals."""
import hashlib,json,time
import numpy as np
from scipy.optimize import least_squares
from parameter_regimes import ROOT,NAMES,DESIGNS,save
from fast_column_sensitivity import sensitive

OUT=ROOT/'runs_high_affinity_profile'/'pilot_v1'
GRID=ROOT/'runs_regime_grid'/'pilot_v1'
LOCAL=ROOT/'runs_local_regime_uncertainty'/'pilot_v1'
DATA=ROOT/'runs_parameter_regimes'/'pilot_v1'/'high_affinity'
FREE=np.array([0,1,3]);BOUNDS=np.log([[1e-4,5,.005,.15],[.01,200,.4,.75]])

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    files={n:ROOT/n for n in ['high_affinity_profile.py','fast_column_sensitivity.py','parameter_regimes.py','physics_model.py']}
    files.update(grid_result=GRID/'result.json',local_result=LOCAL/'result.json',observations=DATA/'observations.npz')
    provenance={n:hashlib.sha256(p.read_bytes()).hexdigest() for n,p in files.items()}
    if (OUT/'provenance.json').exists():assert json.loads((OUT/'provenance.json').read_text())==provenance
    else:save(OUT/'provenance.json',provenance)
    bases=json.loads((GRID/'result.json').read_text())['results'];local=json.loads((LOCAL/'result.json').read_text())['results'];obs=np.load(DATA/'observations.npz')
    tasks=[]
    for design in DESIGNS:
        base=next(r for r in bases if r['case']=='high_affinity' and r['design']==design);g=next(r for r in local if r['case']=='high_affinity' and r['design']==design)['with_porosity'];x=np.log([base['estimate'][n] for n in NAMES]);std=g['log_parameter_std'][2];cov=np.array(g['covariance'])
        for z in [0,-1,1,-2,2]:
            target=float(np.clip(x[2]+z*std,BOUNDS[0,2],BOUNDS[1,2]));delta=target-x[2]
            predicted=np.clip(x[FREE]+cov[FREE,2]/cov[2,2]*delta,BOUNDS[0,FREE]+1e-9,BOUNDS[1,FREE]-1e-9)
            for start_name,start in [('local_prediction',predicted)]+([('baseline_nuisance',x[FREE])] if abs(z)==2 else []):
                tasks.append(dict(key=f'{design}_{z}_{start_name}',design=design,z_requested=z,target_log_KL=target,target_KL=float(np.exp(target)),boundary_clipped=bool(abs(delta-z*std)>1e-12),actual_local_z=float(delta/std),start_name=start_name,start=start.tolist(),base_log_parameters=x.tolist(),base_cost=base['cost']))
    save(OUT/'protocol.json',dict(case='high_affinity',concentration_index=0,porosity_index=0,fit_cells=512,reference_cells=1024,max_nfev=20,parameter_order=NAMES,tasks=tasks,scope='Post-pilot development diagnostic: KL fixed at fitted log value plus 0,+/-1,+/-2 local standard deviations, clipped to original bounds. Other 3 parameters reoptimized. Extra baseline-nuisance start only at extremes. 14 fits, same observed objective and bounds. No threshold-based confidence interval or coverage claim; sparse points, not a full profile.'))
    rows=json.loads((OUT/'fits.json').read_text()) if (OUT/'fits.json').exists() else []
    for task in tasks:
        if any(r['key']==task['key'] for r in rows):continue
        design=task['design'];cachedx=None;cached=None;tic=time.perf_counter()
        def evaluate(u):
            nonlocal cachedx,cached
            if cachedx is None or not np.array_equal(cachedx,u):
                x=np.array(task['base_log_parameters']);x[2]=task['target_log_KL'];x[FREE]=u
                pairs=[sensitive(x,feed,np.linspace(.005,1,count),512) for feed,count in DESIGNS[design].copy()];pred=np.concatenate([p[0] for p in pairs]);J=np.vstack([p[1] for p in pairs])[:,FREE]/.02
                cached=(np.r_[(pred-obs[design+'_0'])/.02,(x[3]-np.log(obs['porosity_observed'][0]))/.05],np.vstack([J,[0,0,20]]));cachedx=u.copy()
            return cached
        print('PROFILE',task['key'],'KL',task['target_KL'],flush=True)
        try:
            initial=float(.5*np.sum(evaluate(np.array(task['start']))[0]**2));fit=least_squares(lambda u:evaluate(u)[0],task['start'],jac=lambda u:evaluate(u)[1],bounds=BOUNDS[:,FREE],max_nfev=20,ftol=1e-8,xtol=1e-8,gtol=1e-6);x=np.array(task['base_log_parameters']);x[2]=task['target_log_KL'];x[FREE]=fit.x
            row=dict(task,estimate=dict(zip(NAMES,np.exp(x).tolist())),initial_cost=initial,cost=float(fit.cost),success=bool(fit.success),status=int(fit.status),message=fit.message,nfev=fit.nfev,optimality=float(fit.optimality),free_active_bounds=fit.active_mask.tolist())
        except Exception as exc:row=dict(task,success=False,error=repr(exc))
        row['seconds']=time.perf_counter()-tic;rows.append(row);save(OUT/'fits.json',rows);print('SAVED',task['key'],row.get('cost',row.get('error')),flush=True)
    points=[]
    for design in DESIGNS:
        center=[r for r in rows if r['design']==design and r['z_requested']==0 and 'cost' in r]
        if not center:continue
        baseline=min(r['cost'] for r in center)
        for z in [-2,-1,0,1,2]:
            candidates=[r for r in rows if r['design']==design and r['z_requested']==z and 'cost' in r]
            if not candidates:continue
            best=min(candidates,key=lambda r:r['cost']);points.append(dict(best,twice_delta_cost=2*(best['cost']-baseline),local_quadratic_prediction=best['actual_local_z']**2,start_cost_gap=max(r['cost'] for r in candidates)-best['cost'],available_starts=len(candidates)))
    save(OUT/'result.json',dict(fits=rows,points=points));print('DONE',OUT,flush=True)

if __name__=='__main__':main()
