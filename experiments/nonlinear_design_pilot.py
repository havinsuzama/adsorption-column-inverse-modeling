"""Paired nonlinear inverse pilot with separate data and optimizer seeds.

Fit log(k, rho*qmax, KL, eps); independent noisy log-density determines qmax.
Exact forward sensitivities avoid four extra PDE solves per Jacobian.
"""
import argparse,hashlib,json,time,uuid
from datetime import datetime
from dataclasses import replace,asdict
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
from scipy.sparse import lil_matrix
from physics_model import Column,solve
ROOT=Path(__file__).resolve().parent
NAMES=['k_min','density_capacity_product','KL_l_mg','porosity']

def sensitive(logp,feed,times,n=128):
    k,prod,kl,eps=np.exp(logp);p=replace(Column(),k_min=k,qmax_mg_g=prod/.8,KL_l_mg=kl,porosity=eps,inlet_mg_l=feed);A=p.A;beta=p.beta;Pe=p.Pe;rate=p.duration_min*k;a=kl*feed;cap=p.qmax_mg_g/p.qref_mg_g
    def rhs(t,y):
        state=y.reshape(5,2,n);c,s=state[0];E=cap*a*c/(1+a*c);Ec=cap*a/(1+a*c)**2;ds=rate*(E-s)
        flux=np.empty(n+1);flux[0]=1;flux[-1]=c[-1];flux[1:n]=c[:-1]-n/Pe*np.diff(c);dc=-A*n*np.diff(flux)-beta*ds
        result=np.empty_like(state);result[0,0]=dc;result[0,1]=ds
        for j in range(4):
            u,v=state[j+1];sv=rate*(Ec*u-v)
            if j==0:sv=sv+ds
            if j==1:sv=sv+rate*E
            if j==2:sv=sv+rate*E/(1+a*c)
            f=np.empty(n+1);f[0]=0;f[-1]=u[-1];f[1:n]=u[:-1]-n/Pe*np.diff(u)
            if j==3:f[1:n]-=n/Pe*np.diff(c)
            su=-A*n*np.diff(f)-beta*sv
            if j==3:su=su-dc
            result[j+1,0]=su;result[j+1,1]=sv
        return result.ravel()
    # Arrow-block sparsity: each sensitivity depends on the base and itself.
    pattern=lil_matrix((10*n,10*n),dtype=int)
    for block in range(5):
        offset=2*n*block
        for i in range(n):
            neighbors=list(range(max(0,i-1),min(n,i+2)))
            pattern[offset+i,[offset+j for j in neighbors]+[offset+n+i]]=1
            pattern[offset+n+i,[offset+i,offset+n+i]]=1
            if block:
                pattern[offset+i,neighbors+[n+i]]=1;pattern[offset+n+i,[i,n+i]]=1
    sol=solve_ivp(rhs,(0,1),np.zeros(10*n),method='BDF',t_eval=times,rtol=2e-8,atol=1e-10,jac_sparsity=pattern.tocsr())
    if not sol.success:raise RuntimeError(sol.message)
    y=sol.y.reshape(5,2,n,-1);return y[0,0,-1],y[1:,0,-1].T

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--noise-count',type=int,default=3);ap.add_argument('--max-nfev',type=int,default=40);args=ap.parse_args()
    if args.noise_count<1 or args.max_nfev<1:ap.error('positive counts required')
    out=ROOT/'runs_nonlinear_design'/(datetime.now().strftime('%Y%m%d_%H%M%S')+'_'+uuid.uuid4().hex[:6]);out.mkdir(parents=True);print('Output',out,flush=True)
    p=Column();truth=np.array([p.k_min,p.bulk_density_kg_l*p.qmax_mg_g,p.KL_l_mg,p.porosity]);checks=[]
    # Verify augmented solver and all four log sensitivities against independent
    # central differences of the original finite-volume implementation.
    for feed,values in [(100.,truth),(200.,truth*np.array([1.2,.9,1.3,1.1]))]:
        ts=np.linspace(.005,1,20);curve,J=sensitive(np.log(values),feed,ts,64);b=replace(p,k_min=values[0],qmax_mg_g=values[1]/.8,KL_l_mg=values[2],porosity=values[3],inlet_mg_l=feed);direct=solve(b,64,times=ts,rtol=2e-9,atol=1e-11)['c'][:,-1];columns=[]
        for j in range(4):
            plus=values.copy();minus=values.copy();plus[j]*=np.exp(.001);minus[j]*=np.exp(-.001)
            def forward(v):return solve(replace(b,k_min=v[0],qmax_mg_g=v[1]/.8,KL_l_mg=v[2],porosity=v[3]),64,times=ts,rtol=2e-9,atol=1e-11)['c'][:,-1]
            columns.append((forward(plus)-forward(minus))/.002)
        fd=np.stack(columns,1);ce=float(np.max(abs(curve-direct)));je=float(np.linalg.norm(J-fd)/np.linalg.norm(fd));assert ce<1e-5 and je<1e-3,(ce,je);checks.append(dict(feed=feed,curve_max_error=ce,jacobian_relative_error=je));print('CHECK',checks[-1],flush=True)
    (out/'checks.json').write_text(json.dumps(checks,indent=2));designs={'one_condition':[(100.,40)],'two_concentrations':[(100.,20),(200.,20)]};clean={};geometry={};reference=[]
    for name,design in designs.items():
        curves=[];geometry[name]=[]
        for feed,count in design:
            ts=np.linspace(.005,1,count);b=replace(p,inlet_mg_l=feed);r=solve(b,512,times=ts,rtol=2e-9,atol=1e-11);curves.append(r['c'][:,-1]);geometry[name].append((feed,ts));mass=float(np.max(abs((r['c']+b.beta*r['s']).mean(1)-r['integrated_net_flux']))/max(r['integrated_net_flux'].max(),1));assert mass<1e-6;reference.append(dict(design=name,feed=feed,mass_relative_error=mass))
        clean[name]=np.concatenate(curves)
    starts=[np.log([.0007,38.4,.024,.32]),np.log([.0015,57.6,.072,.5])];bounds=np.log([[1e-4,5.,.005,.15],[.01,200.,.4,.75]])
    observations={};rows=[];tic=time.perf_counter()
    for noise_index in range(args.noise_count):
        data_seed=20261026+noise_index;rng=np.random.default_rng(data_seed);noise=rng.normal(size=40);density_z=float(rng.normal());density_observed=p.bulk_density_kg_l*np.exp(.05*density_z)
        for name in designs:
            observed=clean[name]+.02*noise;observations[f'{name}_{noise_index}']=observed
            for start_index,start in enumerate(starts):
                cachex=None;cachevalue=None;calls=0
                def evaluate(x):
                    nonlocal cachex,cachevalue,calls
                    if cachex is None or not np.array_equal(x,cachex):
                        pairs=[sensitive(x,feed,ts,128) for feed,ts in geometry[name]];cachevalue=(np.concatenate([v[0] for v in pairs]),np.vstack([v[1] for v in pairs]));cachex=x.copy();calls+=1
                    return cachevalue
                def residual(x):return (evaluate(x)[0]-observed)/.02
                def jac(x):return evaluate(x)[1]/.02
                print('FIT',noise_index,name,start_index,flush=True);starttime=time.perf_counter();fit=least_squares(residual,start,jac=jac,bounds=bounds,max_nfev=args.max_nfev,xtol=1e-7,ftol=1e-7,gtol=1e-7);values=np.exp(fit.x);fitcurve=evaluate(fit.x)[0];elapsed=time.perf_counter()-starttime
                row=dict(noise_index=noise_index,data_seed=data_seed,start_index=start_index,initial_values=np.exp(start).tolist(),estimate=dict(zip(NAMES,values.tolist())),relative_parameter_errors=(values/truth-1).tolist(),density_observed=density_observed,qmax_estimate=values[1]/density_observed,qmax_relative_error=(values[1]/density_observed)/p.qmax_mg_g-1,cost=float(fit.cost),optimality=float(fit.optimality),success=bool(fit.success),status=int(fit.status),message=fit.message,nfev=int(fit.nfev),sensitivity_evaluations=calls,active_bounds=fit.active_mask.tolist(),noise_free_outlet_rmse=float(np.sqrt(np.mean((fitcurve-clean[name])**2))),observed_outlet_rmse=float(np.sqrt(np.mean((fitcurve-observed)**2))),seconds=elapsed,design=name)
                rows.append(row);(out/'partial_results.json').write_text(json.dumps(rows,indent=2));print('RESULT',json.dumps(row),flush=True)
    data=dict(truth=dict(zip(NAMES,truth.tolist())),assumptions=asdict(p),checks=checks,reference_checks=reference,fit_grid=128,reference_grid=512,noise_sigma=.02,log_density_sigma=.05,noise_count=args.noise_count,starts=2,max_nfev=args.max_nfev,bounds_physical=np.exp(bounds).tolist(),results=rows,seconds=time.perf_counter()-tic,scope='Paired Gaussian standardized noise vectors, distinct fixed starts. No clipping. Identifiable product fitted with density gauge .8; qmax derived from an independent noisy log-density measurement, exact separability for this model. Porosity unknown with explicit bounds. Three noise realizations are pilot evidence, not calibrated confidence intervals or blind final data.')
    (out/'result.json').write_text(json.dumps(data,indent=2),encoding='utf-8');np.savez_compressed(out/'observations.npz',**observations,**{n+'_clean':v for n,v in clean.items()})
    (out/'provenance.json').write_text(json.dumps({n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in ['physics_model.py','nonlinear_design_pilot.py']},indent=2));print('DONE',out,flush=True)
if __name__=='__main__':main()
