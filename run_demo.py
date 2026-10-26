"""Short CPU entry point; default forward demo, optional coarse inverse demo."""
import argparse,json
from dataclasses import replace
from pathlib import Path
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from physics_model import Column,solve

ROOT=Path(__file__).resolve().parent

def forward(out):
    times=np.linspace(0,1,121);fig=Figure(figsize=(8,4),layout='constrained');FigureCanvasAgg(fig);ax=fig.subplots();rows=[]
    for feed in [100.,150.,200.]:
        p=replace(Column(),inlet_mg_l=feed);r=solve(p,128,times=times)
        stored=np.mean(r['c']+p.beta*r['s'],axis=1);flux=r['integrated_net_flux'];mass=float(np.max(abs(stored-flux))/max(np.max(abs(flux)),1e-12))
        ax.plot(times*p.duration_min,r['c'][:,-1],label=f'C0={feed:g} mg/L');rows.append(dict(inlet_mg_l=feed,cells=128,relative_mass_residual=mass,outlet_min=float(r['c'][:,-1].min()),outlet_max=float(r['c'][:,-1].max())))
    ax.set(xlabel='Time (min)',ylabel='Normalized outlet concentration C/C0',title='Synthetic Langmuir–LDF column');ax.legend();fig.savefig(out/'forward_demo.png',dpi=160)
    data=dict(scope='Illustrative synthetic forward run at one grid, not a convergence study or experimental validation.',results=rows)
    (out/'forward_demo.json').write_text(json.dumps(data,indent=2));print(json.dumps(data,indent=2))

def inverse(out,seed):
    from scipy.optimize import least_squares
    from fast_column_sensitivity import sensitive
    times=np.linspace(.005,1,20);truth=np.array([.001,48.,.04,.4]);rng=np.random.default_rng(seed);feeds=[100.,200.]
    clean=np.concatenate([solve(replace(Column(),inlet_mg_l=f),256,times=times)['c'][:,-1] for f in feeds]);observed=clean+.02*rng.normal(size=40);eps_obs=.4*np.exp(.05*rng.normal());lastx=None;cached=None
    def evaluate(x):
        nonlocal lastx,cached
        if lastx is None or not np.array_equal(lastx,x):
            pairs=[sensitive(x,f,times,64) for f in feeds];pred=np.concatenate([p[0] for p in pairs]);J=np.vstack([p[1] for p in pairs])/.02
            cached=(np.r_[(pred-observed)/.02,(x[3]-np.log(eps_obs))/.05],np.vstack([J,[0,0,0,20]]),pred);lastx=x.copy()
        return cached
    fit=least_squares(lambda x:evaluate(x)[0],np.log([.0007,38.4,.024,.32]),jac=lambda x:evaluate(x)[1],bounds=np.log([[1e-4,5,.005,.15],[.01,200,.4,.75]]),max_nfev=15,ftol=1e-8,xtol=1e-8,gtol=1e-6)
    estimate=np.exp(fit.x);prediction=evaluate(fit.x)[2];names=['k_min','density_capacity_product','KL_l_mg','porosity'];data=dict(seed=seed,scope='Coarse synthetic example: 256-cell data, 64-cell fit, one start, 15 objective evaluations maximum. Not the archived benchmark.',truth=dict(zip(names,truth.tolist())),estimate=dict(zip(names,estimate.tolist())),porosity_observed=float(eps_obs),relative_errors=(estimate/truth-1).tolist(),cost=float(fit.cost),success=bool(fit.success),status=int(fit.status),message=fit.message,nfev=fit.nfev,optimality=float(fit.optimality),active_bounds=fit.active_mask.tolist())
    (out/'inverse_demo.json').write_text(json.dumps(data,indent=2));np.savez_compressed(out/'inverse_demo.npz',times=times,observed=observed,clean=clean,prediction=prediction)
    fig=Figure(figsize=(9,4),layout='constrained');FigureCanvasAgg(fig);axes=fig.subplots(1,2)
    for i,(ax,f) in enumerate(zip(axes,feeds)):
        sl=slice(20*i,20*(i+1));ax.plot(times,clean[sl],label='Synthetic reference');ax.scatter(times,observed[sl],s=14,label='Noisy observations');ax.plot(times,prediction[sl],'--',label='Coarse inverse fit');ax.set(title=f'C0={f:g} mg/L',xlabel='Normalized time',ylabel='C/C0');ax.legend(fontsize=8)
    fig.savefig(out/'inverse_demo.png',dpi=160);print(json.dumps(data,indent=2))

def results():
    source=ROOT/'experiments'/'runs_parameter_regimes'/'pilot_v1'
    print('Archived synthetic regimes: median absolute relative errors; NOT independent-cell confidence intervals.')
    for case in ['slow_kinetics','fast_kinetics','low_affinity','high_affinity']:
        d=json.loads((source/case/'result.json').read_text());print('\n'+case)
        for design in ['one_condition','two_concentrations']:
            rows=[r for r in d['selected'] if r['design']==design];err=100*np.median(np.abs([r['relative_errors'] for r in rows]),axis=0);rms=np.median([r['common_150_rmse'] for r in rows]);print(f'  {design}: k={err[0]:.2f}%, P={err[1]:.2f}%, KL={err[2]:.2f}%, eps={err[3]:.2f}%; outlet RMS={rms:.6f}')

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--mode',choices=['forward','inverse','results'],default='forward');parser.add_argument('--seed',type=int,default=20261121);parser.add_argument('--out',type=Path,default=ROOT/'generated');args=parser.parse_args()
    if args.mode=='results':results();return
    args.out.mkdir(parents=True,exist_ok=True)
    if args.mode=='forward':forward(args.out)
    else:inverse(args.out,args.seed)
    print('Saved outputs:',args.out.resolve())

if __name__=='__main__':main()
