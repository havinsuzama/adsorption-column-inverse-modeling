"""Portable CPU experiment: balanced early-window PINN warm-start training.

Training needs no SciPy and no reference labels.
Reference fields are read only by the optional --validate evaluation.
"""
import argparse,json,math,time,uuid
from datetime import datetime
from pathlib import Path
try:
    import numpy as np
    import torch
    from torch import nn
except ImportError as exc:
    raise SystemExit('Missing package. From the repository root: python -m pip install -r requirements-pinn.txt\n'+str(exc))

ROOT=Path(__file__).resolve().parent
L=.06;T=24000.;EPS=.4;RHO=.8;C0=100.;QREF=48.1
K=.001;QMAX=60.;KL=.04;H=.01
V=.157e-6/(math.pi*(.01/2)**2*EPS)
A=T*V/L;PE=V*L/7.496197819628419e-6;BETA=RHO*1000*QREF/(EPS*C0)
BINS=((1e-7,1e-5),(1e-5,1e-4),(1e-4,1e-3),(1e-3,H))

def grad(y,z):
    return torch.autograd.grad(y,z,torch.ones_like(y),create_graph=True)[0]

def equilibrium(c):
    return QMAX/QREF*(KL*C0*c)/(1+KL*C0*c)

def envelope(x,tau):
    v=A*H;d=v/PE;t=tau.clamp_min(1e-12)
    u=(x-v*t)/(2*torch.sqrt(d*t));w=(x+v*t)/(2*torch.sqrt(d*t))
    return .5*torch.erfc(u)+torch.exp(-u*u)*(torch.sqrt(v*v*t/(math.pi*d))-.5*(1+v*x/d+v*v*t/d)*torch.special.erfcx(w))

class Onset(nn.Module):
    """float64 onset model; state keys match start_model.pt."""
    def __init__(self,nquad=24):
        super().__init__()
        self.mnet=nn.Sequential(nn.Linear(2,24),nn.Tanh(),nn.Linear(24,1))
        self.hnet=nn.Sequential(nn.Linear(3,32),nn.Tanh(),nn.Linear(32,32),nn.Tanh(),nn.Linear(32,1))
        self.local=nn.Sequential(nn.Linear(4,32),nn.Tanh(),nn.Linear(32,32),nn.Tanh(),nn.Linear(32,1))
        self.set_quad(nquad)
    def set_quad(self,n):
        q,w=np.polynomial.legendre.leggauss(n)
        self.register_buffer('nodes',torch.tensor((q+1)/2,dtype=torch.float64))
        self.register_buffer('weights',torch.tensor(w/2,dtype=torch.float64));self.nquad=n
    def concentration(self,x,tau):
        g=x-x**8/8;logt=torch.log1p(100*tau)/math.log(101)
        m=torch.nn.functional.softplus(self.mnet(torch.cat([tau,logt],1)))*torch.sqrt(tau.clamp_min(1e-16)/(tau+.001))
        root=torch.sqrt(tau.clamp_min(1e-16));xi=g/(g+torch.sqrt((A*H/PE)*tau.clamp_min(1e-16)))
        h=torch.nn.functional.softplus(self.hnet(torch.cat([g,tau,logt],1))+self.local(torch.cat([g,tau,root,xi],1)))
        base=envelope(g,tau);b0=envelope(torch.zeros_like(g),tau)
        return torch.where(tau>0,base*torch.exp(-m*g-g*g*h)/(1+m*b0/PE),torch.zeros_like(base))
    def forward(self,z):
        x,tau=z[:,:1],z[:,1:]/H;c=self.concentration(x,tau)
        u=tau*self.nodes.reshape(1,-1)
        cu=self.concentration(x.expand_as(u).reshape(-1,1),u.reshape(-1,1)).reshape(len(z),-1)
        kt=T*K*H;s=kt*tau*torch.sum(self.weights*torch.exp(-kt*(tau-u))*equilibrium(cu),dim=1,keepdim=True)
        return torch.cat([c,s],1)

def residuals(model,z,exact=False):
    c,s=model(z).split(1,1);gc=grad(c,z);cxx=grad(gc[:,:1],z)[:,:1]
    analytic=T*K*(equilibrium(c)-s)
    st=grad(s,z)[:,1:] if exact else analytic
    return gc[:,1:]/A+gc[:,:1]-cxx/PE+BETA/A*st,st/(T*K)-(equilibrium(c)-s)

def balanced_points(n,generator):
    """Equal counts/log-time measure in each bin; half uniform x, half front."""
    groups=[]
    for lo,hi in BINS:
        u=torch.rand(n,2,dtype=torch.float64,generator=generator)
        theta=10**(math.log10(lo)+(math.log10(hi)-math.log10(lo))*u[:,1:])
        tau=theta/H;extent=(3*torch.sqrt(A*H/PE*tau)+A*H*tau).clamp(max=1)
        x=u[:,:1].clone();x[n//2:]*=extent[n//2:]
        groups.append(torch.cat([x,theta],1))
    return torch.cat(groups).requires_grad_(True)

def loss_bins(model,z,n):
    r,_=residuals(model,z);return (r.reshape(4,n)**2).mean(1)

def evaluate_bins(model,n):
    """Independent fixed Sobol uniform-x/log-time points, not training points."""
    rows=[]
    for i,(lo,hi) in enumerate(BINS):
        u=torch.quasirandom.SobolEngine(2,scramble=True,seed=320041+i).draw(n).double()
        u[:,1]=10**(math.log10(lo)+(math.log10(hi)-math.log10(lo))*u[:,1])
        sums=[];maxima=[]
        for part in u.split(128):
            z=part.clone().requires_grad_(True);r,k=residuals(model,z,exact=True)
            sums.append(torch.cat([r,k],1).detach());maxima.append(float(r.abs().max().detach()))
        a=torch.cat(sums)
        rows.append(dict(theta_lo=lo,theta_hi=hi,mass_rms=float(a[:,0].square().mean().sqrt()),kinetic_rms=float(a[:,1].square().mean().sqrt()),mass_max_abs=max(maxima)))
    return rows

def predict(model,z):
    with torch.no_grad():return np.concatenate([model(torch.tensor(p,dtype=torch.float64)).numpy() for p in np.array_split(z,max(1,int(np.ceil(len(z)/1024))))])

def validate(model):
    """Earlier early-window gate + explicit dense local check, never joined validation."""
    torch.manual_seed(7654321);z=torch.rand(4096,2,dtype=torch.float64);z[:,1]=1e-6+(H-1e-6)*z[:,1]
    a=[]
    for p in z.split(128):
        r,k=residuals(model,p.clone().requires_grad_(True),True);a.append(torch.cat([r,k],1).detach())
    a=torch.cat(a);m=dict(uniform_mass_rms=float(a[:,0].square().mean().sqrt()),uniform_kinetic_rms=float(a[:,1].square().mean().sqrt()))
    t=torch.linspace(1e-6,H,256,dtype=torch.float64).reshape(-1,1)
    zi=torch.cat([torch.zeros_like(t),t],1).requires_grad_(True);c=model(zi)[:,:1];ri=c-grad(c,zi)[:,:1]/PE-1
    zo=torch.cat([torch.ones_like(t),t],1).requires_grad_(True);ro=grad(model(zo)[:,:1],zo)[:,:1]
    m.update(inlet_rms=float(ri.square().mean().sqrt().detach()),outlet_rms=float(ro.square().mean().sqrt().detach()))
    ref=np.load(ROOT/'reference_early.npz');tt,xx=np.meshgrid(ref['t'],ref['x'],indexing='ij');y=predict(model,np.c_[xx.ravel(),tt.ravel()]).reshape(len(ref['t']),len(ref['x']),2)
    m.update(c_rmse=float(np.sqrt(np.mean((y[:,:,0]-ref['ref_c'])**2))),s_rmse=float(np.sqrt(np.mean((y[:,:,1]-ref['ref_s'])**2))),c_transfer_rmse=float(np.sqrt(np.mean((y[-1,:,0]-ref['ref_c'][-1])**2))))
    gen=torch.Generator().manual_seed(719);q=balanced_points(256,gen).detach()
    with torch.no_grad():
        s24=model(q)[:,1];model.set_quad(48);s48=model(q)[:,1];model.set_quad(24)
    m['quad24_48_max']=float((s24-s48).abs().max())
    m['early_window_gate']=bool(m['uniform_mass_rms']<.05 and m['uniform_kinetic_rms']<.03 and m['inlet_rms']<1e-5 and m['outlet_rms']<1e-5 and m['c_rmse']<.03 and m['s_rmse']<.03 and m['quad24_48_max']<.001)
    dx=np.linspace(0,.15,257);dt=np.geomspace(1e-5,1e-4,41);tt,xx=np.meshgrid(dt,dx,indexing='ij');p=np.c_[xx.ravel(),tt.ravel()];rs=[]
    for part in np.array_split(p,int(np.ceil(len(p)/128))):
        r,_=residuals(model,torch.tensor(part,dtype=torch.float64,requires_grad=True),True);rs.append(r.detach().numpy().ravel())
    r=np.concatenate(rs);m.update(dense_front_mass_rms=float(np.sqrt(np.mean(r*r))),dense_front_mass_max=float(abs(r).max()))
    m['scope']='One synthetic early-window case. Fixed thresholds, dense local diagnostic; no late-window, multi-seed or experimental validation. Quadrature check uses new balanced test points.'
    return m

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=['smoke','quick','standard'],default='quick')
    parser.add_argument('--steps',type=int);parser.add_argument('--points-per-bin',type=int);parser.add_argument('--lbfgs',type=int)
    parser.add_argument('--seed',type=int,default=20261020);parser.add_argument('--threads',type=int,default=2)
    parser.add_argument('--validate',action='store_true');parser.add_argument('--no-plot',action='store_true')
    parser.add_argument('--out',type=Path,default=ROOT/'runs')
    args=parser.parse_args();defaults={'smoke':(2,8,32,0),'quick':(80,16,256,0),'standard':(600,32,1024,80)}
    steps,n,neval,lbfgs=defaults[args.mode];steps=steps if args.steps is None else args.steps;n=n if args.points_per_bin is None else args.points_per_bin;lbfgs=lbfgs if args.lbfgs is None else args.lbfgs
    if steps<0 or n<4 or lbfgs<0 or args.threads<1:parser.error('need steps/lbfgs >= 0, points-per-bin >= 4, threads >= 1')
    torch.set_num_threads(args.threads);torch.manual_seed(args.seed)
    path=args.out.resolve()/(datetime.now().strftime('%Y%m%d_%H%M%S')+'_'+uuid.uuid4().hex[:6]);path.mkdir(parents=True,exist_ok=False)
    start=time.perf_counter();ck=torch.load(ROOT/'start_model.pt',map_location='cpu',weights_only=True)
    model=Onset(ck['nquad']).double();model.load_state_dict(ck['state_dict'])
    print('Output folder:',path,flush=True);print('CPU/float64; mode:',args.mode,'steps:',steps,'points:',4*n,flush=True)
    before=evaluate_bins(model,neval);print('Initial per-bin mass RMS:',[round(v['mass_rms'],5) for v in before],flush=True)
    gen=torch.Generator().manual_seed(args.seed)
    calibration=balanced_points(max(64,n),torch.Generator().manual_seed(args.seed+990))
    scales=loss_bins(model,calibration,max(64,n)).detach().clamp_min(1e-4)
    # Smooth maximum of per-bin relative MSE, frozen scales, temperature=0.25.
    # This is a balancing hypothesis; it does not impose a hard no-regression constraint.
    def objective(z):
        losses=loss_bins(model,z,n);ratios=losses/scales
        return .25*(torch.logsumexp(ratios/.25,dim=0)-math.log(4)),losses
    opt=torch.optim.Adam(model.parameters(),lr=.0003);history=[];train_start=time.perf_counter()
    for step in range(1,steps+1):
        opt.zero_grad();z=balanced_points(n,gen);loss,parts=objective(z);loss.backward(inputs=list(model.parameters()));opt.step()
        if step==1 or step%max(1,steps//5)==0 or step==steps:
            row=dict(step=step,objective=float(loss.detach()),bin_mse=parts.detach().tolist());history.append(row);print('Adam',row,flush=True)
    calls=[0]
    if lbfgs:
        fixed=balanced_points(n,gen);optimizer=torch.optim.LBFGS(model.parameters(),lr=.5,max_iter=lbfgs,max_eval=int(lbfgs*1.5),history_size=30,line_search_fn='strong_wolfe')
        def closure():
            optimizer.zero_grad();fixed.grad=None;loss,_=objective(fixed);loss.backward(inputs=list(model.parameters()));calls[0]+=1
            if calls[0]%20==0:print('LBFGS',calls[0],float(loss.detach()),flush=True)
            return loss
        optimizer.step(closure)
    training_seconds=time.perf_counter()-train_start
    after=evaluate_bins(model,neval);ratios=[b['mass_rms']/max(a['mass_rms'],1e-12) for a,b in zip(before,after)]
    recommend=bool(max(ratios)<=1.05 and sum(ratios)/4<1)
    torch.save(dict(state_dict=model.state_dict(),nquad=model.nquad,architecture='Onset'),path/'candidate.pt')
    validation=validate(model) if args.validate else dict(early_window_gate=None,status='Not evaluated. Quick per-bin diagnostics are not physical validation; use --validate.')
    result=dict(mode=args.mode,seed=args.seed,steps=steps,points_per_bin=n,test_points_per_bin=neval,lbfgs_limit=lbfgs,lbfgs_closures=calls[0],threads=args.threads,
                before=before,after=after,per_bin_rms_ratios=ratios,candidate_meets_development_rule=recommend,
                development_rule='Each bin RMS <=1.05*baseline, mean of four RMS ratios <1. Development-only selection, not blind validation.',
                training_scale_mse=scales.tolist(),loss='Smooth maximum of four relative MSEs; fixed calibration scales floored at 1e-4; temperature .25.',history=history,
                validation=validation,training_seconds=training_seconds,total_seconds=time.perf_counter()-start,
                versions=dict(torch=torch.__version__,numpy=np.__version__),scope='Fixed synthetic parameters, early window only; one warm-start chain; no reference labels in training. start_model.pt is never overwritten.')
    (path/'result.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    if not args.no_plot:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,ax=plt.subplots(figsize=(8,4),layout='constrained');x=np.arange(4)
        ax.plot(x,[b['mass_rms'] for b in before],'o-',label='Start');ax.plot(x,[b['mass_rms'] for b in after],'o-',label='Balanced training')
        ax.set_xticks(x,['0.0024-0.24','0.24-2.4','2.4-24','24-240']);ax.set(xlabel='Time bin (min, log-time sampling)',ylabel='Mass-balance RMS residual',yscale='log');ax.legend()
        fig.savefig(path/'comparison.png',dpi=160);plt.close(fig)
    print('RESULT:',json.dumps(dict(after=[v['mass_rms'] for v in after],rms_ratios=ratios,development_rule_pass=recommend,early_window_gate=validation['early_window_gate'],training_seconds=training_seconds),indent=2),flush=True)
    print('Saved:',path/'result.json',flush=True)

if __name__=='__main__':main()
