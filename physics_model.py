"""Conservative, dimensionally explicit Langmuir-LDF column model.
Synthetic assumptions. No radioactive decay, competing ions, or experimental calibration.
"""
from dataclasses import dataclass, asdict, replace
import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import lil_matrix

@dataclass(frozen=True)
class Column:
    length_m: float = .06
    diameter_m: float = .01
    flow_ml_min: float = .157
    porosity: float = .4
    bulk_density_kg_l: float = .8
    inlet_mg_l: float = 100.
    qref_mg_g: float = 48.1
    duration_min: float = 24000.
    dispersion_m2_min: float = 7.496197819628419e-6
    k_min: float = .001
    qmax_mg_g: float = 60.
    KL_l_mg: float = .04

    @property
    def velocity(self): return self.flow_ml_min*1e-6/(np.pi*(self.diameter_m/2)**2*self.porosity)
    @property
    def A(self): return self.duration_min*self.velocity/self.length_m
    @property
    def Pe(self): return self.velocity*self.length_m/self.dispersion_m2_min
    @property
    def beta(self): return self.bulk_density_kg_l*1000*self.qref_mg_g/(self.porosity*self.inlet_mg_l)

def solve(column=Column(),n=128,times=None,rtol=2e-7,atol=1e-9):
    p=column
    if n<3 or min(p.porosity,1-p.porosity,p.bulk_density_kg_l,p.dispersion_m2_min,p.k_min,p.qmax_mg_g,p.KL_l_mg)<=0:
        raise ValueError('Invalid model/grid parameters')
    times=np.linspace(0,1,201) if times is None else np.asarray(times)
    def rhs(t,y):
        c,s=y[:n],y[n:2*n]
        seq=(p.qmax_mg_g/p.qref_mg_g)*(p.KL_l_mg*p.inlet_mg_l*c)/(1+p.KL_l_mg*p.inlet_mg_l*c)
        ds=p.duration_min*p.k_min*(seq-s)
        flux=np.empty(n+1);flux[0]=1;flux[-1]=c[-1]
        flux[1:n]=c[:-1]-n/p.Pe*(c[1:]-c[:-1])
        return np.r_[-p.A*n*np.diff(flux)-p.beta*ds,ds,p.A*(1-c[-1])]
    jac=lil_matrix((2*n+1,2*n+1),dtype=int)
    for i in range(n):
        for j in range(max(0,i-1),min(n,i+2)): jac[i,j]=1
        jac[i,n+i]=1;jac[n+i,i]=1;jac[n+i,n+i]=1
    jac[-1,n-1]=1
    sol=solve_ivp(rhs,(0,1),np.zeros(2*n+1),method='BDF',t_eval=times,
                  rtol=rtol,atol=atol,jac_sparsity=jac.tocsr())
    if not sol.success: raise RuntimeError(sol.message)
    return dict(t=sol.t,c=sol.y[:n].T,s=sol.y[n:2*n].T,
                integrated_net_flux=sol.y[-1],x=(np.arange(n)+.5)/n,nfev=sol.nfev)
