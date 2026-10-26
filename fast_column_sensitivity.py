"""Langmuir-LDF forward sensitivities with an analytic sparse state Jacobian.

Same upwind finite-volume equations as nonlinear_design_pilot.sensitive.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import bmat,diags,eye,csc_matrix
from physics_model import Column

def system(logp,feed,n):
    k,prod,kl,eps=np.exp(logp);p=Column();A=p.duration_min*p.flow_ml_min*1e-6/(np.pi*(p.diameter_m/2)**2*eps*p.length_m);Pe=A*p.length_m**2/(p.duration_min*p.dispersion_m2_min);beta=1000*.8*p.qref_mg_g/(eps*feed);rate=p.duration_min*k;a=kl*feed;cap=prod/.8/p.qref_mg_g
    adv=A*n;diff=adv*n/Pe;diag=np.full(n,-adv-2*diff);diag[[0,-1]]=-adv-diff
    L=diags([np.full(n-1,adv+diff),diag,np.full(n-1,diff)],[-1,0,1],format='csc');dd=np.full(n,-2*diff);dd[[0,-1]]=-diff;Ld=diags([np.full(n-1,diff),dd,np.full(n-1,diff)],[-1,0,1],format='csc');forcing=np.zeros(n);forcing[0]=adv;I=eye(n,format='csc');zero=csc_matrix((n,n))
    def rhs(t,y):
        state=y.reshape(5,2,n);c,s=state[0];U=state[1:,0];V=state[1:,1];den=1+a*c;E=cap*a*c/den;Ec=cap*a/den**2;ds=rate*(E-s);dc=L@c+forcing-beta*ds
        sv=rate*(U*Ec-V);sv[0]+=ds;sv[1]+=rate*E;sv[2]+=rate*E/den;su=(L@U.T).T-beta*sv;su[3]+=-dc+Ld@c
        result=np.empty_like(state);result[0,0]=dc;result[0,1]=ds;result[1:,0]=su;result[1:,1]=sv;return result.ravel()
    def jac(t,y):
        state=y.reshape(5,2,n);c=state[0,0];U=state[1:,0];den=1+a*c;Ec=cap*a/den**2;Ecc=-2*cap*a*a/den**3;EKc=cap*a*(1-a*c)/den**3;dc=diags(rate*Ec,format='csc');D=bmat([[L-beta*dc,beta*rate*I],[dc,-rate*I]],format='csc');blocks=[[None]*5 for _ in range(5)];blocks[0][0]=D
        for j in range(4):
            base=rate*Ecc*U[j]
            if j in [0,1]:base=base+rate*Ec
            if j==2:base=base+rate*EKc
            vc=diags(base,format='csc');vs=-rate*I if j==0 else zero;uc=-beta*vc;us=-beta*vs
            if j==3:uc=uc-L+beta*dc+Ld;us=us-beta*rate*I
            blocks[j+1][0]=bmat([[uc,us],[vc,vs]],format='csc');blocks[j+1][j+1]=D
        return bmat(blocks,format='csc')
    return rhs,jac

def sensitive(logp,feed,times,n=128):
    rhs,jac=system(logp,feed,n);sol=solve_ivp(rhs,(0,1),np.zeros(10*n),method='BDF',t_eval=times,rtol=2e-8,atol=1e-10,jac=jac)
    if not sol.success:raise RuntimeError(sol.message)
    y=sol.y.reshape(5,2,n,-1);return y[0,0,-1],y[1:,0,-1].T

def checks():
    from dataclasses import replace
    from physics_model import solve
    x=np.log([.001,48.,.04,.4]);times=np.linspace(.005,1,20)
    c,J=sensitive(x,150.,times,32)
    reference=solve(replace(Column(),inlet_mg_l=150.),32,times=times,rtol=2e-9,atol=1e-11)['c'][:,-1]
    curve_error=float(np.max(abs(c-reference)))
    direction=np.array([.2,-.3,.4,.1]);h=1e-4
    plus=sensitive(x+h*direction,150.,times,32)[0]
    minus=sensitive(x-h*direction,150.,times,32)[0]
    fd=(plus-minus)/(2*h);error=float(np.linalg.norm(J@direction-fd)/np.linalg.norm(fd))
    assert curve_error<1e-5 and error<1e-3,(curve_error,error)
    return dict(curve_max_error=curve_error,directional_jacobian_relative_error=error)

if __name__=='__main__':
    import json
    print(json.dumps(checks(),indent=2))
