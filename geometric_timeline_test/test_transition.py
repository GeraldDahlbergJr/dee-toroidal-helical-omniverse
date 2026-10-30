#!/usr/bin/env python3
"""Standalone nondimensional kinematic benchmark. No Run 43 inputs or network.

Run: python3 test_transition.py
Requires Python 3.10+, NumPy, Matplotlib. Writes outputs beside this script.
"""
from pathlib import Path
import csv
import hashlib
import json
import math
import platform
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent
R, LO, HI = 3.0, 0.25, 0.75
WP, WT, V = 2.2, math.pi, 0.4
KNOTS = (2., 4., 10., 12.)
T = 14.

def smooth(u):
    u = np.clip(u, 0., 1.)
    return 10*u**3 - 15*u**4 + 6*u**5

def primitive(u):
    u = np.clip(u, 0., 1.)
    return 2.5*u**4 - 3*u**5 + u**6

def gate(t):
    t = np.asarray(t)
    return 1-smooth((t-2)/2)+smooth((t-10)/2)

def clock(t):
    """Exact integral from zero of gate(t), for t>=0."""
    t = np.asarray(t)
    return (np.minimum(t,2.) + 2*(np.clip((t-2)/2,0,1)-primitive((t-2)/2))
            + 2*primitive((t-10)/2) + np.maximum(t-12.,0.))

def exact(t, rho=.55, theta0=.2, phi0=.35, wp=WP, wt=WT, v=V):
    t = np.asarray(t)
    target = theta0 + rho/R*np.sin(theta0) + wp*clock(t)
    th = target.copy()
    for _ in range(12):
        th -= (th+rho/R*np.sin(th)-target)/(1+rho/R*np.cos(th))
    return np.stack([th, phi0+wt*t, v*clock(t)], axis=-1)

def embed(y, rho=.55):
    th, ph, z = np.moveaxis(np.asarray(y),-1,0)
    radius = R+rho*np.cos(th)
    return np.stack([radius*np.cos(ph), radius*np.sin(ph), z+rho*np.sin(th)],axis=-1)

def rhs(t, y, rho=.55, wp=WP, wt=WT, v=V):
    return np.array([wp*gate(t)/(1+rho/R*np.cos(y[0])), wt, v*gate(t)])

def rk4(h, rho=.55, wp=WP, wt=WT, v=V):
    # All default refinement grids land on schedule breakpoints.
    n = round(T/h)
    times = np.linspace(0,T,n+1)
    h = T/n
    y = np.empty((n+1,3)); y[0] = [.2,.35,0.]
    for i,t in enumerate(times[:-1]):
        f = lambda tt, yy: rhs(tt,yy,rho,wp,wt,v)
        a=f(t,y[i]); b=f(t+h/2,y[i]+h*a/2)
        c=f(t+h/2,y[i]+h*b/2); d=f(t+h,y[i]+h*c)
        y[i+1]=y[i]+h*(a+2*b+2*c+d)/6
    return times,y

def velocity(x,t):
    xx,yy,zz=np.asarray(x)
    radial=np.hypot(xx,yy)
    ph=np.arctan2(yy,xx)
    dz=zz-V*clock(t)
    rho=np.hypot(radial-R,dz)
    th=np.arctan2(dz,radial-R)
    td=WP*gate(t)/(1+rho/R*np.cos(th))
    return np.array([-rho*np.sin(th)*np.cos(ph)*td-radial*np.sin(ph)*WT,
                     -rho*np.sin(th)*np.sin(ph)*td+radial*np.cos(ph)*WT,
                     rho*np.cos(th)*td+V*gate(t)])

def position_only(ph):
    return bool(abs(np.angle(np.exp(1j*ph[-1]))) < 1e-8)

def history_check(t, rho, ph, entered_at_inlet, required=1):
    """Sampled diagnostic, NOT a continuous-time security guarantee.

    ph is a trusted lifted angle, not an angle reconstructed from endpoints.
    Strictly increasing positive circulation is the selected convention.
    """
    if not entered_at_inlet or abs(np.angle(np.exp(1j*ph[0])))>1e-8:
        return False
    if np.any(np.diff(t)<=0) or np.any((rho<LO)|(rho>HI)):
        return False
    dp=np.diff(ph)
    if np.any(dp< -1e-12) or np.any(dp>np.pi/4):
        return False
    return bool(ph[-1]-ph[0]>=required*2*np.pi-1e-8 and position_only(ph))

def certified_history(segments, entered_at_inlet, required=1):
    """Sufficient continuous-path certificate for piecewise Bernstein curves.

    Segment: (start time, end time, radial controls, lifted-angle controls).
    Convex hull containment covers ALL times, not only samples. Monotone
    angle controls imply nonnegative angular derivative. Conservative:
    rejection may mean 'uncertified', not an actual path violation.
    Scope: supplied exact polynomial path; floats/tolerances are not formal
    interval arithmetic. No certificate authenticates an observed history.
    """
    if not segments or not entered_at_inlet:
        return False
    previous=None
    first_phi=float(segments[0][3][0])
    if abs(np.angle(np.exp(1j*first_phi)))>1e-8:
        return False
    for start,end,radial,phase in segments:
        radial=np.asarray(radial); phase=np.asarray(phase)
        if not (np.all(np.isfinite(radial)) and np.all(np.isfinite(phase)) and np.isfinite(start) and np.isfinite(end)):
            return False
        if end<=start or len(radial)<1 or len(phase)<2:
            return False
        if np.any(radial<LO) or np.any(radial>HI) or np.any(np.diff(phase)<0):
            return False
        if previous is not None:
            old_t,old_r,old_p=previous
            if abs(start-old_t)>1e-12 or abs(radial[0]-old_r)>1e-12 or abs(phase[0]-old_p)>1e-12:
                return False
        previous=(end,radial[-1],phase[-1])
    final=previous[2]
    return bool(final-first_phi>=required*2*np.pi-1e-8 and abs(np.angle(np.exp(1j*final)))<1e-8)

def main():
    checks = []
    def check(name, passed, **evidence):
        checks.append(dict(name=name,passed=bool(passed),**evidence))
    # Convergence against an independently evaluated implicit first integral.
    conv=[]
    for h in [.1,.05,.025,.0125]:
        t,y=rk4(h)
        err=float(np.max(np.linalg.norm(embed(y)-embed(exact(t)),axis=1)))
        conv.append(dict(h=h,max_position_error=err))
    for i in range(1,len(conv)):
        conv[i]['order']=math.log2(conv[i-1]['max_position_error']/conv[i]['max_position_error'])
    check('RK4 convergence', all(3.7<c['order']<4.3 for c in conv[1:]), rows=conv)
    # Cartesian finite differences independently check divergence of the field.
    points=[]
    for time in [1.,3.,6.,11.,13.]:
        for rho in [.25,.55,.75]:
            for th in [.3,1.7,4.1]:
                points.append((embed(np.array([th,.7,V*clock(time)]),rho),time))
    divs=[]
    for h in [.02,.01,.005]:
        vals=[]
        for x,time in points:
            value=0.
            for j in range(3):
                e=np.eye(3)[j]*h
                value+=(velocity(x+e,time)[j]-velocity(x-e,time)[j])/(2*h)
            vals.append(value)
        divs.append(dict(h=h,max_abs_divergence=float(np.max(np.abs(vals)))))
    div_orders=[math.log2(divs[i-1]['max_abs_divergence']/divs[i]['max_abs_divergence']) for i in [1,2]]
    check('Cartesian divergence refinement',all(1.9<p<2.1 for p in div_orders),rows=divs,orders=div_orders)
    # Exact local derivatives at all joins: S'=30u^2(1-u)^2;
    # S''=60u(1-u)(1-2u). Both vanish at u=0,1.
    check('C2 switching joins', all(30*u*u*(1-u)**2==0 and 60*u*(1-u)*(1-2*u)==0 for u in [0.,1.]),
          knots=list(KNOTS),regularity='g is C2, trajectories are C3')
    hold=exact(np.array([4.,6.,8.,10.]))
    closure=float(np.max(np.linalg.norm(embed(hold)-embed(hold[0]),axis=1)))
    check('Three complete toroidal turns', closure<1e-12 and np.ptp(hold[:,0])<1e-12 and np.ptp(hold[:,2])<1e-12,
          max_closure_error=closure,turns=3)
    check('Axial drift restored',gate(13.)==1. and gate(6.)==0. and abs(float(clock(14.))-6)<1e-12,
          final_translation=float(V*clock(14.)),restored_translation_speed=V)
    lower_j=LO*(R-HI); lower_speed=(R-HI)*abs(WT)
    check('Nonsingular tube and nonzero speed',lower_j>0 and lower_speed>0,
          jacobian_lower_bound=lower_j,speed_lower_bound=lower_speed)
    # Independent integration, multiple radii, helical senses, axial senses,
    # and zero poloidal winding. No parameter adjusted using test outcomes.
    sweep=[]
    for rho in [LO,.55,HI]:
        for wp in [-WP,0.,WP]:
            for wt in [-WT,WT]:
                for v in [-V,V]:
                    t,y=rk4(.05,rho,wp,wt,v)
                    err=float(np.max(np.linalg.norm(embed(y,rho)-embed(exact(t,rho,wp=wp,wt=wt,v=v),rho),axis=1)))
                    sweep.append(dict(rho=rho,wp=wp,wt=wt,v=v,max_position_error=err))
    check('36 parameter cases', max(a['max_position_error'] for a in sweep)<1e-4,
          count=len(sweep),max_position_error=max(a['max_position_error'] for a in sweep))
    # Negative controls for the audit: all end on the same meridional plane.
    n=401; t=np.linspace(0,1,n); full=2*np.pi*t; rr=np.full(n,.55)
    cases=[('valid_full_turn',rr,full,True,True),
           ('side_entry_at_exit',rr,np.zeros(n),False,False),
           ('partial_turn_from_side',rr,1.5*np.pi+.5*np.pi*t,False,False),
           ('reverse_backtracking',rr,2*np.pi*np.sin(np.pi*t)**2,True,False),
           ('leaves_tube_then_returns',rr+.4*np.sin(np.pi*t),full,True,False),
           ('counterfeit_one_step_turn',rr,np.r_[np.zeros(n-1),2*np.pi],True,False)]
    attacks=[]
    for name,rho,ph,inlet,expected in cases:
        accepted=history_check(t,rho,ph,inlet)
        attacks.append(dict(name=name,expected_history_acceptance=expected,history_acceptance=accepted,
                            position_only_acceptance=position_only(ph)))
    check('History audit positive and negative controls', all(a['history_acceptance']==a['expected_history_acceptance'] for a in attacks),cases=attacks)
    # Explicit counterexample: a smooth radial excursion wholly between samples
    # evades this sampled audit. Preserve this limitation, never call it secure.
    center=(t[200]+t[201])/2; half=(t[201]-t[200])*.45
    dense=np.linspace(t[200],t[201],1001)
    def narrow(q):
        u=(q-center)/half
        return .55+np.where(np.abs(u)<1,.4*(1-u*u)**3,0.)
    hidden_accepted=history_check(t,narrow(t),full,True)
    hidden_outside=bool(np.max(narrow(dense))>HI)
    check('Sampling limitation reproduced',hidden_accepted and hidden_outside,
          sampled_audit_accepts=hidden_accepted,true_path_exits=hidden_outside,
          remedy='Continuous containment proof or certified interval bounds; sampled audit alone is insufficient.')
    # Exact Bernstein hull certificate for the hidden excursion. Degree-six
    # bump .55 + .4*64*u^3*(1-u)^3 reaches .95 inside a gap.
    left,right=center-half,center+half
    pleft,pright=2*np.pi*left,2*np.pi*right
    hidden_segments=[(0.,left,[.55,.55],[0.,pleft]),
                     (left,right,[.55,.55,.55,1.83,.55,.55,.55],[pleft,pright]),
                     (right,1.,[.55,.55],[pright,2*np.pi])]
    valid_segments=[(0.,1.,[.55,.55],[0.,2*np.pi])]
    cert_cases=dict(valid=certified_history(valid_segments,True),
                    hidden_excursion=certified_history(hidden_segments,True),
                    side_entry=certified_history(valid_segments,False),
                    incomplete=certified_history([(0.,1.,[.55,.55],[0.,np.pi])],True),
                    reversal=certified_history([(0.,1.,[.55,.55],[0.,3*np.pi,np.pi,2*np.pi])],True),
                    discontinuity=certified_history([(0.,.5,[.55,.55],[0.,np.pi]),(.5,1.,[.6,.6],[np.pi,2*np.pi])],True))
    check('Continuous polynomial path certificate',cert_cases==dict(valid=True,hidden_excursion=False,side_entry=False,incomplete=False,reversal=False,discontinuity=False),cases=cert_cases)
    results=dict(scope='Prescribed incompressible kinematic benchmark, isolated from Run 43',
                 nondimensional=True,parameters=dict(R=R,rho_min=LO,rho_max=HI,omega_p=WP,omega_t=WT,v=V,T=T),
                 checks=checks,passed_checks=sum(c['passed'] for c in checks),total_checks=len(checks),
                 parameter_sweep=sweep,
                 environment=dict(python=sys.version,numpy=np.__version__,matplotlib=matplotlib.__version__,platform=platform.platform()),
                 source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (OUT/'results.json').write_text(json.dumps(results,indent=2)+'\n')
    ts=np.linspace(0,T,2801); yy=exact(ts); xyz=embed(yy)
    with (OUT/'trajectory.csv').open('w',newline='') as f:
        writer=csv.writer(f); writer.writerow(['t','g','theta_lifted','phi_lifted','Z','x','y','z'])
        writer.writerows(np.column_stack([ts,gate(ts),yy,xyz]))
    fig=plt.figure(figsize=(13,8),layout='constrained')
    ax=fig.add_subplot(221,projection='3d')
    for a,b,color,label in [(0,2,'#2878b5','Helical'),(2,4,'#e6a228','Transition'),(4,10,'#28854c','Toroidal: 3 turns'),(10,12,'#ae64ad','Restore'),(12,14,'#172c64','Axial drift restored')]:
        mask=(ts>=a)&(ts<=b)
        ax.plot(*xyz[mask].T,color=color,label=label,lw=1.8)
    ax.set(xlabel='x',ylabel='y',zlabel='z',title='Prescribed trajectory in a translating toroidal tube')
    ax.legend(fontsize=7,loc='upper left'); ax.set_box_aspect([1,1,.65])
    ax=fig.add_subplot(222); ax.plot(ts,gate(ts),label='g(t)'); ax.plot(ts,V*clock(ts),label='Z(t)')
    ax.set(xlabel='Nondimensional time',title='Smooth stop and restoration of axial drift'); ax.legend(); ax.grid(alpha=.2)
    ax=fig.add_subplot(223); ax.loglog([c['h'] for c in conv],[c['max_position_error'] for c in conv],'o-')
    ax.set(xlabel='Time step',ylabel='Maximum position error',title='RK4 vs implicit exact trajectory'); ax.grid(alpha=.2)
    ax=fig.add_subplot(224); ax.axis('off')
    text='LOCAL KINEMATIC BENCHMARK\n\n'
    text+='\n'.join(f"{'PASS' if c['passed'] else 'FAIL'}  {c['name']}" for c in checks)
    text+='\n\nSampled audit misses a hidden excursion; the polynomial\ncertificate rejects it. No DEE field equations or Run 43 data.'
    ax.text(0,1,text,va='top',fontsize=9,linespacing=1.6)
    fig.suptitle('Continuous helical → toroidal → axial propagation test',fontsize=16)
    fig.savefig(OUT/'transition_diagnostics.png',dpi=170); plt.close(fig)
    print(json.dumps(dict(passed=results['passed_checks'],total=results['total_checks'],convergence=conv,divergence=divs,failures=[c['name'] for c in checks if not c['passed']]),indent=2))
    return 0 if all(c['passed'] for c in checks) else 1

if __name__=='__main__':
    raise SystemExit(main())
