"""Periodic 3-D CMC conformal initial data for smooth helical DEE scalar waves.

This is a selected regular seed of the archived three-scalar action, not
an asserted globally wound phase or the missing toroidal spacetime metric.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np

L=8.; G=.005; R=.0+1.7; SIGMA=.65; M=2; NZ=1
A=.30; B=.20; C=.10; LV=4.; V=1.

def deriv(u,h,axis): return (np.roll(u,-1,axis)-np.roll(u,1,axis))/(2*h)
def lap(u,h):
    return sum((np.roll(u,-1,axis)-2*u+np.roll(u,1,axis))/h**2 for axis in range(3))

def fields(n):
    h=L/n; axis=(np.arange(n)-n//2)*h
    x,y,z=np.meshgrid(axis,axis,axis,indexing='ij')
    radial2=x*x+y*y
    # Smooth Cartesian polynomial for azimuthal mode m=2; no atan2 branch cut.
    ring=np.exp(-((radial2-R*R)/(2*R*SIGMA))**2)
    axial=np.exp(-.5*(L/(np.pi*SIGMA))**2*np.sin(np.pi*z/L)**2)
    envelope=ring*axial
    wave=((x+1j*y)/R)**M*np.exp(1j*2*np.pi*NZ*z/L)
    rho=1-.025*envelope
    theta=.035*envelope*wave.real
    psi=.028*envelope*(wave*np.exp(.4j)).imag
    return np.stack((rho,theta,psi),axis=-1)

def energy(f,conformal,h):
    grad=np.stack([deriv(f,h,axis) for axis in range(3)],axis=-1)
    rho=f[...,0]; kt=1+A*rho*rho; kp=1+B*rho*rho; mix=C*rho*rho
    norms=np.sum(grad*grad,axis=-1)/conformal[...,None]**4
    cross=np.sum(grad[...,1,:]*grad[...,2,:],axis=-1)/conformal**4
    potential=LV/4*(rho*rho-V*V)**2
    return (potential+.5*(norms[...,0]+kt*norms[...,1]+kp*norms[...,2])
            +mix*cross)

def solve(n,tol=5e-13):
    h=L/n; f=fields(n); psi=np.ones((n,n,n));
    freq=np.fft.fftfreq(n)
    kx,ky,kz=np.meshgrid(freq,freq,freq,indexing='ij')
    eig=-4/h**2*(np.sin(np.pi*kx)**2+np.sin(np.pi*ky)**2+np.sin(np.pi*kz)**2)
    eig[0,0,0]=1
    for iteration in range(1000):
        e=energy(f,psi,h)
        k2=24*np.pi*G*np.mean(psi**5*e)/np.mean(psi**5)
        source=psi**5*(k2/12-2*np.pi*G*e)
        spectrum=np.fft.fftn(source); spectrum[0,0,0]=0
        candidate=1+np.fft.ifftn(spectrum/eig).real
        if np.max(np.abs(candidate-psi))<tol:
            psi=candidate; break
        psi=candidate
    else: raise RuntimeError('conformal constraint iteration did not converge')
    e=energy(f,psi,h)
    k2=24*np.pi*G*np.mean(psi**5*e)/np.mean(psi**5)
    k=-np.sqrt(k2)
    # Conformal identity gives an independent discretized H residual.
    hres=-8*psi**-5*lap(psi,h)+2*k2/3-16*np.pi*G*e
    return f,psi,k,dict(n=n,iterations=iteration+1,
        hamiltonian_rms=float(np.sqrt(np.mean(hres*hres))),
        hamiltonian_max=float(np.max(np.abs(hres))),
        psi_min=float(psi.min()),psi_max=float(psi.max()),
        rho_min=float(f[...,0].min()),kinetic_det_min=float(np.min(
          (1+A*f[...,0]**2)*(1+B*f[...,0]**2)-(C*f[...,0]**2)**2)),
        k=float(k),momentum_max=0.)

def main():
    rows=[solve(n)[3] for n in (16,32,64)]
    result=dict(scope='3-D periodic CMC conformally flat helical-scalar ring; initial constraints only',
                parameters=dict(L=L,G=G,R=R,sigma=SIGMA,m=M,axial_periods=NZ),runs=rows)
    Path(__file__).with_name('toroidal_helical_initial_data_checkpoint.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__': main()
