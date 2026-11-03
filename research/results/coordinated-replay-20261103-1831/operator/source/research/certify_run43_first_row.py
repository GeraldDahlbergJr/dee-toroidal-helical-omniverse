"""Full first-row matrix certificate, not a fitted or reference-gradient budget.

Certifies the exact product-rule decomposition of frozen tensor divergence
minus the direct second-order Cartesian Hessian on EVERY input coefficient.
"""
import sys,json
from pathlib import Path
import numpy as np
from scipy import sparse as sp
from run43_second_hessian import SecondHessianGrid


def maximum(m):return float(np.max(abs(m.data))) if m.nnz else 0.


def main():
    out=Path(sys.argv[1])
    if out.exists():raise FileExistsError(out)
    out.mkdir(parents=True);rows=[]
    for n in (8,16):
        g=SecondHessianGrid(n);er,et,ep=g.basis;r,t,q=g.r,g.t,g.q
        a=np.array([er,et/r,ep/q]).transpose(1,0,2)
        da=np.empty((3,3,3,g.size)) # [chart derivative,Cartesian component,chart component,node]
        da[0]=np.array([np.zeros_like(er),-et/r**2,-np.cos(t)*ep/q**2]).transpose(1,0,2)
        da[1]=np.array([et,-er/r,r*np.sin(t)*ep/q**2]).transpose(1,0,2)
        da[2]=np.array([np.cos(t)*ep,-np.sin(t)*ep/r,(-np.cos(t)*er+np.sin(t)*et)/q]).transpose(1,0,2)
        idx=np.arange(g.size).reshape(g.shape)[1].ravel();diag=sp.diags
        # Direct pure chart D2 are recovered algebraically from stored Hessians.
        q2=[g.components[0],diag(r*r)@(g.components[1]-diag(1/r)@g.d[0]),diag(q*q)@(g.components[2]-diag(np.cos(t)/q)@g.d[0]+diag(np.sin(t)/(r*q))@g.d[1])]
        B=[g.d[k]@g.d[k]-q2[k] for k in range(3)]
        K={};principal={};leibniz={}
        for i in range(3):
            for j in range(3):
                principal[i,j]=[diag(a[i,k,idx]*a[j,k,idx])@B[k][idx] for k in range(3)]
                L=sp.csr_matrix((len(idx),g.size))
                for alpha in range(3):
                    for beta in range(3):
                        coef=a[j,beta]
                        E=g.d[alpha][idx]@diag(coef)-diag(coef[idx])@g.d[alpha][idx]-sp.csr_matrix((da[alpha,j,beta,idx],(np.arange(len(idx)),idx)),shape=(len(idx),g.size))
                        L+=diag(a[i,alpha,idx])@E@g.d[beta]
                leibniz[i,j]=L;K[i,j]=sum(principal[i,j])+L
        blocks=[];radial=[];angular=[];product=[];directblocks=[]
        for i in range(3):
            blockrow=[];rr=[];aa=[];ll=[];dd=[]
            for j in range(3):
                hij=sum(diag(c[idx])@op[idx] for c,op in zip([er[i]*er[j],et[i]*et[j],ep[i]*ep[j],er[i]*et[j]+et[i]*er[j],er[i]*ep[j]+ep[i]*er[j],et[i]*ep[j]+ep[i]*et[j]],g.components))
                lap=sum(g.components[:3])[idx]
                H=(lap if i==j else sp.csr_matrix(lap.shape))+hij/3
                T=(sum(g.grad[k][idx]@g.grad[k] for k in range(3)) if i==j else sp.csr_matrix(lap.shape))+g.grad[j][idx]@g.grad[i]-2/3*(g.grad[i][idx]@g.grad[j])
                blockrow.append(T-H);dd.append(H)
                def tensor_part(parts):return (sum(parts[k,k] for k in range(3)) if i==j else sp.csr_matrix(lap.shape))+parts[j,i]-2/3*parts[i,j]
                rr.append(tensor_part({key:val[0] for key,val in principal.items()}))
                aa.append(tensor_part({key:val[1]+val[2] for key,val in principal.items()}))
                ll.append(tensor_part(leibniz))
            blocks.append(blockrow);radial.append(rr);angular.append(aa);product.append(ll);directblocks.append(dd)
        difference=sp.bmat(blocks,format='csr');rad=sp.bmat(radial,format='csr');ang=sp.bmat(angular,format='csr');lei=sp.bmat(product,format='csr')
        err=difference-rad-ang-lei;scale=maximum(difference)
        assert maximum(err)<2e-10*max(1,scale)
        # The radial part at row 1 is exactly (I+er er^T/3) h delta3/4.
        d3=sp.lil_matrix((len(idx),g.size));h=.2/n
        for offset,value in enumerate((-1,3,-3,1)):
            d3[np.arange(len(idx)),idx+(offset-1)*n*n]=value/(4*h*h)
        d3=d3.tocsr()
        formula=sp.bmat([[diag((np.ones(len(idx)) if i==j else np.zeros(len(idx)))+er[i,idx]*er[j,idx]/3)@d3 for j in range(3)] for i in range(3)],format='csr')
        assert maximum(rad-formula)<1e-8
        for name,m in [('difference',difference),('radial',rad),('angular',ang),('leibniz',lei)]:sp.save_npz(out/f'{name}_{n}.npz',m)
        rows.append(dict(n=n,shape=list(difference.shape),max_coefficient=scale,identity_max_abs=maximum(err),identity_relative=maximum(err)/scale,radial_formula_max_abs=maximum(rad-formula),radial_max=maximum(rad),angular_max=maximum(ang),leibniz_max=maximum(lei)))
    (out/'certificate.json').write_text(json.dumps(dict(identity='T-H = radial pure-composition defect + angular pure-composition defect + exact discrete Leibniz defects',rows=rows),indent=2)+'\n')
    print(json.dumps(rows,indent=2))


if __name__=='__main__':main()
