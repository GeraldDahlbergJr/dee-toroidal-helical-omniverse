import os,resource,time,json,traceback
# Enforce a single logical CPU and 4 GiB address-space ceiling before imports.
os.sched_setaffinity(0,{min(os.sched_getaffinity(0))})
resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
import experiment as ex
from evaluator import evaluate,traces
import numpy as np
out={'base':'5a83752a0e8df98e9db5c25041aebcd3e87ed6f9','status':'INCONCLUSIVE','rows':[],'stop':None}
path=ex.ROOT/'results.json'
def save():path.write_text(json.dumps(out,indent=2)+'\n')
try:
    for n in (8,16,32):
        space=ex.Space(n,6);c,rels=ex.solve(space)
        np.savez_compressed(ex.ROOT/f'coefficients_{n}.npz',coefficients=c,n=n)
        row={'n':n,'linear_relative_max':max(rels),'traces':traces(c,n)}
        del space
        for k in (6,8):
            row[f'evaluator_q{k}']=evaluate(c,n,k)
        row['elapsed_s']=time.monotonic()-ex.START;row['peak_rss_KiB']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        out['rows'].append(row);save();print(json.dumps(row),flush=True)
        tr=row['traces']
        if np.any(np.array(tr['radial_max'])>tr['tolerances']) or np.any(np.array(tr['seam_max'])>tr['tolerances']):
            out['status']='FAILED';out['stop']='trace criterion';break
        if len(out['rows'])>1:
            prev=out['rows'][-2];orders={};failure=False
            for region,err in row['evaluator_q8']['residual_rms'].items():
                ec=np.array(prev['evaluator_q8']['residual_rms'][region]);ef=np.array(err)
                uc=2*abs(ec-np.array(prev['evaluator_q6']['residual_rms'][region]));uf=2*abs(ef-np.array(row['evaluator_q6']['residual_rms'][region]))
                lower=np.log2((ec-uc)/(ef+uf));orders[region]=lower.tolist()
                failure |= bool(np.any(lower<1.5) or np.any(~np.isfinite(lower)))
            row['conservative_orders_quadrature_only']=orders
            if failure:out['status']='FAILED';out['stop']='independent convergence criterion';save();break
        save()
    else:out['status']='INCONCLUSIVE';out['stop']='all grids complete; final uncertainty and quadrature review required'
except Exception as e:
    out['stop']=str(e);out['exception']=traceback.format_exc();print(out['exception'],flush=True)
finally:
    out['elapsed_s']=time.monotonic()-ex.START;save()
