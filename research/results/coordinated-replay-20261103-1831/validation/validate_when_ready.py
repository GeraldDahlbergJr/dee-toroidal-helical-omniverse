"""Await completed genuine solves, then independently evaluate all arrays."""
from pathlib import Path
import json, os, subprocess, sys, time
here=Path(__file__).parent
source=Path('/workspace/run43-investigation/solver/source-current')
started=time.monotonic(); results=[]
for variant in ('baseline','fourth','local'):
    fields=source/'research/results'/f'live-20261103-{variant}'
    while not (fields/'report.json').is_file():
        exitfile=source/'live-evidence'/f'{variant}-solve.exit'
        if exitfile.is_file() and int(exitfile.read_text()) != 0:
            raise RuntimeError(f'{variant} solve failed; inspect preserved solver log')
        if time.monotonic()-started>3600:
            raise TimeoutError(f'{variant}: no completed all-n report within one hour')
        time.sleep(5)
    command=[sys.executable,str(here/'verify_fresh.py'),str(source),str(fields),str(here/f'{variant}-independent-validation.json')]
    print(json.dumps({'variant':variant,'command':command}),flush=True)
    with (here/f'{variant}-independent-validation.log').open('w') as logfile:
        result=subprocess.run(command,stdout=logfile,stderr=subprocess.STDOUT,env={**os.environ,'OPENBLAS_NUM_THREADS':'1','PYTHONDONTWRITEBYTECODE':'1'})
    results.append({'variant':variant,'exit_code':result.returncode})
    (here/'validation-exit-status.json').write_text(json.dumps(results,indent=2)+'\n')
    if result.returncode:
        raise RuntimeError(f'{variant} validator execution failed; inspect validation log')
    print(json.dumps({'variant':variant,'status':'completed'}),flush=True)
