"""Rebuild all seven current manuscript figures from package-local source data."""
from pathlib import Path
import argparse,os,sys,subprocess,json,time,hashlib

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output-dir',type=Path,required=True)
    args=ap.parse_args();out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
    here=Path(__file__).resolve().parent;pkg=here.parents[1]
    env=os.environ.copy()
    for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:env[name]='1'
    env['GRADIENT_FIG_INPUT']=str(pkg/'inputs_private/figures')
    env['GRADIENT_FIG_OUTPUT']=str(out)
    env['GRADIENT_FIG_CORTEX']=str(out/'cortical_layers')
    env['MPLCONFIGDIR']=str(out/'mplconfig')
    env['PYTHONIOENCODING']='utf-8'
    rows=[]
    for script in ['cortical_layers.py','figure1_build.py','build_figure2.py','build_figure3.py','build_figure4.py','supplements.py','build_supplement3.py']:
        t=time.perf_counter()
        with (out/f'{Path(script).stem}.log').open('w',encoding='utf-8') as log:
            # Bootstrap DLL lookup before importing plotting/scientific packages.
            proc=subprocess.run([sys.executable,'-c',"import sys,runpy;sys.path.insert(0,sys.argv[1]);sys.path.insert(0,sys.argv[3]);import guarded_runner;sys.addaudithook(guarded_runner.audit);import runtime_init;runpy.run_path(sys.argv[2],run_name='__main__')",str(here),str(here/script),str(pkg)],env=env,stdout=log,stderr=subprocess.STDOUT)
        rows.append({'script':script,'exit_code':proc.returncode,'seconds':time.perf_counter()-t})
        print(script,proc.returncode,round(rows[-1]['seconds'],1),flush=True)
        if proc.returncode:
            print((out/f'{Path(script).stem}.log').read_text(encoding='utf-8')[-5000:],flush=True)
            break
    expected=[f'Figure_{i}_revised' for i in range(1,5)]+[f'Supplementary_Figure_{i}' for i in range(1,4)]
    missing=[f'{name}.{ext}' for name in expected for ext in ['pdf','svg','png'] if not (out/'figures'/f'{name}.{ext}').exists()]
    outputs=[{'file':str(p.relative_to(out)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size} for p in (out/'figures').glob('*')] if (out/'figures').exists() else []
    result={'status':'PASS' if len(rows)==7 and not missing and all(r['exit_code']==0 for r in rows) else 'FAIL',
        'scope':'Figure rendering from numeric display inputs; model reproduction is checked by core/connectivity/external modules.',
        'scripts':rows,'missing':missing,'outputs':outputs}
    (out/'verification.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    return 0 if result['status']=='PASS' else 1

if __name__=='__main__':sys.exit(main())
