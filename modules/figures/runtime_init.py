"""Locate native libraries inside the active Python environment on Windows."""
import os,sys
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
env=Path(sys.executable).parent
os.environ['PATH']=str(env/'Library/bin')+os.pathsep+str(env)+os.pathsep+os.environ.get('PATH','')
DLL_HANDLES=[os.add_dll_directory(str(p)) for p in (env,env/'Library/bin',env/'DLLs') if os.name=='nt' and p.exists()]
