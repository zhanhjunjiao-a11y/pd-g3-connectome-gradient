"""Shared local setup; importing this module writes no files."""
import os
import sys
from pathlib import Path
for _key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:
    os.environ[_key]='1'
DLL_HANDLES=[]
if os.name=='nt':
    _env=Path(sys.executable).parent
    os.environ['PATH']=str(_env/'Library/bin')+os.pathsep+str(_env)+os.pathsep+os.environ.get('PATH','')
    DLL_HANDLES=[os.add_dll_directory(str(p)) for p in [_env,_env/'Library/bin',_env/'DLLs'] if p.exists()]
PKG=Path(__file__).resolve().parents[2]
INPUTS=PKG/'inputs_private/external'
VERIFY=PKG/'verification/external'


def sha(path):
    import hashlib
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def compare_frames(got,expected,keys,tolerance=1e-9):
    import numpy as np
    import pandas as pd
    a=got.set_index(keys).sort_index(); b=expected.set_index(keys).sort_index()
    issues=[];errors={}
    if not a.index.equals(b.index):issues.append('row keys differ')
    if list(got.columns)!=list(expected.columns):issues.append('column order differs')
    if not issues:
        for col in b.columns:
            if pd.api.types.is_numeric_dtype(b[col]) and not pd.api.types.is_bool_dtype(b[col]):
                x=pd.to_numeric(a[col],errors='coerce').to_numpy(float);y=b[col].to_numpy(float)
                if not np.array_equal(np.isnan(x),np.isnan(y)):
                    issues.append(col+': missingness differs');continue
                finite=np.isfinite(y)
                err=float(np.max(np.abs(x[finite]-y[finite]))) if finite.any() else 0.
                errors[col]=err
                if not np.allclose(x,y,rtol=tolerance,atol=tolerance,equal_nan=True):issues.append(col+': numerical mismatch')
            else:
                x=a[col].fillna('').astype(str);y=b[col].fillna('').astype(str)
                if not x.equals(y):issues.append(col+': text/boolean mismatch')
    return {'status':'PASS' if not issues else 'FAIL','rows':len(got),'keys':keys,
            'max_abs_numeric_error':max(errors.values(),default=0),'field_errors':errors,'issues':issues}
