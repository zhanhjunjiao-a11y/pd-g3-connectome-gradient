"""PPMI finalization: fixed regional gradients and conditional effect-map agreement.

Portable package adaptation; statistical operations and random-number configuration retained.
Run --preflight before resampling. No input, historical result or manuscript is changed.
"""
from pathlib import Path
import os, sys, json, hashlib, time, argparse
for name in ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
             'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS']:
    os.environ[name] = '1'
if os.name == 'nt':
    env_root = Path(sys.executable).parent
    os.environ['PATH'] = str(env_root/'Library/bin')+';'+str(env_root)+';'+os.environ.get('PATH','')
    DLL_HANDLES = [os.add_dll_directory(str(dll))
                   for dll in [env_root,env_root/'Library/bin',env_root/'DLLs'] if dll.exists()]
from common import INPUTS, PKG
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
import pandas as pd
import scipy
from scipy import stats
import statsmodels
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
HERE = PKG / 'runs/external/maps'
REG = INPUTS / 'registry.csv'
GRAD = INPUTS / 'gradients.npy'
DISC = INPUTS / 'discovery_effects.csv'
ORDER = INPUTS / 'subject_order.txt'
N_NULL, N_CI, BATCH = 9999, 5000, 256
NULL_SEED, CI_SEED = 20261030, 20261031
MODEL_NAMES = ['age_sex_motion', 'scanner_adjusted_motion', 'site_adjusted_motion', 'site_shared_motion']

def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()

def serial(value):
    return value.item() if isinstance(value, np.generic) else str(value)

def save_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, default=serial), encoding='utf-8')

def bh(values):
    return multipletests(np.asarray(values), method='fdr_bh')[1]

def load_inputs():
    reg = pd.read_csv(REG, dtype={'SITE_final': 'string', 'sub_id': str})
    assert len(reg) == 438 and reg.sub_id.is_unique
    assert reg.external_row_index.notna().all()
    assert sorted(reg.external_row_index.astype(int)) == list(range(438))
    reg = reg.sort_values('external_row_index').reset_index(drop=True)
    ids = [s.strip() for s in ORDER.read_text().splitlines() if s.strip() and not s.startswith('#')]
    assert reg.sub_id.tolist() == ids
    assert reg.group.tolist() == [0]*79 + [1]*359
    required = ['group', 'age_final', 'age_original', 'sex', 'meanFD_final', 'scanner_exact']
    assert reg[required].notna().all().all()
    for name in ['group', 'age_final', 'age_original', 'sex', 'meanFD_final']:
        assert np.isfinite(reg[name].astype(float)).all()
    assert (reg.meanFD_final >= 0).all()
    assert set(reg.sex.astype(float)) <= {0., 1.}
    reg['SITE_final'] = reg['SITE_final'].str.strip().replace('', pd.NA)
    g = np.load(GRAD)
    assert g.shape == (438, 283, 10) and np.isfinite(g).all()
    disc = pd.read_csv(DISC)
    disc = pd.concat([disc.loc[disc.gradient.eq(f'G{k}')].sort_values('roi_id')
                      for k in [1, 2, 3]], ignore_index=True)
    assert len(disc) == 849
    for k in [1, 2, 3]:
        assert disc.loc[disc.gradient.eq(f'G{k}'), 'roi_id'].tolist() == list(range(1, 284))
    assert np.isfinite(disc.estimate).all()
    y = g[:, :, :3].transpose(0, 2, 1).reshape(438, 849)
    return reg, y, disc

def design(reg, model, allowed_site_file=None):
    mask = np.ones(len(reg), dtype=bool)
    if model in {'site_adjusted_motion', 'site_shared_motion'}:
        mask &= reg.SITE_final.notna().to_numpy()
        if model == 'site_shared_motion':
            assert allowed_site_file is not None, 'Shared-site model requires the explicit approved support rule'
        if allowed_site_file is not None:
            # Only an explicit parent-approved fixed support rule may restrict sites.
            allowed = json.loads(Path(allowed_site_file).read_text(encoding='utf-8'))
            assert isinstance(allowed, dict) and allowed.get('approved_rule')
            mask &= reg.SITE_final.isin([str(v) for v in allowed['allowed_sites']]).to_numpy()
    rr = reg.loc[mask].reset_index(drop=True)
    x = rr[['group', 'age_final', 'sex', 'meanFD_final']].astype(float).copy()
    for column in ['age_final', 'meanFD_final']:
        scale = x[column].std(ddof=0)
        assert scale > 0
        x[column] = (x[column] - x[column].mean()) / scale
    covariates = 'age_final + sex + meanFD_final'
    if model == 'scanner_adjusted_motion':
        x = pd.concat([x, pd.get_dummies(rr.scanner_exact.astype(str), prefix='scanner',
                                      drop_first=True, dtype=float)], axis=1)
        covariates += ' + scanner_exact fixed effects'
    elif model in {'site_adjusted_motion', 'site_shared_motion'}:
        x = pd.concat([x, pd.get_dummies(rr.SITE_final.astype(str), prefix='SITE',
                                      drop_first=True, dtype=float)], axis=1)
        covariates += ' + blinded clinical SITE_final fixed effects'
    return sm.add_constant(x, has_constant='add'), rr, mask, covariates

def inspect_design(xd, rr, model, covariates):
    x = xd.to_numpy(float)
    gi = list(xd.columns).index('group')
    z = np.delete(x, gi, axis=1)
    pinv, zpinv = np.linalg.pinv(x), np.linalg.pinv(z)
    h = np.einsum('ij,ji->i', x, pinv)
    hz = np.einsum('ij,ji->i', z, zpinv)
    rank, zrank = np.linalg.matrix_rank(x), np.linalg.matrix_rank(z)
    group_resid = x[:, gi] - z @ (zpinv @ x[:, gi])
    reasons = []
    if rank != x.shape[1]: reasons.append('Full design is rank-deficient')
    if zrank != z.shape[1]: reasons.append('Restricted design is rank-deficient')
    if rank != zrank+1 or float(group_resid @ group_resid) < 1e-10:
        reasons.append('Group effect not identified conditional on covariates')
    if h.max() >= .999 or hz.max() >= .999:
        reasons.append('Leverage >= .999: HC2/HC3 invalid or numerically unstable; no automatic exclusion')
    cat = 'SITE_final' if model in {'site_adjusted_motion', 'site_shared_motion'} else 'scanner_exact'
    counts = pd.crosstab(rr[cat], rr.group).reindex(columns=[0, 1], fill_value=0)
    result = dict(model=model, n=len(x), n_PD=int(rr.group.sum()), n_HC=int((rr.group==0).sum()),
                  parameters=x.shape[1], rank=rank, df=len(x)-rank, null_rank=zrank,
                  group_residual_sum_squares=float(group_resid@group_resid),
                  condition_number=np.linalg.cond(x), max_leverage=h.max(),
                  null_max_leverage=hz.max(), n_leverage_ge0999=int((h>=.999).sum()),
                  n_leverage_gt02=int((h>.2).sum()), n_leverage_gt05=int((h>.5).sum()),
                  categorical_field=cat, category_count=len(counts),
                  singleton_categories=int((counts.sum(axis=1)==1).sum()),
                  categories_with_both_groups=int(((counts[0]>0)&(counts[1]>0)).sum()),
                  covariates=covariates, design_columns=list(xd.columns),
                  passed=not reasons, blocking_reasons=reasons)
    # Individual identities and blinded category counts remain in LOCAL files only.
    counts.to_csv(HERE/f'LOCAL_{model}_category_counts.csv')
    private = rr[['sub_id', 'external_row_index', 'group', 'age_final', 'meanFD_final',
                  'scanner_exact', 'SITE_final']].copy()
    private['leverage_full'] = h
    private['leverage_restricted'] = hz
    private.to_csv(HERE/f'LOCAL_{model}_design_audit.csv', index=False)
    return result

def rank_rho(a, reference_rank):
    ranks = stats.rankdata(a, axis=-1)
    ranks -= ranks.mean(axis=-1, keepdims=True)
    denominator = np.linalg.norm(ranks, axis=-1) * np.linalg.norm(reference_rank, axis=-1)
    assert np.all(denominator > 0)
    return np.einsum('bgr,gr->bg', ranks, reference_rank) / denominator

def resample_rho(weighted_residual, beta_base, ref_rank, n_boot, seed,
                 x, fitted, adjusted_residual, indices, gi, checks, label):
    rng = np.random.default_rng(seed)
    rhos = np.empty((n_boot, 3))
    for start in range(0, n_boot, BATCH):
        count = min(BATCH, n_boot-start)
        w = rng.integers(0, 2, size=(count, len(x)), dtype=np.int8)*2-1
        b = w.astype(float) @ weighted_residual + beta_base
        rhos[start:start+count] = rank_rho(b.reshape(count, 3, 283), ref_rank)
        if start == 0:
            errors = []
            for k in range(3):
                ystar = fitted[:, indices] + adjusted_residual[:, indices]*w[k, :, None]
                explicit = np.array([sm.OLS(ystar[:, j], x).fit().params[gi]
                                     for j in range(len(indices))])
                errors.append(float(np.max(np.abs(explicit-b[k, indices]))))
            checks.append(dict(check=label+'_bootstrap_explicit_refit',
                               max_abs_error=max(errors), passed=max(errors)<1e-10))
            assert max(errors) < 1e-10
        if start % (BATCH*10) == 0:
            print(f'{label}: {min(start+count,n_boot)}/{n_boot}', flush=True)
    assert np.isfinite(rhos).all()
    return rhos

def run_one(model, xd, rr, mask, covariates, y_all, disc, diagnostics):
    y, x = y_all[mask], xd.to_numpy(float)
    p, n = x.shape[1], len(x)
    gi = list(xd.columns).index('group')
    pinv = np.linalg.pinv(x)
    b = pinv@y
    fitted = x@b
    resid = y-fitted
    h = np.einsum('ij,ji->i', x, pinv)
    c, beta = pinv[gi], b[gi]
    se = np.sqrt(np.sum((c[:, None]*resid/(1-h[:, None]))**2, axis=0))
    assert np.isfinite(se).all() and (se>0).all()
    df = n-p
    tvalue = beta/se
    pvalue = 2*stats.t.sf(np.abs(tvalue), df)
    qvalue = np.concatenate([bh(row) for row in pvalue.reshape(3,283)])
    crit = stats.t.ppf(.975, df)
    out = disc[['gradient', 'roi_id', 'roi_name', 'anatomical_description', 'network']].copy()
    out.insert(0, 'model', model)
    for key, val in dict(estimate=beta, se=se, t=tvalue, df=df,
                         ci95_low=beta-crit*se, ci95_high=beta+crit*se,
                         p_two_sided=pvalue, q_BH_283=qvalue, n=n,
                         n_PD=int(rr.group.sum()), n_HC=int((rr.group==0).sum()),
                         standard_errors='HC3', covariates=covariates,
                         effect_direction='PD minus HC', unit='Aligned gradient coordinate',
                         discovery_estimate=disc.estimate.to_numpy(),
                         discovery_q_BH=disc.q_BH.to_numpy()).items():
        out[key] = val
    out['significant_q05'] = out.q_BH_283 < .05
    checks = []
    # Check every regional fit, not just a sample, against independent statsmodels fitting.
    errors = {'beta':0., 'se':0., 'p':0., 'ci':0.}
    for j in range(849):
        fit = sm.OLS(y[:,j], x).fit(cov_type='HC3', use_t=True)
        errors['beta'] = max(errors['beta'], abs(fit.params[gi]-beta[j]))
        errors['se'] = max(errors['se'], abs(fit.bse[gi]-se[j]))
        errors['p'] = max(errors['p'], abs(fit.pvalues[gi]-pvalue[j]))
        errors['ci'] = max(errors['ci'], float(np.max(np.abs(
            fit.conf_int()[gi] - [beta[j]-crit*se[j], beta[j]+crit*se[j]]))))
    assert max(errors.values()) < 1e-9
    checks.append(dict(check=model+'_all849_vs_statsmodels_HC3', max_abs_error=max(errors.values()), passed=True))
    z = np.delete(x, gi, axis=1)
    zpinv = np.linalg.pinv(z)
    zfit = z@(zpinv@y)
    hz = np.einsum('ij,ji->i', z, zpinv)
    assert h.max()<.999 and hz.max()<.999
    rz, rf = (y-zfit)/np.sqrt(1-hz[:,None]), resid/np.sqrt(1-h[:,None])
    assert np.max(np.abs(c@zfit))<1e-10
    ref = disc.estimate.to_numpy().reshape(3,283)
    ref_rank = stats.rankdata(ref, axis=-1)
    ref_rank -= ref_rank.mean(axis=-1, keepdims=True)
    observed = rank_rho(beta.reshape(1,3,283), ref_rank)[0]
    indices = [0,56,200,282,283,400,565,566,790,848]
    null = resample_rho(c[:,None]*rz, np.zeros(849), ref_rank, N_NULL, NULL_SEED,
                        x,zfit,rz,indices,gi,checks,model+'_null')
    boot = resample_rho(c[:,None]*rf, beta, ref_rank, N_CI, CI_SEED,
                        x,fitted,rf,indices,gi,checks,model+'_CI')
    wb_p = (1+(np.abs(null)>=np.abs(observed)[None,:]).sum(axis=0))/(N_NULL+1)
    wb_q = bh(wb_p)
    percentile = np.quantile(boot,[.025,.975],axis=0)
    boot_z = np.arctanh(np.clip(boot,-1+1e-12,1-1e-12))
    bounds = np.quantile(boot_z,[.025,.975],axis=0)
    basic = np.tanh(2*np.arctanh(observed)[None,:]-bounds[::-1])
    agreement, overlap = [], []
    for k in range(3):
        selected = disc.loc[disc.gradient.eq(f'G{k+1}'),'q_BH'].to_numpy()<.05
        bvec = beta.reshape(3,283)[k]
        same = np.sign(bvec)==np.sign(ref[k])
        significant = qvalue.reshape(3,283)[k]<.05
        agreement.append(dict(model=model,gradient=f'G{k+1}',n=n,n_PD=int(rr.group.sum()),
            n_HC=int((rr.group==0).sum()),covariates=covariates,n_regions=283,
            spearman_rho=observed[k],ci95_low=basic[0,k],ci95_high=basic[1,k],
            ci_method='Fisher-z basic unrestricted HC2 wild bootstrap; conditional on discovery map',
            p_wild_two_sided=wb_p[k],q_BH_3=wb_q[k],
            null_hypothesis='No external regional group effects, conditional on fixed discovery map and design',
            n_null=N_NULL,n_ci=N_CI,null_seed=NULL_SEED,ci_seed=CI_SEED,
            null_mean=float(null[:,k].mean()),null_sd=float(null[:,k].std(ddof=1)),
            bootstrap_mean=float(boot[:,k].mean()),
            percentile_ci_low_diagnostic=percentile[0,k],percentile_ci_high_diagnostic=percentile[1,k],
            pearson_r_descriptive=stats.pearsonr(ref[k],bvec).statistic,
            same_sign_regions=int(same.sum()),same_sign_proportion=float(same.mean()),
            unique_discovery_effects=len(np.unique(ref[k])),unique_external_effects=len(np.unique(bvec))))
        overlap.append(dict(model=model,gradient=f'G{k+1}',discovery_significant=int(selected.sum()),
            external_significant=int(significant.sum()),
            discovery_selected_same_direction=int((selected&same).sum()),
            discovery_selected_external_significant=int((selected&significant).sum()),
            discovery_selected_external_significant_same_direction=int((selected&significant&same).sum()),
            scope='Descriptive counts; not a separate inferential family'))
    # Write completed models independently before producing combined files.
    dest = HERE/model
    dest.mkdir(exist_ok=True)
    out.to_csv(dest/'regional_effects_all_849.csv',index=False)
    pd.DataFrame(agreement).to_csv(dest/'map_agreement.csv',index=False)
    pd.DataFrame(overlap).to_csv(dest/'regional_overlap_descriptive.csv',index=False)
    pd.DataFrame(checks).to_csv(dest/'verification.csv',index=False)
    np.savez_compressed(dest/'bootstrap_distributions.npz', null_rho=null, unrestricted_rho=boot)
    save_json(dest/'model_diagnostics.json',diagnostics)
    print(pd.DataFrame(agreement)[['model','gradient','spearman_rho','p_wild_two_sided','q_BH_3']].to_string(index=False),flush=True)
    return out, agreement, overlap, checks

def main():
    global HERE
    parser = argparse.ArgumentParser()
    parser.add_argument('--preflight',action='store_true')
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--models',nargs='+',choices=MODEL_NAMES,default=['age_sex_motion','scanner_adjusted_motion','site_shared_motion'])
    parser.add_argument('--approved-site-support',default=str(INPUTS/'site_support.json'))
    args = parser.parse_args()
    HERE=args.output_dir.resolve(); HERE.mkdir(parents=True,exist_ok=True)
    start = time.time()
    reg,y,disc = load_inputs()
    designs, diagnostics = {}, []
    for model in args.models:
        xd,rr,mask,cov = design(reg,model,args.approved_site_support)
        diagnostic = inspect_design(xd,rr,model,cov)
        diagnostics.append(diagnostic)
        designs[model] = (xd,rr,mask,cov)
    save_json(HERE/'preflight_diagnostics.json',diagnostics)
    print(json.dumps([{k:v for k,v in d.items() if k!='design_columns'} for d in diagnostics],
                     ensure_ascii=False,indent=2,default=serial),flush=True)
    if args.preflight:
        return
    assert all(d['passed'] for d in diagnostics), 'Blocked design: review preflight; no model exclusion/merging is automatic'
    regional, agreement, overlap, checks = [], [], [], [dict(check='frozen_registry_original_gradient_row_order',passed=True)]
    for d in diagnostics:
        model = d['model']
        out,aa,oo,cc = run_one(model,*designs[model],y,disc,d)
        regional.append(out); agreement.extend(aa); overlap.extend(oo); checks.extend(cc)
    pd.concat(regional,ignore_index=True).to_csv(HERE/'regional_effects_all_849.csv',index=False)
    pd.DataFrame(agreement).to_csv(HERE/'map_agreement.csv',index=False)
    pd.DataFrame(overlap).to_csv(HERE/'regional_overlap_descriptive.csv',index=False)
    pd.DataFrame(checks).to_csv(HERE/'verification.csv',index=False)
    assert all(c['passed'] for c in checks)
    paths=[REG,GRAD,DISC,ORDER,Path(__file__)]
    if args.approved_site_support:paths.append(Path(args.approved_site_support))
    manifest=dict(input_sha256={str(p):sha(p) for p in paths},
        numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,statsmodels=statsmodels.__version__,
        python=sys.version,executable=sys.executable,BLAS_threads=1,
        elapsed_seconds=time.time()-start,models=diagnostics,
        n_null=N_NULL,n_ci=N_CI,null_seed=NULL_SEED,ci_seed=CI_SEED,
        approved_site_support=args.approved_site_support,
        limits=['Conditional on fixed discovery effects/reference; discovery uncertainty omitted',
                'No spatial-specificity test or whole-brain cortical spin',
                'SITE_final is a blinded clinical site identifier, not verified physical MRI location',
                'Participant-level resampling assumes independence; no site-cluster bootstrap',
                'Updated covariates do not resolve uncertain MRI medication state'])
    save_json(HERE/'input_manifest.json',manifest)
    print('Finished in',round(time.time()-start,2),'seconds',flush=True)

if __name__=='__main__':
    main()
