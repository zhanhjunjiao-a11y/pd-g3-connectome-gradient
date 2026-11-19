"""Refit current gradient/system/clinical models from package-private inputs.

No imports or reads of the historical project are performed. Expected aggregate
tables are read only after result generation, for numerical verification.
"""
import os,sys
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
ENV=Path(sys.executable).parent
os.environ['PATH']=str(ENV/'Library/bin')+os.pathsep+str(ENV)+os.pathsep+os.environ.get('PATH','')
DLL=[os.add_dll_directory(str(p)) for p in (ENV,ENV/'Library/bin',ENV/'DLLs') if os.name=='nt' and p.exists()]
import argparse,json,hashlib,platform,time
import numpy as np
import pandas as pd
import scipy
from scipy import stats
import statsmodels
import statsmodels.api as sm
from statsmodels.stats.diagnostic import het_breuschpagan

PKG=Path(__file__).resolve().parents[2]
INPUT=PKG/'inputs_private/core';VERIFY=PKG/'verification/core'
PRIMARY=['compression','profile_angle','within_sm_fc','within_cb_fc','between_fc','segregation']
AUX=['sm_g3','cb_g3','profile_pearson_z','profile_rms','gradient_distance_3d','gradient_distance_10d','global_fc']
CLINICAL=['fixed19_g3','compression','gradient_distance_3d','gradient_distance_10d','profile_angle','within_sm_fc','within_cb_fc','between_fc']
GEOMETRY=['compression','gradient_distance_3d','gradient_distance_10d']
FC=['profile_angle','within_sm_fc','within_cb_fc','between_fc']
CONTEXT=[f'{s}_within_fc' for s in ['DMN','SM','CB']]+[f'{a}_{b}_{m}' for a,b in [('DMN','SM'),('DMN','CB'),('SM','CB')] for m in ['g2_gap','profile_angle','between_fc']]
KEYS={2:['gradient','roi_id'],6:['metric'],7:['gradient','roi_id'],8:['model','metric'],9:['metric'],10:['metric'],12:['metric'],13:['model','roi_id'],15:['analysis','gradient','roi_id','metric'],21:['metric'],26:['model','metric']}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(name):return pd.read_csv(INPUT/name,float_precision='round_trip')
def bh(p):
    p=np.asarray(p,float);order=np.argsort(p,kind='stable');v=p[order]*len(p)/np.arange(1,len(p)+1)
    out=np.empty_like(v);out[order]=np.clip(np.minimum.accumulate(v[::-1])[::-1],0,1);return out
def adjust(rows,field='q_BH'):
    for r,q in zip(rows,bh([r['p_two_sided'] for r in rows])):r[field]=float(q);r['significant_q05']=bool(q<.05)
    return rows
def design(frame,terms):return np.column_stack([np.ones(len(frame)),frame[terms].to_numpy(float)])
def fit(y,x,j=1,robust=True):
    y=np.asarray(y,float);x=np.asarray(x,float)
    assert np.isfinite(y).all() and np.isfinite(x).all() and np.linalg.matrix_rank(x)==x.shape[1]
    model=sm.OLS(y,x).fit();est=model.get_robustcov_results(cov_type='HC3',use_t=True) if robust else model
    ci=est.conf_int(alpha=.05)[j]
    r=dict(n=len(y),estimate=float(est.params[j]),se=float(est.bse[j]),ci95_low=float(ci[0]),ci95_high=float(ci[1]),t=float(est.tvalues[j]),df=float(est.df_resid),p_two_sided=float(est.pvalues[j]),standard_errors='HC3' if robust else 'Conventional OLS')
    # Independent NumPy sandwich/conventional-SE numerical check for every fit.
    pinv=np.linalg.pinv(x);b=pinv@y;res=y-x@b;h=np.sum(x*pinv.T,axis=1)
    sem=np.sqrt(np.sum(pinv[j]**2*(res/(1-h))**2)) if robust else np.sqrt(np.sum(res**2)/model.df_resid*np.sum(pinv[j]**2))
    assert abs(r['estimate']-b[j])<1e-8 and abs(r['se']-sem)<1e-8
    return r,model
def diagnostic(model,x,analysis,gradient='',roi_id='',metric=''):
    resid=model.resid;inf=model.get_influence();cook=inf.cooks_distance[0]
    r=dict(analysis=analysis,gradient=gradient,roi_id=roi_id,metric=metric,
        residual_shapiro_p=float(stats.shapiro(resid).pvalue),breusch_pagan_p=float(het_breuschpagan(resid,x)[1]),maximum_cooks_distance=float(cook.max()),n_cooks_above_4_over_n=int((cook>4/len(resid)).sum()),maximum_leverage=float(inf.hat_matrix_diag.max()))
    if analysis=='Regional paired gradients':
        r.update(residual_skew=float(stats.skew(resid)),residual_excess_kurtosis=float(stats.kurtosis(resid)),maximum_absolute_standardized_residual=float(np.max(np.abs(resid/resid.std(ddof=1)))))
    elif analysis=='Acute motor improvement':
        r.update(residual_skew=float(stats.skew(resid,bias=False)),residual_excess_kurtosis=float(stats.kurtosis(resid,bias=False)),n=len(resid),design_rank=x.shape[1],removed_for_influence=0)
    else:r.update(n=len(resid),removed_for_influence=0)
    return r
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output-dir',type=Path,required=True);args=ap.parse_args()
    out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
    assert out!=INPUT.resolve() and INPUT.resolve() not in out.parents
    t0=time.monotonic();manifest=json.loads((INPUT/'input_manifest.json').read_text(encoding='utf-8'))
    for e in manifest['prepared']:assert sha(PKG/e['path'])==e['sha256'],e['path']
    features=read('features.csv');assert features.definition.eq('anatomical_all').all()
    bas=features.loc[features.role.isin(['HC','PD'])].copy()
    assert len(bas)==475 and bas.person_id.nunique()==475
    rd=read('regional_design.csv');g=np.load(INPUT/'baseline_gradients.npy',allow_pickle=False)
    assert g.shape==(475,283,10)
    pair=read('paired_registry236.csv').set_index('sub_id');offg=np.load(INPUT/'paired_off_gradients.npy',allow_pickle=False);ong=np.load(INPUT/'paired_on_gradients.npy',allow_pickle=False)
    assert offg.shape==ong.shape==(236,283,10)
    policy=json.loads((INPUT/'paired_analysis_policy.json').read_text(encoding='utf-8'))
    active_ids=set(policy['included_sub_ids']);assert len(active_ids)==235
    good=pair.index.isin(active_ids);assert good.sum()==235
    assert np.array_equal(good,pair.dfd.notna().to_numpy())
    assert active_ids.issubset(set(features.loc[features.role.eq('PD'),'sub_id']))
    yp=(ong-offg)[good];xp=design(pair.loc[good],['dfd'])
    atlas=read('atlas_networks.csv');smix=np.flatnonzero(atlas.category.eq('Somatomotor'));cbix=np.flatnonzero(atlas.category.eq('Cerebellum'));assert len(smix)==30 and len(cbix)==31
    tables={};diags=[];model_counts={};fits_done=0
    print('Inputs validated: baseline475, active paired235 from preserved source pool236, SMN30/CB31.',flush=True)
    # S02: all regions/components, conventional OLS, separate BH283.
    rr=[];x=design(rd,['group','age','sex','meanFD'])
    for k in range(3):
        family=[]
        for j in range(283):
            r,_=fit(g[:,j,k],x,robust=False);r.update(gradient=f'G{k+1}',roi_id=j+1,n_PD=334,n_HC=141,family_n=283);family.append(r)
        rr+=adjust(family)
    tables[2]=pd.DataFrame(rr);model_counts[2]=849;fits_done+=849
    disease=tables[2].set_index(['gradient','roi_id'])
    # S13: G3 only, archived four sensitivity specifications.
    rr=[]
    for mode in ['linear_meanFD','nonlinear_meanFD','strict_meanFD_lt_0.2','HC3_robust_meanFD']:
        keep=rd.meanFD.lt(.2).to_numpy() if mode=='strict_meanFD_lt_0.2' else np.ones(475,bool)
        d=rd.loc[keep].copy();d['meanFD_squared']=d.meanFD**2
        terms=['group','age','sex','meanFD']+(['meanFD_squared'] if mode=='nonlinear_meanFD' else [])
        x=design(d,terms);family=[]
        for j in range(283):
            r,_=fit(g[keep,j,2],x,robust=mode=='HC3_robust_meanFD');r.update(model=mode,gradient='G3',roi_id=j+1,n_PD=int(d.group.eq(1).sum()),n_HC=int(d.group.eq(0).sum()),family_n=283);family.append(r)
        rr+=adjust(family)
    tables[13]=pd.DataFrame(rr);model_counts[13]=1132;fits_done+=1132
    print('Regional disease and G3 sensitivities fitted.',flush=True)
    # System baseline + S26 motion checks; no outcome redefinition.
    rr=[];bas['meanFD_squared']=bas.meanFD**2
    for mode,d,terms in [('reference_linear_FD',bas,['group','age','sex','meanFD']),('low_FD_lt_0.2',bas.loc[bas.meanFD.lt(.2)],['group','age','sex','meanFD']),('quadratic_FD',bas,['group','age','sex','meanFD','meanFD_squared'])]:
        x=design(d,terms)
        for fam,metrics in [('primary_6',PRIMARY),('auxiliary_7',AUX)]:
            family=[]
            for metric in metrics:
                r,_=fit(d[metric],x)
                groups=[d.loc[d.group.eq(k),metric] for k in [0,1]]
                pooled=np.sqrt(sum((len(z)-1)*z.var(ddof=1) for z in groups)/(len(d)-2))
                r.update(model=mode,metric=metric,family_id=fam,family_n=len(metrics),n_PD=int(d.group.eq(1).sum()),n_HC=int(d.group.eq(0).sum()),standardized_effect=r['estimate']/pooled,standardized_ci_low=r['ci95_low']/pooled,standardized_ci_high=r['ci95_high']/pooled);family.append(r)
            rr+=adjust(family)
    tables[26]=pd.DataFrame(rr);model_counts[26]=39;fits_done+=39
    tables[6]=tables[26].loc[tables[26].model.eq('reference_linear_FD')&tables[26].metric.isin(PRIMARY)].copy();model_counts[6]=6
    tables[12]=tables[26].loc[tables[26].model.eq('reference_linear_FD')&tables[26].metric.isin(AUX)].copy();model_counts[12]=7
    # S07 and its S15 diagnostics, from paired regional arrays.
    rr=[]
    for k in range(3):
        family=[]
        for j in range(283):
            y=yp[:,j,k];r,model=fit(y,xp,j=0);dr=disease.loc[(f'G{k+1}',j+1)]
            r.update(gradient=f'G{k+1}',roi_id=j+1,OFF_raw_mean=float(offg[good,j,k].mean()),OFF_raw_sd=float(offg[good,j,k].std(ddof=1)),ON_raw_mean=float(ong[good,j,k].mean()),ON_raw_sd=float(ong[good,j,k].std(ddof=1)),raw_change_mean=float(y.mean()),raw_change_sd=float(y.std(ddof=1)),disease_effect=float(dr.estimate),disease_q_BH=float(dr.q_BH),original_disease_difference_region=bool(dr.q_BH<.05),opposite_to_disease_direction=bool(r['estimate']*dr.estimate<0),family_n=283)
            family.append(r);diags.append(diagnostic(model,xp,'Regional paired gradients',f'G{k+1}',j+1))
        rr+=adjust(family)
    for r,q in zip(rr,bh([r['p_two_sided'] for r in rr])):r.update(q_BH_849_sensitivity=float(q),significant_BH849=bool(q<.05))
    tables[7]=pd.DataFrame(rr);model_counts[7]=849;fits_done+=849
    # S08: frozen atlas3 system pairs, preserving the full twelve-test family.
    context=read('context_features.csv');assert context.definition.eq('atlas3').all()
    oc=context.loc[context.role.eq('OFF')].set_index('person_id');nc=context.loc[context.role.eq('ON')].set_index('person_id').reindex(oc.index)
    assert len(oc)==len(nc)==236 and oc.sub_id.equals(nc.sub_id)
    delta=nc[CONTEXT]-oc[CONTEXT];delta['dfd']=oc.dfd
    delta=delta.loc[oc.sub_id.isin(active_ids)].copy()
    assert len(delta)==235 and delta.dfd.notna().all()
    rr=[]
    for mode,d,terms in [('motion_adjusted_235',delta,['dfd']),('same_235_raw',delta,[])]:
        family=[];x=design(d,terms)
        for metric in CONTEXT:
            y=d[metric];r,_=fit(y,x,j=0);sd=y.std(ddof=1)
            r.update(model=mode,metric=metric,family_n=12,raw_change_mean=float(y.mean()),raw_change_sd=float(sd),standardized_effect=r['estimate']/sd,standardized_ci_low=r['ci95_low']/sd,standardized_ci_high=r['ci95_high']/sd);family.append(r)
        rr+=adjust(family)
    tables[8]=pd.DataFrame(rr);model_counts[8]=24;fits_done+=24
    # S21: geometry recomputed directly from complete-system gradient coordinates.
    mg_off=offg[:,smix].mean(axis=1)-offg[:,cbix].mean(axis=1);mg_on=ong[:,smix].mean(axis=1)-ong[:,cbix].mean(axis=1)
    geom={'compression':(mg_off[:,2],mg_on[:,2]),'gradient_distance_3d':(np.linalg.norm(mg_off[:,:3],axis=1),np.linalg.norm(mg_on[:,:3],axis=1)),'gradient_distance_10d':(np.linalg.norm(mg_off,axis=1),np.linalg.norm(mg_on,axis=1))}
    rr=[]
    for metric,(a,b) in geom.items():
        y=(b-a)[good];r,_=fit(y,xp,j=0);rawci=stats.t.interval(.95,len(y)-1,loc=y.mean(),scale=stats.sem(y))
        r.update(metric=metric,adjusted_change=r['estimate'],se_HC3=r['se'],ci_low=r['ci95_low'],ci_high=r['ci95_high'],p=r['p_two_sided'],raw_mean_change=float(y.mean()),raw_mean_ci_low=float(rawci[0]),raw_mean_ci_high=float(rawci[1]),raw_sd_change=float(y.std(ddof=1)),raw_descriptive_ci_df=234,family_n=3,two_sided=True,use_t=True,exploratory=True);rr.append(r)
    tables[21]=pd.DataFrame(adjust(rr,'q_BH3'));model_counts[21]=3;fits_done+=3
    print('System disease/paired, all regional paired, and geometry fitted.',flush=True)
    # S09/S10: one imaging predictor per model, separate eight-test families.
    b=features.loc[features.role.eq('PD')].set_index('sub_id').sort_index().copy()
    off=features.loc[features.role.eq('OFF')].set_index('sub_id').sort_index().copy();on=features.loc[features.role.eq('ON')].set_index('sub_id').reindex(off.index).copy()
    scores=read('fixed19_baseline.csv').set_index('sub_id');sp=read('fixed19_paired.csv').set_index('sub_id')
    b['fixed19_g3']=scores.reindex(b.index).comp_full
    off['fixed19_g3']=sp.reindex(off.index).baseline_comp;on['fixed19_g3']=sp.reindex(on.index).baseline_comp+sp.reindex(on.index).dcomp
    frozen=read('paired_clinical235.csv').set_index('sub_id').sort_index()
    clinical=frozen[['age','sex','UPDRS_OFF','UPDRS_ON','dfd','improvement']].copy()
    assert len(b)==334 and len(clinical)==235 and set(clinical.index)==set(pair.index[good])
    for m in CLINICAL:
        clinical['off_'+m]=off.loc[clinical.index,m];clinical['delta_'+m]=on.loc[clinical.index,m]-off.loc[clinical.index,m]
    clinextra=[]
    for snum,frame,outcome in [(9,b,'UPDRS_OFF'),(10,clinical,'improvement')]:
        rr=[]
        for metric in CLINICAL:
            term=metric if snum==9 else 'delta_'+metric
            covs=['age','sex','meanFD'] if snum==9 else ['UPDRS_OFF','off_'+metric,'age','sex','dfd']
            y=frame[outcome];x=design(frame,[term]+covs);r,model=fit(y,x);factor=frame[term].std(ddof=1)/y.std(ddof=1)
            r.update(metric=metric,family_n=8,standardized_beta=r['estimate']*factor,standardized_ci95_low=r['ci95_low']*factor,standardized_ci95_high=r['ci95_high']*factor)
            reduced=sm.OLS(y,design(frame,covs)).fit();r['partial_r2_in_sample']=float((reduced.ssr-model.ssr)/reduced.ssr)
            rr.append(r)
            if metric in FC:diags.append(diagnostic(model,x,'Baseline clinical association' if snum==9 else 'Acute motor improvement',metric=metric))
        tables[snum]=pd.DataFrame(adjust(rr,'q_BH8'));model_counts[snum]=8;fits_done+=8
    tables[15]=pd.DataFrame(diags);model_counts[15]=857
    # Results have now been generated exclusively from individual-level inputs.
    for n,frame in tables.items():frame.to_csv(out/f'S{n:02d}_refit.csv',index=False,float_format='%.17g')
    print('All estimates generated; beginning comparison with frozen published tables.',flush=True)
    for e in manifest['expected_aggregates']:
        assert sha(PKG/e['path'])==e['sha256'], 'Frozen comparison table changed: '+e['path']
    checks=[]
    for n,actual in sorted(tables.items()):
        p=next((VERIFY/'expected').glob(f'Supplementary_Data_S{n:02d}_*.csv'))
        old=pd.read_csv(p,float_precision='round_trip',keep_default_na=False)
        a=actual.copy().fillna('');keys=KEYS[n]
        def kval(v):
            if isinstance(v,(int,float,np.integer,np.floating)) and np.isfinite(v) and float(v).is_integer():return str(int(v))
            return str(v)
        lookup={tuple(kval(r[k]) for k in keys):r for r in a.to_dict('records')}
        missing=[];diffs=[];compared=0;errors={};numeric_columns=[]
        for col in old.columns:
            if col in keys:continue
            values=[v for v in old[col] if v not in ['', 'nan']]
            if not values:continue
            try:[float(v) for v in values]
            except (ValueError,TypeError):continue
            numeric_columns.append(col)
        for r in old.to_dict('records'):
            key=tuple(kval(r[k]) for k in keys);new=lookup.get(key)
            if new is None:missing.append(key);continue
            for col in numeric_columns:
                if r[col] in ['', 'nan']:continue
                if col not in new or new[col]=='':missing.append(key+(col,));continue
                expected=float(r[col]);value=float(new[col]);e=abs(value-expected);errors[col]=max(errors.get(col,0.),e);compared+=1
                # Truncated historical CSVs and BLAS/library variation: pointwise tolerance.
                if not np.isclose(value,expected,rtol=2e-9,atol=2e-10):diffs.append({'key':key,'field':col,'expected':expected,'refit':value,'abs_error':e})
        goodcheck=not missing and not diffs and len(old)==len(actual)
        check={'data':f'S{n:02d}','status':'PASS' if goodcheck else 'DIFFERENCE','expected_rows':len(old),'refit_rows':len(actual),'model_count':model_counts[n],
               'numeric_values_compared':compared,'numeric_fields':numeric_columns,'max_abs_error':max(errors.values(),default=0.),'max_abs_error_by_field':errors,
               'missing_keys_or_fields':missing,'differences':diffs,'reason':'Independently refitted from package individual inputs; numeric columns compared after fitting' if goodcheck else 'See differences; original expectations were not changed',
               'expected_sha256':sha(p),'refit_sha256':sha(out/f'S{n:02d}_refit.csv')}
        checks.append(check)
    report={'status':'PASS' if all(x['status']=='PASS' for x in checks) else 'DIFFERENCE','checks':checks,'unique_inferential_fits':fits_done,'diagnostic_rows':len(diags),
            'input_manifest_sha256':sha(INPUT/'input_manifest.json'),'code_sha256':sha(__file__),'packages':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'statsmodels':statsmodels.__version__},
            'elapsed_seconds':time.monotonic()-t0,'original_files_modified':False,'source_aggregates_used_as_results':False,
            'scope':'Inference from frozen individual gradient/system features; no fMRI preprocessing, FC estimation or gradient embedding in this core module. S15 diagnoses historical four-FC models, not new BH4 inference.'}
    (out/'verification.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    (VERIFY/'verification.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({'status':report['status'],'fits':fits_done,'checks':[{k:x[k] for k in ['data','status','refit_rows','max_abs_error']} for x in checks],'elapsed_seconds':report['elapsed_seconds']},ensure_ascii=False),flush=True)
    raise SystemExit(0 if report['status']=='PASS' else 1)
if __name__=='__main__':main()
