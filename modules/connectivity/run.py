"""Recompute FC and profile analyses from the package's restricted FC inputs.

No legacy scripts are imported. No reference result is an inference input.
The original 283-region families and frozen sample memberships are retained.
"""
from pathlib import Path
import os,sys,json,hashlib,argparse,time,platform
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
if os.name=='nt':
    dll=Path(sys.prefix)/'Library/bin'
    if dll.exists():os.environ['PATH']=str(dll)+os.pathsep+os.environ['PATH'];_dll=os.add_dll_directory(str(dll))
import numpy as np
import pandas as pd
import scipy
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
sys.stdout.reconfigure(encoding='utf-8')
PKG=Path(__file__).resolve().parents[2]
INPUT=PKG/'inputs_private/connectivity'
VERIFY=PKG/'verification/connectivity'
NETS=['Visual','Somatomotor','DorsAttn','VentAttn','Limbic','Frontoparietal','Default','Subcortex','Cerebellum']
EXTERNAL=['Visual','DorsAttn','VentAttn','Frontoparietal','Default','Limbic','Subcortex']
METRICS=['compression','sm_g3','cb_g3','profile_angle','profile_pearson_z','profile_rms','gradient_distance_3d','gradient_distance_10d','within_sm_fc','within_cb_fc','between_fc','segregation','global_fc']
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
def save(d,p):d.to_csv(p,index=False,float_format='%.15g',encoding='utf-8-sig')
def bh(p):return multipletests(np.asarray(p),method='fdr_bh')[1]
def embed(fc):
    """Actual saved implementation: angular affinity, affinity sparsity, alpha=.5,t=1."""
    x=fc.copy();np.fill_diagonal(x,0)
    norms=np.linalg.norm(x,axis=1,keepdims=True);norms[norms==0]=1e-10;x=x/norms
    a=1-np.arccos(np.clip(x@x.T,-1,1))/np.pi
    a[a<np.percentile(a,90,axis=1)[:,None]]=0;a=(a+a.T)/2
    da=a.sum(1)**-.5;m=a*np.outer(da,da);d=m.sum(1)
    p=m/d[:,None];p=p*np.sqrt(d)[:,None]/np.sqrt(d)[None,:]
    val,vec=np.linalg.eigh((p+p.T)/2);ix=np.argsort(val)[::-1][1:11]
    grad=vec[:,ix]/np.sqrt(d[:,None]);val=val[ix]
    grad*=np.where(val>0,val,1)[None,:]
    return grad,val
def align(source,target):
    s=source-source.mean(0);t=target-target.mean(0);u,_,vt=np.linalg.svd(t.T@s)
    return s@(u@vt).T+target.mean(0)
def atlas_from_source():
    a=pd.read_csv(INPUT/'roi_source_metadata.csv').sort_values('Label').reset_index(drop=True)
    assert np.array_equal(a.Label,np.arange(1,284))
    net=pd.read_csv(INPUT/'atlas_network_assignment.csv').sort_values('Label')
    assert np.array_equal(net.category,a.category)
    o=pd.DataFrame(dict(roi_id=a.Label,roi_name=a.roi_name,anatomical_description=a['Modified cyto-architectonic'],network=a.category))
    o['structure']=np.select([a.Label.le(210),a.Label.le(252)],['Cortex','Subcortex and brainstem'],default='Cerebellum')
    o['lobe']=a.Lobe;o['gyrus']=a.Gyrus
    # The anatomical class boundary (210 cortex) is not the atlas-source boundary:
    # Brainnetome contributes 210 cortical + 24 subcortical parcels; CIT168 adds 18.
    o['atlas_source']=np.select([a.Label.le(234),a.Label.le(252)],['Brainnetome','CIT168'],default='SUIT-related')
    o['somatomotor_system']=a.category.eq('Somatomotor');o['cerebellar_system']=a.category.eq('Cerebellum')
    o['external_target']=~(o.somatomotor_system|o.cerebellar_system);o['voxel_count_2mm']=a['Voxel num(2mm)']
    return a,o
def extract_all(fc,grad,reg,atlas,out):
    a=np.flatnonzero(atlas.network.eq('Somatomotor'));b=np.flatnonzero(atlas.network.eq('Cerebellum'));t=np.flatnonzero(atlas.external_target)
    assert (len(a),len(b),len(t))==(30,31,222)
    ii,jj=np.triu_indices(283,1);rows=[];ps=[];qs=[]
    for k in range(len(reg)):
        x=np.array(fc[k]);np.fill_diagonal(x,0);g=grad[k]
        p=x[np.ix_(a,t)].mean(0);q=x[np.ix_(b,t)].mean(0);ps.append(p);qs.append(q)
        ma=g[a].mean(0);mb=g[b].mean(0)
        wa=x[np.ix_(a,a)][np.triu_indices(len(a),1)].mean();wb=x[np.ix_(b,b)][np.triu_indices(len(b),1)].mean();ab=x[np.ix_(a,b)].mean()
        cos=np.dot(p,q)/(np.linalg.norm(p)*np.linalg.norm(q));corr=np.corrcoef(p,q)[0,1]
        d=dict(compression=ma[2]-mb[2],sm_g3=ma[2],cb_g3=mb[2],profile_angle=1-np.arccos(np.clip(cos,-1,1))/np.pi,profile_pearson_z=np.arctanh(np.clip(corr,-.99999999,.99999999)),profile_rms=np.sqrt(np.mean((p-q)**2)),gradient_distance_3d=np.linalg.norm(ma[:3]-mb[:3]),gradient_distance_10d=np.linalg.norm(ma-mb),within_sm_fc=wa,within_cb_fc=wb,between_fc=ab,segregation=(wa+wb)/2-ab,global_fc=x[ii,jj].mean())
        for net in NETS:
            j=atlas.iloc[t].network.eq(net).to_numpy()
            for seed,z in [('SM',p),('CB',q)]:d[f'net_{seed}_{net}']=z[j].mean() if j.any() else np.nan
        rows.append(d)
    d=pd.concat([reg.reset_index(drop=True),pd.DataFrame(rows)],axis=1);d['definition']='anatomical_all'
    save(d,out/'LOCAL_anatomical_all_recomputed.csv')
    np.savez_compressed(out/'LOCAL_profiles.npz',p=np.array(ps),q=np.array(qs),roi_ids=t+1)
    return d,np.array(ps),np.array(qs),t
def vector_models(reg,y,meta,model):
    cov=['age','sex','meanFD']+(['global_fc'] if model=='global_FC_adjustment' else [])
    x=np.column_stack([np.ones(len(reg)),reg.group.to_numpy(),*[(reg[c]-reg[c].mean())/reg[c].std(ddof=1) for c in cov]])
    xp=np.linalg.pinv(x);h=np.einsum('ij,ji->i',x,xp);w=xp[1]/(1-h);n,p=x.shape;df=n-p
    hm=reg.group.eq(0).to_numpy();pm=~hm;result=meta.copy();store={k:np.zeros(y.shape[1]) for k in ['estimate','se','HC_adjusted_mean','PD_adjusted_mean','standardized_effect','HC_raw_mean','HC_raw_sd','PD_raw_mean','PD_raw_sd']}
    for start in range(0,y.shape[1],2048):
        sl=slice(start,min(start+2048,y.shape[1]));z=np.asarray(y[:,sl]);b=xp@z;r=z-x@b
        pooled=np.sqrt(((hm.sum()-1)*z[hm].var(0,ddof=1)+(pm.sum()-1)*z[pm].var(0,ddof=1))/(n-2))
        vals=dict(estimate=b[1],se=np.sqrt(w*w@(r*r)),HC_adjusted_mean=b[0],PD_adjusted_mean=b[0]+b[1],standardized_effect=b[1]/pooled,HC_raw_mean=z[hm].mean(0),HC_raw_sd=z[hm].std(0,ddof=1),PD_raw_mean=z[pm].mean(0),PD_raw_sd=z[pm].std(0,ddof=1))
        for k,v in vals.items():store[k][sl]=v
    for k,v in store.items():result[k]=v
    result['t']=result.estimate/result.se;result['p_two_sided']=2*stats.t.sf(abs(result.t),df)
    crit=stats.t.ppf(.975,df);result['ci95_low']=result.estimate-crit*result.se;result['ci95_high']=result.estimate+crit*result.se
    result['q_BH']=bh(result.p_two_sided);result['family_n']=len(result);result['n']=n;result['n_PD']=int(pm.sum());result['n_HC']=int(hm.sum());result['df']=df
    result['standard_errors']='HC3';result['covariates']=' + '.join(cov);result['model']=model
    result['effect_direction']='PD minus HC';result['unit']='Mean Fisher-z FC';result['significant_q05']=result.q_BH.lt(.05)
    return result
def fc_analyses(fc,reg,rawatlas,atlas,feat,out):
    mask=reg.role.isin(['HC','PD']).to_numpy();r=reg.loc[mask].reset_index(drop=True).copy();r['global_fc']=feat.loc[mask,'global_fc'].to_numpy()
    ii,jj=np.triu_indices(283,1);seeds=rawatlas.loc[rawatlas.selected_G3,'Label'].tolist();assert len(seeds)==19
    em=pd.DataFrame(dict(roi_i=ii+1,roi_j=jj+1,name_i=atlas.roi_name.iloc[ii].to_numpy(),name_j=atlas.roi_name.iloc[jj].to_numpy(),network_i=atlas.network.iloc[ii].to_numpy(),network_j=atlas.network.iloc[jj].to_numpy()))
    em['touches_G3_selected19']=em.roi_i.isin(seeds)|em.roi_j.isin(seeds)
    nm=pd.DataFrame([dict(roi_id=i+1,roi_name=atlas.roi_name.iloc[i],seed_network=atlas.network.iloc[i],target_network=n,n_target_regions=int((atlas.network==n).sum()-(atlas.network.iloc[i]==n)),G3_selected19=i+1 in seeds) for i in range(283) for n in NETS])
    y=np.lib.format.open_memmap(out/'LOCAL_edges.npy',mode='w+',dtype=np.float64,shape=(475,39903));z=np.lib.format.open_memmap(out/'LOCAL_ROI_targets.npy',mode='w+',dtype=np.float64,shape=(475,2547))
    membership=np.column_stack([(atlas.network==n).to_numpy(float) for n in NETS]);counts=nm.n_target_regions.to_numpy().reshape(283,9)
    for k,i in enumerate(np.flatnonzero(mask)):
        x=np.array(fc[i]);np.fill_diagonal(x,0);y[k]=x[ii,jj];z[k]=(x@membership/counts).reshape(-1)
    y.flush();z.flush();sens=[]
    for model in ['primary','global_FC_adjustment','exclude_one_record_linkage_flagged_control']:
        keep=np.ones(len(r),bool) if model!='exclude_one_record_linkage_flagged_control' else ~r.exclude_flagged_control.to_numpy(bool)
        rr=r.loc[keep].reset_index(drop=True)
        for kind,arr,meta in [('edge',y,em),('ROI-to-network',z,nm)]:
            print('Fitting FC',model,kind,len(rr),arr.shape[1],flush=True)
            d=vector_models(rr,arr if keep.all() else arr[keep],meta,model)
            if model=='primary':
                if kind=='edge':
                    save(d,out/'Supplementary_Data_S03_wholebrain_edge_FC.csv');selected=[]
                    for seed in seeds:
                        s=d[(d.roi_i==seed)|(d.roi_j==seed)].copy();s['seed_roi']=seed;s['seed_name']=atlas.roi_name.iloc[seed-1];s['seed_network']=atlas.network.iloc[seed-1]
                        s['target_roi']=np.where(s.roi_i==seed,s.roi_j,s.roi_i);s['target_name']=atlas.roi_name.iloc[s.target_roi.to_numpy()-1].to_numpy();s['target_network']=atlas.network.iloc[s.target_roi.to_numpy()-1].to_numpy();selected.append(s)
                    save(pd.concat(selected,ignore_index=True),out/'Supplementary_Data_S04_G3_selected19_edge_display.csv')
                else:save(d,out/'Supplementary_Data_S05_ROI_to_network_mean_FC.csv')
            else:
                if kind=='edge':d=d[d.touches_G3_selected19].copy()
                d.insert(0,'analysis_level',kind);sens.append(d)
    save(pd.concat(sens,ignore_index=True),out/'Supplementary_Data_S11_edge_and_ROI_mean_sensitivity.csv')
def profile_analyses(reg,feat,atlas,ps,qs,targets,out):
    mask=reg.role.isin(['HC','PD']).to_numpy();r=reg.loc[mask].reset_index(drop=True);feat=feat.loc[mask].reset_index(drop=True);p=ps[mask];q=qs[mask]
    a=atlas.iloc[targets].reset_index(drop=True);u=p/np.linalg.norm(p,axis=1)[:,None];v=q/np.linalg.norm(q,axis=1)[:,None];dist=.5*(u-v)**2
    pc=p-p.mean(1)[:,None];qc=q-q.mean(1)[:,None];uc=pc/np.linalg.norm(pc,axis=1)[:,None];vc=qc/np.linalg.norm(qc,axis=1)[:,None];cdist=.5*(uc-vc)**2
    identity=float(np.max(abs(dist.sum(1)-(1-(u*v).sum(1)))));assert identity<1e-12
    x=np.column_stack([np.ones(len(r)),r[['group','age','sex','meanFD']].to_numpy()]);x0=np.r_[1,0,r[['age','sex','meanFD']].mean().to_numpy()];x1=x0.copy();x1[1]=1
    def fit(y,metric):
        rob=sm.OLS(y,x).fit().get_robustcov_results(cov_type='HC3',use_t=True);ci=rob.conf_int()[1];g0=np.asarray(y)[r.group.eq(0)];g1=np.asarray(y)[r.group.eq(1)]
        pooled=np.sqrt(((len(g0)-1)*g0.var(ddof=1)+(len(g1)-1)*g1.var(ddof=1))/(len(y)-2))
        d=dict(metric=metric,n=len(y),n_HC=len(g0),n_PD=len(g1),estimate=rob.params[1],se=rob.bse[1],ci_low=ci[0],ci_high=ci[1],t=rob.tvalues[1],df=rob.df_resid,p=rob.pvalues[1],standardized_effect=rob.params[1]/pooled,HC_raw_mean=g0.mean(),PD_raw_mean=g1.mean(),HC_adjusted_mean=x0@rob.params,PD_adjusted_mean=x1@rob.params,covariates='age + sex + meanFD',standard_errors='HC3',exploratory=True)
        for tag,xx in [('HC',x0),('PD',x1)]:
            se=np.sqrt(xx@rob.cov_params()@xx);vv=xx@rob.params;crit=stats.t.ppf(.975,rob.df_resid);d[tag+'_adjusted_ci_low']=vv-crit*se;d[tag+'_adjusted_ci_high']=vv+crit*se
        return d
    def family(rows):
        d=pd.DataFrame(rows);d['q_BH']=bh(d.p);d['family_n']=len(d);return d
    rows=[];center=[];raw=[];gaps=[]
    for net in EXTERNAL:
        j=a.network.eq(net).to_numpy();n=int(j.sum())
        d=fit(dist[:,j].sum(1),'distance_'+net);d.update(network=net,n_targets=n,per_target_estimate=d['estimate']/n,per_target_ci_low=d['ci_low']/n,per_target_ci_high=d['ci_high']/n);rows.append(d)
        d=fit(cdist[:,j].sum(1),'centered_distance_'+net);d.update(network=net,n_targets=n);center.append(d)
        for seed,z in [('SM',p),('CB',q)]:
            d=fit(z[:,j].mean(1),'rawFC_'+seed+'_'+net);d.update(seed=seed,network=net,n_targets=n);raw.append(d)
        d=fit((p[:,j]-q[:,j]).mean(1),'raw_SM_minus_CB_'+net);d.update(network=net,n_targets=n);gaps.append(d)
    frames=[]
    for data,label in [(rows,'cosine_distance_BH7'),(center,'centered_distance_BH7'),(raw,'raw_FC_BH14'),(gaps,'raw_SM_minus_CB_BH7')]:
        d=family(data);d.insert(0,'analysis_family',label);frames.append(d)
    save(pd.concat(frames,ignore_index=True),out/'Supplementary_Data_S17_profile_target_networks.csv')
    save(pd.DataFrame([fit(dist.sum(1),'total_cosine_distance'),fit(cdist.sum(1),'total_centered_distance'),fit(1-np.arccos(np.clip((u*v).sum(1),-1,1))/np.pi,'profile_angle')]),out/'total_distance_results.csv')
    roirows=[];rawrows=[]
    for j in range(len(targets)):
        ar=a.iloc[j];meta={k:ar[k] for k in ['roi_id','roi_name','anatomical_description','network','structure']};d=fit(dist[:,j],'ROI_distance_'+str(ar.roi_id));d.update(meta)
        for role,g in [('HC',0),('PD',1)]:
            m=r.group.eq(g).to_numpy();gap=(u-v)[m,j]
            d.update({role+'_u_mean':u[m,j].mean(),role+'_v_mean':v[m,j].mean(),role+'_normalized_gap_mean':gap.mean(),role+'_normalized_gap_variance':gap.var(ddof=0),role+'_mean_gap_component':.5*gap.mean()**2,role+'_gap_variance_component':.5*gap.var(ddof=0)})
        roirows.append(d)
        for seed,z in [('SM',p),('CB',q)]:
            d=fit(z[:,j],'ROI_rawFC_'+seed+'_'+str(ar.roi_id));d.update(meta,seed=seed);rawrows.append(d)
    rd=family(roirows);rf=family(rawrows);save(rf,out/'ROI_raw_FC_results.csv')
    for seed in ['SM','CB']:
        cols=['roi_id','estimate','ci_low','ci_high','p','q_BH','family_n','HC_adjusted_mean','PD_adjusted_mean'];d=rf[rf.seed.eq(seed)][cols]
        rd=rd.merge(d.rename(columns={c:seed+'_rawFC_'+c for c in cols if c!='roi_id'}),on='roi_id',validate='one_to_one')
    save(rd,out/'Supplementary_Data_S18_profile_target_regions.csv')
    assoc=[]
    for role,g in [('HC',0),('PD',1)]:
        m=r.group.eq(g).to_numpy();rr=r.loc[m]
        for name in ['compression','gradient_distance_3d','gradient_distance_10d']:
            y=feat.loc[m,name].to_numpy();y=(y-y.mean())/y.std(ddof=1);pred=dist.sum(1)[m];pred=(pred-pred.mean())/pred.std(ddof=1)
            xx=np.column_stack([np.ones(m.sum()),pred,rr[['age','sex','meanFD']].to_numpy()]);rob=sm.OLS(y,xx).fit().get_robustcov_results(cov_type='HC3',use_t=True);ci=rob.conf_int()[1]
            assoc.append(dict(role=role,outcome='signed_SM_minus_CB_G3_difference' if name=='compression' else name,predictor='cosine_distance',n=int(m.sum()),standardized_beta=rob.params[1],ci_low=ci[0],ci_high=ci[1],p=rob.pvalues[1],covariates='age + sex + meanFD',standard_errors='HC3',interpretation='descriptive correspondence from the same FC input, not causal mediation'))
    save(family(assoc),out/'Supplementary_Data_S19_profile_geometry_correspondence.csv')
    disp=rd.groupby('network')[['HC_mean_gap_component','PD_mean_gap_component','HC_gap_variance_component','PD_gap_variance_component']].sum().reset_index();disp['mean_gap_change']=disp.PD_mean_gap_component-disp.HC_mean_gap_component;disp['dispersion_change']=disp.PD_gap_variance_component-disp.HC_gap_variance_component;disp['interpretation']='unadjusted descriptive decomposition; between-participant dispersion of normalized seed differences; not adjusted causal attribution'
    save(disp,out/'Supplementary_Data_S20_profile_distance_mean_dispersion.csv')
    return dict(cosine_distance_identity_max_error=identity,network_group_coefficient_additivity=float(abs(sum(z['estimate'] for z in rows)-fit(dist.sum(1),'total')['estimate'])))
def compare_outputs(out):
    checks=[]
    for ref in sorted((VERIFY/'reference').glob('Supplementary_Data_*.csv')):
        p=out/ref.name;old=pd.read_csv(ref);new=pd.read_csv(p);missing=[c for c in old if c not in new]
        check=dict(file=ref.name,rows_expected=len(old),rows_new=len(new),missing_columns=missing,columns={})
        if len(old)!=len(new) or missing:check['pass']=False;checks.append(check);continue
        for c in old:
            x=old[c];y=new[c]
            if pd.api.types.is_numeric_dtype(x) and not pd.api.types.is_bool_dtype(x):
                xx=x.to_numpy(float);yy=pd.to_numeric(y,errors='coerce').to_numpy(float);ok=np.isclose(xx,yy,rtol=2e-8,atol=2e-10,equal_nan=True);finite=np.isfinite(xx)&np.isfinite(yy)
                check['columns'][c]=dict(pass_=bool(ok.all()),mismatch_n=int((~ok).sum()),max_abs_error=float(np.max(abs(xx[finite]-yy[finite]))) if finite.any() else 0.)
            else:
                ok=x.fillna('').astype(str).eq(y.fillna('').astype(str));check['columns'][c]=dict(pass_=bool(ok.all()),mismatch_n=int((~ok).sum()))
        check['pass']=all(v['pass_'] for v in check['columns'].values());checks.append(check)
    return checks
def main():
    pa=argparse.ArgumentParser();pa.add_argument('--output-dir',type=Path,required=True);pa.add_argument('--skip-gradient-rebuild',action='store_true',help='Use frozen gradients; default reconstructs all 947 from FC');args=pa.parse_args();out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=True);start=time.time()
    manifest=json.loads((INPUT/'input_manifest.json').read_text(encoding='utf-8'))
    for name,info in manifest['packed_files'].items():assert sha(INPUT/name)==info['sha256'],'Packaged input hash mismatch: '+name
    reg=pd.read_csv(INPUT/'LOCAL_records.csv');fc=np.load(INPUT/'LOCAL_fc_all.npy',mmap_mode='r');saved=np.load(INPUT/'LOCAL_saved_gradients.npy',mmap_mode='r');ref=np.load(INPUT/'reference_gradients.npy')
    rawatlas,atlas=atlas_from_source();save(atlas,out/'Supplementary_Data_S01_atlas_and_networks.csv')
    hcmean=np.asarray(fc[reg.role.eq('HC')]).mean(0);savedmean=np.load(INPUT/'reference_mean_hc_fc.npy');meanerr=float(np.max(abs(hcmean-savedmean)));newref,lambdas=embed(hcmean);referr=float(np.max(abs(align(newref,ref)-ref)))
    gchecks=dict(reference_mean_FC_max_abs_error=meanerr,reference_embedding_after_orthogonal_alignment_max_abs_error=referr,reference_n_HC=int(reg.role.eq('HC').sum()),algorithm='normalized-angle affinity; row90 percentile affinity threshold then symmetric average; alpha0.5,t1; centered orthogonal Procrustes, no scaling',all_recomputed=not args.skip_gradient_rebuild)
    if args.skip_gradient_rebuild:grad=saved
    else:
        grad=np.lib.format.open_memmap(out/'LOCAL_recomputed_gradients.npy',mode='w+',dtype=np.float64,shape=saved.shape);errs=[]
        for i in range(len(reg)):
            grad[i]=align(embed(fc[i])[0],ref);errs.append(float(np.max(abs(grad[i]-saved[i]))))
            if (i+1)%150==0:print('Reconstructed gradients',i+1,'/947',flush=True)
        grad.flush();gchecks['max_abs_error_saved_gradients']=max(errs);gchecks['all_within_1e7']=bool(max(errs)<1e-7);save(pd.DataFrame({'row':np.arange(len(errs)),'max_abs_error':errs}),out/'LOCAL_gradient_reconstruction_errors.csv')
    # Statistical analysis uses frozen saved gradient coordinates, preserving exact legacy values;
    # FC->gradient reconstruction above is separately verified rather than silently replacing inputs.
    feat,ps,qs,t=extract_all(fc,saved,reg,atlas,out)
    original=pd.read_csv(VERIFY/'LOCAL_original_anatomical_all.csv').set_index(['role','person_id']).loc[pd.MultiIndex.from_frame(reg[['role','person_id']])]
    featurechecks={c:float(np.nanmax(abs(feat[c].to_numpy()-original[c].to_numpy()))) for c in METRICS}
    (out/'system_feature_verification.json').write_text(json.dumps(featurechecks,indent=2),encoding='utf-8')
    print('System features ready: 947 records, 13 metrics; max error',max(featurechecks.values()),flush=True)
    fc_analyses(fc,reg,rawatlas,atlas,feat,out);profilechecks=profile_analyses(reg,feat,atlas,ps,qs,t,out)
    comparisons=compare_outputs(out);report=dict(status='PASS' if all(c['pass'] for c in comparisons) and max(featurechecks.values())<1e-10 and gchecks.get('all_within_1e7',True) else 'DISCREPANCY',runtime_seconds=time.time()-start,n_role_records=len(reg),gradient=gchecks,system_features_max_errors=featurechecks,profile_identities=profilechecks,tables=comparisons,input_manifest_sha256=sha(INPUT/'input_manifest.json'),run_script_sha256=sha(__file__),versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,pandas=pd.__version__,statsmodels=sm.__version__),source_files_changed=False,reference_results_used_for_fitting=False)
    (out/'verification.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8');(VERIFY/'verification.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({k:report[k] for k in ['status','runtime_seconds','n_role_records']}));print([(c['file'],c['pass']) for c in comparisons])
if __name__=='__main__':main()
