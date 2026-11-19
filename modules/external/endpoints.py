"""Reconstruct eight endpoints from frozen gradients/FC and refit all eight models."""
from common import INPUTS
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests


def run(output):
    output.mkdir(parents=True,exist_ok=True)
    reg=pd.read_csv(INPUTS/'registry.csv',dtype={'sub_id':str,'SITE_final':str})
    order=[x.strip() for x in (INPUTS/'subject_order.txt').read_text().splitlines() if x.strip() and not x.startswith('#')]
    assert reg.sub_id.tolist()==order and len(reg)==438
    g=np.load(INPUTS/'gradients.npy');assert g.shape==(438,283,10) and np.isfinite(g).all()
    atlas=pd.read_csv(INPUTS/'atlas.csv');assert np.array_equal(atlas.Label,np.arange(1,284))
    si=np.flatnonzero(atlas.category.eq('Somatomotor'));ci=np.flatnonzero(atlas.category.eq('Cerebellum'))
    external=np.setdiff1d(np.arange(283),np.r_[si,ci]);assert (len(si),len(ci),len(external))==(30,31,222)
    d=g[:,si].mean(axis=1)-g[:,ci].mean(axis=1)
    cfg=json.loads((INPUTS/'fixed19_config.json').read_text())
    up=np.array(cfg['cortical_pole_roi_1based'])-1;dn=np.array(cfg['cerebellar_pole_roi_1based'])-1
    assert len(up)==10 and len(dn)==9
    feat=pd.DataFrame({'sub_id':reg.sub_id,'fixed19_g3':g[:,up,2].mean(1)-g[:,dn,2].mean(1),
       'signed_g3_difference':d[:,2],'gradient_distance_3d':np.linalg.norm(d[:,:3],axis=1),'gradient_distance_10d':np.linalg.norm(d,axis=1)})
    edges=np.load(INPUTS/'fc_edges.npy',mmap_mode='r');assert edges.shape==(438,39903)
    ij=np.triu_indices(283,1);us=np.triu_indices(30,1);uc=np.triu_indices(31,1)
    vals=[]
    for ev in edges:
        m=np.zeros((283,283));m[ij]=ev;m[(ij[1],ij[0])]=ev
        a=m[np.ix_(si,external)].mean(0);b=m[np.ix_(ci,external)].mean(0)
        cosine=float(a@b/(np.linalg.norm(a)*np.linalg.norm(b)))
        vals.append([m[np.ix_(si,si)][us].mean(),m[np.ix_(ci,ci)][uc].mean(),m[np.ix_(si,ci)].mean(),1-np.arccos(np.clip(cosine,-1,1))/np.pi])
    feat[['within_sm_fc','within_cb_fc','between_fc','profile_angle']]=np.array(vals)
    feat.to_csv(output/'endpoints.csv',index=False)
    counts=pd.crosstab(reg.SITE_final,reg.group_name)
    sites=counts.loc[(counts.HC>0)&(counts.iPD>0)].index.tolist()
    approved=json.loads((INPUTS/'site_support.json').read_text(encoding='utf-8'))
    assert sites==approved['allowed_sites']
    mask=reg.SITE_final.isin(sites).to_numpy()
    metrics={'fixed19_g3':('Fixed 19-region G3 score','fixed19','gradient coordinate'),
      'signed_g3_difference':('SMN minus cerebellum G3 mean','geometry3','gradient coordinate'),
      'gradient_distance_3d':('3D centroid distance','geometry3','gradient distance'),
      'gradient_distance_10d':('10D centroid distance','geometry3','gradient distance'),
      'within_sm_fc':('SMN within-system FC','FC4','Fisher z'),
      'within_cb_fc':('Cerebellar within-system FC','FC4','Fisher z'),
      'between_fc':('SMN–cerebellar FC','FC4','Fisher z'),
      'profile_angle':('External connectivity-profile similarity','FC4','normalized angular similarity')}
    models={
      'age_sex_motion':(np.ones(len(reg),bool),'age_final',True,None),
      'scanner_adjusted_motion':(np.ones(len(reg),bool),'age_final',True,'scanner_exact'),
      'site_shared_motion':(mask,'age_final',True,'SITE_final'),
      'exclude_ambiguous_age':((~reg.age_ambiguous).to_numpy(),'age_final',True,None),
      'exclude_screen_flag':((~reg.screen_flag).to_numpy(),'age_final',True,None),
      'original_age_motion':(np.ones(len(reg),bool),'age_original',True,None),
      'legacy_age_sex':(np.ones(len(reg),bool),'age_original',False,None),
      'legacy_scanner':(np.ones(len(reg),bool),'age_original',False,'scanner_exact')}
    rows=[];diagnostics=[];manual_errors=[]
    for model,(use,age,fd,cat) in models.items():
        d=reg.loc[use].copy();f=feat.loc[use]
        x=pd.DataFrame({'Intercept':np.ones(len(d)),'PD':d.group.to_numpy(),'age':((d[age]-d[age].mean())/d[age].std(ddof=0)).to_numpy(),'sex':d.sex.to_numpy()})
        if fd:x['meanFD']=((d.meanFD_final-d.meanFD_final.mean())/d.meanFD_final.std(ddof=0)).to_numpy()
        if cat:x=pd.concat([x,pd.get_dummies(d[cat],prefix=cat,drop_first=True,dtype=float).reset_index(drop=True)],axis=1)
        X=x.to_numpy(dtype=float);pinv=np.linalg.pinv(X);rank=np.linalg.matrix_rank(X);h=np.sum(X*pinv.T,axis=1)
        assert rank==X.shape[1] and h.max()<.999999,(model,rank,h.max())
        covariates='age, sex'+(', meanFD' if fd else '')+(', scanner model' if cat=='scanner_exact' else ', blinded clinical SITE' if cat else '')
        diagnostics.append({'model':model,'n':len(d),'n_PD':int(d.group.sum()),'n_HC':int((d.group==0).sum()),'rank':rank,'df_residual':len(d)-rank,'max_leverage':float(h.max()),'condition_number':float(np.linalg.cond(X))})
        for metric,(label,family,unit) in metrics.items():
            y=f[metric].to_numpy();fit=sm.OLS(y,X).fit(cov_type='HC3',use_t=True)
            beta=pinv@y;res=y-X@beta;cov=(pinv*(res/(1-h))**2)@pinv.T
            berr=abs(beta[1]-fit.params[1]);serr=abs(np.sqrt(cov[1,1])-fit.bse[1]);assert berr<1e-12 and serr<1e-10
            manual_errors.append({'model':model,'metric':metric,'beta_error':berr,'HC3_SE_error':serr})
            interval=fit.conf_int(alpha=.05)[1]
            rows.append({'model':model,'metric':metric,'metric_label':label,'family':family,'n':len(d),'n_PD':int(d.group.sum()),'n_HC':int((d.group==0).sum()),'df_residual':len(d)-rank,'estimate':float(fit.params[1]),'se_HC3':float(fit.bse[1]),'ci95_low':float(interval[0]),'ci95_high':float(interval[1]),'t':float(fit.tvalues[1]),'p_two_sided':float(fit.pvalues[1]),'q_BH':np.nan,'multiplicity':'Unadjusted; one fixed score' if family=='fixed19' else 'BH across three geometry metrics' if family=='geometry3' else 'BH across four system FC metrics','covariates':covariates,'effect_direction':'PD minus HC','unit':unit,'report_role':'legacy reproduction only' if model.startswith('legacy') else 'primary' if model=='scanner_adjusted_motion' else 'sensitivity'})
    out=pd.DataFrame(rows)
    for (_,family),idx in out.groupby(['model','family']).groups.items():
        if family!='fixed19':out.loc[idx,'q_BH']=multipletests(out.loc[idx,'p_two_sided'],method='fdr_bh')[1]
    out.to_csv(output/'eight_metric_results.csv',index=False)
    pd.DataFrame(diagnostics).to_csv(output/'model_diagnostics.csv',index=False)
    pd.DataFrame(manual_errors).to_csv(output/'manual_HC3_checks.csv',index=False)
    for family,filename in [('FC4','Supplementary_Data_S14_PPMI_complete_system_external.csv'),('fixed19','Supplementary_Data_S16_PPMI_fixed19_score_external.csv'),('geometry3','Supplementary_Data_S23_PPMI_system_geometry.csv')]:
        out.loc[out.family.eq(family)&~out.model.str.startswith('legacy')].to_csv(output/filename,index=False)
    return feat,out


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output-dir',type=Path,required=True)
    run(parser.parse_args().output_dir.resolve())
