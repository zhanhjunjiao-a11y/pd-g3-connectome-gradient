"""Recompute S22 from frozen individual/metadata inputs; no stored numeric summaries."""
from common import INPUTS
import numpy as np
import pandas as pd


def run(output):
    output.mkdir(parents=True,exist_ok=True)
    reg=pd.read_csv(INPUTS/'registry.csv',dtype={'sub_id':str})
    reg['PATNO']=reg.sub_id.str.removeprefix('sub-')
    history=pd.read_csv(INPUTS/'clinical_history.csv',dtype={'PATNO':str}).set_index('PATNO')
    mds=pd.read_csv(INPUTS/'clinical_mds.csv',dtype={'PATNO':str}).set_index('PATNO')
    scan=pd.read_csv(INPUTS/'scan_mapping.csv',dtype={'PATNO':str}).set_index('PATNO')
    pd_stc=pd.read_csv(INPUTS/'pd_stc_metadata.csv').set_index('subject')
    hc_stc=pd.read_csv(INPUTS/'hc_metadata_audit.csv',dtype={'subject':str})
    hc_stc['subject']='sub-'+hc_stc.subject.str.removeprefix('sub-')
    hc_stc=hc_stc.set_index('subject')
    series={g:pd.read_csv(INPUTS/f'{prefix}_series_metadata.tsv',sep='\t').drop_duplicates('subject',keep='last').set_index('subject')
            for g,prefix in [('PD','pd'),('HC','hc')]}
    definitions=pd.read_csv(INPUTS/'s22_definitions.csv',keep_default_na=False)
    rows=[]
    flagvars={'age_corrected':lambda d:(d.age_final-d.age_original).abs().gt(1e-8),
              'age_ambiguous':lambda d:d.age_ambiguous,
              'blinded_clinical_SITE_available':lambda d:d.SITE_final.notna(),
              'later_screening_failure_label':lambda d:d.screen_flag}
    for definition in definitions.to_dict('records'):
        group,var,level=definition['group'],definition['variable'],definition['level']
        d=reg.loc[reg.group.eq(1 if group=='PD' else 0)];n=len(d)
        s=series[group].loc[d.sub_id]
        r={**definition,'cohort_n':n,'available_n':0,'missing_n':n,'count':np.nan,'percent':np.nan,
           **{k:np.nan for k in ['mean','sd','median','q1','q3','min','max']}}
        if var in flagvars:
            values=flagvars[var](d);count=int(values.sum())
            r.update(available_n=n,missing_n=0,count=count,percent=100*count/n)
        elif var in ['sex','scanner_manufacturer','scanner_model','scan_visit_anchoring','TR_distribution','saved_timepoints_distribution','SliceTimingCorrected']:
            if var=='sex': values=np.where(d.sex.eq(1),'male','female')
            elif var=='scanner_manufacturer':values=d.vendor.to_numpy()
            elif var=='scanner_model':values=d.scanner_exact.to_numpy()
            elif var=='scan_visit_anchoring':
                a=scan.loc[d.PATNO];values=np.where(a.eligible_for_clinical_join.eq(True),a.EVENT_ID,'unresolved')
            elif var=='TR_distribution':values=s.TR.astype(float).map(str).to_numpy()
            elif var=='saved_timepoints_distribution':values=s.n_timepoints.astype(int).map(str).to_numpy()
            else: values=(pd_stc.loc[d.sub_id,'SliceTimingCorrected'] if group=='PD' else hc_stc.loc[d.sub_id,'JSON_STC']).map(str).to_numpy()
            values=pd.Series(values).dropna().astype(str);values=values.loc[~values.isin(['','None'])]
            count=int(values.eq(level).sum());r.update(available_n=len(values),missing_n=n-len(values),count=count,percent=100*count/n)
        else:
            if var=='age_used_in_imaging_models':values=d.age_final
            elif var=='meanFD':values=d.meanFD_final
            elif var=='TR':values=s.TR
            elif var=='saved_timeseries_timepoints':values=s.n_timepoints
            elif var=='nominal_retained_duration':values=s.TR*s.n_timepoints/60
            elif var in ['symptom_duration_years','diagnosis_duration_years','education_years']:values=history.loc[d.PATNO,var]
            elif var in ['untreated_U3','OFF_U3','ON_U3','non_ON_U3','non_ON_HY']:values=mds.loc[d.PATNO,var]
            else:raise ValueError('Unmapped S22 variable '+var)
            values=pd.to_numeric(values,errors='coerce').dropna();values=values.loc[np.isfinite(values)]
            r.update(available_n=len(values),missing_n=n-len(values))
            if len(values):r.update(mean=values.mean(),sd=values.std(ddof=1),median=values.median(),q1=values.quantile(.25),q3=values.quantile(.75),min=values.min(),max=values.max())
        rows.append(r)
    cols=['group','variable','unit','level','cohort_n','available_n','missing_n','count','percent','mean','sd','median','q1','q3','min','max','status','note']
    result=pd.DataFrame(rows)[cols]
    result.to_csv(output/'Supplementary_Data_S22_PPMI_characteristics.csv',index=False,encoding='utf-8-sig')
    return result
