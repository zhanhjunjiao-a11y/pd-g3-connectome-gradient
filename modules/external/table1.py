"""Rebuild Table 1 using shared clinical records and paper-specific frozen rosters/FD."""
from common import PKG,INPUTS,VERIFY,compare_frames
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
from scipy import stats


def run(output):
    output.mkdir(parents=True,exist_ok=True)
    clinical=pd.read_csv(INPUTS/'table1_clinical.csv').set_index('record_key')
    inputs=pd.read_csv(INPUTS/'table1_analysis_inputs.csv')
    inputs=inputs.loc[inputs.paper.eq('gradient')]
    policy=json.loads((PKG/'inputs_private/core/paired_analysis_policy.json').read_text(encoding='utf-8'))
    active_records=set(policy['included_record_keys']);assert len(active_records)==235
    inputs=inputs.loc[~inputs.role.isin(['paired_OFF','paired_ON'])|inputs.record_key.isin(active_records)].copy()
    flags=pd.read_csv(INPUTS/'table1_cohorts.csv').set_index('record_key')
    frames={}
    for role in ['disease_PD','disease_HC','paired_OFF','paired_ON']:
        x=inputs.loc[inputs.role.eq(role)].set_index('record_key')
        assert x.index.is_unique
        d=clinical.loc[x.index].copy()
        d['age']=x.age_frozen;d['sex']=x.sex_frozen;d['meanFD']=x.meanFD_frozen
        d['education']=x.education_curated
        frames[role]=d
    pd0,hc,off,on=(frames[k] for k in ['disease_PD','disease_HC','paired_OFF','paired_ON'])
    assert [len(x) for x in [pd0,hc,off,on]]==[334,141,235,235]
    assert off.index.equals(on.index)
    assert set(pd0.index)==set(flags.index[flags.gradient_disease_PD])
    assert set(hc.index)==set(flags.index[flags.gradient_disease_HC])
    assert set(off.index)==active_records
    assert active_records.issubset(set(flags.index[flags.gradient_clinical_pair]))
    assert active_records==set(flags.index[flags.gradient_adjusted_imaging_pair])
    paired_motion=off.meanFD.notna()&on.meanFD.notna()
    assert paired_motion.sum()==235
    assert (pd0.education<=30).all() # confirmed 67 -> 12 overlay, not a new exclusion
    cols=['Characteristic','PD_334','HC_141','PD_HC_P','Paired_OFF_235','Paired_ON_235','OFF_ON_P']
    rows=[];numeric=[]
    def pformat(p):return '<0.001' if p<.001 else f'{p:.3f}'
    def desc(series,label,group,median=False,dp=1):
        x=pd.to_numeric(series,errors='coerce').dropna()
        rec={'variable':label,'group':group,'n':len(x),'mean':x.mean(),'sd':x.std(ddof=1),'median':x.median(),'q1':x.quantile(.25),'q3':x.quantile(.75)}
        numeric.append(rec)
        return f'{rec["median"]:.{dp}f} [{rec["q1"]:.{dp}f}, {rec["q3"]:.{dp}f}]' if median else f'{rec["mean"]:.{dp}f} ± {rec["sd"]:.{dp}f}'
    sexcounts=[int(x.sex.eq(1).sum()) for x in [pd0,hc,off]]
    sexp=stats.chi2_contingency([[sexcounts[0],len(pd0)-sexcounts[0]],[sexcounts[1],len(hc)-sexcounts[1]]],correction=False).pvalue
    rows.append(['Male sex, n (%)',f'{sexcounts[0]} ({sexcounts[0]/len(pd0)*100:.1f})',f'{sexcounts[1]} ({sexcounts[1]/len(hc)*100:.1f})',pformat(sexp),f'{sexcounts[2]} ({sexcounts[2]/len(off)*100:.1f})','—','—'])
    specs=[('Age, years','age',None,True,False,1),('Education, yearsᵃ','education',None,True,False,1),
       ('Mean framewise\ndisplacement, mmᵇ','meanFD','meanFD',True,False,3),
       ('Disease duration, yearsᶜ','disease_duration',None,False,False,1),
       ('Hoehn–Yahr stageᵈ','HY_OFF','HY_ON',False,True,1),
       ('MDS-UPDRS Part III','UPDRSIII_OFF','UPDRSIII_ON',False,False,1),
       ('PIGD score','PIGD_OFF','PIGD_ON',False,False,1),('TD score','TD_OFF','TD_ON',False,False,1),
       ('LEDD, mg/day','LEDD',None,False,False,1),('MMSEᵉ','MMSE',None,True,False,1),
       ('MoCAᵉ','MoCA',None,True,False,1),('HAMAᵉ','HAMA',None,True,False,1),('HAMD-17ᵉ','HAMD',None,True,False,1)]
    tests=[{'variable':'sex','contrast':'PD-HC','method':'Pearson chi-square without continuity correction','p':float(sexp)}]
    for label,field,onfield,compare_hc,median,dp in specs:
        pv=desc(pd0[field],field,'PD334',median,dp)
        hv=desc(hc[field],field,'HC141',median,dp) if compare_hc else '—'
        ov=desc(off[field],field,'OFF235',median,dp)
        nv=desc(on[onfield],onfield,'ON235',median,dp) if onfield else '—'
        gp='—';pp='—'
        if compare_hc:
            fit=stats.ttest_ind(pd0[field].dropna(),hc[field].dropna(),equal_var=False);gp=pformat(fit.pvalue)
            tests.append({'variable':field,'contrast':'PD-HC','method':'Welch t','p':float(fit.pvalue)})
        if onfield and not median:
            keep=off[field].notna()&on[onfield].notna()
            fit=stats.ttest_rel(on.loc[keep,onfield],off.loc[keep,field]);pp=pformat(fit.pvalue)
            tests.append({'variable':field,'contrast':'ON-OFF','method':'paired t','n':int(keep.sum()),'p':float(fit.pvalue)})
        rows.append([label,pv,hv,gp,ov,nv,pp])
    table=pd.DataFrame(rows,columns=cols)
    table.to_csv(output/'Table_1_aggregate_values.csv',index=False,encoding='utf-8-sig')
    pd.DataFrame(numeric).to_csv(output/'Table_1_numeric_summaries.csv',index=False)
    pd.DataFrame(tests).to_csv(output/'Table_1_tests.csv',index=False)
    expected=pd.read_csv(VERIFY/'reference/Table_1_aggregate_values.csv',keep_default_na=False)
    result=compare_frames(table,expected,['Characteristic'])
    (output/'verification.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    return table,result


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args();_,result=run(args.output_dir.resolve());print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
