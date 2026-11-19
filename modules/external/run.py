"""Local portable PPMI and Table 1 reproduction; references are used only to verify."""
from common import INPUTS,VERIFY,compare_frames,sha
from pathlib import Path
import argparse,json,subprocess,sys,time,traceback
import pandas as pd
import endpoints,descriptives,table1


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args();out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
    started=time.time();checks={}
    try:
        # Model fitting/reconstruction occurs before any expected output is read.
        endpoints.run(out/'endpoints')
        descriptives.run(out/'descriptives')
        subprocess.run([sys.executable,'-X','utf8',str(Path(__file__).resolve().parents[2]/'guarded_runner.py'),str(Path(__file__).with_name('map_validation.py')),'--output-dir',str(out/'maps')],check=True)
        _,table_check=table1.run(out/'table1');checks['Table1']=table_check
        maps=pd.read_csv(out/'maps/map_agreement.csv')
        mapcols=['model','gradient','n','n_PD','n_HC','covariates','n_regions','spearman_rho','p_wild_two_sided','q_BH_3','null_hypothesis','n_null','null_seed']
        maps[mapcols].to_csv(out/'maps/Supplementary_Data_S24_PPMI_effect_map_agreement.csv',index=False)
        regional=pd.read_csv(out/'maps/regional_effects_all_849.csv')
        regional.to_csv(out/'maps/Supplementary_Data_S25_PPMI_regional_gradients.csv',index=False)
        mapping={14:('endpoints',['model','metric']),16:('endpoints',['model','metric']),22:('descriptives',['group','variable','level']),23:('endpoints',['model','metric']),24:('maps',['model','gradient']),25:('maps',['model','gradient','roi_id'])}
        for number,(folder,keys) in mapping.items():
            ref=next((VERIFY/'reference').glob(f'Supplementary_Data_S{number:02d}_*.csv'))
            got=pd.read_csv(out/folder/ref.name);expected=pd.read_csv(ref)
            checks[f'S{number:02d}']=compare_frames(got,expected,keys,tolerance=1e-9)
        for key,rel,reference,keys in [
          ('endpoints','endpoints/endpoints.csv','endpoints.csv',['sub_id']),
          ('all_64_models','endpoints/eight_metric_results.csv','eight_metric_results.csv',['model','metric']),
          ('full_map_bootstrap','maps/map_agreement.csv','map_agreement_full.csv',['model','gradient'])]:
            checks[key]=compare_frames(pd.read_csv(out/rel),pd.read_csv(VERIFY/'reference'/reference),keys,tolerance=1e-9)
        status='PASS' if all(x['status']=='PASS' for x in checks.values()) else 'FAIL'
        summary={'status':status,'scope':'Reconstructs PPMI endpoints from frozen FC/gradients, refits all scalar models and full wild-bootstrap maps, recomputes S22 from individual metadata and all 14 Table1 rows from shared clinical records and paper rosters. Does not rerun raw MRI preprocessing.','checks':checks,'elapsed_seconds':time.time()-started,'input_manifest_sha256':sha(INPUTS/'input_manifest.json')}
    except Exception:
        summary={'status':'FAIL','checks':checks,'error':traceback.format_exc(),'elapsed_seconds':time.time()-started}
    (out/'verification.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(summary,ensure_ascii=False,indent=2),flush=True)
    return 0 if summary['status']=='PASS' else 1


if __name__=='__main__':sys.exit(main())
