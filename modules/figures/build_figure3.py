"""Figure 3: display frozen observed data and existing model estimates.

Run: python code/build_figure3.py
Requirements: numpy, pandas, matplotlib, pillow. No models are refitted.
All imaging effects are ON minus OFF; clinical improvement is OFF minus ON.
"""
from pathlib import Path
import os, json, hashlib
ROOT = Path(os.environ['GRADIENT_FIG_OUTPUT'])
DATA, OUT = Path(os.environ['GRADIENT_FIG_INPUT']), ROOT/'figures'
QA = ROOT/'qa/figure3'
QA.mkdir(parents=True, exist_ok=True)
OUT.mkdir(exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR', str(QA/'mplconfig'))
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

contract = {
 'claim':'Acute challenge is accompanied by changes in complete-system geometry, selected regional G3 values and within-cerebellar FC; these do not establish normalization or explain individual motor improvement.',
 'panels':{'A':'235 observed OFF/ON clinical pairs','B':'Three complete-system gradient geometry estimates, BH3','C':'19 disease-difference regions; state q values from all 283 G3 regions','D':'Four complete-system connectivity estimates, existing BH12'},
 'archetype':'Paired observations and coordinated effect estimates',
 'size_mm':[185,205], 'dpi':600,
 'statistics':'Frozen HC3 OLS intercepts at uncentered DeltaFD=0, n=235, df=233; pointwise 95% CI; all tests two-sided.',
 'clinical':'All 235 raw observations, connected pairs; arithmetic group means and stored unadjusted paired-change Student-t CI, df=234.',
 'transforms':'None to stored estimates or CI. Jitter changes x position only. FC and dimensionless angle similarity use separate axes.',
 'not_shown':'Redundant raw paired FC distributions; original sources and SI S08 retained. Target decomposition remains supplementary.'
}
(QA/'figure_contract.json').write_text(json.dumps(contract,indent=2),encoding='utf-8')

raw=pd.read_csv(DATA/'Figure_3A_anonymous_motor_scores.csv')
summary=pd.read_csv(DATA/'Figure_3A_motor_score_summary.csv').iloc[0]
geo=pd.read_csv(DATA/'Figure_3B_complete_system_geometry.csv')
reg=pd.read_csv(DATA/'Figure_3C_G3_regional_changes.csv')
fc=pd.read_csv(DATA/'Figure_3D_full_system_changes.csv')
assert len(raw)==235 and len(geo)==3 and len(reg)==19 and len(fc)==4
assert np.allclose(raw.UPDRS_ON-raw.UPDRS_OFF,raw.ON_minus_OFF)
assert np.isclose(raw.ON_minus_OFF.mean(),summary.ON_minus_OFF_mean)
assert set(geo.n)==set(reg.n)==set(fc.n)=={235}
assert set(geo.family_n)=={3} and set(fc.family_n)=={12}
assert set(geo.df)=={233} and set(geo.standard_errors)=={'HC3'}
assert (reg.q_BH_283_within_gradient<.05).sum()==4
assert (geo.q_BH3<.05).all() and (fc.q_BH<.05).sum()==1
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
 'font.size':7.2,'axes.labelsize':7.2,'xtick.labelsize':7,'ytick.labelsize':7,
 'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6,
 'xtick.major.width':.6,'ytick.major.width':.6,'xtick.major.size':2.5,'ytick.major.size':2.5,
 'svg.fonttype':'none','pdf.fonttype':42,'savefig.facecolor':'white'})
BLUE,ORANGE,CB,DARK,GREY='#3B68A8','#C24A2F','#8B6BA8','#343B43','#ADB3B9'
fig=plt.figure(figsize=(7.2834645669,8.0708661417))  # 185 x 205 mm
fig.text(.025,.982,'Acute levodopa challenge: clinical and imaging changes',fontsize=10,weight='bold',va='top')
fig.text(.025,.956,'Clinical and motion-adjusted imaging pairs: n = 235',fontsize=7.5,color=DARK)

def header(letter,title,x,y):
 fig.text(x,y,letter,weight='bold',fontsize=10,va='top')
 fig.text(x+.031,y,title,weight='bold',fontsize=8,va='top')
def signed(v,dp=5): return f'{v:.{dp}f}'.replace('-','−')
def qfmt(q): return '< 0.001' if q<.001 else f'= {q:.4f}'
def forest(ax,b,lo,hi,q,color=BLUE):
 ax.axvline(0,color=GREY,lw=.6,zorder=0)
 ax.errorbar(b,0,xerr=[[b-lo],[hi-b]],fmt='o',color=color,ecolor=color,
             mfc=color if q<.05 else 'white',mec=color,ms=4.5,lw=1,capsize=2)
 ax.set_ylim(-.7,.7); ax.set_yticks([])
 for side in ['top','right','left']: ax.spines[side].set_visible(False)

# A: all raw clinical observations, with real subject pairing.
header('A','Motor scores',.025,.919)
fig.text(.056,.894,'Observed OFF / ON pairs',fontsize=7,color=DARK)
ax=fig.add_axes([.095,.638,.235,.238])
jitter=((np.arange(len(raw))*0.61803398875)%1-.5)*.13  # x-only offsets; all measured values preserved
for j,off,on in zip(jitter,raw.UPDRS_OFF,raw.UPDRS_ON):
 ax.plot([j,1+j],[off,on],color=GREY,alpha=.20,lw=.35,zorder=0)
for x,key,color in [(0,'UPDRS_OFF',BLUE),(1,'UPDRS_ON',ORANGE)]:
 vals=raw[key].to_numpy()
 ax.scatter(x+jitter,vals,s=4.5,alpha=.45,c=color,linewidths=0)
 ax.plot([x-.18,x+.18],[vals.mean()]*2,color='black',lw=1.5)
ax.set_xticks([0,1],['OFF','ON']); ax.set_xlim(-.4,1.4)
ax.set_ylabel('MDS-UPDRS-III score (points)'); ax.set_ylim(-2,max(raw.UPDRS_OFF.max(),raw.UPDRS_ON.max())+5)
fig.text(.056,.604,f'Mean ON − OFF: {signed(summary.ON_minus_OFF_mean,1)} points\n95% CI {signed(summary.paired_mean_CI95_low,1)} to {signed(summary.paired_mean_CI95_high,1)} (unadjusted)',fontsize=7.1,va='top',linespacing=1.5)

# B: geometry gets a dedicated, prominent panel. No raw means relabeled as effects.
header('B','Complete-system gradient geometry',.395,.919)
fig.text(.426,.894,'SMN (30 regions) / cerebellum (31); BH: 3 tests',fontsize=7,color=DARK)
labels=['Signed G3 mean difference (SMN − cerebellum)','3D centroid distance','10D centroid distance']
for k,(row,label) in enumerate(zip(geo.itertuples(index=False),labels)):
 y=.864-k*.094
 fig.text(.426,y,label,fontsize=7.2,va='bottom')
 a=fig.add_axes([.43,y-.045,.525,.034])
 forest(a,row.adjusted_change,row.ci_low,row.ci_high,row.q_BH3)
 a.set_xlim(-.010,.010); a.set_xticks([-.01,0,.01])
 if k<2:
  a.set_xticklabels([]); a.tick_params(bottom=False); a.spines['bottom'].set_visible(False)
 else: a.set_xticklabels(['−0.010','0','0.010'])
 fig.text(.426,y-(.077 if k==2 else .062),f'{signed(row.adjusted_change)} [{signed(row.ci_low)}, {signed(row.ci_high)}]; q {qfmt(row.q_BH3)}',fontsize=7,va='top')
fig.text(.69,.565,'Adjusted ON − OFF (gradient units)',ha='center',fontsize=7.2)

# C: all selected disease regions, with q from the full regional state family.
header('C','Regional G3 changes',.025,.526)
fig.text(.056,.501,'19 disease-difference regions; BH across 283 G3 tests',fontsize=7,color=DARK)
names={'A4ul_l':'L A4ul','A4ul_r':'R A4ul','A4t_l':'L A4t','A4t_r':'R A4t',
 'A4ll_l':'L A4ll','A4ll_r':'R A4ll','A1/2/3tru_l':'L A1/2/3tru','A1/2/3tru_r':'R A1/2/3tru',
 'G_r':'R G','rTtha_l':'L rTtha','Left_VI':'L lobule VI','Vermis_VI':'Vermis VI','Right_VI':'R lobule VI',
 'Left_CrusI':'L Crus I','Right_CrusI':'R Crus I','Right_CrusII':'R Crus II',
 'Left_VIIIb':'L lobule VIIIb','Vermis_VIIIb':'Vermis VIIIb','Vermis_IX':'Vermis IX'}
ax=fig.add_axes([.183,.137,.32,.342])
ax.axvline(0,color=GREY,lw=.6)
for i,row in enumerate(reg.itertuples(index=False)):
 b,lo,hi,q=row.adjusted_ON_minus_OFF,row.ci95_low,row.ci95_high,row.q_BH_283_within_gradient
 color=BLUE if q<.05 else '#707982'
 ax.errorbar(b,i,xerr=[[b-lo],[hi-b]],fmt='o',color=color,ecolor=color,mfc=color if q<.05 else 'white',ms=3.5,capsize=1.5,lw=.8)
ax.set_yticks(range(19),[names[v] for v in reg.roi_name]); ax.tick_params(axis='y',length=0)
ax.set_ylim(18.7,-.7); ax.set_xlim(-.012,.019); ax.set_xticks([-.01,0,.01],['−0.01','0','0.01'])
ax.spines['left'].set_visible(False)
ax.set_xlabel('Adjusted ON − OFF (G3 units)')

# D: the three Fisher-z FC changes share an axis; angle similarity has its own axis.
header('D','Complete-system connectivity',.555,.526)
fig.text(.586,.501,'Adjusted ON − OFF; BH across 12 tests',fontsize=7,color=DARK)
fig.text(.586,.473,'Mean FC (Fisher z)',fontsize=7.4,weight='bold')
dlabels=['Within SMN','Within cerebellum','Between systems']
ax_fc=fig.add_axes([.790,.280,.180,.174])
ax_fc.axvline(0,color=GREY,lw=.6,zorder=0)
ax_fc.set_ylim(-.5,2.5)
ax_fc.set_yticks([])
ax_fc.set_xlim(-.025,.065)
ax_fc.set_xticks([-.02,0,.06],['−0.02','0','0.06'])
for side in ['top','right','left']: ax_fc.spines[side].set_visible(False)
for k,(row,label) in enumerate(zip(fc.iloc[:3].itertuples(index=False),dlabels)):
 y_data=2-k
 y_fig=.280+((y_data+.5)/3)*.174
 ax_fc.errorbar(row.beta,y_data,xerr=[[row.beta-row.ci_low],[row.ci_high-row.beta]],
     fmt='o',color=BLUE,ecolor=BLUE,mec=BLUE,mfc=BLUE if row.q_BH<.05 else 'white',
     ms=4.5,lw=1,capsize=2)
 fig.text(.586,y_fig,label,fontsize=7.4,va='center')
 fig.text(.586,y_fig-.021,
     f'{signed(row.beta,4)} [{signed(row.ci_low,4)}, {signed(row.ci_high,4)}]; q {qfmt(row.q_BH)}',
     fontsize=7,va='center',color=DARK)

fig.text(.586,.228,'External-profile similarity',fontsize=7.4,weight='bold')
fig.text(.586,.207,'Normalized-angle similarity (dimensionless)',fontsize=7,color=DARK)
row=fc.iloc[3]
ax_profile=fig.add_axes([.590,.153,.380,.033])
forest(ax_profile,row.beta,row.ci_low,row.ci_high,row.q_BH)
ax_profile.set_xlim(-.015,.055)
ax_profile.set_xticks([-.01,0,.02,.04],['−0.01','0','0.02','0.04'])
fig.text(.586,.112,
 f'{signed(row.beta,4)} [{signed(row.ci_low,4)}, {signed(row.ci_high,4)}]; q {qfmt(row.q_BH)}',
 fontsize=7,va='center',color=DARK)
fig.text(.025,.058,'B–D: adjusted estimates and pointwise 95% CI; n = 235. Filled markers indicate BH q < 0.05.',fontsize=7.1)
fig.text(.025,.036,'Imaging change = ON − OFF. Motor improvement = OFF − ON. Fixed scan order: OFF → ON.',fontsize=7.1)
fig.savefig(OUT/'Figure_3_revised.svg')
fig.savefig(OUT/'Figure_3_revised.pdf')
fig.savefig(OUT/'Figure_3_revised.png',dpi=600)
fig.savefig(OUT/'Figure_3_preview.png',dpi=300)
plt.close(fig)
qa={'input_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in DATA.glob('Figure_3*.csv')},
 'matplotlib_version':mpl.__version__,'numpy_version':np.__version__,'pandas_version':pd.__version__,
 'checks':'235 raw pairs, 235 imaging pairs, 3 BH3 geometry estimates, 19 regional estimates, 4 BH12 connectivity estimates; all source checks passed.'}
(QA/'build_checks.json').write_text(json.dumps(qa,indent=2),encoding='utf-8')
print(qa['checks'])
