"""Figure 2: G2/G3 regional FC, disease geometry, and complete-system FC.

Run with Python, numpy, pandas, scipy, matplotlib, pillow installed.
All points/effects come from bundled saved outputs; no models are fitted here.
"""
from pathlib import Path
import os, json
ROOT = Path(os.environ['GRADIENT_FIG_OUTPUT'])
DATA, OUT, QA = Path(os.environ['GRADIENT_FIG_INPUT']), ROOT / 'figures', ROOT / 'qa' / 'figure2'
OUT.mkdir(exist_ok=True)
QA.mkdir(parents=True, exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR', str(QA / 'mplconfig'))
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

mpl.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 7.5, 'axes.labelsize': 7.5, 'axes.titlesize': 8,
    'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5,
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.linewidth': .65, 'xtick.major.width': .6, 'ytick.major.width': .6,
    'pdf.fonttype': 42, 'svg.fonttype': 'none', 'legend.frameon': False,
    'figure.facecolor': 'white', 'savefig.facecolor': 'white',
})
C = {'HC': '#3B68A8', 'PD': '#C24A2F', 'SMN': '#3C8583', 'CB': '#8563A1', 'effect': '#245E8B'}
DIV = LinearSegmentedColormap.from_list('group_difference', ['#27558A', '#A7C2DB', '#FFFFFF', '#E7B493', '#AF472D'])
obs = pd.read_csv(DATA / 'Figure_2_anonymous_observations.csv')
stats = pd.read_csv(DATA / 'Figure_2_system_effects.csv').set_index('metric')
regional = pd.read_csv(DATA / 'Figure_2A_network_mean_effects.csv')
registry = pd.read_csv(DATA / 'Figure_2_region_registry.csv')
g2 = pd.read_csv(DATA / 'Figure_2B_G2_target_effects.csv')
geo = pd.read_csv(DATA / 'Figure_2C_system_geometry.csv').set_index('metric')
assert obs.groupby('group').size().to_dict() == {'HC': 141, 'PD': 334}
assert len(stats) == 6 and stats.family_n.eq(6).all()
assert len(regional) == 171 and regional.family_n.eq(2547).all()
assert len(g2) == 78 and g2.family_n.eq(2547).all()
assert geo.family_n.eq(7).all()
for frame in [stats, g2, geo]:
    assert frame.n.eq(475).all() and frame.df.eq(470).all()

order = [57,58,59,60,67,68,161,162,164,225,257,258,259,260,261,264,271,272,275]
names = ['L precentral (upper limb)', 'R precentral (upper limb)',
         'L precentral (trunk)', 'R precentral (trunk)', 'L paracentral', 'R paracentral',
         'L postcentral (trunk)', 'R postcentral (trunk)', 'R insula', 'L thalamus',
         'L cerebellum VI', 'Vermis VI', 'R cerebellum VI', 'L Crus I', 'R Crus I',
         'R Crus II', 'L VIIIb', 'Vermis VIIIb', 'Vermis IX']
targets = ['Visual','Somatomotor','DorsAttn','VentAttn','Limbic','Frontoparietal','Default','Subcortex','Cerebellum']
abbr = ['VIS','SMN','DAN','VAN','LIM','FPN','DMN','SUB','CB']
eff = regional.pivot(index='roi', columns='target_network', values='standardized_effect').reindex(index=order, columns=targets)
q = regional.pivot(index='roi', columns='target_network', values='q_BH_full_family').reindex(index=order, columns=targets)
assert registry.Label.tolist() == order and np.isfinite(eff.to_numpy()).all()
lim = float(np.ceil(np.abs(eff.to_numpy()).max()*10)/10)

fig = plt.figure(figsize=(185/25.4, 235/25.4))
def panel(x, y, letter, title):
    fig.text(x, y, letter, fontsize=11, fontweight='bold')
    fig.text(x+.035, y, title, fontsize=8.8)
def signed(v):
    return f'{v:+.4f}'.replace('-', '−')

# A: all original 19 x 9 effects, with an explicit numeric-difference color scale.
panel(.025,.969,'A','G3-associated regional FC')
fig.text(.025,.946,'19 G3-altered regions → nine target systems')
ax_a = fig.add_axes([.267,.511,.290,.383])
im = ax_a.imshow(eff.to_numpy(), aspect='auto', interpolation='nearest', cmap=DIV, vmin=-lim, vmax=lim)
ax_a.set_yticks(np.arange(19), names)
ax_a.set_xticks(np.arange(9), abbr, rotation=55, ha='left', rotation_mode='anchor')
ax_a.xaxis.tick_top()
ax_a.tick_params(length=0, pad=3)
ax_a.spines[['top','right','bottom','left']].set_visible(False)
for boundary in [8.5,9.5]: ax_a.axhline(boundary,color='white',lw=1.5)
for (i,j), v in np.ndenumerate(q.to_numpy()):
    if v < .05:
        ax_a.plot(j,i,'o',ms=2.3,markeredgewidth=0,
                  color='white' if abs(eff.iloc[i,j]) > .65*lim else '#252525')
for tick, category in zip(ax_a.get_yticklabels(),registry.category):
    tick.set_color(C['SMN'] if category == 'Somatomotor' else C['CB'] if category == 'Cerebellum' else '#222222')
cax = fig.add_axes([.267,.489,.210,.008])
cb=fig.colorbar(im,cax=cax,orientation='horizontal',ticks=[-lim,0,lim])
cb.outline.set_visible(False)
cb.ax.tick_params(labelsize=7.5,pad=2,length=2)
fig.text(.267,.455,'Adjusted (PD OFF − HC) / pooled SD')
fig.text(.267,.439,'Dots: q < 0.05 (BH, 2,547 tests)')

# B: all 39 selected G2 regions to two prespecified display targets. Each dot
# is one saved region-to-target effect, never a group-mean or extra test.
panel(.605,.969,'B','G2-associated regional FC')
fig.text(.605,.946,'39 G2-altered regions → SMN / VAN')
ax_b=fig.add_axes([.747,.743,.228,.161])
groups=[('Somatomotor','increased','SMN: G2↑'),('Somatomotor','decreased','SMN: G2↓'),
        ('VentAttn','increased','VAN: G2↑'),('VentAttn','decreased','VAN: G2↓')]
for row,(target,direction,label) in enumerate(groups):
    part=g2.loc[g2.target_network.eq(target)&g2.G2_direction.eq(direction)].sort_values('roi_id')
    jitter=.19*np.sin(np.arange(len(part))*2.399963229728653)
    sig=part.q_BH.to_numpy()<.05
    for yes in [False,True]:
        mask=sig==yes
        ax_b.scatter(part.estimate.to_numpy()[mask],row+jitter[mask],s=13,
            facecolors=C['effect'] if yes else 'white',edgecolors=C['effect'],linewidths=.7,zorder=3 if yes else 2)
ax_b.set_yticks(range(4),[g[2] for g in groups])
ax_b.set_ylim(3.5,-.5)
ax_b.set_xlim(-.18,.18)
assert g2.estimate.between(*ax_b.get_xlim()).all()
ax_b.set_xticks([-.15,0,.15],['−0.15','0','0.15'])
ax_b.axvline(0,color='#999999',ls='--',lw=.65,zorder=0)
ax_b.axhline(1.5,color='#dddddd',lw=.65,zorder=0)
ax_b.tick_params(length=2.5,pad=3)
fig.text(.792,.702,'Adjusted PD OFF − HC (Fisher z)',ha='center')
fig.text(.605,.678,'Each point = one region; filled: q < 0.05')
fig.text(.605,.662,'Open: q ≥ 0.05 (BH, 2,547 tests)')

# C: disease effects on system centroid distances; q is preserved from BH7.
panel(.605,.619,'C','System centroid distance')
fig.text(.605,.600,'Complete SMN and cerebellum')
ax_c=fig.add_axes([.723,.531,.252,.055])
for i,key in enumerate(['gradient_distance_3d','gradient_distance_10d']):
    r=geo.loc[key]
    ax_c.errorbar(r.estimate,i,xerr=[[r.estimate-r.ci95_low],[r.ci95_high-r.estimate]],
        fmt='o',ms=4.5,color=C['effect'],ecolor=C['effect'],elinewidth=1.1,capsize=2.5)
ax_c.axvline(0,color='#999999',ls='--',lw=.65)
ax_c.set_xlim(-.027,.004)
ax_c.set_ylim(1.55,-.55)
ax_c.set_yticks([0,1],['3D','10D'])
ax_c.set_xticks([-.02,-.01,0],['−0.02','−0.01','0'])
ax_c.tick_params(length=2.5,pad=2)
fig.text(.795,.492,'Adjusted PD OFF − HC (gradient units)',ha='center')
for i,(key,label) in enumerate([('gradient_distance_3d','3D'),('gradient_distance_10d','10D')]):
    r=geo.loc[key]
    fig.text(.605,.467-i*.017,f'{label}: {signed(r.estimate)} [{signed(r.ci95_low)}, {signed(r.ci95_high)}]')
fig.text(.605,.430,'95% CI; both q < 0.001 (BH, 7 tests)')

def raw_distribution(ax,metric):
    for j,group in enumerate(['HC','PD']):
        v=obs.loc[obs.group.eq(group),metric].to_numpy()
        violin=ax.violinplot([v],positions=[j],widths=.66,showextrema=False,bw_method='scott',points=160)
        for body in violin['bodies']:
            body.set_facecolor(C[group]);body.set_edgecolor(C[group]);body.set_alpha(.14);body.set_linewidth(.7)
        offset=.19*np.sin(np.arange(len(v))*2.399963229728653)
        ax.scatter(j+offset,v,s=2.9,c=C[group],alpha=.32,linewidths=0,rasterized=True,zorder=2)
        ax.boxplot([v],positions=[j],widths=.21,showfliers=False,patch_artist=True,zorder=3,
            boxprops={'facecolor':'white','edgecolor':C[group],'linewidth':.8},
            medianprops={'color':C[group],'linewidth':1.2},whiskerprops={'color':C[group],'linewidth':.7},
            capprops={'color':C[group],'linewidth':.7})
    ax.set_xlim(-.55,1.55);ax.set_xticks([0,1],['HC','PD']);ax.tick_params(length=2.5,pad=2)
def effect_annotation(metric,x):
    r=stats.loc[metric]
    assert r.q_BH<.001
    lines=[f'Δ = {signed(r.estimate)}','95% CI',f'[{signed(r.ci95_low)}, {signed(r.ci95_high)}]','q < 0.001']
    for i,line in enumerate(lines):fig.text(x,.165-i*.015,line,ha='center',fontsize=7.5)

# D–F: unchanged individual values, four saved BH6 effects with larger type.
panel(.025,.392,'D','Within-system FC')
panel(.489,.392,'E','Between-system FC')
panel(.773,.392,'F','External profiles')
fig.text(.025,.369,'Raw participant values; adjusted PD OFF − HC effects and 95% CIs below (BH, 6 tests)')
metrics=['within_sm_fc','within_cb_fc','between_fc','profile_angle']
axes=[]
for x,w,metric,title in [(.085,.135,metrics[0],'SMN'),(.295,.135,metrics[1],'Cerebellum'),
                        (.555,.160,metrics[2],'SMN ↔ cerebellum'),(.830,.145,metrics[3],'Profile similarity')]:
    ax=fig.add_axes([x,.200,w,.139]);raw_distribution(ax,metric);ax.set_title(title,pad=7)
    axes.append(ax);effect_annotation(metric,x+w/2)
for ax in axes[:2]:
    ax.set_ylim(0,1.10);ax.set_yticks([0,.25,.5,.75,1])
axes[0].set_ylabel('Mean FC (Fisher z)',labelpad=2)
axes[1].tick_params(labelleft=False)
axes[2].set_ylim(-.48,.58);axes[2].set_yticks([-.4,-.2,0,.2,.4]);axes[2].set_ylabel('Mean FC (Fisher z)',labelpad=2)
axes[2].axhline(0,color='#999999',ls='--',lw=.65,zorder=0)
axes[3].set_ylim(0,1);axes[3].set_yticks([0,.25,.5,.75,1]);axes[3].set_ylabel('Normalized-angle similarity',labelpad=2)
for ax,metric in zip(axes,metrics):
    lo,hi=ax.get_ylim();assert obs[metric].min()>lo and obs[metric].max()<hi

# Compact F construction inset: no generated brain image or synthetic profile.
diagram=fig.add_axes([.025,.024,.950,.073]);diagram.axis('off')
diagram.text(0,.97,'F · Construction of external profiles',fontsize=7.5,va='top')
def node(x,y,w,h,label,color):
    diagram.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.008,rounding_size=0.018',
        facecolor='white',edgecolor=color,linewidth=.8))
    diagram.text(x+w/2,y+h/2,label,ha='center',va='center',fontsize=7.5,color=color)
node(.010,.46,.173,.29,'SMN (30 regions)',C['SMN'])
node(.010,.01,.173,.29,'Cerebellum (31)',C['CB'])
node(.340,.175,.220,.40,'Same 222 external regions','#444444')
for y,col in [(.605,C['SMN']),(.155,C['CB'])]:
    diagram.add_patch(FancyArrowPatch((.19,y),(.33,.375),arrowstyle='-|>',mutation_scale=8,color=col,linewidth=.8))
diagram.text(.26,.59,'Mean FC',ha='center',fontsize=7.5)
diagram.add_patch(FancyArrowPatch((.574,.375),(.633,.375),arrowstyle='-|>',mutation_scale=8,color='#555555',linewidth=.8))
diagram.text(.650,.65,'Two 222-dimensional FC profiles',fontsize=7.5)
diagram.text(.650,.37,'Similarity = 1 − arccos(cosine) / π',fontsize=7.5)
diagram.text(.650,.09,'Signed FC; no centering or thresholding',fontsize=7.5)
fig.text(.025,.007,'HC, n = 141; PD, n = 334 · All participant values shown in D–F · Models adjusted for age, sex and meanFD',fontsize=7.5)

fig.canvas.draw()
renderer=fig.canvas.get_renderer()
text_bounds=[]
minfont=100
for text in fig.findobj(mpl.text.Text):
    if not text.get_visible() or not text.get_text():continue
    minfont=min(minfont,text.get_fontsize())
    bb=text.get_window_extent(renderer).transformed(fig.transFigure.inverted())
    text_bounds.append({'text':text.get_text(),'font_pt':text.get_fontsize(),'bounds':list(bb.extents)})
assert minfont>=7.5
for ext in ['pdf','svg','png']:fig.savefig(OUT/f'Figure_2_revised.{ext}',dpi=600)
fig.savefig(OUT/'Figure_2_preview.png',dpi=220)
from matplotlib.transforms import Bbox
for name,box in {'A':(0,.420,.585,.995),'BC':(.590,.420,.998,.995),'DEF':(0,0,.998,.412)}.items():
    bbox=Bbox.from_extents(box[0]*fig.get_figwidth(),box[1]*fig.get_figheight(),box[2]*fig.get_figwidth(),box[3]*fig.get_figheight())
    fig.savefig(QA/f'panel_{name}.png',dpi=250,bbox_inches=bbox)
plt.close(fig)
report={'dimensions_mm':[185,235],'minimum_text_pt':minfont,'n_HC':141,'n_PD':334,
    'models_refitted':False,'heatmap_cells':171,'G2_effect_points':78,'geometric_effects':2,
    'raw_points_per_distribution':475,'BH_families':{'A_B':2547,'C':7,'D_E_F':6},
    'centroid_CI':'pointwise HC3 t-based 95% CI, df=470; not simultaneous intervals',
    'panels':{'A':'G3 region-to-target effects / pooled raw SD','B':'G2 region-to-SMN and VAN adjusted Fisher-z effects',
        'C':'3D and 10D centroid distance, PD OFF minus HC','D':'within-system FC',
        'E':'between-system FC','F':'external-profile similarity and construction'},
    'text_bounds':text_bounds}
(QA/'render_contract.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf8')
print(json.dumps({k:v for k,v in report.items() if k!='text_bounds'},ensure_ascii=False,indent=2))
