"""Rebuild SFigures 1 and 2 using the established display code and all source values."""
from pathlib import Path
import os,json
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle

HERE=Path(__file__).resolve().parent
DATA=Path(os.environ['GRADIENT_FIG_INPUT'])
OUT=Path(os.environ['GRADIENT_FIG_OUTPUT'])/'figures'
OUT.mkdir(parents=True,exist_ok=True)
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
 'font.size':8,'axes.labelsize':8,'axes.titlesize':9,'xtick.labelsize':7,'ytick.labelsize':7,
 'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,
 'axes.linewidth':.6,'xtick.major.width':.6,'ytick.major.width':.6,'legend.frameon':False,
 'figure.facecolor':'white','savefig.facecolor':'white'})
C={'HC':'#3B68A8','PD':'#C24A2F','SM':'#478D89','CB':'#8B6BA8'}
DIV=LinearSegmentedColormap.from_list('direction',['#27558A','#A7C2DB','#FFFFFF','#E7B493','#AF472D'])
reg=pd.read_csv(DATA/'supplement1/selected19_roi_registry.csv')
edges=pd.read_csv(DATA/'supplement1/selected19_to_all_primary.csv')
net=pd.read_csv(DATA/'supplement1/yeo7_283_assignment.csv')
order=reg.Label.tolist()
assert len(order)==19 and len(edges)==19*282 and edges.family_n.eq(39903).all()
names=['L precentral (upper limb)','R precentral (upper limb)','L precentral (trunk)','R precentral (trunk)',
 'L paracentral','R paracentral','L postcentral (trunk)','R postcentral (trunk)','R insula','L thalamus',
 'L cerebellum VI','Vermis VI','R cerebellum VI','L Crus I','R Crus I','R Crus II','L VIIIb','Vermis VIIIb','Vermis IX']
cats=['Visual','Somatomotor','DorsAttn','VentAttn','Limbic','Frontoparietal','Default','Subcortex','Cerebellum']
abbr=['VIS','SMN','DAN','VAN','LIM','FPN','DMN','SUB','CB']
targets=[];blocks=[]
for cat in cats:
    vals=net.loc[net.category.eq(cat),'Label'].tolist()
    blocks.append((len(targets),len(vals)));targets+=vals
assert len(set(targets))==283
eff=edges.pivot(index='seed_roi',columns='target_roi',values='standardized_effect').reindex(index=order,columns=targets)
q=edges.pivot(index='seed_roi',columns='target_roi',values='q_BH_full_family').reindex(index=order,columns=targets)
shown=eff.where(q.lt(.05));limit=float(np.ceil(np.nanmax(np.abs(eff.to_numpy()))*10)/10)
exec(compile((HERE/'supplement1_plot.py').read_text(encoding='utf-8'),'supplement1_plot.py','exec'))
def save(fig,name):
    for ext in ['pdf','svg','png']:fig.savefig(OUT/f'{name}.{ext}',dpi=600)
    plt.close(fig)
save(fig,'Supplementary_Figure_1')
exec(compile((HERE/'supplement2_plot.py').read_text(encoding='utf-8'),'supplement2_plot.py','exec'))
print(json.dumps({'supplement1_rows':len(edges),'supplement2_rows':len(d2),'figures':2}))
