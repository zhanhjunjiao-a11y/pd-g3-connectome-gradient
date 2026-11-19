"""Figure 1 annotation revision using the original anatomical rendering and layout.

Reads locked regional statistics only. No refitting, no regional reselection,
no statistical-value smoothing, and no value clipping. All meshes are real.
"""
from pathlib import Path
import os,sys,json,hashlib,gc,warnings,zipfile
HERE=Path(__file__).resolve().parent
BASE=Path(os.environ['GRADIENT_FIG_OUTPUT'])
OUT=BASE/'figures'
INPUT=Path(os.environ['GRADIENT_FIG_INPUT'])/'figure1'
QA=BASE/'qa'/'figure1'
OUT.mkdir(parents=True,exist_ok=True)
QA.mkdir(parents=True,exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR',str(QA/'fig1_mplcache'))
import numpy as np
import pandas as pd
import nibabel as nib
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.collections import PolyCollection,LineCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.stats import mode
from PIL import Image
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
    'svg.fonttype':'none','pdf.fonttype':42,'font.size':7.5,'axes.linewidth':.5,
    'savefig.facecolor':'white','figure.facecolor':'white'})

# Frozen aggregate data and existing cortical render assets are packaged locally.
# The original cortical renders are already raster layers in the source SVG.
SOURCE=INPUT/'Figure_1_source_data.csv'
source_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in INPUT.iterdir() if p.is_file()}
df=pd.read_csv(SOURCE).sort_values(['gradient','roi_id'])
healthy=df.pivot(index='roi_id',columns='gradient',values='healthy_mean_gradient').reindex(index=range(1,284),columns=['G1','G2','G3']).to_numpy()
assert healthy.shape==(283,3) and np.isfinite(healthy).all()
assert len(df)==849 and set(df.n_records_PD)=={334} and set(df.n_records_HC)=={141}
counts={g:int((d.q_BH_within_gradient<.05).sum()) for g,d in df.groupby('gradient')}
assert counts=={'G1':1,'G2':39,'G3':19}
atlas=nib.load(INPUT/'Atlas283.nii')
atlasdata=atlas.get_fdata().astype(int)
assert np.array_equal(np.unique(atlasdata),np.arange(284))
inv=np.linalg.inv(atlas.affine)
SRES=INPUT

# The stronger diverging palette is fully opaque; its color bar shows the exact
# linear value scale. Intermediate control points strengthen low/mid contrast.
CMAP=LinearSegmentedColormap.from_list('opaque_blue_orange',[
    (0.,'#173864'),(.18,'#3B68A8'),(.34,'#6099C4'),(.49,'#EEF3F5'),
    (.5,'#FFFFFF'),(.51,'#FCF2E7'),(.66,'#DC8B50'),(.82,'#C24A2F'),(1.,'#701B11')],N=256)
CMAP.set_bad('#D6D6D6')

# SUIT flatmap: sample MNI/FSL pial–white coordinates to match the standard-MNI
# volume grid, not SUIT-template coordinates. The source atlas has a standard
# 2-mm MNI grid but no fine-grained template-family metadata; this matching
# assumption is explicit in QA. No spatial registration is performed here.
flat=nib.load(SRES/'FLAT.surf.gii')
fxy=flat.darrays[0].data; ffaces=flat.darrays[1].data
cp=nib.load(SRES/'PIAL_FSL.surf.gii').darrays[0].data
cw=nib.load(SRES/'WHITE_FSL.surf.gii').darrays[0].data
sample=[]
for depth in [0,.2,.4,.6,.8,1]:
    coord=cp*(1-depth)+cw*depth
    vox=np.rint(nib.affines.apply_affine(inv,coord)).astype(int)
    valid=np.all((vox>=0)&(vox<np.array(atlasdata.shape)),axis=1)
    v=np.zeros(len(vox),int); v[valid]=atlasdata[tuple(vox[valid].T)]
    # The four deep nuclei (280..283) are not cerebellar surface parcels.
    v[(v<253)|(v>279)]=0
    sample.append(v)
samples=np.array(sample)
# Ignore empty samples when finding the modal label, as old SUIT rendering
# ignored out-of-atlas/background values. No nearest-neighbour filling of gaps.
cbparcel=np.zeros(samples.shape[1],int)
for i in range(samples.shape[1]):
    nonzero=samples[:,i][samples[:,i]>0]
    if len(nonzero): cbparcel[i]=np.bincount(nonzero).argmax()
cbface=mode(cbparcel[ffaces],axis=1,keepdims=False).mode
borders=np.loadtxt(SRES/'borders.txt',delimiter=',')
cbunique=np.unique(cbparcel[cbparcel>0])
print('Cerebellar parcels projected:',len(cbunique),flush=True)

def render_cortex_row(values,lim,stem):
    # Reuse the exact six existing raster cortical layers. No anatomy is redrawn.
    out=Path(os.environ['GRADIENT_FIG_CORTEX'])/f'{stem}.png'
    assert out.exists()
    return out

canonical=nib.as_closest_canonical(atlas)
cdata=canonical.get_fdata().astype(int)

def draw_sub(ax,values,lim):
    # Same plane as the old figure, confined to actual subcortical labels.
    zi=round((8-canonical.affine[2,3])/canonical.affine[2,2])
    lab=cdata[:,:,zi]; roi=(lab>=211)&(lab<=252)
    x,y=np.where(roi)
    lab=lab[x.min()-3:x.max()+4,y.min()-3:y.max()+4]
    valid=(lab>=211)&(lab<=252)
    base=np.where(valid,.82,np.nan)
    ax.imshow(base.T,origin='lower',cmap='gray',vmin=0,vmax=1,interpolation='nearest')
    data=np.r_[np.nan,values][lab]; data[~valid]=np.nan
    transparent=CMAP.copy(); transparent.set_bad((0,0,0,0))
    ax.imshow(data.T,origin='lower',cmap=transparent,vmin=lim[0],vmax=lim[1],interpolation='nearest')
    ax.set_axis_off()

def draw_cb(ax,values,lim):
    val=np.r_[np.nan,values][cbface]
    color=CMAP(Normalize(*lim)(np.nan_to_num(val)))
    color[np.isnan(val),:3]=.84; color[:,3]=1
    coll=PolyCollection(fxy[ffaces,:2],facecolors=color,edgecolors='none',linewidths=0,antialiased=False)
    coll.set_rasterized(True); ax.add_collection(coll)
    # Original SUIT anatomical fissures: real-coordinate separators, no invented borders.
    ax.plot(borders[:,0],borders[:,1],'.',markersize=.34,color='#242424',rasterized=True)
    ax.set_xlim(-96,97); ax.set_ylim(-85,87); ax.set_aspect('equal'); ax.set_axis_off()

def draw_deep_cb(ax,values,lim):
    yi=round((-58-canonical.affine[1,3])/canonical.affine[1,1])
    lab=cdata[:,yi,:]; roi=(lab>=280)&(lab<=283)
    assert set(np.unique(lab[roi]))==set(range(280,284))
    x,z=np.where(roi)
    lab=lab[x.min()-2:x.max()+3,z.min()-2:z.max()+3]
    valid=(lab>=280)&(lab<=283)
    ax.imshow(np.where(valid,.82,np.nan).T,origin='lower',cmap='gray',vmin=0,vmax=1,interpolation='nearest')
    data=np.r_[np.nan,values][lab]; data[~valid]=np.nan
    transparent=CMAP.copy(); transparent.set_bad((0,0,0,0))
    ax.imshow(data.T,origin='lower',cmap=transparent,vmin=lim[0],vmax=lim[1],interpolation='nearest')
    ax.set_axis_off()


# Contract: the current regional findings in the previous six-row anatomical
# layout. The statistical models and source values are unchanged. Brain-surface
# lighting is generated from real mesh normals; statistical values are never
# smoothed, stretched, clipped, or filtered to improve appearance.
W,H=185.,218.
fig=plt.figure(figsize=(W/25.4,H/25.4))

def axmm(x,y,w,h):
    return fig.add_axes([x/W,1-(y+h)/H,w/W,h/H])

def txt(x,y,t,**kw):
    return fig.text(x/W,1-y/H,t,va='top',fontsize=kw.pop('fontsize',7),**kw)

def header(letter,y,title,subtitle):
    txt(2,y,letter,weight='bold',fontsize=9)
    txt(13,y,title,weight='bold',fontsize=9)
    txt(13,y+6,subtitle,fontsize=7)

def columns(y):
    for x,label in [(27,'L lateral'),(55,'L medial'),(83,'R medial'),(111,'R lateral'),(138,'Subcortex'),(161,'Cerebellum')]:
        txt(x,y,label,ha='center',fontsize=7)

header('A',3,'Healthy-control connectome gradients',
       'Mean aligned individual gradients; 141 healthy controls')
columns(17)
header('B',110,'Regional gradient differences in Parkinson’s disease',
       '334 patients / 141 controls; age, sex and mean framewise displacement adjusted')
columns(124)

def row(y,g,values,lim,section):
    png=render_cortex_row(values,lim,f'fig1_cortex_{section}_G{g}')
    ax=axmm(13,y,114,28);ax.imshow(plt.imread(png));ax.set_axis_off()
    ax=axmm(129,y+5,18,19);draw_sub(ax,values,lim)
    ax=axmm(148,y+1,25,21);draw_cb(ax,values,lim)
    ax=axmm(161,y+22,11,4);draw_deep_cb(ax,values,lim)
    txt(148,y+23,'Nuclei',fontsize=7.5)
    txt(2,y+10,f'G{g}',weight='bold',fontsize=8)
    if section=='PD_HC':
        d=df[df.gradient.eq(f'G{g}')]
        n=int((d.q_BH_within_gradient<.05).sum())
        txt(2,y+16,f'{n} ROI'+('s' if n!=1 else ''),fontsize=7.5)
    cax=axmm(174,y+6,1.5,16)
    cb=fig.colorbar(ScalarMappable(norm=Normalize(*lim),cmap=CMAP),cax=cax)
    cb.set_ticks([lim[0],0,lim[1]])
    cb.ax.tick_params(labelsize=7.5,length=2,pad=1)
    cb.outline.set_linewidth(.4)

for i,y in enumerate([22,50,78]):
    print('Rendering HC',i+1,flush=True)
    limit=float(np.ceil(np.abs(healthy[:,i]).max()*100)/100)
    row(y,i+1,healthy[:,i],(-limit,limit),'HC')
for i,y in enumerate([129,157,185]):
    print('Rendering disease',i+1,flush=True)
    d=df[df.gradient.eq(f'G{i+1}')].sort_values('roi_id')
    values=d.t_group.to_numpy().copy()
    values[d.q_BH_within_gradient.to_numpy()>=.05]=np.nan
    row(y,i+1,values,(-8,8),'PD_HC')

txt(174,17,'Value',fontsize=7.5)
txt(183.5,120,'t value (PD OFF − HC)',fontsize=7.5,ha='right')
txt(13,214,'B: coloured regions, BH q < 0.05 within each gradient; grey, non-significant regions.',fontsize=7)
mat=df.pivot(index='gradient',columns='roi_id',values='t_group').reindex(['G1','G2','G3']).to_numpy()
assert mat.shape==(3,283) and np.isfinite(mat).all() and np.abs(mat).max()<8
# All formats now share this rendering, including editable labels and vector
# colour bars. The anatomical image layers and region values are unchanged.
fig.canvas.draw()
renderer=fig.canvas.get_renderer()
text_records=[]
for text in fig.findobj(mpl.text.Text):
    if text.get_visible() and text.get_text():
        box=text.get_window_extent(renderer)
        assert text.get_fontsize()>=7, (text.get_text(),text.get_fontsize())
        assert box.x0>=-0.5 and box.x1<=fig.bbox.width+0.5, (text.get_text(),list(box.extents))
        assert box.y0>=-0.5 and box.y1<=fig.bbox.height+0.5, (text.get_text(),list(box.extents))
        text_records.append({'text':text.get_text(),'font_pt':text.get_fontsize(),'bounds_px':list(box.extents)})
fig.savefig(OUT/'Figure_1_revised.svg',dpi=600,facecolor='white')
fig.savefig(OUT/'Figure_1_revised.pdf',dpi=600,facecolor='white')
fig.savefig(OUT/'Figure_1_revised.png',dpi=600,facecolor='white')
fig.savefig(OUT/'Figure_1_revised.tiff',dpi=600,facecolor='white',pil_kwargs={'compression':'tiff_lzw'})
plt.close(fig)
im=Image.open(OUT/'Figure_1_revised.png')
preview=im.copy();preview.thumbnail((1480,2000))
preview.save(OUT/'Figure_1_revised_preview.png')
assert {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in INPUT.iterdir() if p.is_file()}==source_hashes
qa=dict(core_conclusion='Whole-brain healthy-control organization and regional PD-HC gradient differences.',
    change='Explicit PD OFF minus HC colour-bar direction; 7.5-point tick, nuclei and ROI-count labels; 7-point footnote. Healthy-control mean title retained.',
    archetype='anatomical image plate with quantitative colour scales',
    backend='Python; original Matplotlib layout',width_mm=W,height_mm=H,dpi=600,
    source_rows=len(df),counts=counts,no_models_refit=True,no_statistics_recomputed=True,
    n_clipped_values=0,all_values_preserved=True,
    statistical_values_unchanged=True,source_asset_hashes=source_hashes,
    minimum_text_font_pt=min(r['font_pt'] for r in text_records),
    text_bounds_within_canvas=True,text_records=text_records,
    rendering='Original six cortical raster layers, original subcortical/SUIT routines, vector labels and colour bars.',
    exports={})
for ext in ['svg','pdf','png','tiff']:
    p=OUT/f'Figure_1_revised.{ext}'
    qa['exports'][ext]={'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
qa['raster_pixels']=list(im.size)
qa['raster_dpi']=[float(v) for v in im.info.get('dpi',())]
(QA/'Figure_1_annotation_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'counts':counts,'statistical_values_unchanged':True,'size_mm':[W,H],'pixels':list(im.size),'minimum_text_pt':qa['minimum_text_font_pt']}))
