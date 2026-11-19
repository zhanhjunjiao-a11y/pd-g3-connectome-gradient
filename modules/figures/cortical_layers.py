"""Re-render Figure 1 cortex from the atlas and local Conte69 meshes.

Projection, mesh lighting and camera code preserve the archived fig1_restore.py.
No network download, registration, statistical smoothing or interpolation of ROI values.
"""
from pathlib import Path
import gc, json, os
import numpy as np
import pandas as pd
import nibabel as nib
from scipy.stats import mode
from nilearn import surface
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from PIL import Image

def render(data, out):
    out.mkdir(parents=True,exist_ok=True)
    atlas=nib.load(data/'Atlas283.nii')
    table=pd.read_csv(data/'Figure_1_source_data.csv').sort_values(['gradient','roi_id'])
    meshes={}
    for hemi,short in [('left','lh'),('right','rh')]:
        path=data/'meshes'/f'conte69_32k_{short}.gii'
        surf=surface.load_surf_mesh(path)
        xyz,faces=surf.coordinates,surf.faces
        projected=surface.vol_to_surf(atlas,(xyz,faces),interpolation='nearest_most_frequent',radius=2.0)
        parcel=np.nan_to_num(projected,nan=0).astype(int)
        parcel[parcel>210]=0
        faceparcel=mode(parcel[faces],axis=1,keepdims=False).mode
        tri=xyz[faces]
        normals=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
        normals/=np.maximum(np.linalg.norm(normals,axis=1,keepdims=True),1e-12)
        vnorm=np.zeros_like(xyz,dtype=float)
        for corner in range(3): np.add.at(vnorm,faces[:,corner],normals)
        vnorm/=np.maximum(np.linalg.norm(vnorm,axis=1,keepdims=True),1e-12)
        fnorm=vnorm[faces].mean(1)
        fnorm/=np.maximum(np.linalg.norm(fnorm,axis=1,keepdims=True),1e-12)
        meshes[hemi]=dict(xyz=xyz,faces=faces,faceparcel=faceparcel,normals=fnorm)
    cmap=LinearSegmentedColormap.from_list('opaque_blue_orange',[
        (0.,'#173864'),(.18,'#3B68A8'),(.34,'#6099C4'),(.49,'#EEF3F5'),
        (.5,'#FFFFFF'),(.51,'#FCF2E7'),(.66,'#DC8B50'),(.82,'#C24A2F'),(1.,'#701B11')],N=256)
    reports=[]
    for section in ['HC','PD_HC']:
        for g in ['G1','G2','G3']:
            d=table[table.gradient.eq(g)].sort_values('roi_id')
            if section=='HC':
                values=d.healthy_mean_gradient.to_numpy()
                limit=float(np.ceil(np.abs(values).max()*100)/100)
            else:
                values=d.t_group.to_numpy().copy()
                values[d.q_BH_within_gradient.to_numpy()>=.05]=np.nan
                limit=8.
            norm=Normalize(-limit,limit)
            fig=plt.figure(figsize=(4.7,1.16))
            for n,(hemi,view) in enumerate([('left','lateral'),('left','medial'),('right','medial'),('right','lateral')]):
                m=meshes[hemi];xyz=m['xyz'];faces=m['faces']
                val=np.r_[np.nan,values][m['faceparcel']]
                colors=cmap(norm(np.nan_to_num(val)))
                colors[np.isnan(val),:3]=.90; colors[:,3]=1
                cam=-1 if ((hemi=='left')==(view=='lateral')) else 1
                light=np.array([cam*.65,-.4,.8]);light/=np.linalg.norm(light)
                bright=.68+.32*np.maximum(m['normals']@light,0)
                colors[:,:3]*=bright[:,None]
                ax=fig.add_axes([n*.25+.004,0,.242,1],projection='3d')
                ax.add_collection3d(Poly3DCollection(xyz[faces],facecolors=colors,
                    edgecolors='none',linewidths=0,antialiased=False,zsort='average'))
                span=np.ptp(xyz,axis=0);mid=(xyz.min(0)+xyz.max(0))/2
                ax.set_xlim(mid[0]-span[0]/2,mid[0]+span[0]/2)
                ax.set_ylim(mid[1]-span[1]/2,mid[1]+span[1]/2)
                ax.set_zlim(mid[2]-span[2]/2,mid[2]+span[2]/2)
                ax.set_box_aspect(span,zoom=1.28)
                ax.view_init(elev=0,azim=0 if cam==1 else 180)
                ax.set_proj_type('ortho');ax.set_axis_off()
            name=f'fig1_cortex_{section}_{g}.png'
            fig.savefig(out/name,dpi=600,facecolor='white')
            plt.close(fig);gc.collect()
            a=np.asarray(Image.open(data/name).convert('RGB')).astype(float)
            b=np.asarray(Image.open(out/name).convert('RGB')).astype(float)
            reports.append({'layer':name,'same_shape':a.shape==b.shape,
                'mean_abs_pixel_difference':float(np.abs(a-b).mean()) if a.shape==b.shape else None,
                'exact_pixels':bool(np.array_equal(a,b)),
                'statistical_source':'Figure_1_source_data.csv','cortical_regions':210})
            print('Rendered',name,flush=True)
    (out/'render_verification.json').write_text(json.dumps(reports,indent=2),encoding='utf-8')
    return reports

if __name__=='__main__':
    render(Path(os.environ['GRADIENT_FIG_INPUT'])/'figure1',Path(os.environ['GRADIENT_FIG_CORTEX']))
