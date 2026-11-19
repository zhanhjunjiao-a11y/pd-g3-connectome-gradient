d2=pd.read_csv(DATA/'supplement2/all_roi_results.csv')
fig,axs=plt.subplots(3,1,figsize=(185/25.4,180/25.4),sharex=True)
fig.subplots_adjust(left=.13,right=.97,bottom=.11,top=.92,hspace=.32)
fig.suptitle('Whole-brain regional changes after acute levodopa challenge',fontsize=10,fontweight='bold',y=.97)
for k,(ax,g) in enumerate(zip(axs,['G1','G2','G3'])):
 f=d2.loc[d2.gradient.eq(g)].sort_values('roi_id');assert len(f)==283
 ax.axhline(0,color='#999999',lw=.6);ax.axvspan(210.5,252.5,color='#F4F4F4',zorder=0)
 ax.scatter(f.roi_id,f.t_HC3,color='#AAAAAA',s=9,linewidths=0,rasterized=True)
 sig=f.loc[f.q_BH_283_within_gradient.lt(.05)]
 ax.scatter(sig.roi_id,sig.t_HC3,c=np.where(sig.t_HC3>0,'#C24A2F','#3B68A8'),s=19,linewidths=0,rasterized=True)
 original=f.loc[f.original_disease_roi.astype(str).str.lower().eq('true')]
 ax.scatter(original.roi_id,original.t_HC3,facecolors='none',edgecolors='#222222',s=34,linewidths=.6)
 ax.set_ylabel(f'{g} change\nt statistic');ax.set_xlim(.5,283.5)
 ax.text(-.105,1.05,'ABC'[k],transform=ax.transAxes,fontweight='bold',fontsize=9)
 ax.text(.99,1.025,f'{len(sig)} regions at q < 0.05',ha='right',va='bottom',transform=ax.transAxes,fontsize=7)
 for x in [210.5,252.5]:ax.axvline(x,color='#CCCCCC',lw=.6)
axs[-1].set_xticks([1,50,100,150,210,252,283]);axs[-1].set_xlabel('Region index: cortex 1–210, subcortex 211–252, cerebellum 253–283')
fig.text(.13,.035,'Coloured: q < 0.05 within gradient; outlined: original disease-difference region',fontsize=7)
save(fig,'Supplementary_Figure_2')
