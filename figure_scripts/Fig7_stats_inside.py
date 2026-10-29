import os, sys
sys.path.insert(0, os.getcwd())      # run from the repository root
from config import DATA_DIR, OUT_DIR
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.stats import spearmanr
IN = DATA_DIR + '/'; OUT = OUT_DIR + '/'  # folder with promoter_classes.csv; figures are saved here
d=pd.read_csv(IN+'promoter_classes.csv'); mp=d[(d.GeneType=='Multi')&d.cls.isin(['Major','Minor'])]
C={'Major':'#E64B35','Minor':'#4DBBD5'}
xlim=np.ceil(mp.meanGene.max()*1.02); ylim=np.ceil(mp.meanA.max())+0.5
fig,ax=plt.subplots(figsize=(8.2,5.2)); lines=[]
for k in ['Minor','Major']:                      # Minor first so Major sits on top
    s=mp[mp.cls==k]; rho,p=spearmanr(s.meanGene,s.meanA)
    ax.scatter(s.meanGene,s.meanA,c=C[k],s=6,alpha=0.5,edgecolors='none',label=k)
    lines.append((k,rho,p,len(s)))
ax.plot([0,ylim],[0,ylim],ls='--',color='#555555',lw=1.3)
txt='\n'.join(f"{k}: Spearman ρ = {r:.3f}, {'p < 2.2e-16' if p<2.2e-16 else f'p = {p:.1e}'}, n = {n:,}" for k,r,p,n in sorted(lines))
ax.text(0.98,0.97,txt.replace(', p','\np').replace(', n','\nn').replace('\nMinor','\n\nMinor'),transform=ax.transAxes,va='top',ha='right',multialignment='left',fontsize=9,zorder=5,bbox=dict(boxstyle='round',fc='white',ec='grey',lw=0.8,alpha=0.95))
leg=ax.legend(title='Promoter class',loc='upper left',bbox_to_anchor=(1.02,1.0),markerscale=3,frameon=True,fontsize=9)
for h in leg.legend_handles: h.set_alpha(1)
ax.set_xlim(-0.5,xlim); ax.set_ylim(-0.3,ylim); ax.grid(alpha=0.3,lw=0.5)
ax.set_xlabel('Total Gene Expression',fontsize=12,fontweight='bold'); ax.set_ylabel('Promoter Activity',fontsize=12,fontweight='bold')
fig.tight_layout()
for e in ['png','pdf']: fig.savefig(OUT+'Fig7_correlation_major_minor_oneplot_statsinside.'+e,dpi=600,bbox_inches='tight')
