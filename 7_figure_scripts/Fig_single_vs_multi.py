"""Major promoter activity vs total gene expression: single- vs multi-promoter genes.
Input: promoter_classes.csv (same folder). Output: PNG (600 dpi) + PDF."""
import os, sys
sys.path.insert(0, os.getcwd())      # run from the repository root
from config import DATA_DIR, OUT_DIR
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.stats import spearmanr
IN = DATA_DIR + '/'; OUT = OUT_DIR + '/'
d = pd.read_csv(IN + 'promoter_classes.csv')
maj = d[d.cls == 'Major']                       # one (most active) promoter per expressed gene
C = {'Single': '#3C5488', 'Multi': '#F39B7F'}
LAB = {'Single': 'Single-promoter genes', 'Multi': 'Multi-promoter genes'}
xlim = np.ceil(maj.meanGene.max() * 1.02); ylim = np.ceil(maj.meanA.max()) + 0.5
fig, ax = plt.subplots(figsize=(8.2, 5.2)); stats = []
for k in ['Single', 'Multi']:                   # Single first, Multi on top
    s = maj[maj.GeneType == k]; rho, p = spearmanr(s.meanGene, s.meanA)
    ax.scatter(s.meanGene, s.meanA, c=C[k], s=6, alpha=0.5, edgecolors='none', label=LAB[k], zorder=2)
    stats.append((k, rho, p, len(s)))
ax.plot([0, ylim], [0, ylim], ls='--', color='#555555', lw=1.3, zorder=3)
txt = '\n\n'.join(f"{LAB[k].split('-')[0]}-promoter: Spearman ρ = {r:.3f}\n"
                  f"{'p < 2.2e-16' if p < 2.2e-16 else f'p = {p:.1e}'}\nn = {n:,} genes" for k, r, p, n in stats)
ax.text(0.98, 0.97, txt, transform=ax.transAxes, va='top', ha='right', multialignment='left', fontsize=9, zorder=5,
        bbox=dict(boxstyle='round', fc='white', ec='grey', lw=0.8, alpha=0.95))
leg = ax.legend(title='Gene type', loc='upper left', bbox_to_anchor=(1.02, 1.0), markerscale=3, frameon=True, fontsize=9)
for h in leg.legend_handles: h.set_alpha(1)
ax.set_xlim(-0.5, xlim); ax.set_ylim(-0.3, ylim); ax.grid(alpha=0.3, lw=0.5)
ax.set_xlabel('Total Gene Expression', fontsize=12, fontweight='bold')
ax.set_ylabel('Major Promoter Activity', fontsize=12, fontweight='bold')
fig.tight_layout()
for e in ['png', 'pdf']:
    fig.savefig(OUT + 'Fig_major_promoter_vs_gene_single_vs_multi.' + e, dpi=600, bbox_inches='tight')
print(stats)
