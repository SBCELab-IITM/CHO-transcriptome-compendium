"""Slide 1 panel d: per-sample median log2(TMM-CPM + 1), grouped by cell line / tissue.
Plot box = 2.70 x 1.81 in, the same size as the plot box of panel c on slide 1."""
import os, sys
sys.path.insert(0, os.getcwd())      # run from the repository root
from config import DATA_DIR, OUT_DIR
import pandas as pd, numpy as np, matplotlib
from scipy.stats import mannwhitneyu
from fixax import fig_box, save
t = pd.read_csv(os.path.join(OUT_DIR, 'median_expression_per_sample.csv'))
ORDER = ['CHO-K1', 'CHO-S', 'CHO-DG44', 'CHO-DXB11', 'Tissue']
COL = {'CHO-K1': '#0072B2', 'CHO-S': '#009E73', 'CHO-DG44': '#E69F00', 'CHO-DXB11': '#CC79A7', 'Tissue': '#D55E00'}
y = 'median_log2_TMM_CPM'
p = mannwhitneyu(t.loc[t.group == 'Tissue', y], t.loc[t.group != 'Tissue', y]).pvalue
fig, ax = fig_box(2.70, 1.81, left=0.55, right=0.08, bottom=0.50, top=0.08)
data = [t.loc[t.group == g, y].values for g in ORDER]
bp = ax.boxplot(data, widths=0.6, patch_artist=True, showfliers=False,
                medianprops=dict(color='black', lw=1), whiskerprops=dict(lw=0.7), capprops=dict(lw=0.7))
rng = np.random.default_rng(1)
for i, (g, b) in enumerate(zip(ORDER, bp['boxes'])):
    b.set_facecolor(matplotlib.colors.to_rgba(COL[g], 0.35)); b.set_edgecolor('black'); b.set_linewidth(0.7)   # transparent fill, solid black outline
    ax.scatter(i + 1 + rng.uniform(-0.18, 0.18, len(data[i])), data[i], s=5, c=COL[g], alpha=0.45,
               edgecolors='black', linewidths=0.15, zorder=3)
ax.set_xticks(range(1, 6)); ax.set_xticklabels(ORDER, rotation=40, ha='right', fontsize=6.5)
ax.set_ylabel('Median expression\n(log$_2$(CPM + 1))', fontsize=7, labelpad=6)
ax.tick_params(axis='y', labelsize=6.5); ax.tick_params(length=2.5, width=0.6)
ax.grid(axis='y', alpha=0.3, lw=0.4); ax.set_axisbelow(True)
for sp in ax.spines.values(): sp.set_color('#4D4D4D'); sp.set_linewidth(0.6)   # dark-grey outer frame
ax.tick_params(color='#4D4D4D')
save(fig, 'Fig1f_median_log2CPM')
print('Tissue vs CHO Mann-Whitney p =', p, t.groupby('group')[y].median().round(2).to_dict())
