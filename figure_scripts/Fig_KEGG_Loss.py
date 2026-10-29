"""KEGG Loss dot plot (slide 2), rebuilt from CHO_KEGG_enrichment_with_genes_loss.csv at 600 dpi with wrapped labels."""
import os, sys
sys.path.insert(0, os.getcwd())      # run from the repository root
from config import DATA_DIR, OUT_DIR
import pandas as pd, numpy as np, textwrap, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
WRAP = 24
CMAP = LinearSegmentedColormap.from_list('padj', ['#3B7DC4', '#8E6FA6', '#E06666'])  # low p.adj blue -> high red (as original)
plt.rcParams.update({'font.family': 'DejaVu Sans'})

d = pd.read_csv(os.path.join(DATA_DIR, 'KEGG_Loss.csv'))
d['GR'] = d.GeneRatio.map(lambda s: int(s.split('/')[0]) / int(s.split('/')[1]))
d = d.sort_values('p.adjust', kind='stable').iloc[::-1].reset_index(drop=True)   # most significant plotted at the top
size = lambda c: 4 + 6 * np.asarray(c, float)                             # marker area ~ Count

fig, ax = plt.subplots(figsize=(3.36, 4.32))
sc = ax.scatter(d.GR, d.index, s=size(d.Count), c=d['p.adjust'], cmap=CMAP, vmin=0, vmax=0.05,
                edgecolors='black', linewidths=0.5, zorder=3)
ax.set_yticks(d.index); ax.set_yticklabels([textwrap.fill(t, WRAP) for t in d.Description], fontsize=4.6, linespacing=0.95)
ax.set_ylim(-0.6, len(d) - 0.4); ax.set_xlim(0.008, 0.113); ax.set_xticks([0.025, 0.05, 0.075, 0.1])
ax.tick_params(axis='both', labelsize=4.6, length=2, width=0.5)
ax.set_xlabel('GeneRatio', fontsize=5.5, fontweight='bold', labelpad=2)
ax.grid(axis='x', ls='--', lw=0.4, color='#BDBDBD'); ax.grid(axis='y', lw=0.3, color='#E6E6E6'); ax.set_axisbelow(True)
for s in ax.spines.values(): s.set_linewidth(0.6)

cax = fig.add_axes([0.835, 0.50, 0.035, 0.22])
cb = fig.colorbar(sc, cax=cax, ticks=[0, 0.01, 0.02, 0.03, 0.04]); cb.outline.set_linewidth(0)
cb.ax.tick_params(labelsize=4.4, length=2, width=0.4); cax.set_title('p.adj', fontsize=5, fontweight='bold', loc='left', pad=3)
h = [Line2D([], [], ls='', marker='o', mfc='black', mec='black', ms=np.sqrt(size(c)), label=str(c)) for c in (5, 10, 15)]
fig.legend(handles=h, title='Count', loc='upper left', bbox_to_anchor=(0.815, 0.45), frameon=False, fontsize=4.4,
           title_fontproperties={'size': 5, 'weight': 'bold'}, labelspacing=1.3, borderpad=0.2, handletextpad=0.4, alignment='left')
# plot box identical to the KEGG Reduced panel: 1.367 in wide x 3.66 in tall (frame 0.396-0.803 of 3.36 in, 0.077-0.983 of 4.32 in)
BOXW, H, LEG = 0.407 * 3.36, 4.32, 0.197 * 3.36 + 0.02
fig.canvas.draw(); r = fig.canvas.get_renderer()
labw = (ax.get_window_extent(r).x0 - min(t.get_window_extent(r).x0 for t in ax.get_yticklabels())) / fig.dpi
W = labw + 0.05 + BOXW + LEG
fig.set_size_inches(W, H)
ax.set_position([(labw + 0.05) / W, 0.077, BOXW / W, 0.906])
x1 = (labw + 0.05 + BOXW) / W
cax.set_position([x1 + 0.035 * 3.36 / W, 0.50, 0.035 * 3.36 / W, 0.22])
fig.legends[0].set_bbox_to_anchor((x1 + 0.012 * 3.36 / W, 0.45))
for e in ['png', 'pdf']: fig.savefig(os.path.join(OUT_DIR, f'Fig_KEGG_Loss.{e}'), dpi=600)
print(d[['Description', 'GR', 'p.adjust', 'Count']].iloc[::-1].to_string())
