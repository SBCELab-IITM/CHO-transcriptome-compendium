"""Per category, two plots: (1) activity of ALL promoters of the gene in one grouped box plot
(CHO-K1 vs CHO-DG44 side by side for each promoter; key promoter label in green), (2) gene expression."""
import os, sys
sys.path.insert(0, os.getcwd())      # run from the repository root
from config import DATA_DIR, OUT_DIR
import pickle, pandas as pd, numpy as np, matplotlib.pyplot as plt
from matplotlib.patches import Patch
import fixax
d = pickle.load(open(os.path.join(DATA_DIR, "d.pkl"), 'rb'))
meta = pd.read_csv(os.path.join(DATA_DIR, "64d11fcf-Metadata_Finalversion.csv"), sep=None, engine='python')
meta = meta[meta.Cell_line == 'CHO']; meta['L'] = meta.Cell.replace({'K1 (GS)': 'K1'})
S = {c: [r for r in meta.Run[meta.L == c] if r + '.lognorm' in d.columns] for c in ['K1', 'DG44']}
C = {'K1': '#9B8BC9', 'DG44': '#4DB8A0'}; GREEN = '#00A651'
CASES = [('1_concordant', 'Bgn_1', 2015, 'Concordant (promoter up, gene up)'),
         ('2_discordant', 'Lamtor4_1', 27972, 'Discordant (promoter down, gene up)'),
         ('3_promoter_only', 'LOC113837214', 25951, 'Differential promoter, gene unchanged')]
BP = dict(patch_artist=True, flierprops=dict(marker='o', ms=1.6, mfc='none', mec='#555555', mew=0.4),
          medianprops=dict(color='#222222', lw=0.8), whiskerprops=dict(lw=0.5, color='black'),
          capprops=dict(lw=0.5, color='black'), boxprops=dict(lw=0.5, color='black'))
def style(ax):
    ax.tick_params(labelsize=5.2, length=2, width=0.5)
    for sp in ax.spines.values(): sp.set_color(fixax.FRAME); sp.set_linewidth(0.6)
def prom_plot(ax, allp, pid):
    for i, (_, r) in enumerate(allp.iterrows()):
        for c, dx in [('K1', -0.2), ('DG44', 0.2)]:
            b = ax.boxplot([r[[x + '.lognorm' for x in S[c]]].astype(float).values], positions=[i + dx], widths=0.34, **BP)
            b['boxes'][0].set_facecolor(C[c])
    ax.set_xticks(range(len(allp))); ax.set_xticklabels([f'P{r.promoterPosition}\n({r.promoterId})' for _, r in allp.iterrows()], fontsize=5.2)
    for t, (_, r) in zip(ax.get_xticklabels(), allp.iterrows()):
        if r.promoterId == pid: t.set_color(GREEN); t.set_fontweight('bold')
    ax.set_xlim(-0.6, len(allp) - 0.4); ax.set_ylabel('Absolute promoter activity', fontsize=5.8, labelpad=2)
    ax.set_title('Promoter activity', fontsize=5.8, pad=3); style(ax)
def gene_plot(ax, r):
    b = ax.boxplot([r[[x + '.gene' for x in S[c]]].astype(float).values for c in ['K1', 'DG44']], widths=0.5, **BP)
    for bx, c in zip(b['boxes'], ['K1', 'DG44']): bx.set_facecolor(C[c])
    ax.set_xticks([1, 2]); ax.set_xticklabels(['CHO-K1', 'CHO-DG44'], fontsize=5.2)
    ax.set_ylabel('Gene expression', fontsize=5.8, labelpad=2); ax.set_title('Gene expression', fontsize=5.8, pad=3); style(ax)
H_BOX = 1.35; L, G, R, B, T = 0.45, 0.45, 0.10, 0.36, 0.62
for tag, g, pid, cat in CASES:
    allp = d[d.geneId == g].sort_values('promoterPosition'); gname = g.split('_')[0]; r0 = allp.iloc[0]
    wp = 0.55 * len(allp) + 0.25; wg = 1.0
    W = L + wp + G + wg + R; H = B + H_BOX + T
    fig = plt.figure(figsize=(W, H)); fig.text(L / W, 1 - 0.10 / H, f'{gname}: {cat}', fontsize=6.5, fontweight='bold', va='top')
    a1 = fig.add_axes([L / W, B / H, wp / W, H_BOX / H]); prom_plot(a1, allp, pid)
    a2 = fig.add_axes([(L + wp + G) / W, B / H, wg / W, H_BOX / H]); gene_plot(a2, r0)
    fig.legend(handles=[Patch(fc=C['K1'], ec='#444444', lw=0.5, label=f'CHO-K1 (n = {len(S["K1"])})'),
                        Patch(fc=C['DG44'], ec='#444444', lw=0.5, label=f'CHO-DG44 (n = {len(S["DG44"])})')],
               loc='upper left', bbox_to_anchor=(L / W, 1 - 0.24 / H), ncol=2, frameon=False, fontsize=5.2, handlelength=1, handleheight=0.8, columnspacing=1.2)
    fixax.save(fig, f'Cat{tag}_{gname}_two_plots'); plt.close(fig)
    for k, fn, w in [('promoters', lambda ax: prom_plot(ax, allp, pid), wp), ('gene', lambda ax: gene_plot(ax, r0), wg)]:
        fig, ax = fixax.fig_box(w, H_BOX, left=L, right=R, bottom=B, top=0.22); fn(ax); fixax.save(fig, f'Cat{tag}_{gname}_{k}'); plt.close(fig)
    print('done', gname)
