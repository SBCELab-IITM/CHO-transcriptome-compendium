"""Panel d: per-sample median expression, grouped by CHO cell line / tissue.
Inputs (same folder): Count_raw_all.csv (featureCounts raw counts, genes x samples),
                      Metadata_Finalversion.csv.
Outputs: PNG (600 dpi), PDF, per-sample table, statistics."""
import os, sys
sys.path.insert(0, os.getcwd())      # run from the repository root
from config import DATA_DIR, OUT_DIR
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
IN = DATA_DIR + '/'; OUT = OUT_DIR + '/'
cnt = pd.read_csv(IN + 'Count_raw_all.csv', index_col=0)
meta = pd.read_csv(IN + 'Metadata_Finalversion.csv', sep=None, engine='python')
meta = meta[meta.Run.isin(cnt.columns)].set_index('Run')
cnt = cnt[meta.index]                                   # 350 samples in metadata
cnt = cnt[cnt.sum(axis=1) > 0]                          # drop genes with zero counts in all samples

# ---- TMM normalisation (edgeR algorithm) -> log2 CPM
def tmm_factors(X, logratio_trim=0.3, sum_trim=0.05, a_cutoff=-1e10):
    lib = X.sum(0)
    uq = X.quantile(0.75) / lib
    ref = X.columns[np.argmin(np.abs(uq - uq.mean()))]
    f = {}
    for s in X.columns:
        obs, refc = X[s].values.astype(float), X[ref].values.astype(float)
        nO, nR = lib[s], lib[ref]
        ok = (obs > 0) & (refc > 0)
        o, r = obs[ok] / nO, refc[ok] / nR
        M = np.log2(o / r); A = 0.5 * np.log2(o * r)
        v = (nO - obs[ok]) / nO / obs[ok] + (nR - refc[ok]) / nR / refc[ok]
        keep = A > a_cutoff; M, A, v = M[keep], A[keep], v[keep]
        n = len(M)
        loM, hiM = np.floor(n * logratio_trim) + 1, n + 1 - np.floor(n * logratio_trim)
        loA, hiA = np.floor(n * sum_trim) + 1, n + 1 - np.floor(n * sum_trim)
        rM = pd.Series(M).rank().values; rA = pd.Series(A).rank().values
        k = (rM >= loM) & (rM <= hiM) & (rA >= loA) & (rA <= hiA)
        f[s] = 2 ** (np.nansum(M[k] / v[k]) / np.nansum(1 / v[k]))
    f = pd.Series(f); return f / np.exp(np.mean(np.log(f)))
nf = tmm_factors(cnt)
eff = cnt.sum(0) * nf
logcpm = np.log2(cnt / eff * 1e6 + 1)
log2c = np.log2(cnt + 1)

# ---- groups
grp = meta.apply(lambda r: 'Tissue' if r.Cell_line != 'CHO' else 'CHO-' + r.Cell.replace('K1 (GS)', 'K1'), axis=1)  # K1 (GS) merged into K1
ORDER = ['CHO-K1', 'CHO-S', 'CHO-DG44', 'CHO-DXB11', 'Tissue']
STY = {'CHO-K1': ('#0072B2', '^'), 'CHO-K1 (GS)': ('#56B4E9', 'o'), 'CHO-S': ('#009E73', 's'),
       'CHO-DG44': ('#E69F00', 'v'), 'CHO-DXB11': ('#CC79A7', 'D'), 'Tissue': ('#D55E00', '*')}  # Okabe-Ito colour-blind safe
tab = pd.DataFrame({'group': grp,
                    'median_log2_counts': log2c.median(),
                    'median_log2_TMM_CPM': logcpm.median(),
                    'IQR_log2_TMM_CPM': logcpm.quantile(.75) - logcpm.quantile(.25),
                    'genes_detected': (cnt > 0).sum(),
                    'library_size_M': cnt.sum() / 1e6,
                    'TMM_factor': nf})
tab.to_csv(OUT + 'median_expression_per_sample.csv')

# ---- statistics: tissue vs all CHO
rows = []
for col in ['median_log2_counts', 'median_log2_TMM_CPM', 'IQR_log2_TMM_CPM', 'genes_detected']:
    t = tab.loc[tab.group == 'Tissue', col]; c = tab.loc[tab.group != 'Tissue', col]
    U, p = mannwhitneyu(t, c, alternative='two-sided')
    rows.append(dict(measure=col, tissue_median=t.median(), CHO_median=c.median(), U=U, p=p))
st = pd.DataFrame(rows); st.to_csv(OUT + 'statistics_tissue_vs_CHO.csv', index=False); print(st.round(4))
print(tab.groupby('group')[['median_log2_counts', 'median_log2_TMM_CPM', 'IQR_log2_TMM_CPM', 'genes_detected']].median().round(2))

# ---- plot
def panel(ax, col, ylab, p):
    data = [tab.loc[tab.group == g, col].values for g in ORDER]
    bp = ax.boxplot(data, widths=0.6, patch_artist=True, showfliers=False,
                    medianprops=dict(color='black', lw=1.2), whiskerprops=dict(lw=0.9), capprops=dict(lw=0.9))
    rng = np.random.default_rng(1)
    for i, (g, b) in enumerate(zip(ORDER, bp['boxes'])):
        c, mk = STY[g]; b.set_facecolor(c); b.set_alpha(0.35); b.set_edgecolor('black'); b.set_linewidth(0.9)
        y = data[i]; x = i + 1 + rng.uniform(-0.18, 0.18, len(y))
        ax.scatter(x, y, s=10, marker='o', c=c, alpha=0.45, edgecolors='black', linewidths=0.2, zorder=3)
    ax.set_xticks(range(1, len(ORDER) + 1)); ax.set_xticklabels(ORDER, rotation=40, ha='right', fontsize=6.5)
    ax.set_ylabel(ylab, fontsize=7.5, labelpad=10); ax.tick_params(axis='y', labelsize=6.5)
    ax.grid(axis='y', alpha=0.3, lw=0.5); ax.set_axisbelow(True)
    if False: ax.text(0.02, 0.97, f'Mann–Whitney p = {p:.1e}' if p >= 2.2e-16 else 'Mann–Whitney p < 2.2e-16',
            transform=ax.transAxes, va='top', fontsize=6)
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4))
panel(axes[0], 'median_log2_counts', 'Median expression\n(log2 counts + 1)', st.p[0])
panel(axes[1], 'median_log2_TMM_CPM', 'Median expression\n(log2 TMM-CPM + 1)', st.p[1])
fig.tight_layout()
for e in ['png', 'pdf']: fig.savefig(OUT + 'Fig_median_expression_boxplot.' + e, dpi=600, bbox_inches='tight')

# ---- single panel (TMM-normalised) for slide 1 panel d
fig, ax = plt.subplots(figsize=(3.6, 3.4))
panel(ax, 'median_log2_TMM_CPM', 'Median expression\n(log2 TMM-CPM + 1)', st.p[1])
fig.tight_layout()
for e in ['png', 'pdf']: fig.savefig(OUT + 'Fig_median_expression_TMM.' + e, dpi=600, bbox_inches='tight')
