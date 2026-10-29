"""IGV-style tracks for the three promoter-category genes (CriGri-PICR coordinates).
Tracks: RefSeq transcripts (GTF), promoters (key promoter boxed in green), H3K4me3 and H3K27ac (CHO-K1, mean T5-T17,
depth-normalised; read from the PICRH-aligned bigWigs at the matched position), and RNA-seq coverage for
CHO-K1 / CHO-DG44 when the bigWigs are supplied in RNA (left empty until then)."""
import os, sys
sys.path.insert(0, os.getcwd())      # run from the repository root
from config import DATA_DIR, OUT_DIR
import os, numpy as np, pandas as pd, pyBigWig, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import fixax
EX = pd.read_csv(os.path.join(DATA_DIR, 'exons_3genes.tsv'), sep='\t', header=None, names=['gene', 'tx', 'chr', 's', 'e', 'strand'])
SF = pd.read_csv(os.path.join(DATA_DIR, "Fig9_depth_scale_factors.csv"))
BW = os.path.join(DATA_DIR, "{m}_{t}_treat_pileup.bw"); TPS = ['T5', 'T9', 'T13', 'T17']
RNA = {}   # e.g. RNA['Bgn'] = {'CHO-K1': 'SRR10572657.bw', 'CHO-DG44': 'SRR10493877.bw'}
GENES = {  # name: gene_id, PICR chr, region, PICR->PICRH: picrh_pos = C - picr_pos on chrH, promoters, key promoter
    'Bgn': ('Bgn_1', 'NW_020822680.1', (196800, 210500), 'NC_048604.1', 2450434,
            {'P1 (2015)': 197795, 'P2 (2016)': 202907}, 'P1 (2015)', 'Concordant (promoter up, gene up)'),
    'Lamtor4': ('Lamtor4_1', 'NW_020822534.1', (11181000, 11186000), 'NC_048597.1', 127922188,
                {'P1': 11181861, 'P2 (27972)': 11182007, 'P3': 11182404}, 'P2 (27972)', 'Discordant (promoter down, gene up)'),
    'LOC113837214': ('LOC113837214', 'NW_020822643.1', (3826800, 3839000), 'NC_048602.1', 10414092,
                     {'P1': 3827579, 'P2': 3827819, 'P3 (25951)': 3827976}, 'P3 (25951)', 'Differential promoter, gene unchanged')}
COL = {'H3K4me3': '#3B5BA5', 'H3K27ac': '#B5525C', 'CHO-K1': '#9B8BC9', 'CHO-DG44': '#4DB8A0'}
def chip(mark, chrH, C, a, b, nb):
    vals = []
    for t in TPS:
        f = SF[(SF['mark'] == mark) & (SF.timepoint == t)].scale_factor.iloc[0]
        bw = pyBigWig.open(BW.format(m=mark, t=t))
        v = np.nan_to_num(np.array(bw.values(chrH, C - b, C - a)))[::-1] * f; vals.append(v); bw.close()
    v = np.mean(vals, 0); return v[:len(v) // nb * nb].reshape(-1, len(v) // nb).mean(1) if False else v.reshape(nb, -1).mean(1)
def rna(path, chr_, a, b, nb):
    bw = pyBigWig.open(path); v = np.nan_to_num(np.array(bw.values(chr_, a, b))); bw.close(); return v.reshape(nb, -1).mean(1)
for name, (gid, chr_, (a, b), chrH, C, proms, key, cat) in GENES.items():
    bs = max(5, (b - a) // 400); b = a + (b - a) // bs * bs; nb = (b - a) // bs; x = np.linspace(a, b, nb, endpoint=False)
    sig = [(f'{m}\n(CHO-K1)', chip(m, chrH, C, a, b, nb), COL[m]) for m in ['H3K4me3', 'H3K27ac']]
    rn = [(f'RNA-seq\n{k}', rna(p, chr_, a, b, nb), COL[k]) for k, p in RNA.get(name, {}).items() if os.path.exists(p)]
    tx = EX[EX.gene == gid]; ntx = tx.tx.nunique()
    tracks = rn + sig; H = 0.75 + 0.12 * ntx + 0.25 + 0.42 * len(tracks)
    fig = plt.figure(figsize=(4.2, H)); L, W = 0.75 / 4.2, 3.3 / 4.2; y = H - 0.30
    fig.text(L, (H - 0.08) / H, f'{name}: {cat}', fontsize=6.5, fontweight='bold', va='top')
    # transcripts
    h = 0.12 * ntx + 0.05; y -= h; ax = fig.add_axes([L, y / H, W, h / H]); ax.set_xlim(a, b); ax.set_ylim(-0.5, ntx - 0.5); ax.axis('off')
    for i, (t, e) in enumerate(tx.groupby('tx')):
        yy = ntx - 1 - i; ax.plot([e.s.min(), e.e.max()], [yy, yy], color='#1F3A93', lw=0.5)
        for _, r in e.iterrows(): ax.add_patch(Rectangle((r.s, yy - 0.28), r.e - r.s, 0.56, color='#1F3A93', lw=0))
        ax.text(a - (b - a) * 0.01, yy, t, ha='right', va='center', fontsize=4.2)
    # promoters
    y -= 0.22; ax = fig.add_axes([L, y / H, W, 0.18 / H]); ax.set_xlim(a, b); ax.set_ylim(0, 1); ax.axis('off')
    for j, (lab, pos) in enumerate(proms.items()):
        k = lab == key; ax.plot([pos, pos], [0.1, 0.9], color='#00A651' if k else '#888888', lw=1.2 if k else 0.7)
        ax.text(pos + (j - (len(proms) - 1) / 2) * (b - a) * 0.018 * (min(np.diff(sorted(proms.values()))) < (b - a) * 0.03), 1.0, lab.split(' ')[0], fontsize=4.2, ha='center', va='bottom', color='#00A651' if k else '#666666')
    fig.text(L - 0.66 / 4.2, (y + 0.09) / H, 'Promoters', fontsize=5, va='center')
    # signal tracks
    for lab, v, c in tracks:
        y -= 0.42; ax = fig.add_axes([L, (y + 0.04) / H, W, 0.34 / H]); ax.fill_between(x, v, step='post', color=c, lw=0)
        ax.set_xlim(a, b); ax.set_ylim(0, max(v.max() * 1.1, 1)); ax.set_yticks([0, round(max(v.max(), 1))]); ax.tick_params(labelsize=4.2, length=1.5, width=0.4)
        ax.set_xticks([]); [s.set_visible(False) for k, s in ax.spines.items() if k in ('top', 'right', 'bottom')]; ax.spines['left'].set_linewidth(0.4)
        fig.text(L - 0.66 / 4.2, (y + 0.21) / H, lab, fontsize=5, va='center', color=c)
    # coordinate axis + green box around the key promoter spanning all tracks
    ax0 = fig.add_axes([L, 0.30 / H, W, (H - 0.60) / H], facecolor='none'); ax0.set_xlim(a, b); ax0.set_ylim(0, 1)
    for s in ['top', 'right', 'left']: ax0.spines[s].set_visible(False)
    ax0.spines['bottom'].set_linewidth(0.5); ax0.set_yticks([]); ax0.tick_params(labelsize=4.5, length=2, width=0.4)
    ax0.set_xticks(np.linspace(a, b, 5)); ax0.set_xticklabels([f'{v/1000:,.1f} kb' for v in np.linspace(a, b, 5)])
    ax0.set_xlabel(chr_, fontsize=5, labelpad=1)
    pk = proms[key]; wbox = (b - a) * 0.025
    ax0.add_patch(Rectangle((pk - wbox, 0.0), 2 * wbox, 1.0, fill=False, ec='#00A651', lw=0.9))
    fixax.save(fig, f'IGV_{name}'); plt.close(fig); print('done', name, 'RNA tracks:', len(rn))
