"""Helper: figure whose plot box (axes) has an exact size in inches, so panels line up across slides."""
import os, sys
sys.path.insert(0, os.getcwd())      # run from the repository root
from config import DATA_DIR, OUT_DIR
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
FRAME = '#4D4D4D'   # dark-grey frame, same as Fig 1f
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica'],
                     'pdf.fonttype': 42, 'svg.fonttype': 'none', 'axes.linewidth': 0.6, 'axes.edgecolor': FRAME, 'xtick.color': FRAME, 'ytick.color': FRAME,
                     'xtick.labelcolor': 'black', 'ytick.labelcolor': 'black'})
def fig_box(w, h, left=0.5, right=0.1, bottom=0.45, top=0.1):
    W, H = left + w + right, bottom + h + top
    fig = plt.figure(figsize=(W, H)); ax = fig.add_axes([left / W, bottom / H, w / W, h / H])
    for s in ax.spines.values(): s.set_visible(True); s.set_linewidth(0.6); s.set_color(FRAME)
    return fig, ax
def save(fig, path):
    for e in ['png', 'pdf', 'svg']: fig.savefig(f'{path}.{e}', dpi=600)
