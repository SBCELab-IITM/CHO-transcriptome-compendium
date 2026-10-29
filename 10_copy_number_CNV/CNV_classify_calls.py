# CNV class per CHO line from the ratio to the CHO-K1 average copy number.
# written to reproduce CNV_DESeq2_Integrated_Analysis_Updated_latest-30:3:26.xlsx
import os, sys
sys.path.insert(0, os.getcwd())      # run from the repository root
from config import DATA_DIR, OUT_DIR
import pandas as pd, numpy as np

def cnv_label(r):
    if not np.isfinite(r):
        # K1 average 0 and line > 0 -> infinite ratio, called as the highest amplification class; 0/0 -> Normal
        return "Amplification_8x" if r == np.inf else "Normal"
    if 0.5 <= r <= 2: return "Normal"
    if r > 2: return "Amplification_%dx" % min(2*np.floor(r/2), 8)       # classes 2x,4x,6x,8x (8x = 8 or higher)
    inv = np.inf if r == 0 else 1.0/r
    if inv >= 16: return "Deletion_0.0625x"
    p = 2**np.floor(np.log2(inv))
    return "Deletion_%sx" % {2:"0.5",4:"0.25",8:"0.125"}[int(p)]

def run(norm_csv):
    n = pd.read_csv(norm_csv).set_index("Normalized")
    avg = (n["CHO-K1 ATCC"] + n["CHO-K1 SF"] + n["CHO-K1 PF"]) / 3   # CHO_K1_avg = (PF+SF+ATCC)/3; summation order matters for ratios exactly at 2.0 (e.g. Cd74)
    out = pd.DataFrame({"CHO_K1_avg": avg})
    for name, col in [("CHO_S","CHO-S"),("CHO_DXB11","CHO DXB11"),("CHO_DG44","CHO DG44")]:
        ratio = n[col]/avg
        out[name+"_ratio"] = ratio
        out[name+"_CNV"] = ratio.map(cnv_label)
    return out

if __name__ == "__main__":
    import sys
    run(sys.argv[1] if len(sys.argv) > 1 else os.path.join(DATA_DIR, "Copynumber_norm.csv")).to_csv(os.path.join(OUT_DIR, "CNV_calls.csv"))
