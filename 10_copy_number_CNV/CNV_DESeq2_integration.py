"""Integrate per-gene copy-number calls with DESeq2 results (cell line vs CHO-K1).

Reconstructed from, and checked against, CNV_DESeq2_Integrated_Analysis_Updated_latest-30:3:26.xlsx
(the lab's original join was done interactively and no script survives).

Steps
 1. CNV calls: CNV_classify_calls.run() on Copynumber_norm.csv  (ratio = line / mean(K1 PF, K1 SF, K1 ATCC))
 2. DESeq2 results: <line>_vs_K1_All_DESeq_genes.csv (gene ids carry a '_<n>' suffix; base symbol = suffix stripped)
 3. Join on base symbol (inner join: genes present in both CNV table and DESeq2 table)
 4. Expression_Change : Upregulated  padj < 0.01 & log2FC >  2 ; Downregulated padj < 0.01 & log2FC < -2 ; else Not_Significant
 5. Concordance      : Normal CNV -> No_CNV;  Amplification/Deletion x Up/Down/NS ->
                       CNV_Driven_Up | Compensated_Amp | Silent_Amp  /  Compensated_Del | CNV_Driven_Down | Silent_Del
Inputs (data/): Copynumber_norm.csv, S_vs_K1_All_DESeq_genes.csv, DXB11_vs_K1_All_DESeq_genes.csv, DG44_vs_K1_All_DESeq_genes.csv
Output (results/): CNV_DESeq2_<line>.csv and CNV_DESeq2_summary.csv
"""
import os, sys
sys.path.insert(0, os.getcwd())
from config import DATA_DIR, OUT_DIR
import importlib.util
import pandas as pd, numpy as np

_here = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("CNV_classify_calls", os.path.join(_here, "CNV_classify_calls.py"))
_cc = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(_cc)

PADJ, LFC = 0.01, 2.0
LINES = {"CHO_S": "S", "CHO_DXB11": "DXB11", "CHO_DG44": "DG44"}

def concordance(cnv, expr):
    if cnv == "Normal": return "No_CNV"
    amp = cnv.startswith("Amplification")
    if expr == "Upregulated":   return "CNV_Driven_Up" if amp else "Compensated_Del"
    if expr == "Downregulated": return "Compensated_Amp" if amp else "CNV_Driven_Down"
    return "Silent_Amp" if amp else "Silent_Del"

def integrate(calls, line, prefix):
    d = pd.read_csv(os.path.join(DATA_DIR, f"{prefix}_vs_K1_All_DESeq_genes.csv"))
    d["base"] = d["gene"].str.replace(r"_\d+$", "", regex=True)
    m = calls.reset_index().rename(columns={"Normalized": "Gene_norm"}).merge(d, left_on="Gene_norm", right_on="base", how="inner")
    up = (m.padj < 0.01) & (m.log2FoldChange > LFC); dn = (m.padj < 0.01) & (m.log2FoldChange < -LFC)
    m["Expression_Change"] = np.where(up, "Upregulated", np.where(dn, "Downregulated", "Not_Significant"))
    m["Concordance"] = [concordance(c, e) for c, e in zip(m[f"{line}_CNV"], m.Expression_Change)]
    return m

if __name__ == "__main__":
    calls = _cc.run(os.path.join(DATA_DIR, "Copynumber_norm.csv"))
    summ = []
    for line, prefix in LINES.items():
        m = integrate(calls, line, prefix); m.to_csv(os.path.join(OUT_DIR, f"CNV_DESeq2_{line}.csv"), index=False)
        vc = m.Concordance.value_counts()
        summ.append(dict(Cell_line=line, Total_genes=len(m),
                         CNV_driven=vc.get("CNV_Driven_Up", 0) + vc.get("CNV_Driven_Down", 0),
                         Compensated=vc.get("Compensated_Amp", 0) + vc.get("Compensated_Del", 0),
                         Silent=vc.get("Silent_Amp", 0) + vc.get("Silent_Del", 0)))
    pd.DataFrame(summ).to_csv(os.path.join(OUT_DIR, "CNV_DESeq2_summary.csv"), index=False)
    print(pd.DataFrame(summ))
