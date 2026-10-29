"""SecRecon enrichment of gene sets (lost / reduced in CHO vs tissues) by System and Subsystem.

Method (reconstructed from, and verified against, Loss_Reduced_System_Subsystem_Analysis.xlsx):
  universe M = unique CHO genes in the SecRecon mapping (Sys&Subsys.csv; 965 genes)
  n = genes of the query set that map to SecRecon;  k = of these, genes in the (sub)system;  K = genes of the (sub)system in SecRecon
  p  = one-sided hypergeometric (= Fisher exact, alternative 'greater') P(X >= k)   -> 'p-val' (RAW, not adjusted)
  FC = (k/n) / (K/M)                                                                  -> 'FC' (enrichment ratio)
Also writes a Benjamini-Hochberg adjusted column (p_BH) per table, because the manuscript legend states BH adjustment.

Input : data/Sys&Subsys.csv (CHO GENE SYMBOL, System, Subsystem)
        a gene-set table with a 'CHO GENE SYMBOL' column (the query genes, mapped or unmapped), via --genes
Usage : python SecRecon_overlap_enrichment.py --genes reduced_common_genes.csv --out results/reduced_intersection
"""
import os, sys, argparse
sys.path.insert(0, os.getcwd())
from config import DATA_DIR, OUT_DIR
import pandas as pd, numpy as np
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

def enrich(genes, ref, level):
    ref = ref.drop_duplicates(["CHO GENE SYMBOL", level])
    M = ref["CHO GENE SYMBOL"].nunique()
    q = ref[ref["CHO GENE SYMBOL"].isin(set(genes))]
    n = q["CHO GENE SYMBOL"].nunique()
    rows = []
    for term, g in q.groupby(level):
        k = g["CHO GENE SYMBOL"].nunique(); K = (ref[level] == term).sum()
        rows.append(dict(**{level: term}, k=k, K=K, n=n, M=M,
                         p=hypergeom.sf(k - 1, M, K, n), FC=(k / n) / (K / M),
                         genes=";".join(sorted(g["CHO GENE SYMBOL"]))))
    r = pd.DataFrame(rows).sort_values("p")
    r["p_BH"] = multipletests(r.p, method="fdr_bh")[1]
    return r

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--genes", required=True); ap.add_argument("--out", default=os.path.join(OUT_DIR, "secrecon"))
    a = ap.parse_args()
    ref = pd.read_csv(os.path.join(DATA_DIR, "Sys&Subsys.csv"))
    genes = pd.read_csv(a.genes)["CHO GENE SYMBOL"].dropna().unique()
    for level in ["System", "Subsystem"]:
        enrich(genes, ref, level).to_csv(f"{a.out}_{level}.csv", index=False)
