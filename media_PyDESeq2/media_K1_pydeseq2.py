# Media vs CD-CHO (Gibco) differential expression, K1 samples only, 9 media.
# Tool: PyDESeq2 (Python port of DESeq2), design ~Media. Genes with total count < 10 removed.
# Significant: padj < 0.05 and |log2FC| > 1  (no other filter).
# Inputs (same folder): Count_raw_all.csv (genes x samples), Metadata_Media.csv (Run, Study, Producer, Cell_line, Media, Brand, Library)
# Outputs: DESeq2_K1_media_FULL.xlsx (ALL / UP / DOWN sheet per medium), DESeq2_K1_media_summary.csv
# Needs: pip install pydeseq2 openpyxl pandas
import os, sys
sys.path.insert(0, os.getcwd())      # run from the repository root
from config import DATA_DIR, OUT_DIR
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

counts = pd.read_csv(os.path.join(DATA_DIR, "Count_raw_all.csv"), index_col=0)
counts.columns = [c.replace(".bam", "") for c in counts.columns]

meta = pd.read_csv(os.path.join(DATA_DIR, "Metadata_Media.csv"))
meta = meta[(meta.Cell_line == "K1") & meta.Run.isin(counts.columns)].set_index("Run")

counts = counts[meta.index].round().astype(int).T      # samples x genes
counts = counts.loc[:, counts.sum(axis=0) >= 10]
print("samples:", counts.shape[0], "genes:", counts.shape[1])

ref = "CD-CHO (Gibco)"
dds = DeseqDataSet(counts=counts, metadata=meta[["Media"]], design="~Media", refit_cooks=True, quiet=True)
dds.deseq2()

short = {"ActiPro (Hyclone)": "ActiPro", "CDM4CHO (HyClone)": "CDM4CHO", "CHO-S-SFM II (Gibco)": "CHO-S-SFM II",
         "DA1 (Gibco)": "DA1", "EX-CELL Advanced CHO (SAFC)": "EX-CELL Adv", "EX-CELL CD CHO (SAFC)": "EX-CELL CD CHO",
         "FortiCHO (Gibco)": "FortiCHO", "G13 (Fujifilm)": "G13", "SFM4CHO (Hyclone)": "SFM4CHO"}

summary = []
with pd.ExcelWriter(os.path.join(OUT_DIR, "DESeq2_K1_media_FULL.xlsx")) as writer:
    for medium, name in short.items():
        stats = DeseqStats(dds, contrast=["Media", medium, ref], quiet=True)
        stats.summary()
        res = stats.results_df
        res.index.name = "Gene"
        up = res[(res.padj < 0.05) & (res.log2FoldChange > 1)]
        down = res[(res.padj < 0.05) & (res.log2FoldChange < -1)]
        res, up, down = res.round(6), up.round(6), down.round(6)   # filter on unrounded values, round only for writing
        res.reset_index().to_excel(writer, sheet_name=name + "_ALL", index=False)
        up.reset_index().to_excel(writer, sheet_name=name + "_UP", index=False)
        down.reset_index().to_excel(writer, sheet_name=name + "_DOWN", index=False)
        summary.append({"Medium": name, "Up": len(up), "Down": len(down)})
        print(name, len(up), len(down))

pd.DataFrame(summary).to_csv(os.path.join(OUT_DIR, "DESeq2_K1_media_summary.csv"), index=False)
