# CHO Transcriptome Compendium

Code for the CHO transcriptome compendium paper (Lakshmanan lab, SBCE, IIT Madras).

Paper: TBD | Data: TBD (SRA BioProjects / Zenodo) | Licence: TBD

## File structure

```
.
├── config.R, config.py              # input (data/) and output (results/) paths, used by every script
├── environment/                     # requirements.txt, install_R_packages.R
├── tissue_vs_CHO_gain_loss_reduced/
│   ├── Gain_loss_run_all_subsets.R  # gained / lost / reduced genes, 4 CHO subsets vs tissues
│   └── boxplots.R                   # tissue vs CHO boxplots with Mann-Whitney tests
├── cellline_DEG_vs_K1/
│   └── Cellline_combined.R          # S, DXB11, DG44 vs K1 (K1 + K1 (GS) merged), DESeq2
├── promoter_activity_proActiv/
│   ├── promoter_gene_DESeq2.R       # promoter- and gene-level DESeq2 vs K1 (raw proActiv counts)
│   └── final_from_Prahalaad_2026-10-27/
│       ├── proactiv_run.R               # proActiv on the BAM files
│       ├── proactiv_sample_merging.R    # merge per-sample results
│       └── Promoter_analysis_final.ipynb  # alternative-promoter analysis
├── enrichment_GO_KEGG/
│   ├── Enrichment_clusterProfiler_main.R  # GO and KEGG enrichment
│   └── Enrichment_plots.R
├── PCA_MDS_UMAP/                  # PCA, MDS, UMAP
├── figure_scripts/                # figure panels (Python)
├── media_PyDESeq2/
│   └── media_K1_pydeseq2.py         # media vs CD-CHO in CHO-K1
├── secretory_SecRecon/
│   ├── SecRecon_overlap_enrichment.py        # hypergeometric test, system / subsystem
│   ├── SecRecon_enrichment_bubble_plots.R
│   └── SecRecon_system_subsystem_sunburst.py
└── copy_number_CNV/
    ├── CNV_classify_calls.py        # CNV class from ratio to K1 average
    ├── CNV_DESeq2_integration.py    # CNV x DESeq2 join and concordance classes
    └── CNV_venn.R                   # CNV vs DEG overlap
```

## Running

1. `pip install -r environment/requirements.txt` and `Rscript environment/install_R_packages.R`
2. Put the input files in `data/` (or set `CHO_DATA_DIR`). Results go to `results/` (or `CHO_OUT_DIR`).
3. Run each script from the repository root, for example `Rscript cellline_DEG_vs_K1/Cellline_combined.R`
4. Run `tissue_vs_CHO_gain_loss_reduced`, `cellline_DEG_vs_K1`, `promoter_activity_proActiv`, `PCA_MDS_UMAP`, `media_PyDESeq2` and `copy_number_CNV` first, then `enrichment_GO_KEGG`, `figure_scripts` and `secretory_SecRecon`, which use their outputs.

The count matrix (`Count_raw_all.csv`) is generated from raw reads as described in the Methods. Seeds are fixed to 42 for UMAP and the promoter notebook.
