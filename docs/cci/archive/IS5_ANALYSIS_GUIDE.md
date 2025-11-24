# IS5 Whole Dataset CCI Analysis Guide

This guide describes how to run CCI analysis on the complete IS5 dataset for all receiver clusters.

## Overview

The analysis pipeline:
1. Tests with heavily downsampled data (max 10 cells per patient per cluster)
2. Runs full analysis on IS5 original data for all receiver clusters
3. Saves results to `sobj/cci/whole/`

## Prerequisites

- IS5 Seurat object: `/data/user3/sobj/IS_scvi_251107.qs` (or similar)
- NicheNet data: `/data/user3/git_repo/human`
- Consensus DEG files (optional): `/data/user3/sobj/consensus/final/deg_consensus_ds_*_final_result_*.qs`

## Step 1: Test with Downsampled Data

Test the pipeline with a heavily downsampled dataset:

```bash
cd /home/user3/data_user3/git_repo/_wt/cci
Rscript scripts/cci/test_is5_downsample.R
```

This will:
- Load IS5 data
- Downsample to max 10 cells per patient per cluster
- Save downsampled data to `/data/user3/sobj/IS5_ds10_per_patient_cluster.qs`
- Run CCI analysis on one receiver cluster
- Save test results to `/data/user3/sobj/cci/test_ds10/`

## Step 2: Run Full Analysis

Run CCI analysis for all receiver clusters:

```bash
cd /home/user3/data_user3/git_repo/_wt/cci
Rscript scripts/cci/run_is5_whole_analysis.R
```

This will:
- Load IS5 original data
- For each receiver cluster:
  - Try to load pre-computed DEGs from consensus results
  - If not available, compute DEGs using `FindMarkers`
  - Run CCI analysis
  - Save results to `sobj/cci/whole/receiver_<cluster_name>/`
- Generate summary report

## Output Structure

```
sobj/cci/whole/
├── receiver_<cluster1>/
│   ├── runX/
│   │   ├── nichenet_results.qs
│   │   ├── NicheNet_Circos_LR_*.pdf
│   │   ├── NicheNet_Circos_LR_*.png
│   │   ├── NicheNet_Ligand_Target_Heatmap.png
│   │   └── ...
│   └── ...
├── receiver_<cluster2>/
│   └── ...
└── analysis_summary.qs
```

## DEG Computation Methods

The pipeline supports two methods for DEG computation:

### Method 1: Direct FindMarkers

```r
deg_df <- FindMarkers(
  subset(is5, anno3.scvi == "Platelets / Megakaryocytes"),
  group.by = "g3",
  ident.1 = "2"
) %>%
  marker_filter() %>%
  mutate(cluster = "Platelets / Megakaryocytes")
```

### Method 2: Pre-computed Consensus DEGs

```r
deg_plt <- qs::qread("/data/user3/sobj/consensus/final/deg_consensus_ds_*_final_result_plt.qs")
deg_df <- deg_plt$consensus_scores %>%
  mutate(
    p_val = meta_p,
    p_val_adj = meta_p_adj,
    avg_log2FC = mean_beta,
    cluster = "Platelets / Megakaryocytes"
  )
```

The pipeline automatically tries to load consensus DEGs first, then falls back to FindMarkers if not available.

## Parameters

### DEG Filtering Logic

**Important**: If `top_n_targets_per_ligand` is specified (not NULL), cutoff filters (`p_val_adj_cutoff`, `logfc_cutoff`) are **ignored**. Top N DEGs are selected by:
1. logFC (descending) - highest fold-change first
2. p-value (ascending) - most significant second

This ensures that when `top_n_targets_per_ligand` is set, you get the top N DEGs by effect size regardless of p-value thresholds. This is particularly useful when:
- `p_val_adj_cutoff = 1.0` (all DEGs included, no p-value filter)
- You want to focus on the most strongly upregulated genes

If `top_n_targets_per_ligand` is NULL, cutoff filters are applied as normal.

### Key Parameters

For downsampled test:
- `p_val_adj_cutoff = 1.0` (effectively no filter when using top_n)
- `logfc_cutoff = 0.1` (ignored when using top_n)
- `top_n_targets_per_ligand = 50` (determines top N DEGs by logFC)

For full analysis, default parameters are used (can be adjusted in the script).

## Troubleshooting

1. **No DEGs after filtering**: Lower `p_val_adj_cutoff` and `logfc_cutoff` in test script
2. **SCT assay errors**: Script automatically falls back to RNA assay
3. **Missing consensus DEGs**: Script automatically computes using FindMarkers
4. **Empty ligand-target heatmap**: Script automatically retries with lower cutoff

## Notes

- Analysis time varies by cluster size and DEG count
- Large DEG sets (>3000) automatically use `top_n_targets_per_ligand = 50` for performance
- Circos plots are saved in multiple versions (with/without legend, legend only)
- All results are saved as `.qs` files for easy reloading

