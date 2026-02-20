#!/usr/bin/env Rscript
# ==============================================================================
# Phase 4 Part 2: CD8 trajectory + Gene dynamics fix + Plots
# Picks up from part 1 which completed Mono and CD4 trajectories
# ==============================================================================

# Setup library paths FIRST
.libPaths(c(
  "/home/user3/R/x86_64-pc-linux-gnu-library/4.3",
  "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu",
  .libPaths()
))

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(slingshot)
  library(monocle3)
  library(S4Vectors)
  library(SummarizedExperiment)
  # Load Bioconductor packages BEFORE dplyr to ensure dplyr masks correctly
  library(Seurat)
  library(SeuratObject)
  library(qs)
  library(ggplot2)
  library(patchwork)
  library(tibble)
  library(mgcv)
  # dplyr LAST to ensure first() etc. are available
  library(dplyr)
})

# Source pseudotime functions
source("/data/user3/git_repo/mylit/myR/R/analysis/pseudotime.R")

cat("=== Phase 4 Part 2: CD8 + Gene Dynamics + Plots ===\n")
cat("Started:", format(Sys.time()), "\n\n")

BASE <- "/data/user3/sobj/stroke_hc_v8_2"
TRAJ_DIR <- file.path(BASE, "trajectory")

TARGET_GENES <- c(
  "S100A8", "S100A9", "S100A12", "VCAN", "FCN1", "CD14", "FCGR3A",
  "IL1B", "CXCL8", "CCL2", "CCL3", "TNF", "IL6",
  "ISG15", "IFI6", "IFIT1", "MX1", "OAS1",
  "GZMB", "GZMK", "PRF1", "NKG7", "GNLY", "IFNG",
  "FOXP3", "IL2RA", "CTLA4", "TIGIT", "LAG3", "PDCD1",
  "DDIT4", "TXNIP", "HLA-B", "HLA-C", "CCL4", "HIF1A"
)

# ============================================================
# PART A: CD8 trajectory (the missing piece)
# ============================================================
cat("=" , paste(rep("=", 50), collapse = ""), "\n")
cat("PART A: CD8+ T cell trajectory\n")
cat(paste(rep("=", 50), collapse = ""), "\n\n")

cd8_dir <- file.path(TRAJ_DIR, "cd8")
dir.create(cd8_dir, showWarnings = FALSE, recursive = TRUE)

# Load clean object
cat("[1] Loading clean Seurat object...\n")
sobj <- qread(file.path(BASE, "5_strokev8_clean.qs"))
cat("  Loaded:", ncol(sobj), "cells\n")

# CD8 subset
cat("[2] Creating CD8+ T cell subset...\n")
cd8_anno1 <- c("CD8+ T_Cytotoxic", "CD8+ Trm", "MAIT")
cells_cd8 <- colnames(sobj)[sobj$anno1 %in% cd8_anno1]
cat("  Cells selected:", length(cells_cd8), "\n")
tbl <- table(sobj$anno1[sobj$anno1 %in% cd8_anno1])
for (nm in names(tbl)) cat("    ", nm, ":", tbl[nm], "\n")

sub <- subset(sobj, cells = cells_cd8)
sub <- NormalizeData(sub, verbose = FALSE)
sub <- FindVariableFeatures(sub, nfeatures = 2000, verbose = FALSE)
sub <- ScaleData(sub, verbose = FALSE)
sub <- RunPCA(sub, npcs = 30, verbose = FALSE)
sub <- RunUMAP(sub, dims = 1:20, reduction = "pca",
               reduction.name = "umap_subset", verbose = FALSE)
sub <- FindNeighbors(sub, dims = 1:20, verbose = FALSE)
sub <- FindClusters(sub, resolution = 0.5, verbose = FALSE)
qsave(sub, file.path(cd8_dir, "subset.qs"))
cat("  Saved CD8 subset\n")

# Slingshot
cat("[3] Running Slingshot...\n")
meta <- sub@meta.data
meta_start <- meta[meta$anno1 == "CD8+ T_Cytotoxic", ]
start_tbl <- as.data.frame(table(cluster = meta_start$seurat_clusters),
                            stringsAsFactors = FALSE)
start_cl <- as.character(start_tbl$cluster[which.max(start_tbl$Freq)])
cat("  Start cluster:", start_cl, "\n")

sce <- tryCatch({
  run_slingshot_from_seurat(
    seurat_obj = sub,
    cluster_col = "seurat_clusters",
    reduced_dim_name = "umap_subset",
    start_cluster = start_cl,
    counts_assay_name = "RNA",
    approx_points = 150
  )
}, error = function(e) {
  cat("  Slingshot ERROR:", e$message, "\n")
  NULL
})
if (!is.null(sce)) {
  saveRDS(sce, file.path(cd8_dir, "slingshot_sce.rds"))
  cat("  Saved Slingshot SCE\n")
}

# Monocle3
cat("[4] Running Monocle3...\n")
start_cells <- colnames(sub)[sub$anno1 == "CD8+ T_Cytotoxic"]
if (length(start_cells) > 500) {
  set.seed(42)
  start_cells <- sample(start_cells, 500)
}
cds <- tryCatch({
  run_monocle3_from_seurat(
    seurat_obj = sub,
    counts_assay_name = "RNA",
    reduction_name = "umap_subset",
    root_cells = start_cells,
    preprocess_method = "PCA",
    num_dimensions = 20,
    reduction_method = "UMAP",
    use_partition = FALSE,
    resolution = 1e-4,
    verbose = TRUE
  )
}, error = function(e) {
  cat("  Monocle3 ERROR:", e$message, "\n")
  NULL
})
if (!is.null(cds)) {
  saveRDS(cds, file.path(cd8_dir, "monocle3_cds.rds"))
  pt <- tryCatch(monocle3::pseudotime(cds), error = function(e) NULL)
  if (!is.null(pt)) {
    pt_finite <- pt[is.finite(pt)]
    cat("  Pseudotime range:", round(min(pt_finite), 2), "-",
        round(max(pt_finite), 2), "\n")
    cat("  Inf cells:", sum(is.infinite(pt)), "/", length(pt), "\n")
  }
  cat("  Saved Monocle3 CDS\n")
}

# Free memory
rm(sobj); gc(verbose = FALSE)

# ============================================================
# PART B: Generate UMAP plots for all compartments
# ============================================================
cat("\n", paste(rep("=", 50), collapse = ""), "\n")
cat("PART B: Generating UMAP plots\n")
cat(paste(rep("=", 50), collapse = ""), "\n\n")

for (comp_id in c("mono", "cd4", "cd8")) {
  comp_dir <- file.path(TRAJ_DIR, comp_id)
  cat("Plotting:", comp_id, "...\n")

  # Load subset
  sub_file <- file.path(comp_dir, "subset.qs")
  if (!file.exists(sub_file)) {
    cat("  No subset file. Skipping.\n")
    next
  }
  sub <- qread(sub_file)

  # Load Monocle3 CDS
  cds_file <- file.path(comp_dir, "monocle3_cds.rds")
  if (!file.exists(cds_file)) {
    cat("  No CDS file. Skipping.\n")
    next
  }
  cds <- readRDS(cds_file)

  # Get pseudotime
  pt <- tryCatch(monocle3::pseudotime(cds), error = function(e) NULL)
  if (is.null(pt)) {
    cat("  No pseudotime. Skipping.\n")
    next
  }

  # Build plot data from subset UMAP
  umap_coords <- sub@reductions$umap_subset@cell.embeddings
  common_cells <- intersect(rownames(umap_coords), names(pt))
  plot_df <- data.frame(
    UMAP1 = umap_coords[common_cells, 1],
    UMAP2 = umap_coords[common_cells, 2],
    anno1 = sub$anno1[common_cells],
    pseudotime = pt[common_cells],
    stringsAsFactors = FALSE
  )

  # Add cohort and g3 if available
  if ("cohort" %in% colnames(sub@meta.data)) {
    plot_df$cohort <- sub$cohort[common_cells]
  }
  if ("g3" %in% colnames(sub@meta.data)) {
    plot_df$g3 <- sub$g3[common_cells]
  }

  # Plot 1: anno1
  p1 <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = anno1)) +
    geom_point(size = 0.3, alpha = 0.5) +
    theme_minimal() +
    ggtitle(paste(comp_id, "- Cell Type")) +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(override.aes = list(size = 2)))

  # Plot 2: pseudotime
  p2 <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = pseudotime)) +
    geom_point(size = 0.3, alpha = 0.5) +
    scale_color_viridis_c(na.value = "grey85") +
    theme_minimal() +
    ggtitle(paste(comp_id, "- Monocle3 Pseudotime"))

  # Plot 3: cohort
  p3 <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = cohort)) +
    geom_point(size = 0.3, alpha = 0.3) +
    theme_minimal() +
    ggtitle(paste(comp_id, "- Cohort (HC vs Stroke)")) +
    theme(legend.position = "bottom")

  # Plot 4: g3
  p4 <- ggplot(plot_df[!is.na(plot_df$g3), ],
               aes(x = UMAP1, y = UMAP2, color = factor(g3))) +
    geom_point(size = 0.3, alpha = 0.3) +
    theme_minimal() +
    ggtitle(paste(comp_id, "- g3 (IS only)")) +
    theme(legend.position = "bottom") +
    scale_color_manual(values = c("1" = "#2166AC", "2" = "#B2182B"),
                       name = "g3", labels = c("Good", "Bad"))

  combined <- (p1 + p2) / (p3 + p4)
  ggsave(file.path(comp_dir, "trajectory_overview.pdf"),
         combined, width = 14, height = 12)
  cat("  Saved trajectory_overview.pdf\n")

  # Pseudotime density by g3
  if ("g3" %in% colnames(plot_df)) {
    pt_g3 <- plot_df[!is.na(plot_df$g3), ]
    if (nrow(pt_g3) > 100) {
      p_density <- ggplot(pt_g3,
                          aes(x = pseudotime, fill = factor(g3))) +
        geom_density(alpha = 0.5) +
        scale_fill_manual(values = c("1" = "#2166AC", "2" = "#B2182B"),
                          name = "g3", labels = c("Good", "Bad")) +
        theme_minimal() +
        ggtitle(paste(comp_id, "- Pseudotime density by g3")) +
        xlab("Pseudotime") + ylab("Density")

      p_box <- ggplot(pt_g3,
                      aes(x = factor(g3), y = pseudotime, fill = factor(g3))) +
        geom_boxplot(outlier.size = 0.3) +
        scale_fill_manual(values = c("1" = "#2166AC", "2" = "#B2182B"),
                          name = "g3", labels = c("Good", "Bad")) +
        theme_minimal() +
        ggtitle(paste(comp_id, "- Pseudotime by g3")) +
        xlab("g3 (IAT outcome)") + ylab("Pseudotime")

      combined_g3 <- p_density + p_box
      ggsave(file.path(comp_dir, "pseudotime_g3_comparison.pdf"),
             combined_g3, width = 12, height = 5)
      cat("  Saved pseudotime_g3_comparison.pdf\n")

      # Wilcoxon test
      wt <- wilcox.test(pseudotime ~ factor(g3), data = pt_g3)
      cat("  Wilcoxon test (g3==1 vs g3==2):",
          "W =", wt$statistic,
          ", p =", format.pval(wt$p.value, digits = 3), "\n")
    }
  }

  # Pseudotime by anno1 boxplot
  p_anno_pt <- ggplot(plot_df,
                      aes(x = reorder(anno1, pseudotime, median),
                          y = pseudotime, fill = anno1)) +
    geom_boxplot(outlier.size = 0.3) +
    theme_minimal() +
    ggtitle(paste(comp_id, "- Pseudotime by cell type")) +
    xlab("") + ylab("Pseudotime") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  ggsave(file.path(comp_dir, "pseudotime_by_celltype.pdf"),
         p_anno_pt, width = 8, height = 5)
  cat("  Saved pseudotime_by_celltype.pdf\n")

  rm(sub, cds); gc(verbose = FALSE)
}

# ============================================================
# PART C: Gene dynamics with dplyr::first() fix
# ============================================================
cat("\n", paste(rep("=", 50), collapse = ""), "\n")
cat("PART C: Gene dynamics (with fix)\n")
cat(paste(rep("=", 50), collapse = ""), "\n\n")

# Override the problematic function: replace first() with dplyr::first()
# by monkey-patching analyze_gene_dynamics
# Actually easier: just ensure dplyr is loaded LAST and check
cat("dplyr::first is:", class(dplyr::first), "\n")
cat("first() resolves to:", environmentName(environment(first)), "\n")

# Run simplified gene dynamics for each compartment
for (comp_id in c("mono", "cd4", "cd8")) {
  comp_dir <- file.path(TRAJ_DIR, comp_id)
  cat("\nGene dynamics for:", comp_id, "...\n")

  cds_file <- file.path(comp_dir, "monocle3_cds.rds")
  if (!file.exists(cds_file)) {
    cat("  No CDS. Skipping.\n")
    next
  }
  cds <- readRDS(cds_file)

  # Check pseudotime
  pt <- tryCatch(monocle3::pseudotime(cds), error = function(e) NULL)
  if (is.null(pt) || all(is.infinite(pt))) {
    cat("  No valid pseudotime. Skipping.\n")
    next
  }

  # Check g3
  cell_meta <- as.data.frame(colData(cds))
  if ("g3" %in% colnames(cell_meta) && sum(!is.na(cell_meta$g3)) > 100) {
    condition_col <- "g3"
    cat("  Condition: g3 (", sum(!is.na(cell_meta$g3)), "cells)\n")
  } else {
    cat("  No g3. Skipping dynamics.\n")
    next
  }

  dyn_dir <- file.path(comp_dir, "gene_dynamics")
  dir.create(dyn_dir, showWarnings = FALSE, recursive = TRUE)

  available_genes <- intersect(TARGET_GENES, rownames(cds))
  cat("  Genes available:", length(available_genes), "/",
      length(TARGET_GENES), "\n")

  results <- list()
  for (gene in available_genes) {
    cat("  ", gene, "...")
    res <- tryCatch({
      analyze_gene_dynamics(
        gene_id = gene,
        cds_obj = cds,
        condition_col_name = condition_col,
        output_dir = dyn_dir,
        k_val = 6,
        min_cells_for_fit = 30,
        plot_split = FALSE,
        scale_DR = TRUE
      )
    }, error = function(e) {
      cat(" err:", substr(e$message, 1, 60), "\n")
      list(gene = gene, status = "error", metrics = NULL)
    })
    if (!is.null(res$metrics)) {
      cat(" OK\n")
    }
    results[[gene]] <- res
  }

  # Summary
  summary_df <- do.call(rbind, lapply(results, function(r) {
    if (is.null(r$metrics)) {
      data.frame(gene = r$gene, status = r$status,
                 TV = NA, DR = NA, interaction_pval = NA,
                 stringsAsFactors = FALSE)
    } else {
      data.frame(
        gene = r$gene, status = r$status,
        TV = ifelse(!is.null(r$metrics$TV), r$metrics$TV, NA),
        DR = ifelse(!is.null(r$metrics$DR), r$metrics$DR, NA),
        interaction_pval = ifelse(!is.null(r$metrics$interaction_pval),
                                  r$metrics$interaction_pval, NA),
        stringsAsFactors = FALSE
      )
    }
  }))
  write.csv(summary_df, file.path(dyn_dir, "gene_dynamics_summary.csv"),
            row.names = FALSE)
  cat("  Saved summary:", nrow(summary_df), "genes\n")

  rm(cds); gc(verbose = FALSE)
}

# ============================================================
# Summary
# ============================================================
cat("\n", paste(rep("=", 50), collapse = ""), "\n")
cat("TRAJECTORY ANALYSIS COMPLETE\n")
cat(paste(rep("=", 50), collapse = ""), "\n")
cat("Finished:", format(Sys.time()), "\n\n")

for (comp_id in c("mono", "cd4", "cd8")) {
  comp_dir <- file.path(TRAJ_DIR, comp_id)
  if (dir.exists(comp_dir)) {
    files <- list.files(comp_dir, recursive = TRUE)
    cat(comp_id, ":\n")
    for (f in files) cat("  ", f, "\n")
    cat("\n")
  }
}

cat("Done.\n")
