#!/usr/bin/env Rscript
# ==============================================================================
# Phase 4 v2: Corrected Trajectory Analysis
# FIX: Use original umap.scvi embedding (batch-corrected) instead of recomputing
# FIX: PNG for scatter plots
# ADD: HC vs Stroke (Layer 1) comparison alongside g3 (Layer 2)
# ==============================================================================

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
  library(Seurat)
  library(SeuratObject)
  library(qs)
  library(ggplot2)
  library(patchwork)
  library(tibble)
  library(mgcv)
  library(dplyr)
})

source("/data/user3/git_repo/mylit/myR/R/analysis/pseudotime.R")

cat("=== Phase 4 v2: Trajectory (scVI UMAP, PNG, HC vs Stroke) ===\n")
cat("Started:", format(Sys.time()), "\n\n")

BASE <- "/data/user3/sobj/stroke_hc_v8_2"
TRAJ_DIR <- file.path(BASE, "trajectory_v2")
dir.create(TRAJ_DIR, showWarnings = FALSE, recursive = TRUE)

TARGET_GENES <- c(
  "S100A8", "S100A9", "S100A12", "VCAN", "FCN1", "CD14", "FCGR3A",
  "IL1B", "CXCL8", "CCL2", "CCL3", "TNF", "IL6",
  "ISG15", "IFI6", "IFIT1", "MX1", "OAS1",
  "GZMB", "GZMK", "PRF1", "NKG7", "GNLY", "IFNG",
  "FOXP3", "IL2RA", "CTLA4", "TIGIT", "LAG3", "PDCD1",
  "DDIT4", "TXNIP", "HLA-B", "HLA-C", "CCL4", "HIF1A"
)

# ---- Compartment definitions ----
compartments <- list(
  mono = list(
    name = "Monocyte",
    anno1_include = c("CD14+ Monocyte", "Inflammatory Monocyte",
                      "CD16+ Monocyte", "ISG+ Myeloid"),
    start_hint = "CD14+ Monocyte",
    description = "Classical -> Inflammatory -> Non-classical"
  ),
  cd4 = list(
    name = "CD4+ T cell",
    anno1_include = c("CD4+ T_Naive/Memory", "Treg",
                      "CD4_S100A8_CD14", "ISG+ T_cell"),
    start_hint = "CD4+ T_Naive/Memory",
    description = "Naive/Memory -> Treg / stress-activated"
  ),
  cd8 = list(
    name = "CD8+ T cell",
    anno1_include = c("CD8+ T_Cytotoxic", "CD8+ Trm", "MAIT"),
    start_hint = "CD8+ T_Cytotoxic",
    description = "Cytotoxic -> Trm"
  )
)

# ---- Load data ----
cat("[1] Loading clean Seurat object...\n")
sobj <- qread(file.path(BASE, "5_strokev8_clean.qs"))
cat("  Cells:", ncol(sobj), "| Reductions:", paste(names(sobj@reductions), collapse=", "), "\n")

# Get the ORIGINAL scVI UMAP coordinates (batch-corrected)
umap_scvi <- Embeddings(sobj, "umap.scvi")
cat("  umap.scvi dims:", ncol(umap_scvi), "\n\n")

# ---- Per compartment ----
for (comp_id in names(compartments)) {
  comp <- compartments[[comp_id]]
  cat("\n", paste(rep("=", 60), collapse=""), "\n")
  cat(" ", comp$name, "-", comp$description, "\n")
  cat(paste(rep("=", 60), collapse=""), "\n")

  comp_dir <- file.path(TRAJ_DIR, comp_id)
  dir.create(comp_dir, showWarnings = FALSE, recursive = TRUE)
  dyn_dir <- file.path(comp_dir, "gene_dynamics")
  dir.create(dyn_dir, showWarnings = FALSE, recursive = TRUE)

  # ---- Subset ----
  cells_keep <- colnames(sobj)[sobj$anno1 %in% comp$anno1_include]
  cat("  Cells:", length(cells_keep), "\n")
  tbl <- table(sobj$anno1[sobj$anno1 %in% comp$anno1_include])
  for (nm in names(tbl)) cat("    ", nm, ":", tbl[nm], "\n")

  if (length(cells_keep) < 500) {
    cat("  Too few cells. Skipping.\n"); next
  }

  sub <- subset(sobj, cells = cells_keep)

  # ---- KEY FIX: Use original scVI UMAP, not recomputed ----
  # The umap.scvi reduction is already batch-corrected and respects HC vs Stroke integration
  # We just subset it to current cells
  scvi_umap <- umap_scvi[cells_keep, ]
  cat("  Using original umap.scvi (", nrow(scvi_umap), " cells, ", ncol(scvi_umap), " dims)\n")

  # We still need PCA within subset for Slingshot (it needs cluster distances)
  # But we'll use scVI UMAP for all visualization
  sub <- NormalizeData(sub, verbose=FALSE)
  sub <- FindVariableFeatures(sub, nfeatures=2000, verbose=FALSE)
  sub <- ScaleData(sub, verbose=FALSE)
  sub <- RunPCA(sub, npcs=30, verbose=FALSE)
  sub <- FindNeighbors(sub, dims=1:20, verbose=FALSE)
  sub <- FindClusters(sub, resolution=0.5, verbose=FALSE)

  # Store scVI UMAP as a reduction in subset
  sub[["umap_scvi"]] <- CreateDimReducObject(
    embeddings = scvi_umap,
    key = "UMAPSCVI_",
    assay = "RNA"
  )

  qsave(sub, file.path(comp_dir, "subset.qs"))

  # ---- Slingshot (on umap_scvi 2D) ----
  cat("  Slingshot...\n")
  meta <- sub@meta.data
  meta_start <- meta[meta$anno1 == comp$start_hint, ]
  start_tbl <- as.data.frame(table(cluster=meta_start$seurat_clusters), stringsAsFactors=FALSE)
  start_cl <- as.character(start_tbl$cluster[which.max(start_tbl$Freq)])
  cat("    Start cluster:", start_cl, "\n")

  sce <- tryCatch({
    run_slingshot_from_seurat(
      seurat_obj = sub,
      cluster_col = "seurat_clusters",
      reduced_dim_name = "umap_scvi",
      start_cluster = start_cl,
      counts_assay_name = "RNA",
      approx_points = 150
    )
  }, error = function(e) {
    cat("    Slingshot ERROR:", e$message, "\n"); NULL
  })
  if (!is.null(sce)) {
    saveRDS(sce, file.path(comp_dir, "slingshot_sce.rds"))
    cat("    Saved Slingshot SCE\n")
  }

  # ---- Monocle3 (use scVI UMAP) ----
  cat("  Monocle3...\n")
  start_cells <- colnames(sub)[sub$anno1 == comp$start_hint]
  if (length(start_cells) > 500) {
    set.seed(42); start_cells <- sample(start_cells, 500)
  }

  cds <- tryCatch({
    run_monocle3_from_seurat(
      seurat_obj = sub,
      counts_assay_name = "RNA",
      reduction_name = "umap_scvi",
      root_cells = start_cells,
      preprocess_method = "PCA",
      num_dimensions = 20,
      reduction_method = "UMAP",
      use_partition = FALSE,
      resolution = 1e-4,
      verbose = FALSE
    )
  }, error = function(e) {
    cat("    Monocle3 ERROR:", e$message, "\n"); NULL
  })

  if (!is.null(cds)) {
    saveRDS(cds, file.path(comp_dir, "monocle3_cds.rds"))
    pt <- tryCatch(monocle3::pseudotime(cds), error=function(e) NULL)
    if (!is.null(pt)) {
      pt_f <- pt[is.finite(pt)]
      cat("    Pseudotime:", round(min(pt_f),2), "-", round(max(pt_f),2),
          "| Inf:", sum(is.infinite(pt)), "\n")
    }
  }

  # ---- Plots (all PNG, using scVI UMAP) ----
  cat("  Generating plots (PNG)...\n")

  if (!is.null(cds)) {
    pt <- tryCatch(monocle3::pseudotime(cds), error=function(e) NULL)
  } else {
    pt <- NULL
  }

  plot_df <- data.frame(
    UMAP1 = scvi_umap[cells_keep, 1],
    UMAP2 = scvi_umap[cells_keep, 2],
    anno1 = sub$anno1,
    cohort = sub$cohort,
    stringsAsFactors = FALSE
  )
  if ("g3" %in% colnames(sub@meta.data)) plot_df$g3 <- sub$g3
  if ("index_injury_hand" %in% colnames(sub@meta.data)) plot_df$injury <- sub$index_injury_hand
  if (!is.null(pt)) {
    common <- intersect(rownames(plot_df), names(pt))
    plot_df$pseudotime <- NA_real_
    plot_df[common, "pseudotime"] <- pt[common]
  }

  # P1: Cell type
  p1 <- ggplot(plot_df, aes(UMAP1, UMAP2, color=anno1)) +
    geom_point(size=0.1, alpha=0.3) + theme_minimal() +
    ggtitle(paste(comp$name, "- Cell type (scVI UMAP)")) +
    theme(legend.position="bottom") +
    guides(color=guide_legend(override.aes=list(size=2, alpha=1)))

  # P2: Pseudotime
  if (!is.null(pt)) {
    p2 <- ggplot(plot_df, aes(UMAP1, UMAP2, color=pseudotime)) +
      geom_point(size=0.1, alpha=0.3) +
      scale_color_viridis_c(na.value="grey85") + theme_minimal() +
      ggtitle("Monocle3 Pseudotime")
  } else {
    p2 <- ggplot() + theme_void() + ggtitle("No pseudotime")
  }

  # P3: Cohort (HC vs Stroke)
  p3 <- ggplot(plot_df, aes(UMAP1, UMAP2, color=cohort)) +
    geom_point(size=0.1, alpha=0.2) + theme_minimal() +
    ggtitle("HC vs Stroke") +
    scale_color_manual(values=c("HC"="#2166AC", "Stroke"="#B2182B")) +
    theme(legend.position="bottom")

  # P4: g3 (IS only)
  g3_df <- plot_df[!is.na(plot_df$g3), ]
  if (nrow(g3_df) > 50) {
    p4 <- ggplot(g3_df, aes(UMAP1, UMAP2, color=factor(g3))) +
      geom_point(size=0.1, alpha=0.2) + theme_minimal() +
      ggtitle("g3 (IS only: 1=Good, 2=Bad)") +
      scale_color_manual(values=c("1"="#2166AC","2"="#B2182B"), name="g3") +
      theme(legend.position="bottom")
  } else {
    p4 <- ggplot() + theme_void() + ggtitle("g3 N/A")
  }

  combined <- (p1 + p2) / (p3 + p4)
  ggsave(file.path(comp_dir, "trajectory_overview_scvi.png"),
         combined, width=14, height=12, dpi=150)
  cat("    Saved trajectory_overview_scvi.png\n")

  # ---- Pseudotime density: HC vs Stroke (Layer 1) ----
  if (!is.null(pt)) {
    pt_cohort <- plot_df[is.finite(plot_df$pseudotime), ]
    if (nrow(pt_cohort) > 100) {
      pd1 <- ggplot(pt_cohort, aes(x=pseudotime, fill=cohort)) +
        geom_density(alpha=0.5) +
        scale_fill_manual(values=c("HC"="#2166AC","Stroke"="#B2182B")) +
        theme_minimal() + ggtitle(paste(comp$name, "- Pseudotime: HC vs Stroke")) +
        xlab("Pseudotime") + ylab("Density")
      pb1 <- ggplot(pt_cohort, aes(x=cohort, y=pseudotime, fill=cohort)) +
        geom_boxplot(outlier.size=0.3) +
        scale_fill_manual(values=c("HC"="#2166AC","Stroke"="#B2182B")) +
        theme_minimal() + ggtitle("HC vs Stroke")

      wt_cohort <- wilcox.test(pseudotime ~ cohort, data=pt_cohort)
      cat("    HC vs Stroke Wilcoxon: p =", format.pval(wt_cohort$p.value, digits=3), "\n")

      # g3 density (Layer 2)
      pt_g3 <- pt_cohort[!is.na(pt_cohort$g3), ]
      if (nrow(pt_g3) > 100) {
        pd2 <- ggplot(pt_g3, aes(x=pseudotime, fill=factor(g3))) +
          geom_density(alpha=0.5) +
          scale_fill_manual(values=c("1"="#2166AC","2"="#B2182B"),
                            name="g3", labels=c("Good","Bad")) +
          theme_minimal() + ggtitle("Pseudotime: g3 (IS only)") +
          xlab("Pseudotime") + ylab("Density")
        pb2 <- ggplot(pt_g3, aes(x=factor(g3), y=pseudotime, fill=factor(g3))) +
          geom_boxplot(outlier.size=0.3) +
          scale_fill_manual(values=c("1"="#2166AC","2"="#B2182B"),
                            name="g3", labels=c("Good","Bad")) +
          theme_minimal() + ggtitle("g3 (IS only)")

        wt_g3 <- wilcox.test(pseudotime ~ factor(g3), data=pt_g3)
        cat("    g3 Wilcoxon: p =", format.pval(wt_g3$p.value, digits=3), "\n")

        combined_pt <- (pd1 + pd2) / (pb1 + pb2)
      } else {
        combined_pt <- pd1 + pb1
      }

      ggsave(file.path(comp_dir, "pseudotime_comparison.png"),
             combined_pt, width=12, height=8, dpi=150)
      cat("    Saved pseudotime_comparison.png\n")
    }

    # Pseudotime by cell type
    p_ct <- ggplot(plot_df[is.finite(plot_df$pseudotime), ],
                   aes(x=reorder(anno1, pseudotime, median),
                       y=pseudotime, fill=anno1)) +
      geom_boxplot(outlier.size=0.3) + theme_minimal() +
      ggtitle(paste(comp$name, "- Pseudotime by cell type")) +
      xlab("") + ylab("Pseudotime") +
      theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="none")
    ggsave(file.path(comp_dir, "pseudotime_by_celltype.png"),
           p_ct, width=8, height=5, dpi=150)
    cat("    Saved pseudotime_by_celltype.png\n")
  }

  # ---- Gene dynamics (g3 comparison) ----
  if (!is.null(cds)) {
    cat("  Gene dynamics...\n")
    cell_meta <- as.data.frame(colData(cds))
    has_g3 <- "g3" %in% colnames(cell_meta) && sum(!is.na(cell_meta$g3)) > 100

    if (has_g3) {
      condition_col <- "g3"
      cat("    Condition: g3 (", sum(!is.na(cell_meta$g3)), "cells)\n")
    } else {
      condition_col <- NULL
    }

    if (!is.null(condition_col)) {
      available <- intersect(TARGET_GENES, rownames(cds))
      cat("    Genes:", length(available), "/", length(TARGET_GENES), "\n")

      results <- list()
      for (gene in available) {
        cat("    ", gene, "...")
        res <- tryCatch({
          analyze_gene_dynamics(
            gene_id = gene, cds_obj = cds,
            condition_col_name = condition_col,
            output_dir = dyn_dir, k_val = 6,
            min_cells_for_fit = 30, plot_split = FALSE, scale_DR = TRUE
          )
        }, error = function(e) {
          cat(" err\n")
          list(gene=gene, status="error", metrics=NULL)
        })
        if (!is.null(res$metrics)) cat(" OK\n")
        results[[gene]] <- res
      }

      summary_df <- do.call(rbind, lapply(results, function(r) {
        if (is.null(r$metrics)) {
          data.frame(gene=r$gene, status=r$status, TV=NA, DR=NA, interaction_pval=NA,
                     stringsAsFactors=FALSE)
        } else {
          data.frame(gene=r$gene, status=r$status,
            TV=ifelse(!is.null(r$metrics$TV), r$metrics$TV, NA),
            DR=ifelse(!is.null(r$metrics$DR), r$metrics$DR, NA),
            interaction_pval=ifelse(!is.null(r$metrics$interaction_pval),
                                    r$metrics$interaction_pval, NA),
            stringsAsFactors=FALSE)
        }
      }))
      write.csv(summary_df, file.path(dyn_dir, "gene_dynamics_summary.csv"), row.names=FALSE)
      cat("    Saved summary:", nrow(summary_df), "genes\n")
    }
  }

  cat("  ", comp$name, "complete.\n")
  rm(sub, sce, cds); gc(verbose=FALSE)
}

# ---- Summary ----
cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("TRAJECTORY v2 COMPLETE\n")
cat("Finished:", format(Sys.time()), "\n\n")
for (comp_id in names(compartments)) {
  comp_dir <- file.path(TRAJ_DIR, comp_id)
  if (dir.exists(comp_dir)) {
    files <- list.files(comp_dir, recursive=TRUE)
    cat(comp_id, ":\n")
    for (f in files) cat("  ", f, "\n")
  }
}
cat("Done.\n")
