#!/usr/bin/env Rscript
# HC-only: Validate published aging signatures in Korean PBMC
# Cross-references our FGS results with literature gene sets
#
# Usage: Rscript scripts/hc/validate_aging_signatures.R
# Output: /data/user3/sobj/hc_only_v1/aging_signatures/

.libPaths(c(
  "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu",
  "/home/user3/R/x86_64-pc-linux-gnu-library/4.3",
  .libPaths()
))

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(data.table)
})

cat("=== HC Aging Signature Validation ===\n")
cat(sprintf("Time: %s\n\n", Sys.time()))

out_dir <- "/data/user3/sobj/hc_only_v1/aging_signatures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Define literature gene sets ----
# From major publications
aging_signatures <- list(
  # Luo 2022 (Nat Aging): age-declining naive markers
  naive_decline = c("CCR7", "LEF1", "TCF7", "SELL", "IL7R", "CD27", "CD28",
                     "BACH2", "FOXP1", "SOX4"),

  # Luo 2022: age-increasing effector/cytotoxic
  effector_increase = c("GZMA", "GZMB", "GZMH", "GZMK", "PRF1", "GNLY",
                         "NKG7", "KLRD1", "KLRG1", "FGFBP2"),

  # Connolly 2024 (Aging Cell): loss of identity genes
  identity_loss = c("RPL5", "RPL11", "RPL13", "RPS6", "RPS14", "RPS19",
                     "EEF1A1", "EEF2", "EIF3E", "EIF4A1"),

  # Filippov 2024 consensus aging markers
  filippov_consensus = c("GNLY", "GZMA", "GZMB", "GZMH", "NKG7", "PRF1",
                          "KLRD1", "FGFBP2", "CCL5", "CST7"),

  # Inflammatory / alarmin (age-increasing)
  alarmin_inflammatory = c("S100A8", "S100A9", "S100A12", "S100A4",
                            "CLIC1", "PTGER2", "NEAT1", "MALAT1"),

  # ISG/IFN signaling (age-declining)
  isg_ifn = c("ISG15", "ISG20", "MX1", "MX2", "OAS1", "OAS2", "OAS3",
               "IFIT1", "IFIT2", "IFIT3", "IFI44L", "IFI6", "STAT1",
               "IRF7", "IRF9"),

  # NK senescence (age-increasing)
  nk_senescence = c("CDKN2A", "IGFBP3", "LAG3", "MSC", "FCGR3A",
                     "TIGIT", "PDCD1"),

  # X-escape genes (sex-dimorphic)
  x_escape = c("XIST", "RPS4X", "DDX3X", "EIF2S3", "JPX", "CD99",
                "KDM6A", "KDM5C"),

  # AIDA Korean-relevant
  korean_treg = c("FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2", "RTKN2"),

  # Our FGS top genes (from Treg × age)
  our_fgs_treg_top = c("PLAC8", "EPSTI1", "JAML", "ISG15", "MSC", "DGKA",
                        "ENSG00000274460", "TRAPPC8", "AGMAT", "GPA33"),

  # T cell aging (SOX4-led)
  tcell_aging = c("SOX4", "BLK", "ROBO1", "CASC15", "BCL11B", "SATB1",
                   "THEMIS", "IL6ST", "TOX", "TOX2"),

  # Myeloid aging (SPI1-led)
  myeloid_aging = c("SPI1", "EMP3", "APH1B", "ITGAX", "FCGR3A", "CSF1R",
                     "CD14", "VCAN", "LYZ", "FCN1")
)

# ---- Load Seurat object ----
cat("Loading Seurat object...\n")
sobj <- qread("/data/user3/sobj/hc_only_v1/2_hc_annotated.qs")
cat(sprintf("  %d cells, %d patients\n", ncol(sobj), length(unique(sobj$name))))

# Ensure age_group
if (!"age_group" %in% colnames(sobj@meta.data)) {
  sobj$age_group <- ifelse(sobj$age < 35, "Young",
                           ifelse(sobj$age <= 50, "Middle", "Old"))
}

# ---- Check gene availability ----
all_genes <- rownames(sobj)
cat("\nGene availability in dataset:\n")
sig_summary <- data.frame()
for (sig_name in names(aging_signatures)) {
  genes <- aging_signatures[[sig_name]]
  found <- genes[genes %in% all_genes]
  missing <- genes[!genes %in% all_genes]
  cat(sprintf("  %-25s %d/%d found", sig_name, length(found), length(genes)))
  if (length(missing) > 0 && length(missing) <= 3) {
    cat(sprintf(" (missing: %s)", paste(missing, collapse = ", ")))
  }
  cat("\n")
  sig_summary <- rbind(sig_summary, data.frame(
    signature = sig_name,
    n_total = length(genes),
    n_found = length(found),
    n_missing = length(missing),
    stringsAsFactors = FALSE
  ))
}
write.csv(sig_summary, file.path(out_dir, "signature_availability.csv"),
          row.names = FALSE)

# ---- Score each patient using module scores ----
cat("\nScoring patients with AddModuleScore...\n")
for (sig_name in names(aging_signatures)) {
  genes <- aging_signatures[[sig_name]]
  found <- genes[genes %in% all_genes]
  if (length(found) >= 3) {
    sobj <- AddModuleScore(sobj, features = list(found),
                           name = paste0("sig_", sig_name),
                           ctrl = min(100, length(found) * 5))
  }
}

# Rename (AddModuleScore appends "1")
sig_cols <- grep("^sig_.*1$", colnames(sobj@meta.data), value = TRUE)
for (col in sig_cols) {
  new_name <- sub("1$", "", col)
  sobj@meta.data[[new_name]] <- sobj@meta.data[[col]]
}

# ---- Pseudobulk scoring: patient-level ----
cat("Computing patient-level pseudobulk scores...\n")
sig_cols_clean <- grep("^sig_", colnames(sobj@meta.data), value = TRUE)
sig_cols_clean <- sig_cols_clean[!grepl("1$", sig_cols_clean)]

patient_scores <- sobj@meta.data %>%
  group_by(name) %>%
  summarise(
    age = first(age),
    sex = first(sex),
    age_group = first(age_group),
    n_cells = n(),
    across(all_of(sig_cols_clean), mean, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(patient_scores, file.path(out_dir, "patient_signature_scores.csv"),
          row.names = FALSE)

# ---- Spearman correlation: signature score vs age ----
cat("\nSpearman correlation: signature vs age\n")
cor_results <- data.frame()
for (sig_col in sig_cols_clean) {
  cor_test <- cor.test(patient_scores$age, patient_scores[[sig_col]],
                       method = "spearman", exact = FALSE)
  cor_results <- rbind(cor_results, data.frame(
    signature = sub("sig_", "", sig_col),
    rho = cor_test$estimate,
    p_value = cor_test$p.value,
    direction = ifelse(cor_test$estimate > 0, "increase", "decrease"),
    stringsAsFactors = FALSE
  ))
}
cor_results$padj <- p.adjust(cor_results$p_value, method = "BH")
cor_results <- cor_results[order(cor_results$p_value), ]

write.csv(cor_results, file.path(out_dir, "signature_age_correlation.csv"),
          row.names = FALSE)

cat(sprintf("\n  Results:\n"))
for (i in seq_len(nrow(cor_results))) {
  r <- cor_results[i, ]
  sig_mark <- ifelse(r$padj < 0.05, "***",
                     ifelse(r$padj < 0.1, "**",
                            ifelse(r$p_value < 0.05, "*", "")))
  cat(sprintf("    %-25s rho=%+.3f p=%.2e padj=%.2e %s %s\n",
              r$signature, r$rho, r$p_value, r$padj, r$direction, sig_mark))
}

# ---- Sex stratified ----
cat("\nSex-stratified signature vs age:\n")
sex_cor <- data.frame()
for (sig_col in sig_cols_clean) {
  for (sx in c("F", "M")) {
    sub <- patient_scores[patient_scores$sex == sx, ]
    cor_test <- cor.test(sub$age, sub[[sig_col]],
                         method = "spearman", exact = FALSE)
    sex_cor <- rbind(sex_cor, data.frame(
      signature = sub("sig_", "", sig_col),
      sex = sx,
      rho = cor_test$estimate,
      p_value = cor_test$p.value,
      n = nrow(sub),
      stringsAsFactors = FALSE
    ))
  }
}
write.csv(sex_cor, file.path(out_dir, "signature_age_by_sex.csv"),
          row.names = FALSE)

# ---- Per-compartment scoring ----
cat("\nPer-compartment signature scores vs age:\n")
compartment_cor <- data.frame()
for (comp in sort(unique(sobj$anno2))) {
  comp_meta <- sobj@meta.data[sobj$anno2 == comp, ]
  comp_patient <- comp_meta %>%
    group_by(name) %>%
    summarise(
      age = first(age),
      n_cells = n(),
      across(all_of(sig_cols_clean), mean, na.rm = TRUE),
      .groups = "drop"
    )

  if (nrow(comp_patient) < 20) next

  for (sig_col in sig_cols_clean) {
    cor_test <- cor.test(comp_patient$age, comp_patient[[sig_col]],
                         method = "spearman", exact = FALSE)
    compartment_cor <- rbind(compartment_cor, data.frame(
      compartment = comp,
      signature = sub("sig_", "", sig_col),
      rho = cor_test$estimate,
      p_value = cor_test$p.value,
      n_patients = nrow(comp_patient),
      stringsAsFactors = FALSE
    ))
  }
}
compartment_cor$padj <- p.adjust(compartment_cor$p_value, method = "BH")
compartment_cor <- compartment_cor[order(compartment_cor$p_value), ]
write.csv(compartment_cor, file.path(out_dir, "signature_age_by_compartment.csv"),
          row.names = FALSE)

cat(sprintf("  %d tests, %d sig (padj<0.05)\n",
            nrow(compartment_cor),
            sum(compartment_cor$padj < 0.05, na.rm = TRUE)))

# ---- FGS overlap with literature ----
cat("\nOverlap: Our FGS genes vs literature signatures\n")
fgs_file <- "/data/user3/sobj/hc_only_v1/fgs_continuous_treg_v2/consensus_genes_3plus_methods.csv"
if (file.exists(fgs_file)) {
  fgs_genes <- fread(fgs_file)
  our_genes <- fgs_genes$gene
  cat(sprintf("  Our FGS consensus (3+ methods): %d genes\n", length(our_genes)))

  overlap_results <- data.frame()
  for (sig_name in names(aging_signatures)) {
    lit_genes <- aging_signatures[[sig_name]]
    lit_found <- lit_genes[lit_genes %in% all_genes]
    overlap <- intersect(our_genes, lit_found)

    # Fisher's exact test
    a <- length(overlap)
    b <- length(our_genes) - a
    c <- length(lit_found) - a
    d <- length(all_genes) - a - b - c
    fisher_p <- fisher.test(matrix(c(a, b, c, d), nrow = 2))$p.value

    overlap_results <- rbind(overlap_results, data.frame(
      signature = sig_name,
      n_literature = length(lit_found),
      n_our_fgs = length(our_genes),
      n_overlap = a,
      overlap_genes = paste(overlap, collapse = ", "),
      fisher_p = fisher_p,
      stringsAsFactors = FALSE
    ))
  }
  overlap_results <- overlap_results[order(overlap_results$fisher_p), ]
  write.csv(overlap_results, file.path(out_dir, "fgs_vs_literature_overlap.csv"),
            row.names = FALSE)

  cat("  Overlap results:\n")
  for (i in seq_len(nrow(overlap_results))) {
    r <- overlap_results[i, ]
    cat(sprintf("    %-25s %d/%d overlap (Fisher p=%.2e) %s\n",
                r$signature, r$n_overlap, r$n_literature,
                r$fisher_p,
                ifelse(nchar(r$overlap_genes) > 0,
                       paste0("[", r$overlap_genes, "]"), "")))
  }
}

# ---- Plots ----
cat("\nGenerating plots...\n")

# 1. Signature correlation heatmap
tryCatch({
  # Reshape for heatmap
  hm_data <- cor_results[, c("signature", "rho", "padj")]
  hm_data$sig_label <- ifelse(hm_data$padj < 0.05, "***",
                               ifelse(hm_data$padj < 0.1, "**",
                                      ifelse(hm_data$p_value < 0.05, "*", "")))
  hm_data$signature <- factor(hm_data$signature,
                               levels = hm_data$signature[order(hm_data$rho)])

  p <- ggplot(hm_data, aes(x = "Age", y = signature, fill = rho)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f%s", rho, sig_label)), size = 3) +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                         midpoint = 0, limits = c(-0.6, 0.6)) +
    theme_minimal(base_size = 11) +
    theme(panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA)) +
    labs(title = "Aging Signature Scores vs Age (Spearman rho)",
         x = "", y = "", fill = "rho")
  ggsave(file.path(out_dir, "signature_age_correlation_heatmap.png"),
         p, width = 6, height = 8, dpi = 200, bg = "white")
  cat("  Saved: signature_age_correlation_heatmap.png\n")
}, error = function(e) cat(sprintf("  Heatmap failed: %s\n", e$message)))

# 2. Scatter: top signatures vs age
tryCatch({
  top_sigs <- head(cor_results$signature, 6)
  plots <- list()
  for (sig in top_sigs) {
    sig_col <- paste0("sig_", sig)
    if (sig_col %in% colnames(patient_scores)) {
      p <- ggplot(patient_scores, aes(x = age, y = .data[[sig_col]], color = sex)) +
        geom_point(alpha = 0.7) +
        geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.5) +
        scale_color_manual(values = c("F" = "coral", "M" = "steelblue")) +
        theme_minimal(base_size = 9) +
        theme(panel.background = element_rect(fill = "white", color = NA)) +
        labs(title = sig, x = "Age", y = "Module Score")
      plots[[sig]] <- p
    }
  }
  if (length(plots) > 0) {
    combined <- patchwork::wrap_plots(plots, ncol = 3)
    ggsave(file.path(out_dir, "top_signatures_vs_age.png"),
           combined, width = 15, height = 10, dpi = 200, bg = "white")
    cat("  Saved: top_signatures_vs_age.png\n")
  }
}, error = function(e) cat(sprintf("  Scatter failed: %s\n", e$message)))

# 3. Compartment heatmap
tryCatch({
  if (nrow(compartment_cor) > 0) {
    comp_hm <- compartment_cor %>%
      mutate(sig_star = ifelse(padj < 0.05, "***",
                               ifelse(padj < 0.1, "**",
                                      ifelse(p_value < 0.05, "*", "")))) %>%
      mutate(label = sprintf("%.2f%s", rho, sig_star))

    p <- ggplot(comp_hm, aes(x = compartment, y = signature, fill = rho)) +
      geom_tile(color = "white") +
      geom_text(aes(label = label), size = 2) +
      scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                           midpoint = 0) +
      theme_minimal(base_size = 9) +
      theme(panel.background = element_rect(fill = "white", color = NA),
            plot.background = element_rect(fill = "white", color = NA),
            axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Signature × Compartment × Age",
           x = "", y = "", fill = "rho")
    ggsave(file.path(out_dir, "compartment_signature_heatmap.png"),
           p, width = 12, height = 10, dpi = 200, bg = "white")
    cat("  Saved: compartment_signature_heatmap.png\n")
  }
}, error = function(e) cat(sprintf("  Compartment heatmap failed: %s\n", e$message)))

cat(sprintf("\n=== Aging Signature Validation Complete: %s ===\n", Sys.time()))
cat(sprintf("Output: %s\n", out_dir))
