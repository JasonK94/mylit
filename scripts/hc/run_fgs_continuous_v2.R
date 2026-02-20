#!/usr/bin/env Rscript
# ==============================================================================
# Continuous Variable FGS Pipeline v2 — Top-200, 6 Methods
# ==============================================================================
# Extension of v1: top-200 per method + 2 additional methods (GAM, Elastic Net)
#
# Methods:
#   1. Spearman correlation (pseudobulk per patient)
#   2. LASSO regression (glmnet, family="gaussian")
#   3. Random Forest regression (ranger, importance)
#   4. limma-voom (pseudobulk design ~ age + sex)
#   5. GAM (mgcv: expr ~ s(age) + sex, per gene)
#   6. Elastic Net (alpha=0.5, otherwise same as LASSO)
#
# Reuses pseudobulk from v1 run.
# ==============================================================================

cat("
================================================================================
  Continuous FGS Pipeline v2 — CD4 Treg × Age (Top-200, 6 methods)
  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "
================================================================================
\n")

.libPaths(c(
  "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu",
  "/home/user3/R/x86_64-pc-linux-gnu-library/4.3",
  .libPaths()
))

suppressPackageStartupMessages({
  library(qs)
  library(Matrix)
  library(glmnet)
  library(ranger)
  library(limma)
  library(edgeR)
  library(mgcv)
  library(ggplot2)
  library(dplyr)
  library(pheatmap)
})

set.seed(42)

N_TOP <- 200
OUTPUT_DIR <- "/data/user3/sobj/hc_only_v1/fgs_continuous_treg_v2"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. Load cached pseudobulk
# ==============================================================================
message("[1/6] Loading cached pseudobulk data...")
pb_data <- readRDS("/data/user3/sobj/hc_only_v1/fgs_continuous_treg/pseudobulk_data.rds")
pb_counts  <- pb_data$pb_counts
pb_logcpm  <- pb_data$pb_logcpm
patient_meta <- pb_data$patient_meta
dge <- pb_data$dge

n_genes <- nrow(pb_logcpm)
n_patients <- ncol(pb_logcpm)
message(sprintf("  Pseudobulk: %d genes × %d patients", n_genes, n_patients))
message(sprintf("  Age: %.0f–%.0f (median %.0f)", min(patient_meta$age),
                max(patient_meta$age), median(patient_meta$age)))

age_vec <- patient_meta$age
names(age_vec) <- rownames(patient_meta)
sex_vec <- factor(patient_meta$sex)
names(sex_vec) <- rownames(patient_meta)

all_genes <- rownames(pb_logcpm)

# ==============================================================================
# 2. Method 1: Spearman (recompute for full ranking — fast)
# ==============================================================================
message("\n[2/6] Method 1/6: Spearman Correlation...")
t_m <- Sys.time()

spearman_res <- data.frame(
  gene = all_genes, rho = NA_real_, pvalue = NA_real_, stringsAsFactors = FALSE
)
for (i in seq_len(n_genes)) {
  ct <- cor.test(pb_logcpm[i, ], age_vec, method = "spearman", exact = FALSE)
  spearman_res$rho[i] <- ct$estimate
  spearman_res$pvalue[i] <- ct$p.value
}
spearman_res$padj <- p.adjust(spearman_res$pvalue, method = "BH")
spearman_res$abs_rho <- abs(spearman_res$rho)
spearman_res$rank_p <- rank(spearman_res$pvalue)
spearman_res$rank_rho <- rank(-spearman_res$abs_rho)
spearman_res$combined_rank <- spearman_res$rank_p + spearman_res$rank_rho
spearman_res <- spearman_res[order(spearman_res$combined_rank), ]
rownames(spearman_res) <- spearman_res$gene
spearman_top <- head(spearman_res$gene, N_TOP)

message(sprintf("  Top: %s (rho=%.3f), sig padj<0.05: %d (%.1fs)",
                spearman_res$gene[1], spearman_res$rho[1],
                sum(spearman_res$padj < 0.05, na.rm = TRUE),
                as.numeric(Sys.time() - t_m, units = "secs")))

# ==============================================================================
# 3. Method 2: LASSO (alpha=1)
# ==============================================================================
message("\n[3/6] Method 2/6: LASSO (alpha=1)...")
t_m <- Sys.time()

X_base <- t(pb_logcpm)
sex_dummy <- model.matrix(~ sex_vec - 1)[, -1, drop = FALSE]
colnames(sex_dummy) <- "sexM"
X_with_cov <- cbind(X_base, sex_dummy)
penalty_factor <- c(rep(1, ncol(X_base)), rep(0, ncol(sex_dummy)))

cv_lasso <- cv.glmnet(x = X_with_cov, y = age_vec, family = "gaussian",
                       alpha = 1, penalty.factor = penalty_factor,
                       nfolds = min(10, n_patients), type.measure = "mse")

coef_lasso <- as.numeric(coef(cv_lasso, s = "lambda.1se")[2:(n_genes + 1)])
names(coef_lasso) <- all_genes
lasso_order <- names(sort(abs(coef_lasso), decreasing = TRUE))
lasso_top <- head(lasso_order, N_TOP)

pred_lasso <- predict(cv_lasso, newx = X_with_cov, s = "lambda.1se")
r2_lasso <- cor(pred_lasso, age_vec)^2
n_nonzero_lasso <- sum(coef_lasso != 0)
message(sprintf("  Non-zero: %d, R2=%.3f (%.1fs)", n_nonzero_lasso, r2_lasso,
                as.numeric(Sys.time() - t_m, units = "secs")))

# ==============================================================================
# 4. Method 3: Random Forest
# ==============================================================================
message("\n[4/6] Method 3/6: Random Forest...")
t_m <- Sys.time()

rf_x <- as.data.frame(X_base)
rf_x$sex <- as.numeric(sex_vec == "M")

rf_model <- ranger(x = rf_x, y = age_vec, num.trees = 1000,
                    importance = "permutation", seed = 42)

rf_imp <- ranger::importance(rf_model)
rf_gene_imp <- rf_imp[setdiff(names(rf_imp), "sex")]

# Map mangled names back
mangled_to_orig <- setNames(all_genes, make.names(all_genes))
mapped <- mangled_to_orig[names(rf_gene_imp)]
valid <- !is.na(mapped)
rf_gene_imp <- rf_gene_imp[valid]
names(rf_gene_imp) <- mapped[valid]

rf_order <- names(sort(rf_gene_imp, decreasing = TRUE))
rf_top <- head(rf_order, N_TOP)

message(sprintf("  OOB R2=%.3f, top: %s (%.1fs)", rf_model$r.squared, rf_order[1],
                as.numeric(Sys.time() - t_m, units = "secs")))

# ==============================================================================
# 5. Method 4: limma-voom
# ==============================================================================
message("\n[5/6] Method 4/6: limma-voom...")
t_m <- Sys.time()

design_limma <- model.matrix(~ age + sex, data = patient_meta)
v <- voom(dge, design_limma, plot = FALSE)
fit <- lmFit(v, design_limma)
fit <- eBayes(fit)
tt_age <- topTable(fit, coef = "age", number = Inf, sort.by = "p")
tt_age$gene <- rownames(tt_age)
limma_top <- head(tt_age$gene, N_TOP)

message(sprintf("  Sig (adj.P<0.05): %d, top: %s (%.1fs)",
                sum(tt_age$adj.P.Val < 0.05), limma_top[1],
                as.numeric(Sys.time() - t_m, units = "secs")))

# ==============================================================================
# 6. Method 5: GAM (mgcv) — nonlinear age effect
# ==============================================================================
message("\n  Method 5/6: GAM (mgcv, s(age) + sex)...")
t_m <- Sys.time()

gam_results <- data.frame(gene = all_genes, deviance_explained = NA_real_,
                           p_smooth = NA_real_, edf = NA_real_,
                           stringsAsFactors = FALSE)

for (i in seq_len(n_genes)) {
  g <- all_genes[i]
  expr <- pb_logcpm[g, ]
  gam_df <- data.frame(expr = expr, age = age_vec, sex = sex_vec)
  tryCatch({
    m <- gam(expr ~ s(age, k = 5) + sex, data = gam_df, method = "REML")
    s_table <- summary(m)$s.table
    gam_results$deviance_explained[i] <- summary(m)$dev.expl
    gam_results$p_smooth[i] <- s_table[1, "p-value"]
    gam_results$edf[i] <- s_table[1, "edf"]
  }, error = function(e) {})

  if (i %% 5000 == 0) message(sprintf("    ... %d/%d genes", i, n_genes))
}

gam_results$padj <- p.adjust(gam_results$p_smooth, method = "BH")
gam_results <- gam_results[order(gam_results$p_smooth), ]
rownames(gam_results) <- gam_results$gene
gam_top <- head(gam_results$gene[!is.na(gam_results$p_smooth)], N_TOP)

message(sprintf("  Sig (padj<0.05): %d, top: %s, deviance=%.3f (%.1fs)",
                sum(gam_results$padj < 0.05, na.rm = TRUE),
                gam_top[1], gam_results$deviance_explained[1],
                as.numeric(Sys.time() - t_m, units = "secs")))

# ==============================================================================
# 7. Method 6: Elastic Net (alpha=0.5)
# ==============================================================================
message("\n  Method 6/6: Elastic Net (alpha=0.5)...")
t_m <- Sys.time()

cv_enet <- cv.glmnet(x = X_with_cov, y = age_vec, family = "gaussian",
                      alpha = 0.5, penalty.factor = penalty_factor,
                      nfolds = min(10, n_patients), type.measure = "mse")

coef_enet <- as.numeric(coef(cv_enet, s = "lambda.1se")[2:(n_genes + 1)])
names(coef_enet) <- all_genes
enet_order <- names(sort(abs(coef_enet), decreasing = TRUE))
enet_top <- head(enet_order, N_TOP)

pred_enet <- predict(cv_enet, newx = X_with_cov, s = "lambda.1se")
r2_enet <- cor(pred_enet, age_vec)^2
n_nonzero_enet <- sum(coef_enet != 0)
message(sprintf("  Non-zero: %d, R2=%.3f (%.1fs)", n_nonzero_enet, r2_enet,
                as.numeric(Sys.time() - t_m, units = "secs")))

# ==============================================================================
# 8. Consensus
# ==============================================================================
message("\n[6/6] Computing consensus (top-%d across 6 methods)...", N_TOP)

method_names <- c("spearman", "lasso", "random_forest", "limma_voom", "gam", "elastic_net")
top_lists <- list(
  spearman      = spearman_top,
  lasso         = lasso_top,
  random_forest = rf_top,
  limma_voom    = limma_top,
  gam           = gam_top,
  elastic_net   = enet_top
)

# Build rank matrix
rank_matrix <- matrix(NA, nrow = n_genes, ncol = 6,
                       dimnames = list(all_genes, method_names))

# Spearman
sp_ranks <- setNames(seq_len(nrow(spearman_res)), spearman_res$gene)
rank_matrix[names(sp_ranks), "spearman"] <- sp_ranks

# LASSO
la_ranks <- setNames(seq_along(lasso_order), lasso_order)
rank_matrix[names(la_ranks), "lasso"] <- la_ranks

# Random Forest
rf_ranks <- setNames(seq_along(rf_order), rf_order)
matched_rf <- intersect(names(rf_ranks), all_genes)
rank_matrix[matched_rf, "random_forest"] <- rf_ranks[matched_rf]

# limma-voom
lm_ranks <- setNames(seq_along(tt_age$gene), tt_age$gene)
rank_matrix[names(lm_ranks), "limma_voom"] <- lm_ranks

# GAM
gam_ordered <- gam_results$gene[!is.na(gam_results$p_smooth)]
gam_ranks <- setNames(seq_along(gam_ordered), gam_ordered)
rank_matrix[names(gam_ranks), "gam"] <- gam_ranks

# Elastic Net
en_ranks <- setNames(seq_along(enet_order), enet_order)
rank_matrix[names(en_ranks), "elastic_net"] <- en_ranks

# Consensus metrics
mean_ranks <- rowMeans(rank_matrix, na.rm = TRUE)
n_methods_in_top <- rowSums(rank_matrix <= N_TOP, na.rm = TRUE)

consensus_df <- data.frame(
  gene = all_genes,
  mean_rank = mean_ranks,
  n_methods_top200 = n_methods_in_top,
  rank_spearman = rank_matrix[, "spearman"],
  rank_lasso = rank_matrix[, "lasso"],
  rank_rf = rank_matrix[, "random_forest"],
  rank_limma = rank_matrix[, "limma_voom"],
  rank_gam = rank_matrix[, "gam"],
  rank_enet = rank_matrix[, "elastic_net"],
  stringsAsFactors = FALSE
)

# Add effect sizes
consensus_df$spearman_rho <- spearman_res[consensus_df$gene, "rho"]
limma_lfc <- setNames(tt_age$logFC, tt_age$gene)
consensus_df$limma_logFC <- limma_lfc[consensus_df$gene]
gam_de <- setNames(gam_results$deviance_explained, gam_results$gene)
consensus_df$gam_deviance <- gam_de[consensus_df$gene]

consensus_df <- consensus_df[order(-consensus_df$n_methods_top200, consensus_df$mean_rank), ]
rownames(consensus_df) <- consensus_df$gene

# Summary
for (threshold in c(2, 3, 4, 5, 6)) {
  n <- sum(consensus_df$n_methods_top200 >= threshold)
  message(sprintf("  Genes in top-%d of >= %d methods: %d", N_TOP, threshold, n))
}

# Top consensus genes
cat("\n  Top 30 consensus genes (top-200 in >= 3 methods):\n")
top30 <- head(consensus_df[consensus_df$n_methods_top200 >= 3, ], 30)
cat(sprintf("  %-3s %-18s %6s %8s %8s %8s %8s %8s %8s\n",
            "#", "Gene", "#Meth", "Spear", "LASSO", "RF", "limma", "GAM", "ENet"))
cat(paste(rep("-", 95), collapse=""), "\n")
for (i in seq_len(nrow(top30))) {
  r <- top30[i, ]
  cat(sprintf("  %-3d %-18s %6d %8.0f %8.0f %8.0f %8.0f %8.0f %8.0f\n",
              i, r$gene, r$n_methods_top200,
              r$rank_spearman, r$rank_lasso, r$rank_rf,
              r$rank_limma, r$rank_gam, r$rank_enet))
}

# ==============================================================================
# 9. Enrichment test (Fisher's exact / Poisson)
# ==============================================================================
message("\n  === Enrichment Test ===")
p_single <- N_TOP / n_genes
cat(sprintf("  P(gene in top-%d of 1 method) = %.4f\n", N_TOP, p_single))

for (k in 2:6) {
  p_k <- sum(sapply(k:6, function(j) choose(6, j) * p_single^j * (1 - p_single)^(6 - j)))
  expected <- n_genes * p_k
  observed <- sum(consensus_df$n_methods_top200 >= k)
  fold <- ifelse(expected > 0, observed / expected, Inf)
  pval <- ppois(observed - 1, lambda = expected, lower.tail = FALSE)
  cat(sprintf("  >= %d methods: expected=%.1f, observed=%d, fold=%.1fx, p=%.2e\n",
              k, expected, observed, fold, pval))
}

# ==============================================================================
# 10. Save outputs
# ==============================================================================
message("\n  Saving outputs...")

write.csv(consensus_df, file.path(OUTPUT_DIR, "consensus_all_genes.csv"), row.names = FALSE)

for (threshold in c(2, 3, 4)) {
  sub <- consensus_df[consensus_df$n_methods_top200 >= threshold, ]
  fname <- sprintf("consensus_genes_%dplus_methods.csv", threshold)
  write.csv(sub, file.path(OUTPUT_DIR, fname), row.names = FALSE)
  message(sprintf("  Saved: %s (%d genes)", fname, nrow(sub)))
}

# Per-method top lists
for (mn in names(top_lists)) {
  write.csv(data.frame(rank = seq_along(top_lists[[mn]]), gene = top_lists[[mn]]),
            file.path(OUTPUT_DIR, sprintf("top%d_%s.csv", N_TOP, mn)),
            row.names = FALSE)
}

# Pairwise overlap matrix
overlap_mat <- matrix(0, 6, 6, dimnames = list(method_names, method_names))
for (i in seq_along(method_names)) {
  for (j in seq_along(method_names)) {
    overlap_mat[i, j] <- length(intersect(top_lists[[method_names[i]]],
                                           top_lists[[method_names[j]]]))
  }
}
write.csv(overlap_mat, file.path(OUTPUT_DIR, "pairwise_overlap_top200.csv"))
cat("\n  Pairwise overlap (top-200):\n")
print(overlap_mat)

# Full results
saveRDS(list(
  consensus_df = consensus_df,
  rank_matrix = rank_matrix,
  top_lists = top_lists,
  method_stats = list(
    spearman = list(n_sig_padj05 = sum(spearman_res$padj < 0.05, na.rm = TRUE)),
    lasso = list(n_nonzero = n_nonzero_lasso, R2 = r2_lasso),
    rf = list(OOB_R2 = rf_model$r.squared),
    limma = list(n_sig = sum(tt_age$adj.P.Val < 0.05)),
    gam = list(n_sig = sum(gam_results$padj < 0.05, na.rm = TRUE)),
    enet = list(n_nonzero = n_nonzero_enet, R2 = r2_enet)
  )
), file.path(OUTPUT_DIR, "fgs_continuous_v2_results.rds"))

# ==============================================================================
# 11. Plots
# ==============================================================================
message("\n  Generating plots...")

# --- Upset-style bar: genes by n_methods ---
n_meth_counts <- table(factor(consensus_df$n_methods_top200[consensus_df$n_methods_top200 >= 1],
                               levels = 1:6))
n_meth_df <- data.frame(n_methods = as.integer(names(n_meth_counts)),
                          count = as.integer(n_meth_counts))

p_bar <- ggplot(n_meth_df, aes(x = factor(n_methods), y = count)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
  geom_text(aes(label = count), vjust = -0.3) +
  theme_minimal(base_size = 12) +
  labs(title = sprintf("FGS Consensus: genes in top-%d of N methods (6 methods, %d genes)", N_TOP, n_genes),
       x = "Number of methods in top-200", y = "Number of genes")
ggsave(file.path(OUTPUT_DIR, "barplot_n_methods.png"), p_bar, width = 8, height = 5, dpi = 200)

# --- Heatmap: top consensus genes ---
heatmap_genes <- head(consensus_df$gene[consensus_df$n_methods_top200 >= 3], 60)

if (length(heatmap_genes) >= 5) {
  age_order <- order(patient_meta$age)
  mat <- pb_logcpm[heatmap_genes, age_order]
  mat_z <- t(scale(t(mat)))
  mat_z[mat_z > 3] <- 3; mat_z[mat_z < -3] <- -3

  anno_col <- data.frame(
    Age = patient_meta$age[age_order],
    Sex = patient_meta$sex[age_order],
    row.names = rownames(patient_meta)[age_order]
  )
  anno_row <- data.frame(
    Direction = ifelse(consensus_df[heatmap_genes, "spearman_rho"] > 0, "Up", "Down"),
    N_methods = as.character(consensus_df[heatmap_genes, "n_methods_top200"]),
    row.names = heatmap_genes
  )
  ann_colors <- list(
    Age = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
    Sex = c(F = "#D6604D", M = "#4393C3"),
    Direction = c(Up = "#D73027", Down = "#4575B4"),
    N_methods = c("3" = "#FEE08B", "4" = "#FDAE61", "5" = "#F46D43", "6" = "#D73027")
  )

  png(file.path(OUTPUT_DIR, "heatmap_consensus_top60.png"),
      width = 1600, height = max(800, length(heatmap_genes) * 14 + 200), res = 150)
  pheatmap(mat_z, cluster_cols = FALSE, cluster_rows = TRUE,
           annotation_col = anno_col, annotation_row = anno_row,
           annotation_colors = ann_colors,
           show_colnames = FALSE, fontsize_row = 7,
           main = sprintf("Top %d Consensus Age-Associated Genes (CD4 Treg, 6 methods, top-%d)",
                          length(heatmap_genes), N_TOP),
           color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
           breaks = seq(-3, 3, length.out = 101), border_color = NA)
  dev.off()
  message("  Saved: heatmap_consensus_top60.png")
}

# --- Scatter: top 12 consensus genes ---
scatter_genes <- head(consensus_df$gene[consensus_df$n_methods_top200 >= 3], 12)
if (length(scatter_genes) >= 1) {
  plot_data <- do.call(rbind, lapply(scatter_genes, function(g) {
    data.frame(gene = g, expression = pb_logcpm[g, ],
               age = patient_meta$age, sex = patient_meta$sex)
  }))
  # Labels
  sp_info <- spearman_res[scatter_genes, ]
  label_map <- setNames(
    sprintf("%s\nrho=%+.3f p=%.1e", sp_info$gene, sp_info$rho, sp_info$pvalue),
    sp_info$gene
  )
  plot_data$gene_label <- factor(label_map[plot_data$gene], levels = label_map)

  p_sc <- ggplot(plot_data, aes(x = age, y = expression, color = sex)) +
    geom_point(size = 1.5, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.7) +
    facet_wrap(~ gene_label, scales = "free_y", ncol = 4) +
    scale_color_manual(values = c(F = "#D6604D", M = "#4393C3")) +
    theme_minimal(base_size = 10) +
    labs(title = "Top Consensus Genes × Age (CD4 Treg, 6 methods, top-200)",
         x = "Age (years)", y = "logCPM", color = "Sex") +
    theme(strip.text = element_text(face = "bold", size = 8), legend.position = "bottom")
  ggsave(file.path(OUTPUT_DIR, "scatter_top12_consensus.png"), p_sc,
         width = 16, height = 3.5 * ceiling(length(scatter_genes) / 4), dpi = 200)
  message("  Saved: scatter_top12_consensus.png")
}

cat(sprintf("\n\nDone! Time: %s\n", Sys.time()))
cat(sprintf("Output: %s\n", OUTPUT_DIR))
