#!/usr/bin/env Rscript
# HC-only CellChat Interpretation: Age (Young vs Old) and Sex (M vs F)
# Strategy: Load per-sample CellChat objects, aggregate networks per condition,
# then compare aggregated CellChat objects.
#
# The "merged" qs files from run_cellchat_cli.R contain lists of sample-level
# CellChat objects (not single merged objects). We must aggregate manually.
#
# Usage: Rscript scripts/hc/interpret_cellchat_age.R
# Output: /data/user3/sobj/hc_only_v1/cellchat_interpretation/

.libPaths(c(
  "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu",
  "/home/user3/R/x86_64-pc-linux-gnu-library/4.3",
  .libPaths()
))

suppressPackageStartupMessages({
  library(qs)
  library(CellChat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(tidyr)
})

cat("=== HC CellChat Interpretation ===\n")
cat(sprintf("Time: %s\n\n", Sys.time()))

out_dir <- "/data/user3/sobj/hc_only_v1/cellchat_interpretation"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Helper: aggregate sample-level CellChat into condition summary ----
# Approach: average net$count and net$weight across samples, extract pathway stats
aggregate_cellchat_stats <- function(sample_dir, condition_label) {
  cat(sprintf("  Loading samples for %s...\n", condition_label))

  sample_dirs <- list.dirs(sample_dir, recursive = FALSE)

  all_pathways <- list()
  all_counts <- list()
  all_weights <- list()
  cell_types <- NULL
  n_loaded <- 0

  for (sd in sample_dirs) {
    cc_file <- file.path(sd, "cellchat.qs")
    if (!file.exists(cc_file)) next

    cc <- tryCatch(qread(cc_file), error = function(e) NULL)
    if (is.null(cc)) next

    n_loaded <- n_loaded + 1

    # Get cell type names
    if (is.null(cell_types)) {
      cell_types <- levels(cc@idents)
    }

    # Net count/weight (cell_type x cell_type matrices)
    if (!is.null(cc@net$count) && nrow(cc@net$count) > 0) {
      all_counts[[basename(sd)]] <- cc@net$count
      all_weights[[basename(sd)]] <- cc@net$weight
    }

    # Pathway-level stats
    if (!is.null(cc@netP$pathways) && length(cc@netP$pathways) > 0) {
      for (pw in cc@netP$pathways) {
        prob <- tryCatch(cc@netP$prob[,,pw], error = function(e) NULL)
        if (!is.null(prob)) {
          all_pathways[[length(all_pathways) + 1]] <- data.frame(
            sample = basename(sd),
            condition = condition_label,
            pathway = pw,
            total_prob = sum(prob, na.rm = TRUE),
            n_interactions = sum(prob > 0, na.rm = TRUE),
            max_prob = max(prob, na.rm = TRUE),
            stringsAsFactors = FALSE
          )
        }
      }
    }

    rm(cc); gc(verbose = FALSE)
  }

  cat(sprintf("    Loaded %d samples\n", n_loaded))

  pw_df <- if (length(all_pathways) > 0) bind_rows(all_pathways) else NULL

  # Average count/weight matrices
  avg_count <- NULL
  avg_weight <- NULL
  if (length(all_counts) > 0) {
    # Ensure consistent dimensions
    ct_names <- rownames(all_counts[[1]])
    n_ct <- length(ct_names)

    count_arr <- array(0, dim = c(n_ct, n_ct, length(all_counts)))
    weight_arr <- array(0, dim = c(n_ct, n_ct, length(all_counts)))

    for (i in seq_along(all_counts)) {
      # Match cell types
      if (identical(rownames(all_counts[[i]]), ct_names)) {
        count_arr[,,i] <- all_counts[[i]]
        weight_arr[,,i] <- all_weights[[i]]
      }
    }

    avg_count <- apply(count_arr, c(1, 2), mean)
    avg_weight <- apply(weight_arr, c(1, 2), mean)
    rownames(avg_count) <- colnames(avg_count) <- ct_names
    rownames(avg_weight) <- colnames(avg_weight) <- ct_names
  }

  list(
    pathway_df = pw_df,
    avg_count = avg_count,
    avg_weight = avg_weight,
    cell_types = cell_types,
    n_samples = n_loaded
  )
}

# ---- Helper: compare two conditions ----
compare_conditions <- function(stats1, stats2, name1, name2, label, out_path) {
  cat(sprintf("\n--- Comparing %s: %s vs %s ---\n", label, name1, name2))
  dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

  # 1. Compare total interactions and strength
  if (!is.null(stats1$avg_count) && !is.null(stats2$avg_count)) {
    cat("  Comparing aggregate networks...\n")

    # Ensure same cell types
    ct_common <- intersect(rownames(stats1$avg_count), rownames(stats2$avg_count))
    cat(sprintf("    Common cell types: %d\n", length(ct_common)))

    c1 <- stats1$avg_count[ct_common, ct_common]
    c2 <- stats2$avg_count[ct_common, ct_common]
    w1 <- stats1$avg_weight[ct_common, ct_common]
    w2 <- stats2$avg_weight[ct_common, ct_common]

    # Summary
    cat(sprintf("    %s: avg interactions=%.1f, avg weight=%.3f\n",
                name1, sum(c1), sum(w1)))
    cat(sprintf("    %s: avg interactions=%.1f, avg weight=%.3f\n",
                name2, sum(c2), sum(w2)))
    cat(sprintf("    Diff: interactions=%.1f (%.1f%%), weight=%.3f (%.1f%%)\n",
                sum(c2) - sum(c1), (sum(c2) - sum(c1)) / sum(c1) * 100,
                sum(w2) - sum(w1), (sum(w2) - sum(w1)) / sum(w1) * 100))

    # Differential interaction heatmap
    diff_count <- c2 - c1
    diff_weight <- w2 - w1

    tryCatch({
      # Outgoing (row sums) and incoming (col sums) differences
      out_diff <- data.frame(
        cell_type = ct_common,
        outgoing_count_diff = rowSums(diff_count),
        incoming_count_diff = colSums(diff_count),
        outgoing_weight_diff = rowSums(diff_weight),
        incoming_weight_diff = colSums(diff_weight),
        total_count_1 = rowSums(c1) + colSums(c1),
        total_count_2 = rowSums(c2) + colSums(c2),
        total_weight_1 = rowSums(w1) + colSums(w1),
        total_weight_2 = rowSums(w2) + colSums(w2)
      ) %>% arrange(desc(abs(outgoing_weight_diff)))

      write.csv(out_diff, file.path(out_path, sprintf("%s_celltype_interaction_diff.csv", label)),
                row.names = FALSE)

      # Plot: outgoing/incoming difference barplot
      plot_df <- out_diff %>%
        select(cell_type, outgoing_weight_diff, incoming_weight_diff) %>%
        pivot_longer(-cell_type, names_to = "direction", values_to = "diff") %>%
        mutate(direction = gsub("_weight_diff", "", direction),
               cell_type = factor(cell_type, levels = rev(out_diff$cell_type)))

      p <- ggplot(plot_df, aes(x = diff, y = cell_type, fill = direction)) +
        geom_col(position = "dodge", width = 0.7) +
        geom_vline(xintercept = 0, linetype = "dashed") +
        scale_fill_manual(values = c(outgoing = "#E41A1C", incoming = "#377EB8")) +
        labs(title = sprintf("CellChat: %s vs %s (Interaction Strength Diff)", name2, name1),
             x = sprintf("Weight diff (%s - %s)", name2, name1), y = NULL) +
        theme_minimal(base_size = 10) +
        theme(legend.position = "top")

      ggsave(file.path(out_path, sprintf("%s_interaction_diff_barplot.png", label)),
             p, width = 10, height = 8, dpi = 200, bg = "white")
      cat("  Saved: interaction_diff_barplot.png\n")

    }, error = function(e) cat(sprintf("  Barplot failed: %s\n", e$message)))

    # Differential heatmap
    tryCatch({
      diff_long <- expand.grid(
        sender = ct_common, receiver = ct_common, stringsAsFactors = FALSE
      )
      diff_long$diff_weight <- as.vector(diff_weight)

      # Only show top interactions
      top_senders <- head(out_diff$cell_type[order(-abs(out_diff$outgoing_weight_diff))], 15)
      top_receivers <- head(out_diff$cell_type[order(-abs(out_diff$incoming_weight_diff))], 15)
      top_ct <- unique(c(top_senders, top_receivers))

      diff_sub <- diff_long %>%
        filter(sender %in% top_ct, receiver %in% top_ct)

      p2 <- ggplot(diff_sub, aes(x = receiver, y = sender, fill = diff_weight)) +
        geom_tile() +
        scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                             name = sprintf("Weight\n(%s-%s)", name2, name1)) +
        labs(title = sprintf("CellChat Diff Heatmap: %s", label)) +
        theme_minimal(base_size = 9) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

      ggsave(file.path(out_path, sprintf("%s_diff_heatmap.png", label)),
             p2, width = 12, height = 10, dpi = 200, bg = "white")
      cat("  Saved: diff_heatmap.png\n")

    }, error = function(e) cat(sprintf("  Heatmap failed: %s\n", e$message)))
  }

  # 2. Pathway-level comparison
  if (!is.null(stats1$pathway_df) && !is.null(stats2$pathway_df)) {
    cat("  Comparing pathway-level statistics...\n")

    # Aggregate per pathway per condition (mean across samples)
    pw_all <- bind_rows(stats1$pathway_df, stats2$pathway_df)
    write.csv(pw_all, file.path(out_path, sprintf("%s_all_pathway_stats.csv", label)),
              row.names = FALSE)

    pw_summary <- pw_all %>%
      group_by(condition, pathway) %>%
      summarise(
        mean_prob = mean(total_prob, na.rm = TRUE),
        sd_prob = sd(total_prob, na.rm = TRUE),
        mean_n_int = mean(n_interactions, na.rm = TRUE),
        n_samples = n(),
        .groups = "drop"
      )

    # Wide format for comparison
    pw_wide <- pw_summary %>%
      select(condition, pathway, mean_prob) %>%
      pivot_wider(names_from = condition, values_from = mean_prob, values_fill = 0)

    if (name1 %in% colnames(pw_wide) && name2 %in% colnames(pw_wide)) {
      pw_wide$diff <- pw_wide[[name2]] - pw_wide[[name1]]
      pw_wide$log2fc <- log2((pw_wide[[name2]] + 0.001) / (pw_wide[[name1]] + 0.001))

      # Wilcoxon test per pathway
      pw_pvals <- pw_all %>%
        group_by(pathway) %>%
        summarise(
          p_value = tryCatch(
            wilcox.test(total_prob[condition == name1], total_prob[condition == name2])$p.value,
            error = function(e) NA
          ),
          .groups = "drop"
        )

      pw_wide <- left_join(pw_wide, pw_pvals, by = "pathway")
      pw_wide$p_adj <- p.adjust(pw_wide$p_value, method = "BH")
      pw_wide <- pw_wide[order(-abs(pw_wide$diff)), ]

      write.csv(pw_wide, file.path(out_path, sprintf("%s_differential_pathways.csv", label)),
                row.names = FALSE)

      cat("  Top differential pathways:\n")
      top_pw <- head(pw_wide, 15)
      for (j in seq_len(nrow(top_pw))) {
        sig <- ifelse(!is.na(top_pw$p_adj[j]) && top_pw$p_adj[j] < 0.05, "*", "")
        cat(sprintf("    %-20s %s=%.3f %s=%.3f diff=%.3f log2fc=%.2f p_adj=%.3f %s\n",
                    top_pw$pathway[j],
                    name1, top_pw[[name1]][j],
                    name2, top_pw[[name2]][j],
                    top_pw$diff[j], top_pw$log2fc[j],
                    ifelse(is.na(top_pw$p_adj[j]), 1, top_pw$p_adj[j]),
                    sig))
      }

      # Pathway bar plot
      tryCatch({
        top20 <- head(pw_wide, 20)
        top20$pathway <- factor(top20$pathway, levels = rev(top20$pathway))
        top20$sig <- ifelse(!is.na(top20$p_adj) & top20$p_adj < 0.05, "Sig", "NS")

        p <- ggplot(top20, aes(x = diff, y = pathway, fill = sig)) +
          geom_col() +
          geom_vline(xintercept = 0, linetype = "dashed") +
          scale_fill_manual(values = c(Sig = "#E41A1C", NS = "grey60")) +
          labs(title = sprintf("Pathway Diff: %s vs %s", name2, name1),
               x = sprintf("Mean prob diff (%s - %s)", name2, name1), y = NULL) +
          theme_minimal(base_size = 10)

        ggsave(file.path(out_path, sprintf("%s_pathway_diff_barplot.png", label)),
               p, width = 8, height = 8, dpi = 200, bg = "white")
        cat("  Saved: pathway_diff_barplot.png\n")
      }, error = function(e) cat(sprintf("  Pathway barplot failed: %s\n", e$message)))

      # Scatter plot: condition 1 vs 2 pathway strengths
      tryCatch({
        p2 <- ggplot(pw_wide, aes_string(x = name1, y = name2)) +
          geom_point(aes(color = ifelse(!is.na(p_adj) & p_adj < 0.05, "Sig", "NS")), size = 2) +
          geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
          geom_text(data = pw_wide %>% filter(!is.na(p_adj) & p_adj < 0.05),
                    aes(label = pathway), size = 2.5, hjust = -0.1, vjust = -0.3) +
          scale_color_manual(values = c(Sig = "red", NS = "grey50"), name = "FDR<0.05") +
          labs(title = sprintf("Pathway Strength: %s vs %s", name1, name2),
               x = sprintf("%s (mean prob)", name1),
               y = sprintf("%s (mean prob)", name2)) +
          theme_minimal(base_size = 10) +
          coord_equal()

        ggsave(file.path(out_path, sprintf("%s_pathway_scatter.png", label)),
               p2, width = 8, height = 8, dpi = 200, bg = "white")
        cat("  Saved: pathway_scatter.png\n")
      }, error = function(e) cat(sprintf("  Scatter failed: %s\n", e$message)))
    }
  }
}

# ---- Age comparison ----
cat("=== Age Comparison (Young vs Old) ===\n")
age_sample_dir <- "/data/user3/sobj/hc_only_v1/cellchat/age_v1/samples"

# Get age group mapping
meta <- read.csv("/data/user3/sobj/hc_only_v1/meld/metadata.csv", row.names = 1)
patient_age <- meta %>%
  distinct(name, age, age_group) %>%
  filter(!is.na(age_group))

young_patients <- patient_age$name[patient_age$age_group == "Young"]
old_patients <- patient_age$name[patient_age$age_group == "Old"]
cat(sprintf("  Young: %d, Old: %d patients\n", length(young_patients), length(old_patients)))

# Check which patients have CellChat
all_cc_samples <- basename(list.dirs(age_sample_dir, recursive = FALSE))
young_cc <- intersect(young_patients, all_cc_samples)
old_cc <- intersect(old_patients, all_cc_samples)
cat(sprintf("  With CellChat: Young=%d, Old=%d\n", length(young_cc), length(old_cc)))

# Create temporary dirs with symlinks for age groups
tmp_young <- tempfile("young_"); dir.create(tmp_young)
tmp_old <- tempfile("old_"); dir.create(tmp_old)
for (p in young_cc) file.symlink(file.path(age_sample_dir, p), file.path(tmp_young, p))
for (p in old_cc) file.symlink(file.path(age_sample_dir, p), file.path(tmp_old, p))

stats_young <- aggregate_cellchat_stats(tmp_young, "Young")
stats_old <- aggregate_cellchat_stats(tmp_old, "Old")

age_out <- file.path(out_dir, "age")
compare_conditions(stats_young, stats_old, "Young", "Old", "age", age_out)

# Clean up
unlink(tmp_young, recursive = TRUE)
unlink(tmp_old, recursive = TRUE)

# ---- Sex comparison ----
cat("\n=== Sex Comparison (Female vs Male) ===\n")
sex_sample_dir <- "/data/user3/sobj/hc_only_v1/cellchat/sex_v1/samples"

# Check for sex sample dir or use age samples with sex annotation
if (dir.exists(sex_sample_dir)) {
  female_patients <- patient_age$name[meta %>%
    distinct(name, sex) %>%
    filter(sex == "F") %>%
    pull(name) %>%
    unique()]
  # Actually get from meta
  sex_map <- meta %>% distinct(name, sex) %>% filter(!is.na(sex))
  female_names <- sex_map$name[sex_map$sex == "F"]
  male_names <- sex_map$name[sex_map$sex == "M"]

  all_sex_samples <- basename(list.dirs(sex_sample_dir, recursive = FALSE))
  female_cc <- intersect(female_names, all_sex_samples)
  male_cc <- intersect(male_names, all_sex_samples)

  cat(sprintf("  Female: %d, Male: %d with CellChat\n", length(female_cc), length(male_cc)))

  tmp_f <- tempfile("female_"); dir.create(tmp_f)
  tmp_m <- tempfile("male_"); dir.create(tmp_m)
  for (p in female_cc) file.symlink(file.path(sex_sample_dir, p), file.path(tmp_f, p))
  for (p in male_cc) file.symlink(file.path(sex_sample_dir, p), file.path(tmp_m, p))

  stats_f <- aggregate_cellchat_stats(tmp_f, "Female")
  stats_m <- aggregate_cellchat_stats(tmp_m, "Male")

  sex_out <- file.path(out_dir, "sex")
  compare_conditions(stats_f, stats_m, "Female", "Male", "sex", sex_out)

  unlink(tmp_f, recursive = TRUE)
  unlink(tmp_m, recursive = TRUE)
} else {
  cat("  Sex CellChat sample directory not found, using age samples with sex metadata...\n")

  sex_map <- meta %>% distinct(name, sex) %>% filter(!is.na(sex))
  female_names <- sex_map$name[sex_map$sex == "F"]
  male_names <- sex_map$name[sex_map$sex == "M"]

  female_cc <- intersect(female_names, all_cc_samples)
  male_cc <- intersect(male_names, all_cc_samples)

  cat(sprintf("  Female: %d, Male: %d with CellChat\n", length(female_cc), length(male_cc)))

  tmp_f <- tempfile("female_"); dir.create(tmp_f)
  tmp_m <- tempfile("male_"); dir.create(tmp_m)
  for (p in female_cc) file.symlink(file.path(age_sample_dir, p), file.path(tmp_f, p))
  for (p in male_cc) file.symlink(file.path(age_sample_dir, p), file.path(tmp_m, p))

  stats_f <- aggregate_cellchat_stats(tmp_f, "Female")
  stats_m <- aggregate_cellchat_stats(tmp_m, "Male")

  sex_out <- file.path(out_dir, "sex")
  compare_conditions(stats_f, stats_m, "Female", "Male", "sex", sex_out)

  unlink(tmp_f, recursive = TRUE)
  unlink(tmp_m, recursive = TRUE)
}

cat(sprintf("\n=== CellChat Interpretation Complete: %s ===\n", Sys.time()))
cat(sprintf("Output: %s\n", out_dir))
