#!/usr/bin/env Rscript
# Generate Comparative Circos Plots - Final (Production)
# - Solved Duplication & Interaction Count Instability
# - Perfect Inner Hint Track
# - Global Layout & Metadata
# - Deterministic Results

suppressPackageStartupMessages({
    library(dplyr)
    library(circlize)
    library(RColorBrewer)
    library(qs)
    library(optparse)
})

# Argument Parsing
option_list <- list(
    make_option(c("-f", "--file"),
        type = "character", default = "/data/user3/sobj/cci/mnn/multinichenet_IS2_IS3_permissive/multinichenet_results.qs",
        help = "Path to MultiNicheNet results file (.qs) [default= %default]", metavar = "character"
    ),
    make_option(c("-o", "--outdir"),
        type = "character", default = "/data/user3/sobj/cci/mnn/multinichenet_IS2_IS3_permissive",
        help = "Output directory [default= %default]", metavar = "character"
    ),
    make_option(c("-n", "--top_n"),
        type = "integer", default = 50,
        help = "Number of top interactions to plot [default= %default]", metavar = "integer"
    ),
    make_option(c("-s", "--sort_by"),
        type = "character", default = "score",
        help = "Sorting criterion: 'score' (default), 'p_val_adj', 'logFC', 'activity'", metavar = "character"
    ),
    make_option(c("--p_adj_cutoff"),
        type = "numeric", default = 0.05,
        help = "Adjusted p-value cutoff. Set to 1.0 to disable. [default= %default]", metavar = "numeric"
    ),
    make_option(c("--max_per_sender"),
        type = "integer", default = 0,
        help = "Maximum interactions per sender cell type. 0 to disable. Use this to balance dominant cell types. [default= %default]", metavar = "integer"
    ),
    make_option(c("--max_per_receiver"),
        type = "integer", default = 0,
        help = "Maximum interactions per receiver cell type. 0 to disable. Use this to balance dominant cell types. [default= %default]", metavar = "integer"
    ),
    make_option(c("--logfc_threshold"),
        type = "numeric", default = 0.05,
        help = "LogFC threshold for filtering significant senders [default= %default]", metavar = "numeric"
    ),
    make_option(c("--p_val_threshold"),
        type = "numeric", default = 0.10,
        help = "P-value threshold for filtering significant senders [default= %default]", metavar = "numeric"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

cat("════════════════════════════════════════════════════════════\n")
cat("Generating Comparative Circos Plot\n")
cat("════════════════════════════════════════════════════════════\n")
cat("Options:\n")
cat("  File:         ", opt$file, "\n")
cat("  Top N:        ", opt$top_n, "\n")
cat("  Sort By:      ", opt$sort_by, "\n")
cat("  P-adj Cutoff: ", opt$p_adj_cutoff, "\n")
cat("  Max/Sender:   ", opt$max_per_sender, "\n")
cat("  Max/Receiver: ", opt$max_per_receiver, "\n")
cat("  LogFC Thresh: ", opt$logfc_threshold, "\n")
cat("  P-val Thresh: ", opt$p_val_threshold, "\n")
cat("════════════════════════════════════════════════════════════\n\n")

# 0. Reproducibility
set.seed(1234)

# 1. Load Data
results_path <- opt$file
if (!file.exists(results_path)) stop("Results file not found: ", results_path)
results <- qs::qread(results_path)

meta_info <- "Covariates: Unknown"
if (!is.null(results$.__metadata__) && !is.null(results$.__metadata__$params$covariates)) {
    covs <- results$.__metadata__$params$covariates
    if (length(covs) > 0 && covs != "NA") {
        meta_info <- paste0("Covariates: ", paste(covs, collapse = ", "))
    } else {
        meta_info <- "Covariates: None"
    }
}

# 2. Process Data
la_tbl <- results$prioritization_tables$ligand_activities_target_de_tbl
lr_network <- readRDS("/data/user3/git_repo/human/lr_network_human_21122021.rds") %>%
    distinct(from, to) %>%
    rename(ligand = from, receptor = to)

celltype_de <- results$celltype_de

sig_senders <- celltype_de %>%
    filter(abs(logFC) > opt$logfc_threshold, p_val < opt$p_val_threshold) %>%
    select(gene, cluster_id, logFC) %>%
    rename(ligand = gene, sender = cluster_id, sender_logFC = logFC)

# Base Table Construction
tbl_base <- la_tbl %>%
    mutate(score = activity * abs(logFC)) %>%
    # Ensure logFC is absolute for sorting if needed
    mutate(abs_logFC = abs(logFC)) %>%
    inner_join(lr_network, by = "ligand") %>%
    inner_join(sig_senders, by = "ligand", relationship = "many-to-many") %>%
    mutate(concordant = sign(logFC) == sign(sender_logFC)) %>%
    filter(concordant)

# Filtering
if (opt$p_adj_cutoff < 1.0) {
    cat(sprintf("Filtering by p_val_adj < %.2f\n", opt$p_adj_cutoff))
    tbl_base <- tbl_base %>% filter(p_val_adj < opt$p_adj_cutoff)
}

# Sorting Logic
sort_col <- opt$sort_by
cat(sprintf("Sorting by: %s\n", sort_col))

if (sort_col == "p_val_adj") {
    # Ascending for p-value (lower is better)
    tbl_base <- tbl_base %>% arrange(p_val_adj, desc(abs_logFC))
} else if (sort_col == "logFC" || sort_col == "abs_logFC") {
    tbl_base <- tbl_base %>% arrange(desc(abs(logFC)), p_val_adj)
} else if (sort_col == "activity") {
    tbl_base <- tbl_base %>% arrange(desc(activity), p_val_adj)
} else {
    # Default: Score
    tbl_base <- tbl_base %>% arrange(desc(score), p_val_adj)
}

# Deduplication (Keep best based on sort order)
top_interactions <- tbl_base %>%
    # Critical Fix: Deduplicate AFTER sort
    distinct(ligand, receptor, sender, receiver, .keep_all = TRUE)

# Greedy Selection with Constraints
# This ensures we try to fill 'top_n' by going deeper in the list if necessary,
# rather than just cutting the top list.

candidates <- top_interactions # The full sorted, deduplicated list

if (opt$max_per_sender > 0 || opt$max_per_receiver > 0) {
    cat("Applying greedy selection with constraints...\n")

    keep_indices <- c()
    n_selected <- 0

    # Initialize counts
    # Use standard vectors for speed
    s_counts <- setNames(rep(0, length(unique(candidates$sender))), unique(candidates$sender))
    r_counts <- setNames(rep(0, length(unique(candidates$receiver))), unique(candidates$receiver))

    max_s <- if (opt$max_per_sender > 0) opt$max_per_sender else Inf
    max_r <- if (opt$max_per_receiver > 0) opt$max_per_receiver else Inf

    # Loop through candidates
    for (i in seq_len(nrow(candidates))) {
        if (n_selected >= opt$top_n) break

        s <- as.character(candidates$sender[i])
        r <- as.character(candidates$receiver[i])

        # Check current counts (handle potential new levels safely if any, though uniq covers it)
        curr_s <- if (s %in% names(s_counts)) s_counts[[s]] else 0
        curr_r <- if (r %in% names(r_counts)) r_counts[[r]] else 0

        if (curr_s < max_s && curr_r < max_r) {
            keep_indices <- c(keep_indices, i)

            s_counts[[s]] <- curr_s + 1
            r_counts[[r]] <- curr_r + 1
            n_selected <- n_selected + 1
        }
    }

    top_interactions <- candidates[keep_indices, ]
} else {
    top_interactions <- head(candidates, opt$top_n)
}

# Add direction column which was missing in greedy path
top_interactions <- top_interactions %>%
    mutate(direction = ifelse(logFC > 0, "X2", "X1"))

cat("Interactions Selected (After Dedup):\n")
print(table(top_interactions$direction))

# 3. Layout Setup
senders <- top_interactions %>%
    select(sender, ligand) %>%
    distinct()
receivers <- top_interactions %>%
    select(receiver, receptor) %>%
    distinct()

senders$node <- paste0(senders$ligand, "|", senders$sender)
receivers$node <- paste0(receivers$receptor, "|", receivers$receiver)
senders$type <- "Sender"
receivers$type <- "Receiver"

all_nodes <- bind_rows(
    senders %>% rename(celltype = sender, molecule = ligand),
    receivers %>% rename(celltype = receiver, molecule = receptor)
) %>%
    arrange(desc(type), celltype, molecule)

sector_order <- all_nodes$node

# Colors
all_ct <- unique(all_nodes$celltype)
n_ct <- length(all_ct)
base_cols <- c(brewer.pal(8, "Set1"), brewer.pal(8, "Dark2"), "red", "blue", "green")
if (n_ct > length(base_cols)) {
    ct_pal <- colorRampPalette(base_cols)(n_ct)
} else {
    ct_pal <- base_cols[1:n_ct]
}
names(ct_pal) <- all_ct

# 4. Link Calculation Helper
calculate_link_positions <- function(data_subset) {
    df <- data_subset %>%
        mutate(
            from_node = paste0(ligand, "|", sender),
            to_node = paste0(receptor, "|", receiver)
        ) %>%
        filter(from_node %in% sector_order, to_node %in% sector_order) %>%
        arrange(from_node, receiver, receptor)

    if (nrow(df) == 0) {
        return(df)
    }

    assign_slots <- function(nodes_vec) {
        counts <- table(nodes_vec)
        slots <- list()
        for (node in names(counts)) {
            n <- counts[[node]]
            width <- 1 / n
            intervals <- list()
            for (i in 1:n) {
                intervals[[i]] <- c((i - 1) * width, i * width)
            }
            slots[[node]] <- intervals
        }
        return(slots)
    }

    from_slots <- assign_slots(df$from_node)
    to_slots <- assign_slots(df$to_node)

    df$from_start <- NA
    df$from_end <- NA
    df$to_start <- NA
    df$to_end <- NA

    from_consumed <- setNames(rep(1, length(from_slots)), names(from_slots))
    to_consumed <- setNames(rep(1, length(to_slots)), names(to_slots))

    for (i in 1:nrow(df)) {
        f <- df$from_node[i]
        t <- df$to_node[i]

        idx_f <- from_consumed[[f]]
        range_f <- from_slots[[f]][[idx_f]]
        df$from_start[i] <- range_f[1]
        df$from_end[i] <- range_f[2]
        from_consumed[[f]] <- idx_f + 1

        idx_t <- to_consumed[[t]]
        range_t <- to_slots[[t]][[idx_t]]
        df$to_start[i] <- range_t[1]
        df$to_end[i] <- range_t[2]
        to_consumed[[t]] <- idx_t + 1
    }
    return(df)
}

# 5. Plotting Function
plot_circos_final <- function(data_subset, title_text) {
    if (nrow(data_subset) > 0) {
        links_pos <- calculate_link_positions(data_subset)
    } else {
        links_pos <- data.frame()
    }

    circos.clear()

    gaps <- rep(0.5, length(sector_order))
    names(gaps) <- sector_order
    for (i in 2:length(sector_order)) {
        curr <- all_nodes[i, ]
        prev <- all_nodes[i - 1, ]
        if (curr$type != prev$type) {
            gaps[i - 1] <- 30
        } else if (curr$celltype != prev$celltype) gaps[i - 1] <- 5
    }
    gaps[length(gaps)] <- 30

    circos.par(
        gap.degree = gaps, start.degree = 90,
        canvas.xlim = c(-1.2, 1.2), canvas.ylim = c(-1.2, 1.2),
        cell.padding = c(0, 0, 0, 0)
    )

    # Ensure unique sector names and clean factor levels
    sector_order <- unique(sector_order)
    circos.initialize(factors = factor(sector_order, levels = sector_order), xlim = c(0, 1))

    # Track 1 (Outer): Gene Blocks + Labels
    circos.track(ylim = c(0, 1), track.height = 0.05, panel.fun = function(x, y) {
        sec <- CELL_META$sector.index
        parts <- strsplit(sec, "\\|")[[1]]
        mol <- parts[1]
        ct <- parts[2]

        circos.rect(0, 0, 1, 1, col = ct_pal[ct], border = "black", lwd = 1)
        circos.text(0.5, 2.5, mol, facing = "clockwise", niceFacing = TRUE, cex = 1.0, adj = c(0, 0.5))
    }, bg.border = NA)

    # Track 2 (Inner): Destination Hint Bar
    circos.track(ylim = c(0, 1), track.height = 0.02, bg.border = NA, panel.fun = function(x, y) {
        sec <- CELL_META$sector.index

        # Use Base R to avoid scoping issues with dplyr inside panel.fun
        my_links <- data.frame()
        if (nrow(links_pos) > 0 && "from_node" %in% colnames(links_pos)) {
            my_links <- links_pos[links_pos$from_node == sec, ]
        }

        if (nrow(my_links) > 0) {
            for (k in 1:nrow(my_links)) {
                link_rec_col <- ct_pal[my_links$receiver[k]]
                circos.rect(my_links$from_start[k], 0, my_links$from_end[k], 1,
                    col = link_rec_col, border = NA
                )
            }
        }
    })

    # Draw Links
    if (nrow(links_pos) > 0) {
        for (i in 1:nrow(links_pos)) {
            row <- links_pos[i, ]
            if (row$from_node %in% sector_order && row$to_node %in% sector_order) {
                # Add transparency correctly
                col_val <- adjustcolor(ct_pal[row$sender], alpha.f = 0.8) # 0.8 ~ CC
                circos.link(
                    row$from_node, c(row$from_start, row$from_end),
                    row$to_node, c(row$to_start, row$to_end),
                    col = col_val,
                    border = "black", lwd = 0.1,
                    directional = 1,
                    arr.length = 0.1,
                    arr.width = 0.1,
                    arr.type = "big.arrow"
                )
            }
        }
    }

    title(title_text, cex.main = 1.5)
}

# 6. Output PDF
output_file <- file.path(opt$outdir, "circos_comparison_final.pdf")
if (file.exists(output_file)) {
    i <- 1
    output_file <- file.path(opt$outdir, paste0("circos_comparison_final_", i, ".pdf"))
    while (file.exists(output_file)) {
        i <- i + 1
        output_file <- file.path(opt$outdir, paste0("circos_comparison_final_", i, ".pdf"))
    }
}
dir.create(dirname(output_file), recursive = TRUE)
cat("4. Saving to:", output_file, "\n")

pdf(output_file, width = 18, height = 9)

# Global Title
par(oma = c(0, 0, 3, 0))
layout(matrix(c(1, 2, 3), ncol = 3), widths = c(2, 2, 0.8))

par(mar = c(6, 1, 2, 1), xpd = TRUE)

links_X1 <- top_interactions %>% filter(direction == "X1")
links_X2 <- top_interactions %>% filter(direction == "X2")

# Plot 1
plot_circos_final(links_X1, "Control (X1)")
mtext(
    side = 1, line = 4, cex = 0.8, adj = 0,
    text = paste0(
        "Model: muscat (edgeR)\n",
        meta_info, "\n",
        "Path: ", results_path, "\n",
        "N_interactions: ", nrow(links_X1)
    )
)

# Plot 2
plot_circos_final(links_X2, "Treatment (X2)")
mtext(
    side = 1, line = 4, cex = 0.8, adj = 0,
    text = paste0(
        "Contrast: X2 vs X1\n",
        "Significance: p < ", opt$p_val_threshold, ", |logFC| > ", opt$logfc_threshold, "\n",
        "P-adj Cutoff: < ", opt$p_adj_cutoff, "\n",
        "N_interactions: ", nrow(links_X2)
    )
)

# Global Title
mtext("MultiNicheNetR Analysis: IS2 vs IS3 (Group X2 vs X1)", outer = TRUE, cex = 1.8, font = 2, line = 1)

# Legend
plot.new()
par(mar = c(0, 0, 0, 0))
text(0.05, 0.6, "Cell Types", cex = 1.5, font = 2, adj = c(0, 0))
legend(
    x = 0.05, y = 0.55,
    legend = names(ct_pal),
    fill = ct_pal,
    bty = "n",
    cex = 1.2,
    pt.cex = 2,
    xjust = 0, yjust = 1
)

dev.off()

cat("\n✅ Final Circos Plot Generated!\n")
