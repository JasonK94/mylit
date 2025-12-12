#' Run MultiNicheNet Analysis
#'
#' Performs MultiNicheNet analysis to infer cell-cell communications across conditions.
#'
#' @param sobj Seurat object
#' @param sample_id Character string, metadata column for sample ID (e.g., patient + condition)
#' @param group_id Character string, metadata column for group ID (e.g., condition)
#' @param celltype_id Character string, metadata column for cell type
#' @param covariates Character vector, metadata columns for covariates (default: NA)
#' @param batches Character string, metadata column for batch (default: NA)
#' @param contrasts_oi Character vector, contrasts for DE analysis (e.g., c("GroupA-GroupB"))
#' @param senders_oi Character vector, sender cell types (default: NULL, uses all)
#' @param receivers_oi Character vector, receiver cell types (default: NULL, uses all)
#' @param min_cells Numeric, minimum number of cells per cell type per sample (default: 10)
#' @param species Character string, "human" or "mouse" (default: "human")
#' @param nichenet_data_dir Character string, directory for NicheNet data (default: NULL)
#' @param output_dir Character string, output directory (default: NULL)
#' @param verbose Logical, whether to print progress (default: TRUE)
#' @param cores Numeric, number of cores to use (default: 1)
#'
#' @return List containing MultiNicheNet results
#' @export
#' @importFrom Seurat as.SingleCellExperiment
#' @importFrom SummarizedExperiment colData
#' @importFrom multinichenetr get_abundance_expression_info get_DE_info combine_sender_receiver_de get_ligand_activities_targets_DEgenes generate_prioritization_tables
#' @importFrom dplyr %>% filter select distinct rename mutate arrange top_n
#' @importFrom tibble tibble as_tibble
run_multinichenet_analysis <- function(sobj,
                                       sample_id,
                                       group_id,
                                       celltype_id,
                                       covariates = NA,
                                       batches = NA,
                                       contrasts_oi,
                                       senders_oi = NULL,
                                       receivers_oi = NULL,
                                       min_cells = 10,
                                       species = "human",
                                       nichenet_data_dir = NULL,
                                       output_dir = NULL,
                                       verbose = TRUE,
                                       cores = 1,
                                       # Threshold parameters
                                       min_pct = 0.05,
                                       logfc_thresh = 0.10,
                                       p_val_thresh = 0.05,
                                       p_val_adj = FALSE) {
    if (verbose) message("Starting MultiNicheNet analysis...")

    # Ensure required packages are loaded
    if (!requireNamespace("multinichenetr", quietly = TRUE)) {
        stop("Package 'multinichenetr' is required but not installed.")
    }

    # Convert Seurat to SingleCellExperiment
    if (verbose) message("Converting Seurat object to SingleCellExperiment...")
    sce <- Seurat::as.SingleCellExperiment(sobj, assay = "RNA") # Use RNA assay for DE

    # Ensure metadata columns are characters or factors
    SummarizedExperiment::colData(sce)[[sample_id]] <- make.names(SummarizedExperiment::colData(sce)[[sample_id]])
    # Group ID must be a factor for contrasts to work properly
    SummarizedExperiment::colData(sce)[[group_id]] <- as.factor(make.names(SummarizedExperiment::colData(sce)[[group_id]]))
    SummarizedExperiment::colData(sce)[[celltype_id]] <- make.names(SummarizedExperiment::colData(sce)[[celltype_id]])

    # Define cell types if not provided
    if (is.null(senders_oi)) {
        senders_oi <- unique(SummarizedExperiment::colData(sce)[[celltype_id]])
    }
    if (is.null(receivers_oi)) {
        receivers_oi <- unique(SummarizedExperiment::colData(sce)[[celltype_id]])
    }

    # Load NicheNet data
    if (verbose) message("Loading NicheNet data...")

    base_url <- "https://zenodo.org/records/7074291/files/"
    files_to_load_spec <- list(
        human = list(
            lr_network = "lr_network_human_21122021.rds",
            ligand_target_matrix = "ligand_target_matrix_nsga2r_final.rds"
        ),
        mouse = list(
            lr_network = "lr_network_mouse_21122021.rds",
            ligand_target_matrix = "ligand_target_matrix_nsga2r_final_mouse.rds"
        )
    )

    required_files <- files_to_load_spec[[species]]

    # Helper to load data
    load_nn_data <- function(file_name, url_base, dir_path) {
        # Check for .qs version first
        qs_file_name <- sub("\\.rds$", ".qs", file_name)
        qs_path <- file.path(dir_path, qs_file_name)
        rds_path <- file.path(dir_path, file_name)

        if (file.exists(qs_path)) {
            if (verbose) message("Loading ", qs_file_name, " (qs)...")
            return(qs::qread(qs_path))
        } else if (file.exists(rds_path)) {
            if (verbose) message("Loading ", file_name, " (rds)...")
            return(readRDS(rds_path))
        } else {
            # Download if neither exists
            if (verbose) message("Downloading ", file_name, "...")
            download.file(paste0(url_base, file_name), rds_path, mode = "wb", quiet = !verbose)
            return(readRDS(rds_path))
        }
    }

    if (is.null(nichenet_data_dir)) {
        # Check shared directory first
        shared_dir <- file.path("/data/user3/git_repo", species)
        if (dir.exists(shared_dir)) {
            nichenet_data_dir <- shared_dir
            if (verbose) message("Using shared NicheNet data directory: ", nichenet_data_dir)
        } else {
            nichenet_data_dir <- file.path(getwd(), paste0("nichenet_data_", species))
            if (!dir.exists(nichenet_data_dir)) dir.create(nichenet_data_dir, recursive = TRUE)
            if (verbose) message("Using local NicheNet data directory: ", nichenet_data_dir)
        }
    }

    lr_network <- load_nn_data(required_files$lr_network, base_url, nichenet_data_dir)
    ligand_target_matrix <- load_nn_data(required_files$ligand_target_matrix, base_url, nichenet_data_dir)

    lr_network <- lr_network %>% dplyr::distinct(from, to)

    # Prepare contrast table
    # Assuming contrasts_oi is a vector of strings like "GroupA-GroupB"
    # We need to extract the group name (GroupA) for the contrast table
    if (verbose) message("Preparing contrast table...")

    # Simple parsing: take the part before the first "-" as the group
    # This assumes standard contrast naming "GroupA-GroupB"
    groups <- sapply(contrasts_oi, function(x) {
        parts <- strsplit(x, "-")[[1]]
        if (length(parts) >= 1) {
            return(parts[1])
        } else {
            return(NA)
        }
    })

    contrast_tbl <- tibble::tibble(
        contrast = contrasts_oi,
        group = groups
    )

    if (verbose) {
        message("Contrast table:")
        print(contrast_tbl)
    }

    # Run MultiNicheNet Analysis
    if (verbose) message("Running multi_nichenet_analysis wrapper...")
    if (verbose) {
        message(sprintf(
            "  Params: min_pct=%.2f, logFC=%.2f, p_val=%.2f, adj=%s",
            min_pct, logfc_thresh, p_val_thresh, p_val_adj
        ))
    }

    results <- multinichenetr::multi_nichenet_analysis(
        sce = sce,
        celltype_id = celltype_id,
        sample_id = sample_id,
        group_id = group_id,
        batches = batches,
        covariates = covariates,
        lr_network = lr_network,
        ligand_target_matrix = ligand_target_matrix,
        contrasts_oi = contrasts_oi,
        contrast_tbl = contrast_tbl,
        senders_oi = senders_oi,
        receivers_oi = receivers_oi,
        min_cells = min_cells,
        verbose = verbose,
        n.cores = cores,
        # Pass dynamic parameters
        fraction_cutoff = min_pct,
        logFC_threshold = logfc_thresh,
        p_val_threshold = p_val_thresh,
        p_val_adj = p_val_adj
    )

    if (verbose) message("MultiNicheNet analysis completed.")

    if (!is.null(output_dir)) {
        if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

        # Save with qs (preferred) and RDS (fallback)
        qs_path <- file.path(output_dir, "multinichenet_results.qs")
        rds_path <- file.path(output_dir, "multinichenet_results.rds")

        tryCatch(
            {
                qs::qsave(results, qs_path)
                if (verbose) message("Results saved to: ", qs_path)
            },
            error = function(e) {
                if (verbose) message("Warning: qs::qsave failed, falling back to saveRDS")
                saveRDS(results, rds_path)
                if (verbose) message("Results saved to: ", rds_path)
            }
        )
    }

    return(results)
}

#' Plot MultiNicheNet Circos
#'
#' Generates Circos plots for MultiNicheNet results, separated by group.
#'
#' @param mnn_results List, output from run_multinichenet_analysis
#' @param group_oi Character string, group to visualize (e.g., "Inflamed")
#' @param top_n Numeric, number of top interactions to show (default: 50)
#' @param output_file Character string, path to save the plot (default: NULL)
#' @export
#' @importFrom circlize circos.clear chordDiagram circos.track circos.text uy mm_h CELL_META
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#' @importFrom grDevices pdf dev.off
plot_multinichenet_circos <- function(mnn_results, group_oi, top_n = 50, output_file = NULL) {
    if (!requireNamespace("circlize", quietly = TRUE)) {
        stop("Package 'circlize' is required for Circos plots.")
    }

    prioritization_tbl <- mnn_results$prioritization_tables$group_prioritization_tbl

    if (is.null(prioritization_tbl)) {
        stop("Prioritization table not found in results.")
    }

    # Filter for the group of interest
    # Note: The structure depends on how generate_prioritization_tables was called.
    # If it has a 'group' column, we filter.

    if ("group" %in% colnames(prioritization_tbl)) {
        tbl_filtered <- prioritization_tbl %>% dplyr::filter(group == group_oi)
    } else {
        warning("Group column not found in prioritization table. Using all data.")
        tbl_filtered <- prioritization_tbl
    }

    if (nrow(tbl_filtered) == 0) {
        warning("No interactions found for group: ", group_oi)
        return(NULL)
    }

    # Select top interactions
    # Assuming there is a prioritization score column.
    # multinichenetr usually produces a 'prioritization_score' or similar.
    # Let's check for common score columns.
    score_col <- "prioritization_score"
    if (!score_col %in% colnames(tbl_filtered)) {
        # Try to find a score column
        numeric_cols <- names(which(sapply(tbl_filtered, is.numeric)))
        if (length(numeric_cols) > 0) {
            score_col <- numeric_cols[length(numeric_cols)] # Assume last numeric is score?
            message("Using column '", score_col, "' as prioritization score.")
        } else {
            stop("No numeric score column found in prioritization table.")
        }
    }

    top_tbl <- tbl_filtered %>%
        dplyr::arrange(dplyr::desc(!!rlang::sym(score_col))) %>%
        head(top_n)

    # Prepare data for Circos
    # We want Sender -> Receiver chords, or Ligand -> Receptor.
    # Standard MultiNicheNet Circos: Ligand (colored by Sender) -> Receptor (colored by Receiver or just Receptor)

    # Create links dataframe: from = ligand, to = receptor, value = weight
    links_df <- top_tbl %>%
        dplyr::select(ligand, receptor, !!rlang::sym(score_col)) %>%
        dplyr::rename(from = ligand, to = receptor, weight = !!rlang::sym(score_col))

    # Define colors
    # Sender colors for ligands
    senders <- unique(top_tbl$sender)
    receivers <- unique(top_tbl$receiver)

    # Generate palette for senders
    sender_pal <- RColorBrewer::brewer.pal(max(3, length(senders)), "Set1")[1:length(senders)]
    names(sender_pal) <- senders

    # Map ligands to senders
    ligand_sender_map <- top_tbl %>%
        dplyr::select(ligand, sender) %>%
        dplyr::distinct()
    # Handle ligands appearing in multiple senders (take first or most prominent)
    ligand_sender_map <- ligand_sender_map %>%
        dplyr::group_by(ligand) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup()

    # Map receptors to receivers (optional, or just one color for receptors)
    receptor_receiver_map <- top_tbl %>%
        dplyr::select(receptor, receiver) %>%
        dplyr::distinct()

    # Create grid colors
    grid_cols <- c()

    # Ligands colored by sender
    for (i in 1:nrow(ligand_sender_map)) {
        lig <- ligand_sender_map$ligand[i]
        snd <- ligand_sender_map$sender[i]
        grid_cols[lig] <- sender_pal[snd]
    }

    # Receptors colored by a distinct color or by receiver
    receptor_color <- "grey70"
    for (rec in unique(links_df$to)) {
        grid_cols[rec] <- receptor_color
    }

    # Plot
    if (!is.null(output_file)) {
        pdf(output_file, width = 10, height = 10)
    }

    circlize::circos.clear()
    circlize::chordDiagram(links_df,
        grid.col = grid_cols,
        transparency = 0.25,
        directional = 1,
        direction.type = c("diffHeight", "arrows"),
        link.arr.type = "big.arrow",
        annotationTrack = "grid",
        preAllocateTracks = list(track.height = 0.1)
    )

    circlize::circos.track(track.index = 1, panel.fun = function(x, y) {
        circlize::circos.text(circlize::CELL_META$xcenter,
            circlize::CELL_META$ylim[1],
            circlize::CELL_META$sector.index,
            facing = "clockwise",
            niceFacing = TRUE,
            adj = c(0, 0.5),
            cex = 0.8
        )
    }, bg.border = NA)

    title(paste("MultiNicheNet Interactions -", group_oi))

    # Legend
    legend("topright", legend = names(sender_pal), fill = sender_pal, title = "Sender", bty = "n")

    if (!is.null(output_file)) {
        dev.off()
        message("Circos plot saved to ", output_file)
    }

    return(invisible(links_df))
}
