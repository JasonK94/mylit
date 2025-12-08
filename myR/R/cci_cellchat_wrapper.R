
#' Run CellChat Analysis
#'
#' This function wraps the CellChat analysis pipeline.
#' It takes a Seurat object and runs the standard CellChat workflow.
#'
#' @param sobj A Seurat object.
#' @param group.by Character string. The metadata column to use for cell labels (e.g., cell types).
#' @param species Character string. "human" or "mouse". Default: "human".
#' @param db.use Character vector. Subset of CellChatDB to use. 
#'        Options: "Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact". 
#'        If NULL, uses the entire database. Default: NULL.
#' @param assay_name Character string. Assay to use. Default: "SCT".
#' @param output_dir Character string. Path to save results. Default: NULL.
#' @param n_cores Integer. Number of cores for parallel processing. Default: 4.
#' @param do_parallel Logical. Whether to use parallel processing. Default: TRUE.
#' @param verbose Logical. Default: TRUE.
#'
#' @return A CellChat object with results.
#' @export
#' @importFrom CellChat createCellChat setIdent subsetData identifyOverExpressedGenes identifyOverExpressedInteractions computeCommunProb filterCommunication computeCommunProbPathway aggregateNet netVisual_circle
#' @importFrom Seurat DefaultAssay
#' @importFrom future plan multiprocess
run_cellchat_analysis <- function(sobj,
                                  group.by,
                                  species = c("human", "mouse"),
                                  db.use = NULL,
                                  assay_name = "SCT",
                                  output_dir = NULL,
                                  n_cores = 4,
                                  do_parallel = TRUE,
                                  verbose = TRUE) {
  
  species <- match.arg(species)
  
  if (verbose) message("=== Starting CellChat Analysis ===")
  
  # Check packages
  if (!requireNamespace("CellChat", quietly = TRUE)) {
    stop("Package 'CellChat' is required but not installed.")
  }
  
  # Create output directory
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    if (verbose) message("Output directory: ", output_dir)
  }
  
  # Prepare Seurat object
  if (verbose) message("Preparing Seurat object...")
  if (!group.by %in% colnames(sobj@meta.data)) {
    stop("Column '", group.by, "' not found in metadata.")
  }
  
  # Create CellChat object
  if (verbose) message("Creating CellChat object...")
  # Use the data slot from the specified assay
  # CellChat expects normalized data
  cellchat <- CellChat::createCellChat(object = sobj, group.by = group.by, assay = assay_name)
  
  # Set database
  if (verbose) message("Setting CellChatDB for ", species, "...")
  if (species == "human") {
    CellChatDB <- CellChat::CellChatDB.human
  } else {
    CellChatDB <- CellChat::CellChatDB.mouse
  }
  
  if (!is.null(db.use)) {
    CellChatDB.use <- CellChat::subsetDB(CellChatDB, search = db.use)
    if (verbose) message("Using subset of DB: ", paste(db.use, collapse = ", "))
  } else {
    CellChatDB.use <- CellChatDB
    if (verbose) message("Using full DB")
  }
  cellchat@DB <- CellChatDB.use
  
  # Pre-processing
  if (verbose) message("Preprocessing expression data...")
  cellchat <- CellChat::subsetData(cellchat) # Subset to genes in DB
  
  if (do_parallel) {
    if (verbose) message("Setting up parallel processing with ", n_cores, " cores...")
    future::plan("multisession", workers = n_cores)
  }
  
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  
  # Compute communication probability
  if (verbose) message("Computing communication probabilities...")
  # type = "triMean" is standard
  cellchat <- CellChat::computeCommunProb(cellchat, type = "triMean")
  
  # Filter communication
  if (verbose) message("Filtering communication...")
  cellchat <- CellChat::filterCommunication(cellchat, min.cells = 10)
  
  # Compute pathway probabilities
  if (verbose) message("Computing pathway probabilities...")
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  
  # Aggregate network
  if (verbose) message("Aggregating network...")
  cellchat <- CellChat::aggregateNet(cellchat)
  
  # Visualization and Saving
  if (!is.null(output_dir)) {
    if (verbose) message("Generating and saving plots...")
    
    # 1. Circle plot of aggregated network
    pdf(file.path(output_dir, "net_visual_circle.pdf"), width = 10, height = 10)
    groupSize <- as.numeric(table(cellchat@idents))
    CellChat::netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
    CellChat::netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
    dev.off()
    
    # 2. Bubble plot for all pathways (top 20)
    # Need to check if there are significant pathways
    if (length(cellchat@netP$pathways) > 0) {
        pdf(file.path(output_dir, "net_visual_bubble.pdf"), width = 12, height = 8)
        tryCatch({
            print(CellChat::netVisual_bubble(cellchat, sources.use = NULL, targets.use = NULL, remove.isolate = FALSE))
        }, error = function(e) message("Error in bubble plot: ", e$message))
        dev.off()
    }
    
    # Save CellChat object
    if (verbose) message("Saving CellChat object...")
    saveRDS(cellchat, file.path(output_dir, "cellchat_object.rds"))
    if (requireNamespace("qs", quietly = TRUE)) {
        qs::qsave(cellchat, file.path(output_dir, "cellchat_object.qs"))
    }
  }
  
  if (verbose) message("=== CellChat Analysis Complete ===")
  return(cellchat)
}
