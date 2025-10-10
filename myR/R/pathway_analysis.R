# GO/GSEA Analysis Functions for Stroke PBMC scRNA-seq Data
# Author: Analysis Pipeline
# Purpose: Modular functions for pathway enrichment analysis

# NOTE: Package dependencies should be declared in DESCRIPTION, not with library() calls
# Required: clusterProfiler, org.Hs.eg.db, DOSE, enrichplot, msigdbr, fgsea, dplyr, tibble

#' Convert Gene IDs Between Different Formats
#'
#' This function converts gene identifiers between different formats using
#' the org.Hs.eg.db annotation package.
#'
#' @param genes Character vector of gene identifiers to convert
#' @param from Character string specifying the input ID type (default: "SYMBOL")
#' @param to Character string specifying the output ID type (default: "ENTREZID")
#' @return Data frame with mapping between input and output gene IDs
#' @examples
#' \dontrun{
#' genes <- c("TP53", "BRCA1", "EGFR")
#' entrez_ids <- convert_gene_ids(genes, from = "SYMBOL", to = "ENTREZID")
#' }
#' @export
convert_gene_ids <- function(genes, from = "SYMBOL", to = "ENTREZID") {
  # Convert gene symbols to other formats (ENTREZID, ENSEMBL, etc.)
  converted <- bitr(genes, 
                    fromType = from,
                    toType = to,
                    OrgDb = org.Hs.eg.db,
                    drop = TRUE)
  
  # Remove duplicates and return mapping
  converted <- converted[!duplicated(converted[[from]]), ]
  return(converted)
}

#' Prepare Gene Lists for Pathway Analysis
#'
#' This function processes differential expression results to create gene lists
#' suitable for various pathway enrichment analyses.
#'
#' @param deg_df Data frame containing differential expression results with columns:
#'   gene, avg_log2FC, p_val, p_val_adj
#' @param fc_threshold Numeric threshold for log2 fold change (default: 0.25)
#' @param p_use Character string specifying which p-value to use: "p_val_adj" or "p_val" (default: "p_val_adj")
#' @param pval_threshold Numeric threshold for p-value significance (default: 0.05)
#' @return List containing ranked genes, up/down-regulated genes, and background
#' @examples
#' \dontrun{
#' gene_lists <- prepare_gene_lists(deg_df, fc_threshold = 0.5, pval_threshold = 0.01)
#' }
#' @export
prepare_gene_lists <- function(deg_df, fc_threshold = 0.25, p_use = "p_val_adj", pval_threshold = 0.05) {
  # Validate p_use parameter
  if (!p_use %in% c("p_val_adj", "p_val")) {
    stop("p_use must be either 'p_val_adj' or 'p_val'")
  }
  
  # Select p-value column
  p_col <- deg_df[[p_use]]
  
  # Filter significant genes
  sig_genes <- deg_df %>%
    filter(abs(avg_log2FC) >= fc_threshold & p_col < pval_threshold)
  
  # Create ranked gene list for GSEA (all genes)
  ranked_genes <- deg_df %>%
    arrange(desc(avg_log2FC)) %>%
    pull(avg_log2FC, name = gene)
  
  # Significant gene lists for GO
  up_genes <- sig_genes %>%
    filter(avg_log2FC > 0) %>%
    pull(gene)
  
  down_genes <- sig_genes %>%
    filter(avg_log2FC < 0) %>%
    pull(gene)
  
  all_sig_genes <- sig_genes$gene
  
  return(list(
    ranked = ranked_genes,
    up = up_genes,
    down = down_genes,
    all_sig = all_sig_genes,
    background = deg_df$gene
  ))
}

#' Get Available Pathway Database Sets
#'
#' Returns available pathway database identifiers for enrichment analysis.
#'
#' @param pathway_set Character vector of specific pathway sets to return, 
#'   or NULL to return all available sets
#' @param species Character string specifying species (default: "Homo sapiens")
#' @return Character vector of pathway set identifiers
#' @examples
#' \dontrun{
#' all_sets <- get_pathway_sets()
#' hallmark_only <- get_pathway_sets("H")
#' }
#' @export
get_pathway_sets <- function(pathway_set = NULL, species = "Homo sapiens") {
  # Define available pathway sets
  available_sets <- list(
    "GO_BP" = "GO Biological Process",
    "GO_MF" = "GO Molecular Function", 
    "GO_CC" = "GO Cellular Component",
    "KEGG" = "KEGG",
    "REACTOME" = "Reactome",
    "HALLMARK" = "Hallmark",
    "C2_KEGG" = "MSigDB C2 KEGG",
    "C2_REACTOME" = "MSigDB C2 Reactome",
    "C5_BP" = "MSigDB C5 GO BP",
    "C5_MF" = "MSigDB C5 GO MF",
    "C5_CC" = "MSigDB C5 GO CC",
    "H" = "MSigDB Hallmark"
  )
  
  if (is.null(pathway_set)) {
    return(names(available_sets))
  } else {
    return(pathway_set)
  }
}

#' Run GO Enrichment Analysis
#'
#' Performs Gene Ontology enrichment analysis using clusterProfiler.
#'
#' @param genes Character vector of gene symbols for analysis
#' @param gene_type Character string describing gene set type (for logging)
#' @param ont Character string specifying GO ontology: "BP", "MF", or "CC"
#' @param background Character vector of background genes (optional)
#' @param pval_cutoff Numeric p-value cutoff for significance (default: 0.05)
#' @return enrichResult object from clusterProfiler, or NULL if no results
#' @examples
#' \dontrun{
#' go_result <- run_go_analysis(c("TP53", "BRCA1"), ont = "BP")
#' }
#' @export
run_go_analysis <- function(genes, gene_type = "all", ont = "BP", 
                            background = NULL, pval_cutoff = 0.05) {
  
  # Convert gene symbols to ENTREZID
  gene_ids <- convert_gene_ids(genes, from = "SYMBOL", to = "ENTREZID")
  
  if (nrow(gene_ids) == 0) {
    warning("No genes could be converted to ENTREZID")
    return(NULL)
  }
  
  # Convert background if provided
  if (!is.null(background)) {
    bg_ids <- convert_gene_ids(background, from = "SYMBOL", to = "ENTREZID")
    universe <- bg_ids$ENTREZID
  } else {
    universe <- NULL
  }
  
  # Run GO enrichment
  go_result <- enrichGO(gene = gene_ids$ENTREZID,
                        OrgDb = org.Hs.eg.db,
                        ont = ont,
                        pAdjustMethod = "BH",
                        pvalueCutoff = pval_cutoff,
                        qvalueCutoff = 0.2,
                        universe = universe,
                        readable = TRUE)
  
  return(go_result)
}

#' Run KEGG Pathway Analysis
#'
#' Performs KEGG pathway enrichment analysis using clusterProfiler.
#'
#' @param genes Character vector of gene symbols for analysis
#' @param background Character vector of background genes (optional)
#' @param pval_cutoff Numeric p-value cutoff for significance (default: 0.05)
#' @return enrichResult object from clusterProfiler, or NULL if no results
#' @examples
#' \dontrun{
#' kegg_result <- run_kegg_analysis(c("TP53", "BRCA1"))
#' }
#' @export
run_kegg_analysis <- function(genes, background = NULL, pval_cutoff = 0.05) {
  
  # Convert to ENTREZID
  gene_ids <- convert_gene_ids(genes, from = "SYMBOL", to = "ENTREZID")
  
  if (nrow(gene_ids) == 0) {
    warning("No genes could be converted to ENTREZID")
    return(NULL)
  }
  
  # Convert background if provided
  if (!is.null(background)) {
    bg_ids <- convert_gene_ids(background, from = "SYMBOL", to = "ENTREZID")
    universe <- bg_ids$ENTREZID
  } else {
    universe <- NULL
  }
  
  # Run KEGG enrichment
  kegg_result <- enrichKEGG(gene = gene_ids$ENTREZID,
                            organism = 'hsa',
                            pvalueCutoff = pval_cutoff,
                            pAdjustMethod = "BH",
                            universe = universe)
  
  # Convert back to readable gene symbols
  if (!is.null(kegg_result) && nrow(kegg_result@result) > 0) {
    kegg_result <- setReadable(kegg_result, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  }
  
  return(kegg_result)
}

#' Run Gene Set Enrichment Analysis (GSEA)
#'
#' Performs GSEA using fgsea with MSigDB gene sets.
#'
#' @param ranked_genes Named numeric vector of genes ranked by fold change
#' @param pathway_set Character string specifying pathway database (default: "H")
#' @param min_size Minimum gene set size (default: 15)
#' @param max_size Maximum gene set size (default: 500)
#' @return Data frame with GSEA results
#' @examples
#' \dontrun{
#' ranked_genes <- c(2.5, 1.8, -1.2, -2.1)
#' names(ranked_genes) <- c("TP53", "BRCA1", "EGFR", "MYC")
#' gsea_result <- run_gsea_analysis(ranked_genes, pathway_set = "H")
#' }
#' @export
run_gsea_analysis <- function(ranked_genes, pathway_set = "H", 
                              min_size = 15, max_size = 500) {
  
  # Get MSigDB gene sets
  if (pathway_set == "H") {
    gene_sets <- msigdbr(species = "Homo sapiens", category = "H")
  } else if (pathway_set == "C2_KEGG") {
    gene_sets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
  } else if (pathway_set == "C2_REACTOME") {
    gene_sets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
  } else if (pathway_set == "C5_BP") {
    gene_sets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
  } else if (pathway_set == "C5_MF") {
    gene_sets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF")
  } else if (pathway_set == "C5_CC") {
    gene_sets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:CC")
  } else {
    gene_sets <- msigdbr(species = "Homo sapiens", category = "H")  # Default to Hallmark
  }
  
  # Convert to list format for fgsea
  gene_sets_list <- split(gene_sets$gene_symbol, gene_sets$gs_name)
  
  # Remove genes not in our dataset
  gene_sets_list <- lapply(gene_sets_list, function(x) x[x %in% names(ranked_genes)])
  gene_sets_list <- gene_sets_list[lengths(gene_sets_list) >= min_size & 
                                     lengths(gene_sets_list) <= max_size]
  
  # Run GSEA using fgseaMultilevel (recommended)
  gsea_result <- fgseaMultilevel(pathways = gene_sets_list,
                                 stats = ranked_genes,
                                 minSize = min_size,
                                 maxSize = max_size)
  
  # Convert to data frame and fix data.table issues
  gsea_result <- as.data.frame(gsea_result)
  
  return(gsea_result)
}

#' Format Analysis Results
#'
#' Standardizes the output format from different pathway analysis methods.
#'
#' @param result Results object from pathway analysis
#' @param analysis_type Character string specifying analysis type: "GO", "KEGG", or "GSEA"
#' @return Data frame with standardized columns: pathway, description, p_value, 
#'   adj_p_value, effect_size, gene_count, genes
#' @examples
#' \dontrun{
#' formatted <- format_results(go_result, "GO")
#' }
#' @export
format_results <- function(result, analysis_type = "GO") {
  
  if (is.null(result)) {
    return(data.frame(
      pathway = character(0),
      description = character(0),
      p_value = numeric(0),
      adj_p_value = numeric(0),
      effect_size = numeric(0),
      gene_count = integer(0),
      genes = character(0)
    ))
  }
  
  if (analysis_type == "GSEA") {
    # Format GSEA results
    formatted <- result %>%
      dplyr::arrange(padj) %>%
      dplyr::select(pathway, pval, padj, ES, NES, size, leadingEdge) %>%
      dplyr::rename(
        p_value = pval,
        adj_p_value = padj,
        effect_size = NES,
        gene_count = size,
        genes = leadingEdge
      ) %>%
      dplyr::mutate(
        description = pathway,
        genes = sapply(genes, function(x) paste(x, collapse = "/"))
      )
  } else {
    # Format GO/KEGG results
    res_df <- as.data.frame(result)
    if (nrow(res_df) == 0) {
      return(data.frame(
        pathway = character(0),
        description = character(0),
        p_value = numeric(0),
        adj_p_value = numeric(0),
        effect_size = numeric(0),
        gene_count = integer(0),
        genes = character(0)
      ))
    }
    
    formatted <- res_df %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::select(ID, Description, pvalue, p.adjust, Count, geneID) %>%
      dplyr::rename(
        pathway = ID,
        description = Description,
        p_value = pvalue,
        adj_p_value = p.adjust,
        gene_count = Count,
        genes = geneID
      ) %>%
      dplyr::mutate(
        effect_size = -log10(p_value)  # Use -log10(p-value) as effect size for enrichment
      )
  }
  
  return(formatted)
}

#' Comprehensive GO and GSEA Analysis
#'
#' Main function to perform comprehensive pathway enrichment analysis on 
#' differential expression results from scRNA-seq data.
#'
#' @param DEG Data frame containing differential expression results with required columns:
#'   \code{gene}, \code{avg_log2FC}, \code{p_val}, \code{p_val_adj}
#' @param seurat_obj Seurat object (currently not used, for future compatibility)
#' @param pathway_set Character vector specifying which pathway databases to use.
#'   Options include: "H", "C2_KEGG", "C2_REACTOME", "C5_BP", "C5_MF", "C5_CC".
#'   If NULL, uses default set.
#' @param pathway Character string to filter results for specific pathway names (optional)
#' @param analysis_type Character string specifying analysis type: "GO", "KEGG", "GSEA", or "ALL"
#' @param go_ontology Character vector specifying GO ontologies: "BP", "MF", "CC"
#' @param fc_threshold Numeric threshold for log2 fold change significance (default: 0.25)
#' @param p_use Character string specifying which p-value column to use: "p_val_adj" or "p_val" (default: "p_val_adj")
#' @param pval_threshold Numeric threshold for p-value significance (default: 0.05)
#' @param gsea_min_size Minimum gene set size for GSEA (default: 15)
#' @param gsea_max_size Maximum gene set size for GSEA (default: 500)
#' @param return_plots Logical indicating whether to return plots (default: FALSE, not implemented)
#' @return Named list containing analysis results for each method. Each element is a 
#'   data frame with columns: pathway, description, p_value, adj_p_value, 
#'   effect_size, gene_count, genes
#' @examples
#' \dontrun{
#' # Basic usage
#' results <- myGO(DEG = deg_dataframe, analysis_type = "ALL")
#' 
#' # Use raw p-values instead of adjusted
#' results <- myGO(DEG = deg_dataframe, p_use = "p_val", analysis_type = "GSEA")
#' 
#' # Focus on specific pathways
#' stroke_results <- myGO(DEG = deg_dataframe, pathway = "INFLAMMATORY")
#' 
#' # Access specific results
#' hallmark_gsea <- results$GSEA_H
#' go_bp_up <- results$GO_BP_UP
#' }
#' @export
myGO <- function(DEG = NULL, 
                 seurat_obj = NULL,
                 pathway_set = NULL, 
                 pathway = NULL,
                 analysis_type = c("GO", "KEGG", "GSEA", "ALL"),
                 go_ontology = c("BP", "MF", "CC"),
                 fc_threshold = 0.25,
                 p_use = "p_val_adj",
                 pval_threshold = 0.05,
                 gsea_min_size = 15,
                 gsea_max_size = 500,
                 return_plots = FALSE) {
  
  # Input validation
  analysis_type <- match.arg(analysis_type)
  go_ontology <- match.arg(go_ontology, several.ok = TRUE)
  
  # Extract DEG dataframe if Seurat object provided
  if (!is.null(seurat_obj) && is.null(DEG)) {
    stop("Please provide the DEG dataframe. If you have a Seurat object, extract DEGs first using FindMarkers() or similar.")
  }
  
  if (is.null(DEG)) {
    stop("Please provide either a DEG dataframe or Seurat object")
  }
  
  # Validate DEG dataframe columns
  required_cols <- c("gene", "avg_log2FC", "p_val", "p_val_adj")
  if (!all(required_cols %in% colnames(DEG))) {
    stop(paste("DEG dataframe must contain columns:", paste(required_cols, collapse = ", ")))
  }
  
  # Validate p_use parameter
  if (!p_use %in% c("p_val_adj", "p_val")) {
    stop("p_use must be either 'p_val_adj' or 'p_val'")
  }
  
  # Prepare gene lists
  cat("Preparing gene lists...\n")
  gene_lists <- prepare_gene_lists(DEG, fc_threshold, p_use, pval_threshold)
  
  cat(sprintf("Found %d upregulated, %d downregulated, %d total significant genes\n",
              length(gene_lists$up), length(gene_lists$down), length(gene_lists$all_sig)))
  
  # Initialize results list
  results <- list()
  
  # Determine which pathway sets to analyze
  if (is.null(pathway_set)) {
    pathway_sets_to_analyze <- c("H")  # Default to Hallmark for GSEA
  } else {
    pathway_sets_to_analyze <- get_pathway_sets(pathway_set)
  }
  
  # Filter for specific pathway if requested
  if (!is.null(pathway)) {
    cat(sprintf("Analyzing specific pathway: %s\n", pathway))
  }
  
  # Run analyses based on type requested
  if (analysis_type %in% c("GO", "ALL")) {
    cat("Running GO enrichment analysis...\n")
    
    for (ont in go_ontology) {
      cat(sprintf("  - GO %s analysis...\n", ont))
      
      # Upregulated genes
      if (length(gene_lists$up) > 0) {
        go_up <- run_go_analysis(gene_lists$up, "up", ont, gene_lists$background, pval_threshold)
        if (!is.null(go_up)) {
          results[[paste0("GO_", ont, "_UP")]] <- format_results(go_up, "GO")
        }
      }
      
      # Downregulated genes
      if (length(gene_lists$down) > 0) {
        go_down <- run_go_analysis(gene_lists$down, "down", ont, gene_lists$background, pval_threshold)
        if (!is.null(go_down)) {
          results[[paste0("GO_", ont, "_DOWN")]] <- format_results(go_down, "GO")
        }
      }
      
      # All significant genes
      if (length(gene_lists$all_sig) > 0) {
        go_all <- run_go_analysis(gene_lists$all_sig, "all", ont, gene_lists$background, pval_threshold)
        if (!is.null(go_all)) {
          results[[paste0("GO_", ont, "_ALL")]] <- format_results(go_all, "GO")
        }
      }
    }
  }
  
  if (analysis_type %in% c("KEGG", "ALL")) {
    cat("Running KEGG enrichment analysis...\n")
    
    # KEGG for significant genes
    if (length(gene_lists$all_sig) > 0) {
      kegg_result <- run_kegg_analysis(gene_lists$all_sig, gene_lists$background, pval_threshold)
      if (!is.null(kegg_result)) {
        results[["KEGG"]] <- format_results(kegg_result, "KEGG")
      }
    }
  }
  
  if (analysis_type %in% c("GSEA", "ALL")) {
    cat("Running GSEA analysis...\n")
    
    if (length(gene_lists$ranked) > 0) {
      for (pset in pathway_sets_to_analyze) {
        if (pset %in% c("H", "C2_KEGG", "C2_REACTOME", "C5_BP", "C5_MF", "C5_CC")) {
          cat(sprintf("  - GSEA %s analysis...\n", pset))
          gsea_result <- run_gsea_analysis(gene_lists$ranked, pset, 
                                           gsea_min_size, gsea_max_size)
          if (!is.null(gsea_result) && nrow(gsea_result) > 0) {
            results[[paste0("GSEA_", pset)]] <- format_results(gsea_result, "GSEA")
          }
        }
      }
    }
  }
  
  # Filter for specific pathway if requested
  if (!is.null(pathway)) {
    results <- lapply(results, function(x) {
      if (nrow(x) > 0) {
        x[grepl(pathway, x$pathway, ignore.case = TRUE) | 
            grepl(pathway, x$description, ignore.case = TRUE), ]
      } else {
        x
      }
    })
    
    # Remove empty results
    results <- results[sapply(results, nrow) > 0]
  }
  
  # Print summary
  cat("\n=== Analysis Summary ===\n")
  for (name in names(results)) {
    if (nrow(results[[name]]) > 0) {
      cat(sprintf("%s: %d significant pathways\n", name, nrow(results[[name]])))
    }
  }
  
  return(results)
}