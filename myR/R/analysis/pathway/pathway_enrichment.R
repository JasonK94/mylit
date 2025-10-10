#' Pathway Enrichment Analysis
#'
#' This module provides functions for Gene Ontology (GO), KEGG, and Gene Set
#' Enrichment Analysis (GSEA) using clusterProfiler and fgsea.
#'
#' @name pathway_enrichment
NULL

#' Convert Gene IDs Between Formats
#'
#' Converts gene identifiers between different formats (e.g., SYMBOL, ENTREZID, ENSEMBL).
#'
#' @param genes Character vector of gene identifiers
#' @param from Source ID type (default: "SYMBOL")
#' @param to Target ID type (default: "ENTREZID")
#' @param org_db Organism database (default: org.Hs.eg.db for human)
#' @param remove_na Remove genes that couldn't be converted (default: TRUE)
#'
#' @return Named character vector of converted IDs (names are original IDs)
#'
#' @examples
#' \dontrun{
#' entrez_ids <- convert_gene_ids(c("TP53", "EGFR", "KRAS"))
#' ensembl_ids <- convert_gene_ids(c("TP53", "EGFR"), to = "ENSEMBL")
#' }
#'
#' @export
convert_gene_ids <- function(genes,
                             from = "SYMBOL",
                             to = "ENTREZID",
                             org_db = org.Hs.eg.db::org.Hs.eg.db,
                             remove_na = TRUE) {
  
  if (length(genes) == 0) {
    return(character(0))
  }
  
  # Convert
  converted <- AnnotationDbi::mapIds(
    org_db,
    keys = genes,
    column = to,
    keytype = from,
    multiVals = "first"
  )
  
  if (remove_na) {
    converted <- converted[!is.na(converted)]
  }
  
  n_failed <- sum(is.na(converted))
  if (n_failed > 0) {
    message("Could not convert ", n_failed, "/", length(genes), " gene IDs")
  }
  
  return(converted)
}

#' Prepare Gene Lists for Pathway Analysis
#'
#' Prepares gene lists from differential expression results for pathway analysis.
#'
#' @param deg_results Data frame with DE results (must contain logFC and p-value columns)
#' @param gene_col Name of gene column (default: "gene" or rownames)
#' @param logfc_col Name of log fold change column (default: "logFC" or "avg_log2FC")
#' @param pval_col Name of p-value column (default: "pvalue" or "p_val")
#' @param padj_col Name of adjusted p-value column (default: "padj" or "p_val_adj")
#' @param padj_threshold Adjusted p-value threshold (default: 0.05)
#' @param logfc_threshold Log fold change threshold (default: 0)
#' @param convert_ids Convert gene IDs to ENTREZID (default: TRUE)
#'
#' @return List containing:
#'   \item{ranked_genes}{Named vector of all genes ranked by logFC}
#'   \item{up_genes}{Character vector of upregulated genes}
#'   \item{down_genes}{Character vector of downregulated genes}
#'   \item{sig_genes}{Character vector of all significant genes}
#'   \item{background_genes}{Character vector of all genes}
#'
#' @export
prepare_gene_lists <- function(deg_results,
                               gene_col = "gene",
                               logfc_col = NULL,
                               pval_col = NULL,
                               padj_col = NULL,
                               padj_threshold = 0.05,
                               logfc_threshold = 0,
                               convert_ids = TRUE) {
  
  # Identify columns
  if (!gene_col %in% colnames(deg_results) && !is.null(rownames(deg_results))) {
    deg_results[[gene_col]] <- rownames(deg_results)
  }
  
  if (is.null(logfc_col)) {
    logfc_col <- if ("logFC" %in% colnames(deg_results)) "logFC" else "avg_log2FC"
  }
  
  if (is.null(pval_col)) {
    pval_col <- if ("pvalue" %in% colnames(deg_results)) "pvalue" else "p_val"
  }
  
  if (is.null(padj_col)) {
    padj_col <- if ("padj" %in% colnames(deg_results)) "padj" else "p_val_adj"
  }
  
  # Check required columns exist
  required_cols <- c(gene_col, logfc_col, padj_col)
  missing_cols <- setdiff(required_cols, colnames(deg_results))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Convert IDs if requested
  if (convert_ids) {
    gene_map <- convert_gene_ids(deg_results[[gene_col]])
    deg_results$entrez <- gene_map[deg_results[[gene_col]]]
    deg_results <- deg_results[!is.na(deg_results$entrez), ]
    gene_id_col <- "entrez"
  } else {
    gene_id_col <- gene_col
  }
  
  # Ranked gene list (for GSEA)
  ranked_genes <- deg_results[[logfc_col]]
  names(ranked_genes) <- deg_results[[gene_id_col]]
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)
  
  # Significant genes
  sig_genes <- deg_results[
    deg_results[[padj_col]] < padj_threshold & 
    abs(deg_results[[logfc_col]]) > logfc_threshold,
    gene_id_col
  ]
  
  # Up/down regulated genes
  up_genes <- deg_results[
    deg_results[[padj_col]] < padj_threshold & 
    deg_results[[logfc_col]] > logfc_threshold,
    gene_id_col
  ]
  
  down_genes <- deg_results[
    deg_results[[padj_col]] < padj_threshold & 
    deg_results[[logfc_col]] < -logfc_threshold,
    gene_id_col
  ]
  
  # Background genes
  background_genes <- deg_results[[gene_id_col]]
  
  message("Prepared gene lists:")
  message("  Ranked: ", length(ranked_genes), " genes")
  message("  Significant: ", length(sig_genes), " genes")
  message("  Upregulated: ", length(up_genes), " genes")
  message("  Downregulated: ", length(down_genes), " genes")
  
  return(list(
    ranked_genes = ranked_genes,
    up_genes = up_genes,
    down_genes = down_genes,
    sig_genes = sig_genes,
    background_genes = background_genes
  ))
}

#' Get Available Pathway Databases
#'
#' Returns available pathway database identifiers for GSEA.
#'
#' @param species Species code (default: "Homo sapiens")
#' @param category MSigDB category (e.g., "H", "C2", "C5")
#' @param subcategory MSigDB subcategory (e.g., "CP:KEGG", "GO:BP")
#'
#' @return Character vector of available gene set names
#'
#' @export
get_pathway_sets <- function(species = "Homo sapiens",
                             category = NULL,
                             subcategory = NULL) {
  
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    stop("msigdbr package is required. Install with: install.packages('msigdbr')")
  }
  
  # Get all gene sets
  genesets <- msigdbr::msigdbr(species = species)
  
  # Filter by category and subcategory if specified
  if (!is.null(category)) {
    genesets <- genesets[genesets$gs_cat == category, ]
  }
  
  if (!is.null(subcategory)) {
    genesets <- genesets[genesets$gs_subcat == subcategory, ]
  }
  
  unique(genesets$gs_name)
}

#' Run GO Enrichment Analysis
#'
#' Performs Gene Ontology (GO) enrichment analysis using clusterProfiler.
#'
#' @param genes Character vector of gene IDs (ENTREZID format)
#' @param universe Background gene universe (optional)
#' @param ont GO ontology: "BP", "MF", "CC", or "ALL" (default: "BP")
#' @param pval_cutoff P-value cutoff (default: 0.05)
#' @param qval_cutoff Q-value cutoff (default: 0.05)
#' @param org_db Organism database (default: org.Hs.eg.db)
#'
#' @return enrichResult object from clusterProfiler
#'
#' @export
run_go_analysis <- function(genes,
                            universe = NULL,
                            ont = "BP",
                            pval_cutoff = 0.05,
                            qval_cutoff = 0.05,
                            org_db = org.Hs.eg.db::org.Hs.eg.db) {
  
  if (length(genes) == 0) {
    stop("No genes provided")
  }
  
  ego <- clusterProfiler::enrichGO(
    gene = genes,
    universe = universe,
    OrgDb = org_db,
    ont = ont,
    pAdjustMethod = "BH",
    pvalueCutoff = pval_cutoff,
    qvalueCutoff = qval_cutoff,
    readable = TRUE
  )
  
  if (is.null(ego) || nrow(ego@result) == 0) {
    message("No significant GO terms found")
    return(NULL)
  }
  
  message("Found ", nrow(ego@result), " enriched GO terms (", ont, ")")
  
  return(ego)
}

#' Run KEGG Pathway Enrichment Analysis
#'
#' Performs KEGG pathway enrichment analysis using clusterProfiler.
#'
#' @param genes Character vector of gene IDs (ENTREZID format)
#' @param universe Background gene universe (optional)
#' @param organism KEGG organism code (default: "hsa" for human)
#' @param pval_cutoff P-value cutoff (default: 0.05)
#' @param qval_cutoff Q-value cutoff (default: 0.05)
#'
#' @return enrichResult object from clusterProfiler
#'
#' @export
run_kegg_analysis <- function(genes,
                              universe = NULL,
                              organism = "hsa",
                              pval_cutoff = 0.05,
                              qval_cutoff = 0.05) {
  
  if (length(genes) == 0) {
    stop("No genes provided")
  }
  
  ekegg <- clusterProfiler::enrichKEGG(
    gene = genes,
    universe = universe,
    organism = organism,
    pAdjustMethod = "BH",
    pvalueCutoff = pval_cutoff,
    qvalueCutoff = qval_cutoff
  )
  
  if (is.null(ekegg) || nrow(ekegg@result) == 0) {
    message("No significant KEGG pathways found")
    return(NULL)
  }
  
  message("Found ", nrow(ekegg@result), " enriched KEGG pathways")
  
  return(ekegg)
}

#' Run Gene Set Enrichment Analysis (GSEA)
#'
#' Performs GSEA using fgsea with MSigDB gene sets.
#'
#' @param ranked_genes Named numeric vector of genes ranked by score (e.g., logFC)
#' @param species Species for MSigDB (default: "Homo sapiens")
#' @param category MSigDB category (default: "H" for Hallmark)
#' @param subcategory MSigDB subcategory (optional)
#' @param min_size Minimum gene set size (default: 15)
#' @param max_size Maximum gene set size (default: 500)
#' @param nperm Number of permutations (default: 10000)
#'
#' @return Data frame with GSEA results
#'
#' @export
run_gsea_analysis <- function(ranked_genes,
                              species = "Homo sapiens",
                              category = "H",
                              subcategory = NULL,
                              min_size = 15,
                              max_size = 500,
                              nperm = 10000) {
  
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    stop("msigdbr package is required")
  }
  
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("fgsea package is required")
  }
  
  # Get gene sets
  msig_db <- msigdbr::msigdbr(species = species)
  
  if (!is.null(category)) {
    msig_db <- msig_db[msig_db$gs_cat == category, ]
  }
  
  if (!is.null(subcategory)) {
    msig_db <- msig_db[msig_db$gs_subcat == subcategory, ]
  }
  
  # Convert to list format
  pathways <- split(msig_db$entrez_gene, msig_db$gs_name)
  
  # Remove duplicates in ranked genes
  ranked_genes <- ranked_genes[!duplicated(names(ranked_genes))]
  
  # Run GSEA
  message("Running GSEA with ", length(pathways), " gene sets...")
  
  fgsea_results <- fgsea::fgsea(
    pathways = pathways,
    stats = ranked_genes,
    minSize = min_size,
    maxSize = max_size,
    nperm = nperm
  )
  
  # Sort by p-value
  fgsea_results <- fgsea_results[order(fgsea_results$pval), ]
  
  n_sig <- sum(fgsea_results$padj < 0.05, na.rm = TRUE)
  message("Found ", n_sig, " significant gene sets (padj < 0.05)")
  
  return(as.data.frame(fgsea_results))
}

#' Format Pathway Analysis Results
#'
#' Standardizes output from different pathway analysis methods.
#'
#' @param results Results object from GO, KEGG, or GSEA analysis
#' @param method Method used: "GO", "KEGG", or "GSEA"
#'
#' @return Standardized data frame
#'
#' @export
format_results <- function(results, method = c("GO", "KEGG", "GSEA")) {
  
  method <- match.arg(method)
  
  if (is.null(results)) {
    return(data.frame())
  }
  
  if (method %in% c("GO", "KEGG")) {
    # clusterProfiler results
    df <- as.data.frame(results)
    
    std_df <- data.frame(
      pathway = df$ID,
      description = df$Description,
      gene_ratio = df$GeneRatio,
      bg_ratio = df$BgRatio,
      pvalue = df$pvalue,
      padj = df$p.adjust,
      genes = df$geneID,
      count = df$Count,
      stringsAsFactors = FALSE
    )
    
  } else if (method == "GSEA") {
    # fgsea results
    std_df <- data.frame(
      pathway = results$pathway,
      pvalue = results$pval,
      padj = results$padj,
      ES = results$ES,
      NES = results$NES,
      size = results$size,
      leading_edge = sapply(results$leadingEdge, paste, collapse = ","),
      stringsAsFactors = FALSE
    )
  }
  
  return(std_df)
}

#' Comprehensive Pathway Analysis Wrapper
#'
#' Runs GO, KEGG, and GSEA analyses and returns combined results.
#'
#' @param deg_results Data frame with differential expression results
#' @param run_go Perform GO analysis (default: TRUE)
#' @param run_kegg Perform KEGG analysis (default: TRUE)
#' @param run_gsea Perform GSEA (default: TRUE)
#' @param go_ont GO ontology (default: "BP")
#' @param gsea_category MSigDB category for GSEA (default: "H")
#' @param padj_threshold Adjusted p-value threshold for gene selection (default: 0.05)
#' @param logfc_threshold Log fold change threshold (default: 0)
#'
#' @return List containing results from each analysis method
#'
#' @examples
#' \dontrun{
#' pathway_results <- myGO(
#'   deg_results,
#'   run_go = TRUE,
#'   run_kegg = TRUE,
#'   run_gsea = TRUE
#' )
#' }
#'
#' @export
myGO <- function(deg_results,
                 run_go = TRUE,
                 run_kegg = TRUE,
                 run_gsea = TRUE,
                 go_ont = "BP",
                 gsea_category = "H",
                 padj_threshold = 0.05,
                 logfc_threshold = 0) {
  
  message("=== Pathway Analysis Pipeline ===\n")
  
  # Prepare gene lists
  message("1. Preparing gene lists...")
  gene_lists <- prepare_gene_lists(
    deg_results,
    padj_threshold = padj_threshold,
    logfc_threshold = logfc_threshold,
    convert_ids = TRUE
  )
  
  results <- list()
  
  # GO Analysis
  if (run_go && length(gene_lists$sig_genes) > 0) {
    message("\n2. Running GO enrichment...")
    tryCatch({
      results$GO <- run_go_analysis(
        genes = gene_lists$sig_genes,
        universe = gene_lists$background_genes,
        ont = go_ont
      )
      
      if (!is.null(results$GO)) {
        results$GO_formatted <- format_results(results$GO, "GO")
      }
    }, error = function(e) {
      message("GO analysis failed: ", e$message)
      results$GO <<- NULL
    })
  }
  
  # KEGG Analysis
  if (run_kegg && length(gene_lists$sig_genes) > 0) {
    message("\n3. Running KEGG enrichment...")
    tryCatch({
      results$KEGG <- run_kegg_analysis(
        genes = gene_lists$sig_genes,
        universe = gene_lists$background_genes
      )
      
      if (!is.null(results$KEGG)) {
        results$KEGG_formatted <- format_results(results$KEGG, "KEGG")
      }
    }, error = function(e) {
      message("KEGG analysis failed: ", e$message)
      results$KEGG <<- NULL
    })
  }
  
  # GSEA
  if (run_gsea && length(gene_lists$ranked_genes) > 0) {
    message("\n4. Running GSEA...")
    tryCatch({
      results$GSEA <- run_gsea_analysis(
        ranked_genes = gene_lists$ranked_genes,
        category = gsea_category
      )
      
      if (!is.null(results$GSEA) && nrow(results$GSEA) > 0) {
        results$GSEA_formatted <- format_results(results$GSEA, "GSEA")
      }
    }, error = function(e) {
      message("GSEA failed: ", e$message)
      results$GSEA <<- NULL
    })
  }
  
  message("\n=== Analysis Complete ===")
  
  results$gene_lists <- gene_lists
  
  return(results)
}

