# ==============================================================================
# TML Post-processing Utilities
# ==============================================================================

compute_meta_gene_importance <- function(meta_result,
                                         normalize = TRUE,
                                         normalization_method = "max_abs",
                                         target_model = NULL) {
  if (!is.list(meta_result) ||
    !all(c("best_model", "best_model_name", "l1_signatures", "l2_train") %in% names(meta_result))) {
    stop("meta_result must be the object returned by TML7().")
  }

  `%||%` <- function(x, y) if (!is.null(x)) x else y

  # Select model
  if (!is.null(target_model)) {
    if (!"trained_models" %in% names(meta_result)) {
      stop("trained_models not found in meta_result. Cannot use target_model.")
    }
    if (!target_model %in% names(meta_result$trained_models)) {
      available_models <- names(meta_result$trained_models)
      stop(sprintf(
        "target_model '%s' not found in trained_models. Available: %s",
        target_model, paste(available_models, collapse = ", ")
      ))
    }
    model_type <- target_model
    model <- meta_result$trained_models[[target_model]]
    if (is.null(model)) {
      stop(sprintf("Model '%s' in trained_models is NULL.", target_model))
    }
  } else {
    model_type <- meta_result$best_model_name
    model <- meta_result$best_model
  }

  # Get signature names from l2_train (excluding target)
  sig_names <- base::setdiff(colnames(meta_result$l2_train), ".target")

  if (length(sig_names) == 0) {
    stop("No signature features found in l2_train.")
  }

  # Map sig_names back to original l1_signatures names
  l1_sig_names <- names(meta_result$l1_signatures)
  l1_sig_names_made <- make.names(l1_sig_names)
  name_mapping <- setNames(l1_sig_names, l1_sig_names_made)

  # Check if sig_names are already transformed or original
  direct_match <- sig_names %in% l1_sig_names
  if (all(direct_match)) {
    sig_names_for_l1 <- sig_names
  } else {
    original_sig_names <- name_mapping[sig_names]
    missing_idx <- is.na(original_sig_names)
    if (any(missing_idx)) {
      matches_direct <- sig_names[missing_idx] %in% l1_sig_names
      original_sig_names[missing_idx][matches_direct] <- sig_names[missing_idx][matches_direct]
    }

    if (any(is.na(original_sig_names))) {
      warning(sprintf(
        "Some l2_train columns could not be mapped to l1_signatures: %s",
        paste(sig_names[is.na(original_sig_names)], collapse = ", ")
      ))
      original_sig_names <- original_sig_names[!is.na(original_sig_names)]
    }

    if (length(original_sig_names) == 0) {
      stop(sprintf("No matching signatures found between l2_train columns and l1_signatures."))
    }
    sig_names_for_l1 <- as.character(original_sig_names)
  }

  # Helper to standardise signature definition
  standardise_signature <- function(sig) {
    if (is.null(sig)) {
      return(numeric(0))
    }
    if (is.numeric(sig) && !is.null(names(sig))) {
      return(sig)
    }
    if (is.character(sig)) {
      return(structure(rep(1, length(sig)), names = sig))
    }
    if (is.list(sig) && all(c("up", "down") %in% names(sig))) {
      up <- sig$up %||% character()
      down <- sig$down %||% character()
      nm <- c(up, down)
      wt <- c(rep(1, length(up)), rep(-1, length(down)))
      names(wt) <- nm
      return(wt)
    }
    if (is.list(sig) && all(c("genes", "weights") %in% names(sig))) {
      g <- sig$genes
      w <- sig$weights
      names(w) <- g
      return(w)
    }
    if (is.data.frame(sig)) {
      cn <- tolower(colnames(sig))
      gene_col <- which(cn %in% c("gene", "genes", "feature", "features", "symbol", "id"))[1]
      weight_col <- which(cn %in% c("weight", "weights", "w", "score", "scores", "coef", "coefs"))[1]
      if (!is.na(gene_col) && !is.na(weight_col)) {
        g <- as.character(sig[[gene_col]])
        w <- as.numeric(sig[[weight_col]])
        ok <- !is.na(g) & !is.na(w)
        names(w[ok]) <- g[ok]
        return(w[ok])
      }
    }
    stop("Unsupported signature format.")
  }

  signature_weights <- lapply(meta_result$l1_signatures[sig_names_for_l1], standardise_signature)
  names(signature_weights) <- sig_names_for_l1

  # Obtain signature-level importance
  signature_importance <- NULL

  if (identical(model_type, "glm") && inherits(model$finalModel, "glm")) {
    coefs <- stats::coef(model$finalModel)
    coefs <- coefs[names(coefs) != "(Intercept)"]
    available_sigs <- base::intersect(sig_names, names(coefs))
    if (length(available_sigs) == 0) {
      warning("No matching signatures found in glm coefficients. Returning NULL.")
      return(NULL)
    }
    signature_importance <- coefs[available_sigs]

    # Map names
    mapped_imp_names <- character(length(available_sigs))
    for (i in seq_along(available_sigs)) {
      s <- available_sigs[i]
      if (s %in% l1_sig_names) {
        mapped_imp_names[i] <- s
      } else if (s %in% names(name_mapping)) {
        mapped_imp_names[i] <- name_mapping[[s]]
      } else {
        mapped_imp_names[i] <- NA
      }
    }
    valid <- !is.na(mapped_imp_names)
    signature_importance <- signature_importance[valid]
    names(signature_importance) <- mapped_imp_names[valid]
  } else {
    # Generic importance extraction
    imp_vals <- NULL
    if (identical(model_type, "ranger") && inherits(model$finalModel, "ranger")) {
      imp_vals <- try(ranger::importance(model$finalModel), silent = TRUE)
      if (inherits(imp_vals, "try-error")) imp_vals <- NULL
    }

    if (is.null(imp_vals)) {
      vi <- try(caret::varImp(model, scale = FALSE), silent = TRUE)
      if (!inherits(vi, "try-error") && !is.null(vi)) {
        if (is.data.frame(vi$importance)) {
          if ("Overall" %in% colnames(vi$importance)) {
            imp_vals <- vi$importance$Overall
            names(imp_vals) <- rownames(vi$importance)
          } else if (ncol(vi$importance) >= 1) {
            imp_vals <- vi$importance[, 1]
            names(imp_vals) <- rownames(vi$importance)
          }
        }
      }
    }

    if (is.null(imp_vals)) {
      warning(sprintf("Could not extract importance for model '%s'. Returning NULL.", model_type))
      return(NULL)
    }

    # Map importance names
    imp_names <- names(imp_vals)
    matched_l1_names <- character(length(imp_names))
    for (i in seq_along(imp_names)) {
      s <- imp_names[i]
      if (s %in% l1_sig_names) {
        matched_l1_names[i] <- s
      } else if (s %in% names(name_mapping)) {
        matched_l1_names[i] <- name_mapping[[s]]
      } else {
        s_clean <- gsub("`", "", s)
        if (s_clean %in% l1_sig_names) {
          matched_l1_names[i] <- s_clean
        } else if (s_clean %in% names(name_mapping)) {
          matched_l1_names[i] <- name_mapping[[s_clean]]
        } else {
          matched_l1_names[i] <- NA
        }
      }
    }
    valid <- !is.na(matched_l1_names)
    signature_importance <- imp_vals[valid]
    names(signature_importance) <- matched_l1_names[valid]
  }

  if (is.null(signature_importance) || length(signature_importance) == 0) {
    warning(sprintf("No matching signatures found in importance. Returning NULL."))
    return(NULL)
  }

  # === Normalization Logic ===
  if (normalize) {
    vals <- signature_importance

    if (normalization_method == "max_abs") {
      denom <- max(abs(vals), na.rm = TRUE) %||% 1
      if (denom > 0) signature_importance <- vals / denom
    } else if (normalization_method == "min_max") {
      min_val <- min(vals, na.rm = TRUE)
      max_val <- max(vals, na.rm = TRUE)
      if (max_val > min_val) {
        signature_importance <- (vals - min_val) / (max_val - min_val)
      } else {
        signature_importance <- rep(0, length(vals)) # or 0.5?
      }
    } else if (normalization_method == "softmax") {
      # Softmax: exp(x) / sum(exp(x))
      # For numerical stability, subtract max
      exp_vals <- exp(vals - max(vals, na.rm = TRUE))
      signature_importance <- exp_vals / sum(exp_vals, na.rm = TRUE)
    } else if (normalization_method == "rank") {
      signature_importance <- rank(vals, ties.method = "average")
      # Normalize rank to [0, 1]
      signature_importance <- signature_importance / length(signature_importance)
    } else if (normalization_method == "z_score") {
      sd_val <- stats::sd(vals, na.rm = TRUE)
      if (!is.na(sd_val) && sd_val > 0) {
        signature_importance <- (vals - mean(vals, na.rm = TRUE)) / sd_val
      }
    } else {
      warning(sprintf("Unknown normalization method '%s'. Using max_abs.", normalization_method))
      denom <- max(abs(vals), na.rm = TRUE) %||% 1
      if (denom > 0) signature_importance <- vals / denom
    }
  }

  # Compute gene contributions
  gene_tables <- lapply(names(signature_importance), function(sig) {
    if (!sig %in% names(signature_weights)) {
      return(NULL)
    }
    weights <- signature_weights[[sig]]
    if (length(weights) == 0) {
      return(NULL)
    }

    sig_imp <- signature_importance[[sig]]
    if (is.na(sig_imp) || !is.finite(sig_imp)) {
      return(NULL)
    }

    contrib <- sig_imp * weights

    # Normalize gene contributions as well?
    # Usually we want raw contribution sum, but maybe normalize per signature?
    # The original code normalized here too if normalize=TRUE
    # "if (normalize && any(is.finite(contrib))) { denom <- max(abs(contrib)) ... }"
    # This seems redundant if signature_importance is already normalized.
    # But it ensures that within a signature, the max gene contribution is scaled.
    # Let's keep it consistent with original behavior for "max_abs", but maybe skip for others?
    # Actually, if we use softmax for signature importance, we probably want to preserve the relative magnitude of genes.
    # So we should ONLY normalize signature importance, NOT gene weights again, unless necessary.

    # Original code logic:
    # contrib <- sig_imp * weights
    # if (normalize) contrib <- contrib / max(abs(contrib))
    # This makes the max contribution of ANY gene in THIS signature equal to 1 (or -1).
    # This effectively ignores sig_imp magnitude! Wait.
    # If contrib = sig_imp * weights, and we divide by max(abs(contrib)) = |sig_imp| * max(abs(weights))
    # Then contrib becomes weights / max(abs(weights)) * sign(sig_imp).
    # So signature importance magnitude is LOST! Only sign matters.
    # This seems like a BUG in the original code if the intention was to use L2 importance.

    # Let's check original code again.
    # contrib <- sig_imp * weights
    # if (normalize) ... contrib <- contrib / denom
    # Yes, it normalizes the result.

    # If we want to use L2 importance, we should NOT normalize again here, or normalize differently.
    # Let's modify this to respect signature_importance.

    # New logic:
    # We trust signature_importance is already scaled appropriately.
    # We just multiply by gene weights.
    # But gene weights might have different scales across signatures.
    # standardise_signature returns weights. If they are raw coefficients, they might vary.
    # If we want "contribution", it should be (Signature Importance) * (Gene Weight).

    data.frame(
      signature = sig,
      gene = names(contrib),
      contribution = as.numeric(contrib),
      abs_contribution = abs(as.numeric(contrib)),
      stringsAsFactors = FALSE
    )
  })

  gene_importance <- do.call(rbind, gene_tables)
  if (is.null(gene_importance) || nrow(gene_importance) == 0) {
    warning("No gene contributions could be computed. Returning NULL.")
    return(NULL)
  }

  gene_summary <- stats::aggregate(contribution ~ gene, data = gene_importance, sum)
  gene_summary$abs_contribution <- abs(gene_summary$contribution)

  positive_class <- meta_result$positive_class %||%
    tail(levels(meta_result$l2_train$.target), 1)

  list(
    signature_importance = signature_importance,
    gene_importance = gene_importance,
    gene_summary = gene_summary[order(-gene_summary$abs_contribution), ],
    positive_class = positive_class,
    model_type = model_type,
    normalization_method = normalization_method
  )
}
add_meta_signature_score <- function(
  seurat_obj,
  gene_importance_result,
  signature_name = "meta_signature_score",
  assay = NULL,
  layer = "data",
  normalize = TRUE,
  use_abs_contribution = FALSE
) {
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object.")
  }

  if (!is.list(gene_importance_result) ||
    !"gene_summary" %in% names(gene_importance_result)) {
    stop("gene_importance_result must be the output of compute_meta_gene_importance().")
  }

  gene_summary <- gene_importance_result$gene_summary
  if (!all(c("gene", "contribution") %in% colnames(gene_summary))) {
    stop("gene_summary must contain 'gene' and 'contribution' columns.")
  }

  # Get assay
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(seurat_obj)
  }

  # Get expression matrix
  expr_mat <- Seurat::GetAssayData(seurat_obj, assay = assay, layer = layer)
  if (is.null(expr_mat)) {
    stop(sprintf("Layer '%s' not found in assay '%s'.", layer, assay))
  }

  # Prepare gene weights
  if (use_abs_contribution) {
    if (!"abs_contribution" %in% colnames(gene_summary)) {
      gene_summary$abs_contribution <- abs(gene_summary$contribution)
    }
    weights <- gene_summary$abs_contribution
  } else {
    weights <- gene_summary$contribution
  }
  names(weights) <- gene_summary$gene

  # Find overlapping genes
  available_genes <- base::intersect(names(weights), rownames(expr_mat))
  if (length(available_genes) == 0) {
    stop("No overlapping genes found between signature and expression data.")
  }
  if (length(available_genes) < length(weights)) {
    warning(sprintf(
      "Only %d/%d signature genes found in expression data.",
      length(available_genes), length(weights)
    ))
  }

  # Subset weights to available genes
  w <- weights[available_genes]

  # Calculate signature scores
  # Weighted sum: sum(w * expr) / sum(|w|)
  expr_subset <- expr_mat[available_genes, , drop = FALSE]
  if (inherits(expr_subset, "sparseMatrix")) {
    # For sparse matrices, use Matrix operations
    scores <- as.numeric(Matrix::t(expr_subset) %*% w)
  } else {
    # For dense matrices
    scores <- as.numeric(colSums(expr_subset * w))
  }

  # Normalize by sum of absolute weights
  den <- sum(abs(w))
  if (!is.finite(den) || den == 0) {
    warning("Sum of absolute weights is zero or non-finite. Using unnormalized scores.")
    den <- 1
  }
  scores <- scores / den

  # Z-score normalization if requested
  if (normalize) {
    scores_sd <- stats::sd(scores, na.rm = TRUE)
    if (is.finite(scores_sd) && scores_sd > 0) {
      scores <- as.numeric(scale(scores))
    } else {
      warning("Cannot normalize scores (zero variance). Using unnormalized scores.")
    }
  }

  # Add to metadata
  seurat_obj <- Seurat::AddMetaData(
    object = seurat_obj,
    metadata = scores,
    col.name = signature_name
  )

  message(sprintf(
    "Added signature score '%s' to metadata (%d genes used).",
    signature_name, length(available_genes)
  ))

  return(seurat_obj)
}
