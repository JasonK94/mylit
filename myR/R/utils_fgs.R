#' Calculate Jaccard Similarity between Signatures
#'
#' @param sig_list A named list of gene signatures (character vectors).
#' @return A matrix of Jaccard indices.
#' @export
calculate_jaccard_similarity <- function(sig_list) {
    n <- length(sig_list)
    nms <- names(sig_list)
    mat <- matrix(0, nrow = n, ncol = n, dimnames = list(nms, nms))

    for (i in 1:n) {
        for (j in 1:n) {
            if (i == j) {
                mat[i, j] <- 1
            } else {
                s1 <- sig_list[[i]]
                s2 <- sig_list[[j]]
                u <- base::union(s1, s2)
                if (length(u) > 0) {
                    mat[i, j] <- length(base::intersect(s1, s2)) / length(u)
                } else {
                    mat[i, j] <- 0
                }
            }
        }
    }
    return(mat)
}

# compute_meta_gene_importance removed to avoid conflict with tml_utils.R
