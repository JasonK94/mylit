#' Downsample a Seurat object by a specified ratio
#'
#' This function randomly samples cells from a Seurat object to create a smaller subset.
#' The sampling is done without replacement and maintains the original data structure.
#'
#' @param sobj A Seurat object to be downsampled
#' @param ratio Integer indicating the downsampling ratio (1:ratio). Default is 10, meaning
#'              the output will contain 1/10th of the original cells
#' @param seed Integer specifying the random seed for reproducibility. Default is 1234
#'
#' @return A downsampled Seurat object containing a subset of cells from the input object
#'
#' @examples
#' \dontrun{
#' # Downsample a Seurat object to 1/10th of its original size
#' sobj_down <- downsample_sobj(sobj, ratio = 10)
#' 
#' # Downsample with a different ratio and seed
#' sobj_down <- downsample_sobj(sobj, ratio = 5, seed = 42)
#' }
#' @export
downsample_sobj=function(sobj, ratio=10,seed=1234){
  set.seed(seed)
  barcodes=colnames(sobj)
  barcodes_down=sample(barcodes,length(barcodes)%/%10)
  sobj_down=sobj[,barcodes_down]
  return(sobj_down)
}


printmy=function(markers, sign="+",num=100, pseudobulk=FALSE){
  if(!"gene"%in%names(markers)){
    markers$gene=rownames(markers)
  }
  if(!pseudobulk){
    if(sign=="-"){
      print(paste(markers[markers$avg_log2FC<0,][1:num,]$gene,collapse = ", "))
    }else{
      print(paste(markers[markers$avg_log2FC>0,][1:num,]$gene,collapse = ", "))
    } 
  }else{
    if(sign=="-"){
      print(paste(markers[markers$logFC<0,][1:num,]$gene,collapse = ", "))
    }else{
      print(paste(markers[markers$logFC>0,][1:num,]$gene,collapse = ", "))
    } 
  }
}

printMy=function(markers_list, ...){
  args=list(...)
  for(name in names(markers_list)){
    print(name)
    printmy(markers_list[[name]], pseudobulk=args[["pseudobulk"]])
  }
}

#' helper for cmb, acmb
sort_samples <- function(samples) {
  # Helper function to extract numbers from a string
  extract_numbers <- function(x) {
    nums <- as.numeric(strsplit(x, "\\+")[[1]])
    if (length(nums) == 0 || any(is.na(nums))) return(c(Inf, Inf))
    return(nums)
  }
  
  # Check if all samples are non-numeric (no '+' and can't be converted to number)
  all_non_numeric <- all(sapply(samples, function(x) {
    !grepl("\\+", x) && is.na(suppressWarnings(as.numeric(x)))
  }))
  
  if (all_non_numeric) {
    return(sort(samples))  # If all non-numeric, just return alphabetically sorted
  }
  
  # Separate samples into those with and without "+"
  single_samples <- samples[!grepl("\\+", samples)]
  doublet_samples <- samples[grepl("\\+", samples)]
  
  # Sort single samples
  single_sorted <- single_samples[order(sapply(single_samples, function(x) {
    num <- suppressWarnings(as.numeric(x))
    if (is.na(num)) Inf else num
  }))]
  
  # Sort doublet samples
  doublet_sorted <- doublet_samples[order(
    sapply(doublet_samples, function(x) extract_numbers(x)[1]),  # First number
    sapply(doublet_samples, function(x) extract_numbers(x)[2])   # Second number
  )]
  
  # Combine and return
  return(c(single_sorted, doublet_sorted))
}

## Helper: Sort Sample Identifiers
##
## @description
## Sorts character sample identifiers, handling single-number and dual-number strings (e.g., "1+2").
## @param samples Character vector of sample identifiers.
## @return Sorted character vector of sample IDs.
## @keywords internal
sort_samples <- function(samples) {
  extract_numbers <- function(x) {
    parts <- strsplit(x, "\\+")[[1]]
    nums  <- as.numeric(parts)
    if (any(is.na(nums))) return(c(Inf, Inf))
    nums
  }
  
  is_double <- grepl("\\+", samples)
  single  <- samples[!is_double]
  double  <- samples[ is_double]
  
  single_sorted <- single[order(suppressWarnings(as.numeric(single)), single)]
  double_sorted <- double[order(sapply(double, function(x) extract_numbers(x)[1]),
                                sapply(double, function(x) extract_numbers(x)[2]))]
  c(single_sorted, double_sorted)
}


#' Print Top Genes per Set Combination
#'
#' @param gene_list Named list of character vectors. 각 원소가 유전자 벡터입니다.
#' @param num_print Integer. 각 조합별로 출력할 최대 유전자 수 (기본 100).
#' @examples
#' L <- list(
#'   A = c("TP53","EGFR","MYC","BRCA1"),
#'   B = c("EGFR","MYC","PTEN","BRCA2"),
#'   C = c("MYC","PTEN","ALK")
#' )
#' print_gene_combinations(L, num_print = 2)
#' @export
print_gene_combinations <- function(gene_list, num_print = 100) {
  if (!is.list(gene_list) || is.null(names(gene_list))) {
    stop("gene_list는 반드시 이름이 있는 list여야 합니다.")
  }
  set_names <- names(gene_list)
  n_sets <- length(set_names)
  
  cat(sprintf("Printing top %d genes per combination...\n\n", num_print))
  
  # 1) 각 k-조합에 대해서
  for (k in seq_len(n_sets)) {
    combos <- combn(set_names, k, simplify = FALSE)
    for (comb in combos) {
      # 조합 이름
      comb_name <- paste(comb, collapse = "&")
      
      # 2) 조합 내 모든 리스트에 공통인 유전자
      genes_in <- Reduce(intersect, gene_list[comb])
      
      # 3) 그 밖의 리스트에 있는 유전자는 제외
      other_sets <- setdiff(set_names, comb)
      if (length(other_sets) > 0) {
        genes_out <- unique(unlist(gene_list[other_sets], use.names = FALSE))
        genes_in <- setdiff(genes_in, genes_out)
      }
      
      # 4) 정렬 및 샘플링
      genes_in <- sort(genes_in)
      to_print <- head(genes_in, num_print)
      
      # 5) 출력
      cat(comb_name, ":\n")
      if (length(to_print) == 0) {
        cat("  (none)\n\n")
      } else {
        cat(" ", paste(to_print, collapse = ", "), "\n\n")
      }
    }
  }
}