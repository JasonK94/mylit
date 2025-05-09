
#' 가장 probability가 큰 2개 sample을 골라줌.
#' @export
get_best_two <- function(row) {
  sorted_indices <- order(row, decreasing = TRUE)
  best <- sorted_indices[1]
  second_best <- sorted_indices[2]
  return(c(best, second_best))
}

#' demultiplexing의 핵심.
#' @export
get_barcode_mapping = function(demux_data){
  best_two_indices <- t(apply(demux_data[,-1], 1, get_best_two))
  
  # Get the column names (samples) and probabilities
  best_samples <- names(demux_data)[-1][best_two_indices[,1]]
  second_best_samples <- names(demux_data)[-1][best_two_indices[,2]]
  best_probs <- apply(demux_data[,-1], 1, function(x) max(x))
  second_best_probs <- apply(demux_data[,-1], 1, function(x) sort(x, decreasing = TRUE)[2])
  
  # Calculate the ratio of best to second-best probabilities
  prob_ratio <- best_probs / second_best_probs
  
  # Create a data frame with barcode, best sample, best and second-best probabilities, and the ratio
  barcode_mapping <- data.frame(
    Barcode = demux_data$BARCODE,
    Best_Sample = best_samples,
    Best_Probability = best_probs,
    Second_Best_Sample = second_best_samples,
    Second_Best_Probability = second_best_probs,
    Probability_Ratio = prob_ratio
  )
  return(barcode_mapping)
}

#' doublet 판정
#' @export
is_doublet <- function(sample) {
  return(grepl("\\+", sample))
}

#' 1:6과 같이 숫자 벡터를 넣으면, 1+2, 2+3 형태의 조합을 만들어줌.
#' @importFrom utils combn
#' @export
generate_sample_values <- function(start_num, end_num) {
  singlets <- start_num : end_num
  doublets <- combn(singlets, 2, FUN = function(x) paste(x, collapse = "+"))
  return(c(as.character(singlets), doublets))
}

#' c("a","b","c")과 같이 문자 벡터를 넣으면, "a+b", "b+c" 형태의 조합을 만들어줌.
#' @importFrom utils combn
#' @export
generate_sample_names=function(vector){
  singlets=vector
  doublets=combn(vector, 2, FUN=function(x) paste(x,collapse="+"))
  return(c(singlets,doublets))
}