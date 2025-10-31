#' Get the Two Highest Probability Values
#'
#' This function takes a row of probability values and returns the indices of the two
#' highest values.
#'
#' @param row A numeric vector of probability values
#'
#' @return A numeric vector of length 2 containing the indices of the two highest
#'         probability values
#'
#' @examples
#' \dontrun{
#' # Get indices of two highest probabilities
#' probs <- c(0.1, 0.5, 0.3, 0.8, 0.2)
#' best_two <- get_best_two(probs)  # Returns c(4, 2)
#' }
#' @export
get_best_two <- function(row) {
  sorted_indices <- order(row, decreasing = TRUE)
  best <- sorted_indices[1]
  second_best <- sorted_indices[2]
  return(c(best, second_best))
}

#' Create Barcode Mapping from Demultiplexing Data
#'
#' This function processes demultiplexing data to create a mapping between barcodes
#' and their most likely sample assignments, including probability information.
#'
#' @param demux_data A data frame containing demultiplexing results. Must have a 'BARCODE'
#'                  column and additional columns for sample probabilities
#'
#' @return A data frame containing:
#'         - Barcode: The cell barcode
#'         - Best_Sample: The sample with highest probability
#'         - Best_Probability: The highest probability value
#'         - Second_Best_Sample: The sample with second highest probability
#'         - Second_Best_Probability: The second highest probability value
#'         - Probability_Ratio: Ratio of best to second-best probabilities
#'
#' @examples
#' \dontrun{
#' # Create barcode mapping from demultiplexing results
#' mapping <- get_barcode_mapping(demux_data)
#' 
#' # Filter for high confidence assignments
#' high_conf <- mapping[mapping$Probability_Ratio > 2, ]
#' }
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

#' Check if a Sample is a Doublet
#'
#' This function checks if a sample name contains a '+' character, indicating it is
#' a doublet (combination of two samples).
#'
#' @param sample Character string containing the sample name
#'
#' @return Logical value: TRUE if the sample is a doublet, FALSE otherwise
#'
#' @examples
#' \dontrun{
#' # Check if samples are doublets
#' is_doublet("1")     # Returns FALSE
#' is_doublet("1+2")   # Returns TRUE
#' }
#' @export
is_doublet <- function(sample) {
  return(grepl("\\+", sample))
}

#' Generate Sample Values for Demultiplexing
#'
#' This function generates a vector of sample values including both singlets and
#' doublets for demultiplexing analysis.
#'
#' @param start_num Integer specifying the starting sample number
#' @param end_num Integer specifying the ending sample number
#'
#' @return A character vector containing all possible sample combinations:
#'         singlets (e.g., "1", "2", "3") and doublets (e.g., "1+2", "1+3", "2+3")
#'
#' @examples
#' \dontrun{
#' # Generate sample values for samples 1-3
#' samples <- generate_sample_values(1, 3)
#' # Returns: c("1", "2", "3", "1+2", "1+3", "2+3")
#' }
#' @export
generate_sample_values <- function(start_num, end_num) {
  singlets <- start_num : end_num
  doublets <- combn(singlets, 2, FUN = function(x) paste(x, collapse = "+"))
  return(c(as.character(singlets), doublets))
}

#' Generate Sample Names for Demultiplexing
#'
#' This function generates a vector of sample names including both singlets and
#' doublets for demultiplexing analysis.
#'
#' @param vector Character vector containing sample names
#'
#' @return A character vector containing all possible sample combinations:
#'         singlets (original names) and doublets (combinations with '+')
#'
#' @examples
#' \dontrun{
#' # Generate sample names for samples A, B, C
#' samples <- generate_sample_names(c("A", "B", "C"))
#' # Returns: c("A", "B", "C", "A+B", "A+C", "B+C")
#' }
#' @export
generate_sample_names=function(vector){
  singlets=vector
  doublets=combn(vector, 2, FUN=function(x) paste(x,collapse="+"))
  return(c(singlets,doublets))
}