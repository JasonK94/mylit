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