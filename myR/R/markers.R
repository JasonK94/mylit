#' Filter and process marker genes from Seurat's FindMarkers or FindAllMarkers results
#'
#' This function processes marker gene results by filtering based on log fold change direction,
#' adjusted p-value cutoff, and adding percentage difference information.
#'
#' @param markers A data frame from Seurat's FindMarkers or FindAllMarkers results
#' @param sign Character string indicating the direction of log fold change to keep.
#'             Options: NULL (keep all), "+" (positive only), "-" (negative only)
#' @param p_cutoff Numeric value for adjusted p-value cutoff. Default is NULL (no filtering)
#' @param filter Additional filtering criteria (not used in this function)
#'
#' @return A processed data frame with filtered markers and added pct.diff column
#'
#' @examples
#' \dontrun{
#' # Filter markers to keep only positive log fold changes
#' markers_positive <- marker_trim(markers, sign = "+")
#' 
#' # Filter markers with adjusted p-value < 0.05
#' markers_sig <- marker_trim(markers, p_cutoff = 0.05)
#' }
#' @export
marker_trim=function(markers, sign=NULL, p_cutoff=NULL, filter=NULL){
  if("cluster"%in%names(markers)){}else{ #FindAllMarkers object has "cluster" column, FindMarkers not.
    markers$gene=rownames(markers)
  }
  markers$pct.diff=markers$pct.1-markers$pct.2
  if(is.null(sign)){
    
  }else if(sign=="+"){
    markers=markers[markers$avg_log2FC>0,]
  }else if(sign=="-"){
    markers=markers[markers$avg_log2FC<0,]
  }else{
    print("sign is wrong. put NULL or "+" or "-" ")
  }
  
  if(is.null(p_cutoff)){
    
  }else{
    markers=markers[markers$p_val_adj<p_cutoff,]
  }
  
  return(markers)
}

#' Filter out unwanted genes from marker results
#'
#' This function removes genes matching specific patterns (ribosomal, mitochondrial,
#' hemoglobin, etc.) from marker gene results.
#'
#' @param markers A data frame from Seurat's FindMarkers or FindAllMarkers results
#' @param filter Character vector specifying which gene types to filter out.
#'               Options: "rb" (ribosomal), "mt" (mitochondrial), "hb" (hemoglobin),
#'               "AC" (AC/AL genes), "ENSG" (ENSG genes), "LINC" (LINC genes)
#'
#' @return A filtered data frame with unwanted genes removed
#'
#' @examples
#' \dontrun{
#' # Remove ribosomal and mitochondrial genes
#' markers_filtered <- marker_filter(markers, filter = c("rb", "mt"))
#' 
#' # Remove all unwanted gene types
#' markers_clean <- marker_filter(markers, filter = c("rb", "mt", "hb", "AC", "ENSG", "LINC"))
#' }
#' @export
marker_filter=function(markers, filter=c("rb","mt","hb","AC","ENSG","LINC")){
  rb=mt=AC=ENSG=LINC=""
  if("rb"%in%filter){
    rb="^RPL|^RPS"
  }
  if("mt"%in%filter){
    mt="^MT-"
  }
  if("AC"%in%filter){
    AC="^(AC\\d+|AL\\d+)"
  }
  if("ENSG"%in%filter){
    ENSG="^ENSG"
  }
  if("LINC"%in%filter){
    LINC="^LINC"
  }
  filter_pattern=paste(rb,mt,AC,ENSG,LINC,sep="|")
  if("gene"%in%names(markers)){}else{ # FindAllMarkers object generates "gene" column automatically
    markers$gene=rownames(markers)
  }
  markers=markers[!grepl(filter_pattern,markers$gene),]
  
  if("hb"%in%filter){
    globin_genes <- c(
      "HBA1", "HBA2", "HBM", "HBQ1", "HBQ2", "HBZ", "HBZP1", "HBZP2",
      "HBB", "HBD", "HBG1", "HBG2", "HBE1", "HBBS", "HBBP1", "HBBP2"
    )
    markers=markers[!markers$gene%in%globin_genes,]
  }
  
  return(markers)
}

#' Convert FindAllMarkers results to a list organized by cluster
#'
#' This function takes the results from Seurat's FindAllMarkers and organizes them
#' into a list where each element contains markers for a specific cluster.
#'
#' @param markers A data frame from Seurat's FindAllMarkers results
#'
#' @return A list where each element is a data frame containing markers for one cluster.
#'         Returns NULL if the input doesn't have a 'cluster' column.
#'
#' @examples
#' \dontrun{
#' # Convert markers to a list by cluster
#' marker_list <- all_markers_to_list(markers)
#' 
#' # Access markers for a specific cluster
#' cluster_0_markers <- marker_list[["cluster_0"]]
#' }
#' @export
all_markers_to_list=function(markers){
  if(!"cluster"%in%names(markers)){
    print("column <cluster> should be in the dataframe.")
    return(NULL)
  }else{
    marker_list=list()
    for(i in unique(markers$cluster)){
      name=paste0("cluster_",i)
      marker_list[[name]]=markers[markers$cluster==i,]
    }
  }
  return(marker_list)
}

#' Print top marker genes for each cluster
#'
#' This function prints the top N marker genes for each cluster, either from a
#' marker list or a FindAllMarkers data frame.
#'
#' @param markers Either a list of marker data frames (from all_markers_to_list)
#'               or a FindAllMarkers data frame
#' @param n Integer specifying the number of top markers to print per cluster
#' @param cluster_to_print Character vector specifying which clusters to print.
#'                        If NULL, prints all clusters
#'
#' @return NULL (prints results to console)
#'
#' @examples
#' \dontrun{
#' # Print top 50 markers for all clusters
#' marker_print(markers, n = 50)
#' 
#' # Print top 20 markers for specific clusters
#' marker_print(markers, n = 20, cluster_to_print = c("0", "1"))
#' }
#' @export
marker_print=function(markers, n=100, cluster_to_print=NULL){
  number_to_print=n
  if(is.null(cluster_to_print)){}else{
    if(class(markers)=="list"){
      markers=markers[cluster_to_print]
    }else{
      markers=markers[markers$cluster==cluster_to_print,]
    }
  }
  if(class(markers)=="list"){
    for(i in names(markers)){
      print(i)
      print(paste(markers[[i]][markers[[i]]$avg_log2FC>0,][1:number_to_print,]$gene,collapse = ", "))
    }
  }else{
    for(i in unique(markers$cluster)){
      print(i)
      print(paste(markers[markers$cluster==i,][markers$avg_log2FC>0,][1:number_to_print,]$gene,collapse = ", "))
    }
  }
}
