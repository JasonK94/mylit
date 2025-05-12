#' This function make the findmarker object look good
#' 
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

#' This function make the findmarker object look good
#' 
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
  if("cluster"%in%names(markers)){}else{ #FindAllMarkers object has "cluster" column, FindMarkers not.
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

#' Change FindAllMarkers object to list per each cluster
#' 
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


#' This function make the findmarker object look good
#' 
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
      print(paste(marker_list[[i]][marker_list[[i]]$avg_log2FC>0,][1:number_to_print,]$gene,collapse = ", "))
    }
  }else{
    for(i in unique(markers$cluster)){
      print(i)
      print(paste(marker_list[[i]][marker_list[[i]]$avg_log2FC>0,][1:number_to_print,]$gene,collapse = ", "))
    }
  }
}
