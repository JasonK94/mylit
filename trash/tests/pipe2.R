set.seed(1234)

sample_numbers=c(6,6,6,6,8,8,8,8)
barcode_mapping_list=list()

start_num=1
# SNP demultiplexing & object pipeline preprocessing
for(i in seq_along(list_snp)){
  # demultiplexing
  demux_data=read.csv(list_demux[[i]])
  colnames(demux_data)=c("BARCODE",generate_sample_names(start_num:(start_num+sample_numbers[i]-1)))
  start_num=start_num+sample_numbers[i]
  barcode_mapping_list[[i]]=get_barcode_mapping(demux_data)%>%
    mutate(droplet_demulti=ifelse(is_doublet(Best_Sample),"doublet_demulti","singlet_demulti"), GEM=paste0("GEM",i), Barcode = paste0(Barcode, "_",i))%>%
    {rownames(.)=.$Barcode;
    .} %>%
    mutate(join_key=Best_Sample, day=1)
  
  # object pipeline preprocessing
  ol[[i]]=Read10X(list_snp[[i]])
  sprintf("1. number:%d, raw count:%d",i,ncol(ol[[i]]))
  ol[[i]]=CreateSeuratObject(ol[[i]])
  
  ## demultiplexing - doublet removal
  ol[[i]]=AddMetaData(ol[[i]],barcode_mapping_list[[i]])
  
  ol[[i]]=subset(ol[[i]],droplet_demulti=="singlet_demulti")
  
  ## Pre-processing
  sprintf("2. number:%d, raw count:%d",i,ncol(ol[[i]]))
  # ol[[i]]=subset(ol[[i]],subset=percent.mt<5)
  ol[[i]]=SCTransform(ol[[i]],verbose=F)
  DefaultAssay(ol[[i]])="SCT"
  ol[[i]]=FindVariableFeatures(ol[[i]],verbose=F)
  ol[[i]]=RunPCA(ol[[i]],verbose=F)
  ol[[i]]=FindNeighbors(ol[[i]],verbose=F)
  ol[[i]]=FindClusters(ol[[i]],verbose=F)
  
  ## SoupX
  sprintf("sample: %s, step: soupx",i)
  raw_count=Read10X(list_snp_raw[[i]])
  
  ### SoupChannel(tod, toc, metaData=NULL, calcSoupProfile=TRUE)
  ### tod는 table of droplets으로 cell이 없는 droplet도 포함한 raw count matrix이며, toc는 table of counts로 filetered count matrix
  ### gene list(rownames)를 일치시켜줘야 함
  sc=SoupChannel(raw_count[rownames(ol[[i]]),],ol[[i]]@assays$SCT$counts)
  print("SoupChannel Done")
  
  ### [["seurat_clusters"]]와 $seurat_clusters는 다르다. identical(ol[[i]]$"seurat_clusters",ol[[i]][["seurat_clusters"]])
  sc=setClusters(sc,ol[[i]]$"seurat_clusters")
  print("setCluster Done")
  
  sc=autoEstCont(sc, tfidfMin=0.05, soupQuantile=0.8, forceAccept=TRUE)
  print("autoEstCont Done")
  out=adjustCounts(sc)
  print("adjustCounts Done")
  
  ## re-preprocessing
  sl[[i]]=CreateSeuratObject(out, project="LIT_2024_Stroke_SNPDemultiAndDoubletRemoval_SoupXed_scDblFinderDoubletRemoval", meta.data = ol[[i]]@meta.data)
  sl[[i]]=SCTransform(sl[[i]],verbose=F)
  DefaultAssay(sl[[i]])="SCT"
  sl[[i]]=FindVariableFeatures(sl[[i]],verbose=F)
  sl[[i]]=RunPCA(sl[[i]],verbose=F)
  sl[[i]]=FindNeighbors(sl[[i]],verbose=F)
  sl[[i]]=FindClusters(sl[[i]],verbose=F)
  
  ## scDblFinder doublet removal
  sce <- scDblFinder(GetAssayData(sl[[i]], layer= "counts"), clusters = Idents(sl[[i]]))
  sl[[i]]$"droplet"=sce$scDblFinder.class
  sl[[i]]=subset(sl[[i]],droplet=="singlet")
  sprintf("3. number:%d, raw count:%d",i,ncol(sl[[i]]))
  
  
  # saving
  save_seurat_to_h5ad(
    sl[[i]],
    condaenv="/home/jaecheon/miniconda3/envs/scenvi", #default
    save_path="/data/kjc1/projects/#130.stroke/sobj/h5ad/", #default
    save_name=paste0("stroke_GEM",i,"_" ,format(Sys.time(),"%y-%m-%d-%H-%M"))
  )
}

#HTO의 경우, HTO 갯수가 들어있는 layer도 있어서 명시적으로 layer 불러내어 SCTransform 시행 필요.
for(dir in list_hto){
  sample=basename(dir)
  # Extract first number (prefix)
  prefix <- stringr::str_extract(sample, "^\\d+")
  # Extract last number (time point)
  time_point <- stringr::str_extract(sample, "\\d+$")
  day=case_when(
    time_point=="72" ~3,
    time_point=="48" ~2,
    time_point=="24" ~1,
    time_point=="NA" ~NA
  )
  GEM_hto=basename(list_hto_raw[[list_hto_raw_match[[sample]]]])
  GEM=case_when(
    GEM_hto=="Samples1_4"~9,
    GEM_hto=="Samples5_9"~10,
    GEM_hto=="GEM1_2"~11,
    GEM_hto=="GEM2_2"~12
  )
  
  # HTO based demultiplexing already done
  
  # object pipeline preprocessing
  ol[[sample]]=Read10X(list_hto[[sample]])
  sprintf("1. sample:%s, raw count:%d",sample,ncol(ol[[sample]]))
  ol[[sample]]=CreateSeuratObject(ol[[sample]])
  ol[[sample]]=AddMetaData(ol[[sample]],metadata=sample,col.name = "sample_hto")
  ol[[sample]]=AddMetaData(ol[[sample]],metadata=GEM_hto,col.name = "GEM_hto")
  ol[[sample]]=AddMetaData(ol[[sample]],metadata=GEM,col.name = "GEM")
  ol[[sample]]=AddMetaData(ol[[sample]],metadata=time_point,col.name = "ol_time_point")
  ol[[sample]]=AddMetaData(ol[[sample]],metadata=day,col.name = "day")
  ol[[sample]]=AddMetaData(ol[[sample]],metadata=colnames(ol[[sample]]),col.name = "Barcode")
  ol[[sample]]=AddMetaData(ol[[sample]],metadata=paste0(prefix,"_",day),col.name = "join_key")
  
  sprintf("2. sample:%s, raw count:%d",sample,ncol(ol[[sample]]))
  # ol[[sample]]=subset(ol[[sample]],subset=percent.mt<5)
  ol[[sample]] <- SCTransform(
    ol[[sample]],
    assay           = "RNA",
    layer           = "counts.Gene Expression",
    new.assay.name  = "SCT",
    verbose         = FALSE
  )
  DefaultAssay(ol[[sample]])="SCT"
  ol[[sample]]=FindVariableFeatures(ol[[sample]])
  ol[[sample]]=RunPCA(ol[[sample]])
  ol[[sample]]=FindNeighbors(ol[[sample]])
  ol[[sample]]=FindClusters(ol[[sample]])
  
  ## SoupX
  sprintf("sample: %s, step: soupx",sample)
  raw_count=Read10X(list_hto_raw[[list_hto_raw_match[[sample]]]])
  
  ### SoupChannel(tod, toc, metaData=NULL, calcSoupProfile=TRUE)
  ### tod는 table of droplets으로 cell이 없는 droplet도 포함한 raw count matrix이며, toc는 table of counts로 filetered count matrix
  ### gene list(rownames)를 일치시켜줘야 함
  
  
  #sc=SoupChannel(raw_count$`Gene Expression`[rownames(ol[[sample]]),],ol[[sample]]@assays$SCT$counts,metaData = ol[[i]]@meta.data)
  # Error in SoupChannel(raw_count$`Gene Expression`[rownames(ol[[sample]]),  : 
  #   Rownames of metaData must match column names of table of counts.
  # In addition: Warning message:
  # In sort(colnames(toc)) == sort(rownames(metaData)) :
  #   longer object length is not a multiple of shorter object length
  
  sc=SoupChannel(raw_count$`Gene Expression`[rownames(ol[[sample]]),],ol[[sample]]@assays$SCT$counts)
  print("SoupChannel Done")
  
  ### [["seurat_clusters"]]와 $seurat_clusters는 다르다. identical(ol[[i]]$"seurat_clusters",ol[[i]][["seurat_clusters"]])
  sc=setClusters(sc,ol[[sample]]$"seurat_clusters")
  print("setCluster Done")
  
  sc=autoEstCont(sc, tfidfMin=0.05, soupQuantile=0.8, forceAccept=TRUE)
  print("autoEstCont Done")
  out=adjustCounts(sc)
  print("adjustCounts Done")
  
  sl[[sample]]=CreateSeuratObject(out, project="LIT_2023_Stroke_HTO_SoupXed_scDblFinderDoubletRemoval", meta.data = ol[[sample]]@meta.data)
  sl[[sample]]=SCTransform(sl[[sample]],verbose=F)
  DefaultAssay(sl[[sample]])="SCT"
  sl[[sample]]=FindVariableFeatures(sl[[sample]],verbose=F)
  sl[[sample]]=RunPCA(sl[[sample]],verbose=F)
  sl[[sample]]=FindNeighbors(sl[[sample]],verbose=F)
  sl[[sample]]=FindClusters(sl[[sample]],verbose=F)
  
  sce <- scDblFinder(GetAssayData(sl[[sample]], layer= "counts"), clusters = Idents(sl[[sample]]))
  sl[[sample]]$"droplet"=sce$scDblFinder.class
  sprintf("3. sample:%s, raw count:%f",sample,ncol(sl[[sample]]))
  
  
  # Seurat 객체에서 데이터 추출
  save_seurat_to_h5ad(
    sl[[sample]],
    condaenv="/home/jaecheon/miniconda3/envs/scenvi", #default
    save_path="/data/kjc1/projects/#130.stroke/sobj/h5ad/", #default
    save_name=paste0(sample,"_" ,format(Sys.time(),"%y-%m-%d-%H-%M"))
  )
}

# merging only carries raw/data/scale.data and meta.data slots.
# because there's no safe way to merge reductions, graphs, etc.
# set merge.SCT=TRUE

# features <- SelectIntegrationFeatures(object.list = ol, nfeatures = 3000)
# sct_list <- PrepSCTIntegration(object.list = ol, anchor.features = features) #what's this?

features <- SelectIntegrationFeatures(object.list = sl, nfeatures = 3000)
sct_list <- PrepSCTIntegration(object.list = sl, anchor.features = features) #what's this?
integrated_data <- FindIntegrationAnchors(object.list = sct_list, normalization.method = "SCT",
                                          anchor.features = features,reduction="rpca",
                                          k.anchor=20) #it makes anchorset object. but after IntegrateData it becomes seurat object. before rpca, RunPCA is needed
integrated_data <- IntegrateData(anchorset = integrated_data, normalization.method = "SCT")
# Re-SCTranform to merge multiple SCT layers
DefaultAssay(integrated_data)="SCT"
integrated_data=SCTransform(integrated_data)

# Run PCA on integrated data
DefaultAssay(integrated_data)="integrated"
integrated_data <- RunPCA(integrated_data, verbose = FALSE)
# Run UMAP/t-SNE and clustering
integrated_data <- RunUMAP(integrated_data, dims = 1:30) #The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric. To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
integrated_data <- FindNeighbors(integrated_data, dims = 1:30)
# integrated_data <- FindClusters(integrated_data)
integrated_data <- FindClusters(integrated_data,graph.name = "integrated_snn",resolution = 0.8)

DefaultAssay(integrated_data)="SCT"
integrated_data = PrepSCTFindMarkers(integrated_data)


markers=list()
for(cluster in levels(integrated_data$seurat_clusters)){
  markers[[cluster]]=FindMarkers(integrated_data, ident.1=cluster)
  markers[[cluster]]$gene=rownames(markers[[cluster]])
}

time1=format(Sys.time(),"%y-%m-%d-%H-%M")
saveRDS(integrated_data,file=paste0("/data/kjc1/projects/#130.stroke/sobj/stroke_karo_SCT_integrated_QC_",time1,".rds"))
#saveRDS(integrated_data,file=paste0("/data/kjc1/projects/#130.stroke/sobj/stroke_SCT_integrated_beforeQC_",time1,".rds"))

time2=format(Sys.time(),"%y-%m-%d-%H-%M")
saveRDS(markers,file=paste0("/data/kjc1/projects/#130.stroke/sobj/stroke_karo_markers_",time2,".rds"))
#saveRDS(markers,file=paste0("/data/kjc1/projects/#130.stroke/sobj/stroke_markers_",time2,".rds"))

while (TRUE) {
  # result <- 1 + 1
  # print(result)
  print(format(Sys.time(),"%y-%m-%d-%H-%M"))
  Sys.sleep(100)  # 100초 대기
}