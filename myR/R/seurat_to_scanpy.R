#' @import reticulate Seurat
#' @export
save_seurat_to_h5ad <- function(seurat_obj,
                                condaenv = "/home/jaecheon/miniconda3/envs/scenvi",
                                save_path = "/data/kjc1/projects/#130.stroke/sobj/h5ad/",
                                save_name = NULL,
                                save_SCT = FALSE) {
  # Python과 anndata 설정
  reticulate::use_condaenv(condaenv) # 또는 use_python()으로 Python 지정
  anndata <- import("anndata")
  np <- import("numpy")


  # 기본 카운트 매트릭스 (X로 사용할 데이터)
  # 일반적으로 raw counts를 기본으로 사용
  if (save_SCT) {
    DefaultAssay(seurat_obj) <- "SCT"
  } else {
    DefaultAssay(seurat_obj) <- "RNA"
  }
  if (DefaultAssay(seurat_obj) == "SCT" && dim(seurat_obj[["SCT"]]@counts)[1] > 0) {
    main_matrix <- as.matrix(seurat_obj[["SCT"]]@counts)
  } else {
    main_matrix <- as.matrix(seurat_obj[["RNA"]]@layers$counts)
  }

  # 메타데이터 준비
  metadata <- as.data.frame(seurat_obj@meta.data)

  # 유전자 정보 준비
  var_data <- data.frame(
    gene_symbols = rownames(seurat_obj),
    row.names = rownames(seurat_obj)
  )

  # 모든 assay의 layer를 Python 딕셔너리로 준비
  layers_dict <- dict()

  # RNA assay의 layer 추가
  if ("RNA" %in% names(seurat_obj@assays)) {
    for (layer in Layers(seurat_obj@assays$RNA)) {
      layer_name <- paste0("RNA_", layer)
      layers_dict[layer_name] <- np$array(t(as.matrix(seurat_obj[["RNA"]]@layers[[layer]])))
    }
  }
  if (save_SCT) {
    # SCT assay의 layer 추가; SCT assay가 있다면, 하위 layers를 iteration하여 저장
    if ("SCT" %in% names(seurat_obj@assays)) {
      if (!is.null(Layers(seurat_obj[["SCT"]]))) {
        for (layer_name in Layers(seurat_obj[["SCT"]])) {
          # scale.data는 전체 유전자가 아닌 일부만 포함할 수 있으므로 별도 처리
          # 지금 여기서는 scale.data로는 한정되지 않는다. 실은 data, counts 모두 feature 수가 SCT layer가 RNA에 비해 적다.
          # 따라서 일반화하였으나 변수명은 바꾸지 않았음.
          scale_data <- GetAssayData(seurat_obj[["SCT"]], layer = layer_name)
          scale_genes <- rownames(scale_data)
          scale_cells <- colnames(scale_data)

          # 모든 유전자 x 모든 세포 크기의 매트릭스 생성
          full_scale_data <- matrix(0, nrow = nrow(main_matrix), ncol = ncol(main_matrix))
          rownames(full_scale_data) <- rownames(seurat_obj)
          colnames(full_scale_data) <- colnames(seurat_obj)

          # scale.data 값을 적절한 위치에 복사
          gene_idx <- match(scale_genes, rownames(full_scale_data))
          cell_idx <- match(scale_cells, colnames(full_scale_data))
          gene_idx <- gene_idx[!is.na(gene_idx)]
          cell_idx <- cell_idx[!is.na(cell_idx)]

          layer_key <- paste0("SCT_", layer_name)
          if (length(gene_idx) > 0 && length(cell_idx) > 0) {
            full_scale_data[gene_idx, cell_idx] <- as.matrix(scale_data[
              rownames(scale_data) %in% rownames(full_scale_data),
              colnames(scale_data) %in% colnames(full_scale_data)
            ])
            layers_dict[layer_key] <- np$array(t(full_scale_data))
          }


          # 단순하게 하자면 for 문 아래로 이 것만 있으면 된다. 하지만 그러면 h5ad로 저장이 안 됨.
          # layers_dict[layer_key]= np$array(t(as.matrix(GetAssayData(seurat_obj[["SCT"]],layer=layer_name))))
        }
      }
    }
  }

  # 차원 축소 결과 추가를 위한 obsm 딕셔너리 생성
  obsm_dict <- dict()
  if ("pca" %in% names(seurat_obj@reductions)) {
    obsm_dict["X_pca"] <- np$array(as.matrix(seurat_obj[["pca"]]@cell.embeddings))
  }
  if ("umap" %in% names(seurat_obj@reductions)) {
    obsm_dict["X_umap"] <- np$array(as.matrix(seurat_obj[["umap"]]@cell.embeddings))
  }
  if ("tsne" %in% names(seurat_obj@reductions)) {
    obsm_dict["X_tsne"] <- np$array(as.matrix(seurat_obj[["tsne"]]@cell.embeddings))
  }

  # 유전자 차원 축소 결과(varm)가 있다면 추가 -> HVG에 대한 PC만 계산하여, adata 형식과 incompatible함.
  # varm_dict <- dict()
  # if ("pca" %in% names(seurat_obj@reductions) && !is.null(seurat_obj[["pca"]]@feature.loadings)) {
  #   varm_dict["PCs"]= np$array(as.matrix(seurat_obj[["pca"]]@feature.loadings))
  # }

  # uns에 추가 정보 저장 (선택 사항)
  uns_dict <- dict()
  uns_dict["seurat_default_assay"] <- DefaultAssay(seurat_obj)
  uns_dict["seurat_version"] <- packageVersion("Seurat")
  if (!is.null(seurat_obj@project.name)) {
    uns_dict["project_name"] <- seurat_obj@project.name
  }


  # AnnData 객체 생성 (X는 transpose해야 함 - AnnData는 세포 x 유전자 형식)
  adata <- anndata$AnnData(
    X = np$array(t(main_matrix)),
    obs = metadata,
    var = var_data,
    layers = layers_dict, # SCT의 경우 feature 갯수가 줄어드는 문제가 있음.
    obsm = obsm_dict,
    uns = uns_dict,
    # varm = varm_dict
  )

  # h5ad 파일로 저장
  if (is.null(save_name)) {
    save_name <- format(Sys.time(), "%y-%m-%d-%H-%M")
  }
  file_name <- paste0(save_path, save_name, ".h5ad")
  adata$write(file_name)

  print(paste0("Successfully saved ", save_name, " to ", file_name))
}
