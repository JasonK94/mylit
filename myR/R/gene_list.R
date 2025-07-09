#' Azimuth (AZI) reference gene lists
#'
#' A named list containing gene sets used in Azimuth analysis.
#' https://azimuth.hubmapconsortium.org/references/human_pbmc/
#'
#' @docType data
#' @format A named list of character vectors where each element is a vector of gene symbols.
#' @source <https://azimuth.hubmapconsortium.org/references/human_pbmc/>
#' @keywords datasets
#' @rdname gene_lists
#' @export
genes_azi=list(
  azi_bi=c('MS4A1','TNFRSF13B','IGHM','IGHD','AIM2','CD79A','LINC01857','RALGPS2','BANK1','CD79B'),
  azi_bm=c('MS4A1','COCH','AIM2','BANK1','SSPN','CD79A','TEX9','RALGPS2','TNFRSF13C','LINC01781'),
  azi_bn=c('IGHM','IGHD','CD79A','IL4R','MS4A1','CXCR4','BTG1','TCL1A','CD79B','YBX3'),
  azi_pb=c('IGHA2','MZB1','TNFRSF17','DERL3','TXNDC5','TNFRSF13B','POU2AF1','CPNE5','HRASLS2','NT5DC2'),
  azi_4ctl=c('GZMH','CD4','FGFBP2','ITGB1','GZMA','CST7','GNLY','B2M','IL32','NKG7'),
  azi_4n=c('TCF7','CD4','CCR7','IL7R','FHIT','LEF1','MAL','NOSIP','LDHB','PIK3IP1'),
  azi_4pro=c('MKI67','TOP2A','PCLAF','CENPF','TYMS','NUSAP1','ASPM','PTTG1','TPX2','RRM2'),
  azi_4tcm=c('IL7R','TMSB10','CD4','ITGB1','LTB','TRAC','AQP3','LDHB','IL32','MAL'),
  azi_4tem=c('IL7R','CCL5','FYB1','GZMK','IL32','GZMA','KLRB1','TRAC','LTB','AQP3'),
  azi_treg=c('RTKN2','FOXP3','AC133644.2','CD4','IL2RA','TIGIT','CTLA4','FCRL3','LAIR2','IKZF2'),
  azi_8n=c('CD8B','S100B','CCR7','RGS10','NOSIP','LINC02446','LEF1','CRTAM','CD8A','OXNAD1'),
  azi_8pro=c('MKI67','CD8B','TYMS','TRAC','PCLAF','CD3D','CLSPN','CD3G','TK1','RRM2'),
  azi_8tcm=c('CD8B','ANXA1','CD8A','KRT1','LINC02446','YBX3','IL7R','TRAC','NELL2','LDHB'),
  azi_8tem=c('CCL5','GZMH','CD8A','TRAC','KLRD1','NKG7','GZMK','CST7','CD8B','TRGC2'),
  azi_asdc=c('PPP1R14A','LILRA4','AXL','IL3RA','SCT','SCN9A','LGMN','DNASE1L3','CLEC4C','GAS6'),
  azi_cdc1=c('CLEC9A','DNASE1L3','C1orf54','IDO1','CLNK','CADM1','FLT3','ENPP1','XCR1','NDRG2'),
  azi_cdc2=c('FCER1A','HLA-DQA1','CLEC10A','CD1C','ENHO','PLD4','GSN','SLC38A1','NDRG2','AFF3'),
  azi_pdc=c('ITM2C','PLD4','SERPINF1','LILRA4','IL3RA','TPM2','MZB1','SPIB','IRF4','SMPD3'),
  azi_14mo=c('S100A9','CTSS','S100A8','LYZ','VCAN','S100A12','IL1B','CD14','G0S2','FCN1'),
  azi_16mo=c('CDKN1C','FCGR3A','PTPRC','LST1','IER5','MS4A7','RHOC','IFITM3','AIF1','HES4'),
  azi_nk=c('GNLY','TYROBP','NKG7','FCER1G','GZMB','TRDC','PRF1','FGFBP2','SPON2','KLRF1'),
  azi_nkpro=c('MKI67','KLRF1','TYMS','TRDC','TOP2A','FCER1G','PCLAF','CD247','CLSPN','ASPM'),
  azi_nk56=c('XCL2','FCER1G','SPINK2','TRDC','KLRC1','XCL1','SPTSSB','PPP1R9A','NCAM1','TNFRSF11A'),
  azi_eryth=c('HBD','HBM','AHSP','ALAS2','CA1','SLC4A1','IFIT1B','TRIM58','SELENBP1','TMCC2'),
  azi_hspc=c('SPINK2','PRSS57','CYTL1','EGFL7','GATA2','CD34','SMIM24','AVP','MYB','LAPTM4B'),
  azi_ilc=c('KIT','TRDC','TTLL10','LINC01229','SOX4','KLRB1','TNFRSF18','TNFRSF4','IL1R1','HPGDS'),
  azi_plt=c('PPBP','PF4','NRGN','GNG11','CAVIN2','TUBB1','CLU','HIST1H2AC','RGS18','GP9'),
  azi_dnt=c('PTPN3','MIR4422HG','NUCB2','CAV1','DTHD1','GZMA','MYB','FXYD2','GZMK','AC004585.1'),
  azi_gdt=c('TRDC','TRGC1','TRGC2','KLRC1','NKG7','TRDV2','CD7','TRGV9','KLRD1','KLRG1'),
  azi_mait=c('KLRB1','NKG7','GZMK','IL7R','SLC4A10','GZMA','CXCR6','PRSS35','RBM24','NCR3'),
  mast = c('SIGLEC6','KIT', 'TPSAB1', 'CPA3', 'HDC', 'FCER1A', 'MS4A2', 'HPGD', 'CAVIN2', 'IL1RL1', 'CKLF'),
  th1 = c('TBX21', 'IFNG', 'STAT1', 'STAT4', 'IL12RB2', 'CXCR3', 'CCR5', 'IL2', 'TNF', 'CD40LG'),
  th2 = c('GATA3', 'IL4', 'IL5', 'IL13', 'STAT6', 'IL4R', 'CCR4', 'IL10', 'IRF4', 'PTGDR2'),
  th17 = c('RORC', 'IL17A', 'IL17F', 'IL21', 'IL22', 'STAT3', 'CCR6', 'IL23R', 'CTLA4', 'CSF2'),
  mono = c('CD14', 'LYZ', 'FCGR3A', 'MS4A7', 'CD68', 'ITGAM', 'CCR2', 'CSF1R', 'S100A8', 'S100A9'),
  Treg_markers = c(
    "FOXP3",  # Master transcription factor for Tregs
    "IL2RA",  # CD25, high expression in Tregs
    "CTLA4",  # Immune checkpoint molecule
    "TNFRSF18",  # GITR, involved in Treg activation
    "IKZF2",  # Helios, a marker for thymic-derived Tregs
    "IKZF4",  # Eos, involved in Treg function
    "LAG3",  # Inhibitory receptor
    "CD274",  # PD-L1, involved in immune suppression
    "ENTPD1",  # CD39, contributes to adenosine-mediated suppression
    "NT5E",  # CD73, works with CD39 for adenosine production
    "TGFB1",  # Produces TGF-β, a key immunosuppressive cytokine
    "IL10",  # Produces IL-10, another immunosuppressive cytokine
    "CCR4",  # Chemokine receptor guiding Treg migration
    "CCR8",  # Enriched in suppressive Tregs
    "STAT5B",  # Required for Treg stability
    "CD127"  # IL7R (LOW expression in Tregs)
  ),
  mDC_markers=c(
    "CD1C",  # Classical cDC1 marker
    "CLEC10A",  # Classical cDC2 marker
    "FLT3",  # Key growth factor receptor for DCs
    "IRF8",  # Transcription factor for cDC1 differentiation
    "IRF4",  # Transcription factor for cDC2 differentiation
    "BATF3",  # Involved in cross-presentation function
    "HLA-DRA",  # MHC class II molecule for antigen presentation
    "HLA-DRB1",  # Another MHC-II molecule
    "CD80",  # Co-stimulatory molecule
    "CD86",  # Co-stimulatory molecule
    "CCR7",  # Guides migration to lymph nodes
    "XCR1",  # Marker for cross-presenting cDC1
    "LAMP3",  # Migratory mature DC marker
    "CCL17",  # Chemokine involved in T cell attraction
    "CSF1R",  # Myeloid growth factor receptor
    "SIRPA",  # Differentiates cDC2 from cDC1
    "CD209",  # DC-SIGN, involved in pathogen recognition
    "FCER1A",  # Fc epsilon receptor, found in some DC subsets
    "CCR6"  # Found in migratory DCs
  ),
  cDC1_markers=c("XCR1", "CLEC9A", "IRF8", "BATF3", "FLT3", "HLA-DRA", "CD80", "CD86"),
  cDC2_markers=c("CD1C", "CLEC10A", "IRF4", "SIRPA", "CSF1R", "HLA-DRA", "CD80", "CD86")
  
)

# New pbmc markers ----
#' -------------------- 1.2  fresh 2022‑25 marker additions ---------------
#' @rdname gene_lists
#' @export
genes_pbmc_new <- list(
  ## memory / atypical B cells (FCRL5, TBX21) – Su2023; Hao2024
  memB_FCRL5 = c("FCRL5","T‑BET","ITGAX","ZEB2","CXCR3","CD11C"),                 # :contentReference[oaicite:0]{index=0}
  ## TPEX / progenitor‑exhausted CD8 (TCF7, SLAMF6, CXCR5) – anti‑PD‑1 atlas
  CD8_TPEX   = c("TCF7","SLAMF6","PDCD1","CXCR5","LAG3","CR2","IL7R"),             # :contentReference[oaicite:1]{index=1}
  ## GZMKhigh CD8 TEM (pro‑inflammatory) – CRC & airway atlases
  CD8_GZMK   = c("GZMK","RUNX3","ZBTB38","CXCR3","KLRC1","ITGA1","HLA‑DPA1"),      # :contentReference[oaicite:2]{index=2}
  ## STAB1+ foetal‑like M2 TAM – NSCLC onco‑foetal macrophages
  TAM_STAB1  = c("STAB1","KYNU","NAMPT","MERTK","MARCO","TGFB1","C1QC"),           # :contentReference[oaicite:3]{index=3}
  ## pDC maturation (LAMP5, TCL1A) – updated DC primer 2024
  pDC_mature = c("LILRA4","IL3RA","LAMP5","TCL1A","IRF7","TCF4"),                  # :contentReference[oaicite:4]{index=4}
  ## cDC1 refined
  cDC1_2024  = c("XCR1","CLEC9A","BATF3","IRF8","WDFY4","ID2","NEC3"),             # :contentReference[oaicite:5]{index=5}
  ## cDC2 refined
  cDC2_2024  = c("CD1C","CLEC10A","IRF4","SIRPA","FCER1A","CD301A","LAMP3"),       # :contentReference[oaicite:6]{index=6}
  ## STMN1+ cycling DCs – anti‑PD‑1 atlas
  DC_cycle   = c("MKI67","TOP2A","STMN1","NUSAP1","UBE2C","PTTG1"),                # :contentReference[oaicite:7]{index=7}
  ## IFN‑stimulated monocyte state (ISG high)
  ISG_mono   = c("IFI6","ISG15","MX1","OAS1","IFIT1","IFIT3","CXCL10"),            # :contentReference[oaicite:8]{index=8}
  ## Treg stability panel (updated 2024)
  Treg_core  = c("FOXP3","IL2RA","CTLA4","TNFRSF18","IKZF2","ENTPD1","CCR8")       # :contentReference[oaicite:9]{index=9}
)

# Skin genes ----
#' @rdname gene_lists
#' @export
genes_skin = list(
  # Keratinocytes
  Basal_Keratinocytes = c("KRT5", "KRT14", "TP63"),
  Spinous_Keratinocytes = c("KRT1", "KRT10", "IVL"),
  Granular_Keratinocytes = c("FLG", "LOR", "SPINK5"),
  Corneocytes_Keratinocytes = c("FLG", "LOR", "IVL", "HRNR"),
  Keratinocytes = c("KRT1", "KRT5", "KRT10", "KRT14", "IVL", "FLG", "LOR"),
  # Melanocytes
  Melanocytes = c("MITF", "TYR", "TYRP1", "DCT", "PMEL"),
  # Fibroblasts
  Fibroblasts = c("COL1A1", "DCN", "POSTN", "COL3A1", "FAP", "IGFBP5"),
  # T Cells
  T_Cells = c("CD3D", "CD3E", "CD3G", "TRAC"),
  CD4_T_Helpers = c("CD4", "IL7R", "FOXP3"),
  CD8_T_Cytotoxic = c("CD8A", "CD8B", "GZMB"),
  # Dendritic Cells
  Dendritic_Cells = c("CD1C", "CLEC9A", "ITGAX"),
  Langerhans_Cells = c("CD207", "CD1A", "Langerin"),
  # Macrophages
  Macrophages = c("CD68", "CD163", "CD14"),
  # Mast Cells
  Mast_Cells = c("TPSB2", "KIT", "FCER1A"),
  # Endothelial Cells
  Endothelial_Cells = c("PECAM1", "VWF", "CDH5"),
  # Adipocytes
  Adipocytes = c("ADIPOQ", "FABP4", "LEP"),
  # Merkel Cells
  Merkel_Cells = c("KRT20", "KRT8", "PLEKHA1"),
  
  # Hair Follicle Cells
  Hair_Follicle_Stem = c("SOX9", "LEF1", "CD34"),
  Dermal_Papilla = c("AXIN2", "WIF1", "PDGFRB"),
  
  # Sebaceous Gland Cells
  Sebaceous_Gland = c("SOX9", "MUC1", "KRT7", "KRT19"),
  
  # Eccrine Gland Cells
  Eccrine_Gland = c("KRT7", "KRT19", "AQP5", "SCGB2A2"),
  
  # Lymphatic Endothelial Cells
  Lymphatic_Endothelial = c("PROX1", "LYVE1", "PDPN", "VEGFR3"),
  
  # Stem Cells
  Stem_Cells = c("SOX9", "LGR5", "KRT15", "CD34"),
  
  # Monocytes
  Monocytes = c("CD14", "FCGR3A", "LYZ", "S100A8"),
  schwann = c(  "S100B", "S100A4", "S100A6", "PMP22", "MPZ", "PLP1", "NGFR", "SOX10", "GFAP"),
  neuronal_markers = c(  # Synaptic and Adhesion
    "CNTNAP2", "NRXN1", "NRXN3", "CNTN5", "DSCAM", "DLG2", "DLGAP1", 
    "PTPRD", "ROBO1", "ROBO2", "CDH12", "CDH18", "SHANK2"),
  g_skin_HFSm=c("KRT15", "KRT19", "CD34", "ITGA6", "ITGB1", "LGR5", "LGR6", "SOX9", "SOX10", "MSI2", "KRT73", "KRT83", "FGF22"),
  g_skin_DPCm=c("ALPL", "SOX2", "LEF1", "BMP4", "BMP6", "NOG", "WNT5A", "WNT10B", "PDGFRA", "PDGFRB", "FOXC1", "VCAN", "COL1A1", "COL1A2", "DCN", "ALDH1A1"),
  g_skin_SMCm <- c("ACTA2", "CNN1", "MYH11", "TAGLN", "DES", "SMTN", "CALD1", "VIM", "MYLK", "HSPB1", "LMOD1", "PDGFRB", "CAV1")
)

# RCC genes ----
#' @rdname gene_lists
#' @export
geomx_rcc <- list(
  ccRCC_VHL_related_pathways = c("VEGFA", "CA9", "EGLN3", "HIF1A", "EPAS1"),
  ccRCC_lipid_metabolism = c("PLIN2", "FABP7", "HMGCS1"),
  ccRCC_hypoxia_markers = c("GLUT1", "ENO1", "HK2", "LDHA"),
  ccRCC_markers = c("CD10", "CA9"),
  ccRCC_proliferation_markers = c("MKI67", "PCNA"),
  CAF_activation = c("ACTA2", "FAP", "POSTN", "PDGFRB"),
  CAF_matrix_remodeling = c("COL1A1", "COL3A1", "MMP2", "MMP9", "SPARC"),
  CAF_growth_factors = c("TGFB1", "CTGF", "PDGFA", "PDGFB"),
  immune_T_cells = c("CD3D", "CD3E", "CD3G", "CD8A", "GZMB", "PRF1"),
  immune_Tregs = c("FOXP3", "IL2RA"),
  immune_exhausted_T_cells = c("PDCD1", "CTLA4", "LAG3", "HAVCR2"),
  immune_B_cells = c("CD19", "MS4A1", "CD79A", "CD79B"),
  immune_NK_cells = c("NCAM1", "NKG7", "PRF1", "GZMA"),
  immune_M1_macrophages = c("CD68", "IL6", "TNF", "CXCL9", "CXCL10"),
  immune_M2_macrophages = c("CD163", "MRC1", "TGFB1", "IL10", "CCL22"),
  immune_dendritic_cells = c("CD1C", "CLEC9A", "BATF3", "XCR1"),
  immune_MDSCs = c("ARG1", "NOS2", "S100A8", "S100A9"),
  endothelial_cells = c("PECAM1", "VWF", "FLT1", "KDR"),
  stromal_pericytes = c("PDGFRB", "RGS5", "CSPG4"),
  stromal_adipocytes = c("ADIPOQ", "FABP4", "LEP"),
  stromal_MSCs = c("ENG", "THY1", "NT5E")
)

# RCC 2 ----
#' @rdname gene_lists
#' @export
geomx_rcc2 <- list(
  # Clear Cell Renal Cell Carcinoma (ccRCC) Tumor Cells
  ccRCC_markers = c("CA9", "CD10", "PAX8", "VHL", "NDUFA4L2"),
  general_carcinoma_markers = c(
    "CD274",  # PD-L1
    "MUC1",   # Mucin 1
    "CEACAM5", # CEA (Carcinoembryonic Antigen)
    "VIM",    # Vimentin
    "MKI67",  # Ki-67
    "CD44"    # CD44
  )
  ,
  
  # Endothelial Cells
  endothelial_cells = c("PECAM1", "CD34", "VWF", "KDR"),
  
  # Lymphatic Endothelial Cells
  lymphatic_endothelial_cells = c("FLT4", "PROX1", "PDPN", "LYVE1"),
  
  # Cancer-Associated Fibroblasts (CAFs)
  CAF_activation = c("ACTA2", "FAP", "PDGFRB", "POSTN"),
  CAF_matrix_remodeling = c("COL1A1", "COL3A1", "MMP2", "MMP9"),
  
  # Adipocytes
  stromal_adipocytes = c("ADIPOQ", "FABP4", "PLIN1"),
  
  # Mesenchymal Stem Cells (MSCs)
  stromal_MSCs = c("ENG", "THY1", "NT5E", "PDGFRB"),
  
  # Pericytes
  stromal_pericytes = c("PDGFRB", "RGS5", "CSPG4", "ACTA2"),
  
  # Smooth Muscle Cells
  smooth_muscle_cells = c("ACTA2", "MYH11", "CNN1", "TAGLN"),
  
  # Neuronal Cells
  neuronal_cells = c("RBFOX3", "MAP2", "TUBB3", "NEFL"),
  
  # Glomerular Cells
  glomerular_podocytes = c("NPHS1", "PODXL", "WT1", "SYNPO"),
  
  # Tubular Cells
  proximal_tubule_cells = c("SLC34A1", "LRP2", "AQP1", "CUBN"),
  distal_tubule_cells = c("SLC12A3", "CALB1", "TRPM6", "FXYD2"),
  
  # Immune Cells
  immune_T_cells = c("CD3D", "CD3E", "CD3G", "CD4", "CD8A"),
  immune_Tregs = c("FOXP3", "IL2RA", "CTLA4", "TNFRSF4"),
  immune_exhausted_T_cells = c("PDCD1", "LAG3", "HAVCR2", "TIGIT"),
  immune_B_cells = c("CD19", "MS4A1", "CD79A", "CD79B"),
  immune_NK_cells = c("NCAM1", "NKG7", "KLRD1", "GZMA"),
  immune_M1_macrophages = c("CD68", "NOS2", "IL6", "TNF"),
  immune_M2_macrophages = c("CD163", "MRC1", "TGFB1", "IL10"),
  immune_dendritic_cells = c("CD1C", "CLEC9A", "BATF3", "XCR1"),
  immune_MDSCs = c("ARG1", "S100A8", "S100A9", "IL1B")
)

# gut tissue markers ----
#' @rdname gene_lists
#' @export
genes_gut = list(
  # Epithelial Cells
  Epithelial_General = c("EPCAM", "CDH1", "KRT8", "KRT18"),
  Intestinal_Stem_Cells = c("LGR5", "SOX9", "OLFM4", "ASCL2", "SMOC2", "BMI1", "HOPX"),
  Enterocytes_Absorptive = c("VIL1", "FABP1", "APOA1", "APOA4", "RBP2", "SLC2A2"),
  Goblet_Cells = c("MUC2", "TFF3", "SPDEF", "FCGBP", "KLF4"),
  Paneth_Cells = c("LYZ", "DEFA5", "DEFA6", "MMP7", "SOX9", "DLL4"),
  Enteroendocrine_Cells = c("CHGA", "CHGB", "GCG", "PYY", "SST", "NEUROG3", "NEUROD1"),
  Tuft_Cells = c("DCLK1", "TRPM5", "POU2F3", "IL25"),
  
  # Immune Cells - Lamina Propria
  T_Cells_General = c("CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "TRBC2"),
  CD4_T_Helpers = c("CD4", "IL7R", "CD40LG"),
  Th1 = c("IFNG", "TBX21", "STAT4"),
  Th2 = c("IL4", "IL13", "GATA3", "STAT6"),
  Th17 = c("IL17A", "IL17F", "RORC", "STAT3", "IL23R"),
  T_Regulatory_Cells = c("FOXP3", "IL2RA", "CTLA4", "IKZF2"),
  CD8_T_Cytotoxic = c("CD8A", "CD8B", "GZMA", "GZMB", "PRF1"),
  B_Cells = c("CD19", "MS4A1", "CD79A", "CD79B"),
  Plasma_Cells = c("SDC1", "IGHG1", "IGHA1", "JCHAIN", "XBP1", "PRDM1"),
  Macrophages = c("CD68", "CD163", "CSF1R", "CD14", "ADGRE1"), # ADGRE1 is F4/80
  Monocytes = c("CD14", "FCGR3A", "LYZ", "S100A8", "S100A9"),
  Dendritic_Cells_General = c("ITGAX", "CD83", "CD86", "HLA-DRA"), # ITGAX is CD11c
  cDC1 = c("CLEC9A", "XCR1", "BATF3", "IRF8"),
  cDC2 = c("CD1C", "SIRPA", "CLEC10A", "IRF4"),
  Plasmacytoid_DC = c("CLEC4C", "NRP1", "TCF4", "IRF7"),
  Mast_Cells = c("KIT", "TPSAB1", "CPA3", "FCER1A"),
  Innate_Lymphoid_Cells_General = c("ID2", "IL7R"),
  ILC1 = c("IFNG", "TBX21", "NCR1"), # NCR1 is NKp46
  ILC2 = c("GATA3", "RORA", "IL13", "IL5"),
  ILC3 = c("RORC", "AHR", "IL22", "IL17A", "KIT"),
  
  # Immune Cells - Intraepithelial
  Intraepithelial_Lymphocytes = c("CD3D", "CD8A", "ITGAE", "ENTPD1"), # ITGAE is CD103
  IEL_T_CD8aa = c("CD8A", "ITGAE", "CD244", "ZNF683"), # ZNF683 is HOBIT
  IEL_T_gd = c("TRDC", "TRGC1", "TRGC2", "RORA"),
  
  # Stromal and Other Cells
  Fibroblasts = c("COL1A1", "COL1A2", "COL3A1", "DCN", "PDGFRA", "FAP"),
  Villus_Fibroblasts = c("BMP4", "BMP5", "WNT5A"),
  Crypt_Fibroblasts = c("WNT2B", "DKK3", "RSPO3"),
  Endothelial_Cells_Blood = c("PECAM1", "CDH5", "VWF", "CLDN5", "ENG"),
  Endothelial_Cells_Lymphatic = c("LYVE1", "PROX1", "PDPN"),
  Smooth_Muscle_Cells = c("ACTA2", "MYH11", "TAGLN", "DES", "CNN1"),
  Enteric_Neurons = c("ELAVL4", "TUBB3", "CHAT", "NOS1", "VIP"), # ELAVL4 is HuD
  Enteric_Glia = c("S100B", "GFAP", "SOX10", "PLP1")
)