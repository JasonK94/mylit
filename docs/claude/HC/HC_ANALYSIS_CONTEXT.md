# HC-Only Analysis Context Document

> Created: 2026-02-18 | Dataset: hc_only_v1 | Purpose: Age/Sex effects in healthy PBMC

## 1. Dataset Overview

| Item | Value |
|------|-------|
| **Object** | `/data/user3/sobj/hc_only_v1/2_hc_annotated.qs` |
| **Cells** | 549,395 |
| **Genes** | 33,956 |
| **Patients** | 96 (from 5 HC GEMs: HC_GEM18–HC_GEM24) |
| **Sex** | F=57, M=39 |
| **Age** | 19–66, median 38 |
| **Age groups** | Young (<35): 45, Middle (35–50): 28, Old (>50): 23 |
| **Integration** | scVI (n_latent=30, batch=GEM) |
| **Reductions** | `integrated.scvi` (30D), `umap.scvi` (2D) |
| **Annotation** | anno1 (32 types), anno2 (10 compartments) |

### anno2 Compartments

| Compartment | Cells | % |
|---|---|---|
| T cell | 280,903 | 51.1% |
| Myeloid | 115,608 | 21.0% |
| NK | 96,559 | 17.6% |
| B cell | 23,780 | 4.3% |
| DC | 16,741 | 3.0% |
| Platelet | 7,557 | 1.4% |
| ISG-hi | 4,337 | 0.8% |
| Megakaryocyte | 2,119 | 0.4% |
| RBC | 938 | 0.2% |
| HSPC | 853 | 0.2% |

### anno1 (32 cell types)

T cell: CD4 Naive, CD4 Naive(MAL+), CD4 Naive(Quiescent), CD4 TCM, CD4 Th17-Th22, CD4 Treg, CD8 Effector, CD8 Naive, CD8 Naive(minor), CD8 TEM(GZMK+), MAIT, gdT(Vd2)
Myeloid: CD14 Mono Classical, CD14 Mono DOCK4+, CD14 Mono Inflammatory, CD14 Mono Intermediate, CD16 Mono Non-classical
NK: NK Adaptive, NK CD16hi, NK CD56bright, NK CD56dim, NK Transitional
B cell: B Naive, Plasma Cell
DC: cDC1, cDC2, pDC
Other: HSPC, ISG-hi, Megakaryocyte, Platelet, RBC

### Metadata Columns (31)

`orig.ident`, `nCount_RNA`, `nFeature_RNA`, `sample_id`, `name`, `sample_manifest`, `sample_metadata`, `patient_name`, `gem_name`, `batch_name`, `set_name`, `project_name`, `multiplex_method`, `GEM`, `batch_concat`, `n_genes_by_counts`, `total_counts`, `total_counts_mt`, `pct_counts_mt`, `leiden_res0.5`, `leiden_res1.0`, `leiden_res1.5`, `leiden_res2.0`, `solo_prediction`, `solo_doublet_prob`, `hc_id`, `sex`, `age`, `age_group`, `SET`, `anno1`, `anno2`

---

## 2. Analysis Goals

### Primary Questions

1. **Sex effects**: M vs F에서 cell type frequency, gene expression, cell-cell interaction 차이
2. **Age effects**: Continuous age 및 age group (Young/Middle/Old)에 따른 면역 세포 변화
3. **Treg aging**: Treg에서 특히 age에 따른 functional change (exhaustion, suppressive capacity)
4. **Compartment-specific aging**: 어떤 immune compartment가 aging에 가장 민감한가?

### Relation to stroke_hc_v8_2 Project

이 HC-only 분석은 stroke_hc_v8_2 논문의 **baseline characterization**에 해당:
- HC에서의 age/sex normal variation을 이해해야 stroke-specific change를 정확히 분리 가능
- HC FGS top genes가 stroke DEG와 overlap하는지 확인 (confounding vs true signal)
- HC에서의 cell type frequency age-dependency가 MASC 결과 해석에 영향

---

## 3. Completed Analyses

### 3.1 MrVI (Multi-resolution Variational Inference)

| Item | Value |
|------|-------|
| **Output** | `/data/user3/sobj/hc_only_v1/mrvi/` |
| **Model** | `mrvi_model/model.pt` (17 MB) |
| **Latent** | `mrvi_latent.h5ad` (2.0 GB) |

Sample-level variational inference로 환자 간 transcriptomic heterogeneity 분석. Sample PCA에서 age/sex 축 분리 확인.

### 3.2 DEG Consensus — Sex (M vs F)

| Item | Value |
|------|-------|
| **Output** | `/data/user3/sobj/hc_only_v1/deg_consensus_1/` |
| **Methods** | edgeR-LRT, DESeq2-Wald, limma-voom |
| **Cell types** | 32 (all anno1) |
| **Design** | `~ sex + age` (age as covariate) |
| **Results** | `results_consensus.qs` (15.2 GB) |

#### Key Findings

| Cell Type | Autosomal DEGs (consensus FDR<0.05) | Notable Genes |
|---|---|---|
| CD14 Mono Classical | 435 | Most DEGs of any cell type |
| CD4 Treg | 34 | |
| NK CD56dim | ~200 | |

**Pan-immune sex genes** (significant in many cell types):
- RPS4X (30 types), EIF2S3 (25), JPX (22), CD99 (17)
- Sex-linked genes excluded: XIST, RPS4Y1, DDX3Y, EIF1AY, KDM5D, UTY, etc.

### 3.3 MILO — Sex (M vs F)

| Item | Value |
|------|-------|
| **Output** | `/data/user3/sobj/hc_only_v1/milo/` |
| **Cells** | 100K subsampled |
| **Neighbourhoods** | 8,605 |
| **Significant (SpatialFDR<0.1)** | 117 (85 M-enriched, 32 F-enriched) |
| **Bug fixed** | `testNhoods` needed `reduced.dim = "GRAPH"` (no PCA in object) |

#### Results

- Moderate sex-based DA: 117/8605 nhoods (1.4%) significant
- Male-enriched nhoods slightly outnumber female-enriched (85 vs 32)
- Plots: volcano, boxplot by cell type, significance barplot

### 3.4 Age-Proportion Correlation

| Item | Value |
|------|-------|
| **Output** | `/data/user3/sobj/hc_only_v1/age_continuous/` |
| **Methods** | Spearman correlation + GAM (per cell type vs continuous age) |

Cell type frequency vs age analysis at anno1 and anno2 levels, plus sex-stratified and interaction models.

### 3.5 CellChat — Sex & Age

| Item | Value |
|------|-------|
| **Output (sex)** | `/data/user3/sobj/hc_only_v1/cellchat/sex_v1/` |
| **Output (age)** | `/data/user3/sobj/hc_only_v1/cellchat/age_v1/` |
| **Split** | Per-patient, merged by condition |
| **Group** | anno1 cell types |

Sex: M vs F interaction comparison (96 patients split by sex).
Age: Young vs Old interaction comparison (patients pre-grouped).

### 3.6 MultiNicheNet — Sex & Age

| Item | Value |
|------|-------|
| **Output (sex)** | `/data/user3/sobj/hc_only_v1/mnn/sex_v1/` |
| **Output (age)** | `/data/user3/sobj/hc_only_v1/mnn/age_v1/` |
| **Contrasts** | F-M (sex), Old-Young (age) |
| **Results** | ~690K ligand activities (sex), ~890K (age) |

### 3.7 FGS Continuous — Treg × Age (v1: 4 methods, v2: 6 methods)

| Item | v1 | v2 |
|------|----|----|
| **Output** | `fgs_continuous_treg/` | `fgs_continuous_treg_v2/` |
| **Methods** | Spearman, LASSO, RF, limma-voom | + GAM (mgcv), Elastic Net |
| **Top-N** | 50 per method | 200 per method |
| **Consensus (2+)** | 45 genes | 362 genes |
| **Consensus (3+)** | 16 | 84 |
| **Consensus (4+)** | 10 | 32 |
| **All methods** | 10 (4/4) | 4 (6/6) |

#### Enrichment Significance (v1: top-50)

| Tier | Expected | Observed | Fold | p-value |
|---|---|---|---|---|
| 2+ methods | 0.82 | **45** | 54.7x | 5.5e-61 |
| 3+ methods | 0.0015 | **16** | ~10,000x | 3.5e-59 |
| 4 methods | 1e-6 | **10** | ~10⁷x | 4.1e-67 |

#### v2 Top Consensus Genes (all 6 methods)

| Gene | Mean Rank | Direction | Biological Role |
|---|---|---|---|
| **PLAC8** | 1 | Down with age | Immune regulation, IFN-stimulated |
| **EPSTI1** | 2 | Down with age | Interferon-stimulated |
| **ENSG00000274460** | 3 | Up with age | Unannotated |
| **JAML** | 4 | Down with age | Junctional adhesion, T cell migration |

#### Key Observations

- LASSO R²=0.54 (overfit risk; only 10 non-zero coefs at lambda.1se)
- RF OOB R²=0.033 (realistic — age signal is subtle)
- GAM found 2 genes with padj<0.05, confirming weak but real signal
- **ISG genes** (PLAC8, EPSTI1, ISG15) decrease with age → innate immune signaling decline

### 3.8 FGS by Compartment (anno2 × Age)

| Item | Value |
|------|-------|
| **Output** | `/data/user3/sobj/hc_only_v1/fgs_compartment/` |
| **Compartments** | 10 (all anno2) |
| **Methods** | Spearman, LASSO, RF, limma-voom |
| **Top-N** | 50 per method per compartment |

#### Per-Compartment Age Signal Strength

| Compartment | Cells | Patients | Spearman sig | limma sig | Top Gene (rho) |
|---|---|---|---|---|---|
| **T cell** | 280K | 96 | **898** | **935** | SOX4 (−0.68) |
| **Myeloid** | 116K | 96 | **117** | **60** | SPI1 (+0.51) |
| **NK** | 97K | 96 | **111** | **70** | CDKN2A (+0.61) |
| DC | 17K | 96 | 59 | 11 | SEMA5A (−0.49) |
| B cell | 24K | 96 | 9 | 20 | RHEX (−0.46) |
| Platelet | 8K | 96 | 0 | 0 | RHOC (+0.46) |
| HSPC | 853 | 46 | 0 | 0 | (underpowered) |
| ISG-hi | 4K | 48 | 1 | 0 | (underpowered) |
| Megakaryocyte | 2K | 77 | 0 | 0 | (underpowered) |
| RBC | 938 | 34 | 0 | 2 | (underpowered) |

#### Compartment-Specific Aging Signatures

**T cell** (strongest signal, 898 sig genes):
- SOX4 (−0.68), BLK (−0.62), ROBO1 (−0.62), CASC15 (−0.63) — naive/developmental genes **decrease**
- ENSG00000237471 (+0.58), ENSG00000288727 (+0.57) — **increase**

**NK** (senescence/exhaustion signal):
- CDKN2A (+0.61), IGFBP3 (+0.51), MSC (+0.55), LAG3 (+0.46) — **senescence markers increase with age**
- DAPK2 (+0.54), ZBTB38 (+0.46) — maturation markers

**Myeloid** (activation with age):
- SPI1 (+0.51), EMP3 (+0.49), APH1B (+0.48) — myeloid activation/differentiation increase
- SND1 (−0.44), ADAMTS5 (−0.47) — matrix-related decrease

**B cell** (developmental decline):
- RHEX (−0.46), VAV1 (−0.44), TSHR (−0.43) — signaling/development decrease
- PARP6 (+0.38), EMP3 (+0.45) — stress/activation increase

#### Cross-Compartment Patterns (Treg Top Genes)

| Gene | T cell | NK | Myeloid | B cell | Pattern |
|---|---|---|---|---|---|
| **JAML** | **−0.60*** | −0.09 | −0.32* | −0.03 | T cell-dominant decrease |
| **MSC** | **+0.56*** | **+0.55*** | −0.05 | +0.13 | Lymphoid increase |
| **ISG15** | −0.27* | −0.27* | −0.16 | −0.33* | Pan-immune decrease |
| **PLAC8** | −0.36* | −0.32* | +0.11 | −0.19 | Lymphoid decrease |
| **IL15** | **+0.35*** | +0.16 | +0.09 | +0.28* | Homeostatic increase |
| **ABCF2** | **−0.34*** | −0.27* | **−0.33*** | −0.16 | Pan-immune decrease |
| **DGKA** | **−0.36*** | +0.09 | −0.13 | −0.14 | T cell-specific decrease |

**Pan-compartment gene**: ENSG00000237471 — consensus in NK, Platelet, T cell (mean rho=+0.47)

### 3.9 Treg Exhaustion Markers × Age

| Item | Value |
|------|-------|
| **Output** | `/data/user3/sobj/hc_only_v1/treg_exhaustion_age/` |
| **Markers tested** | 23 requested, 15 found in data, 8 missing |
| **Subclustering** | 7 subclusters (scVI, res=0.5) |

#### Missing Genes (8/23)

FCEB2, OSBPL6, IL1B2, XIBP1, INS3, DHBS2, UBE2QL1, CD276 — possible gene symbol typos (IL1B2→IL1B?, INS3→GINS3?)

#### Sparsity Issue

10/15 found markers have >70% zero patients at pseudobulk level — they are **not typically expressed in Treg**. Only 5 markers are analyzable: NR4A3, NR4A1, GEM, EGR2, PLEKHA7 (<60% zero).

#### Results Summary

| Analysis | Finding |
|---|---|
| Spearman (all 15 markers) | **None** reached p<0.05 |
| limma-voom (sex-adjusted) | **ITGB8 only**: logFC/yr=+0.020, p=0.006 (but 77% zero) |
| Age group (KW: Y vs M vs O) | **No markers** significant |
| Sex-stratified | **ITGB8 in males**: rho=+0.404, p=0.011 |
| Subclustering (7 clusters) | **1/54 tests** p<0.05: Cl1 × GEM (rho=−0.21, padj=0.79) |

**Conclusion**: These exhaustion markers show **no meaningful age association** in HC Treg. This marker set likely derives from CD8 exhaustion or disease-context studies, not healthy Treg aging.

#### Treg Subclusters (scVI, res=0.5)

| Cluster | Cells | Top Markers | Character |
|---|---|---|---|
| 0 | 2,856 | FOXP3, RTKN2, IL2RA, CCR6 | Classic Treg |
| 1 | 2,447 | ARMH1, NOSIP, KLRB1, TSHZ2 | Cytotoxic-like |
| 2 | 748 | PLK1, DLGAP5, CDC20, CCNA2 | Proliferating |
| 3 | 621 | IGFBP7, NCAM1, NCR1, KLRF1 | NK-like / innate |
| 4 | 605 | FXYD2, IL10, SESTD1 | IL10+ suppressive |
| 5 | 553 | CD8A, CD8B, CCL5 | CD8 contaminant? |
| 6 | 370 | LINC02694, FTX, IKZF2 | lncRNA-hi |

---

## 4. Data File Locations

```
/data/user3/sobj/hc_only_v1/
├── 1_hc_seurat.qs                    (4.7G)  # Raw Seurat
├── 2_hc_annotated.qs                 (4.7G)  # Final annotated (USE THIS)
├── hc_only_v1_annotation_map.csv              # Cluster → anno1 mapping
├── hc_only_v1_markers_res2.csv                # FindAllMarkers
│
├── mrvi/                                      # MrVI results
│   ├── mrvi_latent.h5ad              (2.0G)
│   ├── mrvi_model/model.pt
│   └── plots/                                 # 11 PNGs
│
├── deg_consensus_1/                           # DEG Sex (3 methods × 32 types)
│   ├── results_consensus.qs          (15.2G)
│   ├── method_results/                        # Raw per-method CSVs
│   └── [32 cluster subdirs]/                  # Per-type consensus + plots
│
├── milo/                                      # MILO DA (sex)
│   ├── milo_sex_da_results.csv
│   └── milo_sex_*.png                         # Volcano, boxplot, barplot
│
├── age_continuous/                            # Cell type proportion × age
│   ├── spearman_age_anno1.csv
│   └── plots/
│
├── fgs_continuous_treg/                       # FGS Treg v1 (4 methods, top-50)
│   ├── pseudobulk_data.rds           (7.1M)  # Reusable pseudobulk
│   ├── consensus_genes_2plus_methods.csv      # 45 genes
│   └── fgs_continuous_results.rds
│
├── fgs_continuous_treg_v2/                    # FGS Treg v2 (6 methods, top-200)
│   ├── consensus_genes_3plus_methods.csv      # 84 genes
│   ├── consensus_genes_4plus_methods.csv      # 32 genes
│   └── fgs_continuous_v2_results.rds
│
├── fgs_compartment/                           # FGS by compartment × age
│   ├── fgs_compartment_results.rds   (14.6M)
│   ├── fgs_top_genes_across_compartments.csv
│   ├── heatmap_fgs_top_genes_across_compartments.png
│   └── [10 compartment subdirs]/              # Per-compartment results
│
├── cellchat/
│   ├── sex_v1/                                # CellChat M vs F
│   │   └── merged/{F,M}/cellchat_merged.qs
│   └── age_v1/                                # CellChat Young vs Old
│       └── merged/{Young,Old}/cellchat_merged.qs
│
├── mnn/
│   ├── sex_v1/multinichenet_results.qs (372M) # MNN F-M
│   └── age_v1/multinichenet_results.qs (338M) # MNN Old-Young
│
├── treg_analysis/                             # Basic Treg × age
│   └── treg_marker_age_results.csv
│
└── treg_exhaustion_age/                       # Exhaustion markers + subclustering
    ├── treg_subclustered.qs           (107M)
    ├── treg_subcluster_markers.csv
    ├── subcluster_age_associations.csv
    └── exhaustion_markers_age_stats.csv
```

---

## 5. Scripts

| Script | Purpose | Location |
|---|---|---|
| FGS continuous v1 | 4-method Treg × age | `/tmp/run_fgs_continuous.R` |
| FGS continuous v2 | 6-method, top-200 | `/tmp/run_fgs_continuous_v2.R` |
| FGS by compartment | All compartments × age | `/tmp/run_fgs_by_compartment.R` |
| Treg exhaustion markers | 23 markers × age deep | `/tmp/treg_exhaustion_markers_age.R` |
| Treg deep analysis | Q1–Q4 combined | `/tmp/treg_exhaustion_deep.R` |
| MILO pipeline fix | `reduced.dim` bug fix | `myR/R/analysis/milo_pipeline.R` |
| DEG consensus CLI | 3-method pseudobulk DEG | `Git_Repo/_wt/deg-consensus/scripts/consensus/run_deg_consensus_cli.R` |
| CellChat CLI | Per-sample split-merge | `Git_Repo/_wt/cellchat/scripts/cellchat/run_cellchat_cli.R` |
| MNN CLI | Pseudobulk LR analysis | `Git_Repo/_wt/cci/scripts/cci/mnn/run_multinichenet.R` |

---

## 6. Key Findings Summary

### Age Effects — Big Picture

1. **T cell compartment shows the strongest age signal** (898 significant genes): naive/developmental genes (SOX4, BLK, ROBO1) decrease; maturation markers increase
2. **NK cells show clear senescence signature** with age: CDKN2A, IGFBP3, LAG3 increase
3. **Myeloid activation increases** with age: SPI1 (master myeloid TF), EMP3
4. **ISG genes decline** pan-immune with age (ISG15, PLAC8, EPSTI1) — innate immune signaling weakening
5. **Treg-specific exhaustion markers show no age signal** — these markers are not relevant to healthy Treg aging
6. **MSC (musculin)** increases strongly in both T cells and NK with age — potential aging biomarker

### Sex Effects — Big Picture

1. **MILO**: 117 significant neighbourhoods (85 M-enriched, 32 F-enriched) — moderate DA
2. **DEG consensus**: CD14 Mono Classical has the most sex-differential genes (435)
3. **Pan-immune sex genes**: RPS4X, EIF2S3, JPX, CD99 (escape X-inactivation genes)

### Limitations & Caveats

- n=96 is moderate power for detecting subtle age effects (especially in small cell types)
- HC cohort age range (19–66) is not extreme; effects may be larger in elderly (>70)
- LASSO R² is inflated by overfitting (p >> n problem); RF OOB R² is more realistic
- Small compartments (HSPC, RBC, Megakaryocyte) are underpowered (<50 patients)
- CCI results (CellChat, MNN) have not yet been interpreted in detail

---

## 7. Pending / Next Steps

| Priority | Task | Dependencies |
|---|---|---|
| 1 | **CCI interpretation**: Parse CellChat + MNN results for age/sex, Treg focus | Results completed, need parsing |
| 2 | **MASC**: Run for anno1 × sex and anno1 × age_group | Need to re-run (no saved output) |
| 3 | **Treg deep FGS**: Run FGS specifically on Treg subclusters (e.g., classic vs IL10+) | `treg_subclustered.qs` available |
| 4 | **Pathway enrichment**: HALLMARK/KEGG/GOBP on compartment FGS consensus genes | FGS results available |
| 5 | **Cross-dataset**: Compare HC age genes with stroke_hc_v8_2 DEG results | Both datasets available |
| 6 | **Descriptive figures**: UMAP, DotPlot, frequency barplots for HC cohort | Ready to generate |

---

## 8. Known Issues & Fixes

| Issue | Fix |
|---|---|
| MILO `testNhoods` "PCA not found" | Added `reduced.dim = "GRAPH"` to both `calcNhoodDistance()` and `testNhoods()` |
| ranger protection stack overflow | Use `ranger(x=, y=)` interface instead of formula with >18K genes |
| `GEM` gene vs metadata column conflict | Exclude GEM from `DotPlot(features=)` in Treg analysis |
| patient_meta missing `age_group` | Compute from age: `<35 → Young, 35–50 → Middle, >50 → Old` |
| Small compartments (HSPC, ISG-hi) NA age | Filter patients with `is.na(age)` before analysis |
| DEG consensus output dir auto-increment | Check both `deg_consensus/` and `deg_consensus_1/` |
