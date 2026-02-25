# Progress Status — stroke_hc_v8_2

> Updated: 2026-02-20 | 전체 분석 진행 현황 및 시각화 결과물 정리

---

## 완료된 분석 (Completed)

### Phase 1: Pre-processing & Integration
- Pipeline (CellBender / Souporcell / Solo / scVI) ✅
- Clustering (res=2, scvi_snn, 36 clusters) ✅
- Annotation (anno1: 21 types, anno2: 8 compartments) ✅
- Annotation QC (doublet/LQ/RBC 제거 → 205,277 clean cells) ✅
- Subset: `5_1_hc_is.qs` (131K), `5_2_is_g3.qs` (54K) ✅

### Phase 2: Differential Abundance (MASC)
- anno2 × cohort (HC vs Stroke) ✅
- anno1 × cohort ✅
- anno2 × g3 (IS only) ✅
- anno1 × g3 ✅
- anno2 × cohort (no GEM) ✅
- anno1 × cohort (no GEM) ✅

### Phase 3: DEG (Consensus & NEBULA)
- **L1 DEG Consensus**: muscat-edgeR + muscat-DESeq2 + NEBULA (3 methods) ✅
  - 15/20 clusters with consensus (5 skipped: insufficient methods)
- **L2 DEG**: NEBULA only (203K results) ✅ — muscat too sparse for g3
- **Cross-layer concordance**: 134 genes significant in both L1 and L2 ✅
- **Pathway enrichment** (FGS genes): HALLMARK, KEGG, GOBP ✅
- **FGS × DEG overlap**: 218/241 FGS genes significant in L1 DEG ✅

### Phase 4: CCI (Cell-Cell Interaction)

> 상세 문서: [`docs/claude/cci/CCI_ANALYSIS_REPORT.md`](cci/CCI_ANALYSIS_REPORT.md)

#### CellChat v1 (Method 2A: sample-wise merged — original)
- L1 (HC vs Stroke, anno2): 56 samples ✅
- L2 (g3, anno2): 32 samples ✅
- Comparison plots + pathway plots 생성됨, 단 **rankNet/bubble 미생성** (merged @data 빔)

#### CellChat v2 (condition-level from Seurat — 재구현) ✅
- **Seurat에서 condition별 직접 CellChat 생성** → @data 정상, 모든 비교 함수 작동
- **rankNet (information flow)**: L1 + L2 모두 성공 ✅ — 핵심 신규 결과
  - L1: HC-dominant (TGFb, IL2, IL16) vs Stroke-dominant (RESISTIN, MIF, BAFF)
  - L2: Bad-dominant for nearly all pathways (IL1, RESISTIN, FLT3, MIF)
- **Bubble comparison**: L-R pair별 조건 비교 ✅
- **Diff circles/heatmaps**: netVisual_diffInteraction + netVisual_heatmap ✅
- **Signaling role heatmaps** (10-11): per-condition outgoing/incoming centrality ✅
- **Signaling scatter** (12): per-celltype signaling changes between conditions ✅
- **Signaling role scatter** (13): outgoing vs incoming per celltype ✅
- **Pathway chord + contribution**: 9 common pathways × 2 conditions = 36 PNG per layer ✅
- **L1**: 16 main + 36 pathway + 7 scatter = **59 PNG**
- **L2**: 18 main + 36 pathway + 5 scatter = **59 PNG**
- Output: `cci/plots/cellchat_L1_cohort_v2/`, `cellchat_L2_g3_v2/`
- Scripts: `scripts/reimplement_cellchat_v2.R`, `scripts/fix_cellchat_v2_final.R`, `scripts/fix_cellchat_heatmaps.R`

#### MNN (MultiNicheNet) — BUG FIXED & COMPLETE ✅
- **이전**: 6회 실행 모두 `group_prioritization_tbl = 0 rows`
  - 원인: contrast string format mismatch (`"Stroke-HC"` vs `"Stroke - HC"`)
  - 오진: "lr_target_prior_cor sample ≥5 불충족" (실제로는 sample 충분)
- **수정**: `contrasts_oi`를 apostrophe로 wrapping → `makeContrasts`가 string literal로 처리
- **v3 결과** (2026-02-19):

| Config | Cell Types | Prioritized Interactions |
|--------|-----------|-------------------------|
| L1 anno2 (HC vs Stroke) | 5 | 8,244 |
| L1 anno1 (HC vs Stroke) | 18 | 148,958 |
| L2 anno2 (g3=1 vs g3=2) | 5 | 7,600 |
| L2 anno1 (g3=1 vs g3=2) | 17 | 91,111 |

- Output: `cci/mnn/{L1_cohort_anno2_v3, L1_cohort_anno1_v3, L2_g3_anno2_v3, L2_g3_anno1_v3}/`
- 상세: `docs/claude/cci/MNN_FAILURE_ANALYSIS.md`
- **MNN anno1 시각화** (intermediates 기반): L1 6 PNG + 1 PDF, L2 6 PNG + 1 PDF ✅

### Phase 5: FGS (Feature Gene Signature)

#### Whole-dataset FGS (n=50)
- 10 methods: 9/10 succeeded (nmf_loadings FAIL) ✅
- **Ranger best**: AUC 0.713 (xgboost 0.591, elastic_net 0.552)
- **TML**: Ranger ROC 0.869, xgbTree 0.863, GLM 0.804 ✅
- **CMGI**: 241 genes, top: MT-CO1 (-9.21), RPS26 (+6.58), HLA-DQA2 (-6.52), IFI44L (+5.84) ✅
- **시각화**: 8 plots ✅

#### FGS sweep
- n=100: 완료 ✅
- n=200: 완료 ✅ (464 → TBD genes)

#### Within-cell-type FGS (n=50) ✅
- **전체 15 cell types 완료**
- Output: `results/fgs/within_celltype/{celltype}_50_cmgi_genes.csv`
- Cell types: B_cell, CD14_Monocyte, CD16_Monocyte, CD4_S100A8_CD14, CD4_S100A8_CSF3R, CD4_T_Naive_Memory, CD8_T_Cytotoxic, CD8_Trm, cDC2, Inflammatory_Monocyte, ISG_Myeloid, MAIT, NK_cell, Plasma_cell, Treg

### Phase 6: Trajectory

#### Trajectory v2 (Slingshot + Monocle3)
- scVI UMAP, 3 compartments: Mono (72K), CD4 (51K), CD8 (35K) ✅

| | Monocyte | CD4+ T | CD8+ T |
|---|---------|--------|--------|
| Pseudotime | 0–2.69 | 0–2.06 | 0–3.34 |
| Lineages | 5 | 6 | 6 |
| Cohort p | <2e-16 | <2e-16 | 5.16e-13 |
| **g3 p** | **<2e-16** | **0.028** | **0.749 (NS)** |

- Gene dynamics v1 (GAM, no batch correction): 36 genes × 3 compartments ✅
- Gene dynamics v2 (GAMM with batch correction): `analyze_gene_dynamics_v2()` 구현 완료 ✅

#### Gene Dynamics v3 (Batch-corrected GAMM) ✅
- **전체 완료**: 3 compartments × 2 conditions (cohort, g3) = 6 runs
- Model: `expr ~ s(pt) + cond + s(pt, by=cond) + offset(log(nCount)) + pct.mt + s(GEM, bs="re")`
- 각 run: 36 genes → 36 PNG + summary CSV
- Output: `results/trajectory_v3/{mono,cd4,cd8}/gene_dynamics_{cohort,g3}/`

#### Trajectory Effect Sizes ✅ (Phases 1-4 완료)
- Phase 1-2: ABC/rABC/RMISE effect sizes + comparison figures → `trajectory_v3/analysis/`
- Phase 3: Lamian validation (chisq method) → `trajectory_v3/lamian/`
- Phase 4: Pseudobulk GAMM sensitivity → `trajectory_v3/pseudobulk/`
- **핵심 발견**: 10 genes confirmed by all 3 methods (mono/cohort); Mono = ONLY compartment with g3 effect
- Docs: `docs/claude/trajectory/` (4 documents)

### Phase 7: MILO (Neighbourhood DA) ✅
- **L1 (HC vs Stroke)**: 11,739 nhoods, 63% DA (p<0.1)
- **L2 (g3 Good vs Bad)**: 4,664 nhoods, 19% DA (p<0.1)
- **3-method concordance (MILO × MASC × scCODA)**: 14/19 cell types DA by ≥2 methods (L1)
- L2: CD14+ Mono only DA (MILO + scCODA concordant)
- Output: `results/milo/{L1_cohort,L2_g3}/`
- Concordance: `results/milo/L1_concordance_3methods.csv`

### Phase 8: Descriptive Figures
- 63 files v1 (UMAP, DotPlot, bars, frequency, heatmaps, pseudotime, FGS) ✅
- 6 files v2 (paper figures: concordance heatmap, patient FGS scores, pseudobulk heatmap, DEG overlap, MILO dotplot) ✅
- Output v2: `results/figures/v2_paper/`

### Phase 9: Compositional Analysis

#### scCODA (Bayesian Compositional) ✅
- **L1 (HC vs Stroke)**: 14/21 cell types credible effect
  - Credible (composition significantly different): CD14+ Mono, CD16+ Mono, CD4+ T, CD4_S100A8_CSF3R, CD8+ T, CD8+ Trm, ISG+ Myeloid, ISG+ T_cell, MAIT, NK_cell, Treg, cDC2, Platelet/PLA, pDC
  - Not credible: B_cell, CD4_S100A8_CD14, Inflammatory Mono, Plasma_cell, Proliferating, cDC1, Mast_cell
  - Fallback (Mann-Whitney): Platelet/PLA most significant (p=9.7e-12, FC=14.2), MAIT (p=1.1e-10, FC=6.1)
- **L2 (g3 Good vs Bad)**: Only **CD14+ Mono** credible effect
  - All other cell types: not credible (FDR=1.0)
  - Fallback: No cell type significant after FDR correction
- **Interpretation**: L1 composition changes robust; L2 composition differences minimal (g3 effect is expression-level, not composition)
- Output: `results/sccoda/`

### Phase 10: External Validation ✅
- **3 independent bulk datasets** validated via ssGSEA scoring
- **FGS_TOP50** (top 50 CMGI genes): consistent across all 3 datasets
  - GSE16561 (Illumina, n=63, whole blood): **AUC=0.725**, p=0.002, d=0.92
  - GSE22255 (Affy, n=40, PBMC): **AUC=0.677**, p=0.056, d=0.62
  - GSE58294 (Affy, n=92, CE stroke blood): **AUC=0.689**, p=0.007, d=0.71
- **FGS_TOP25_DOWN** (top 25 downregulated genes): strongest discriminator
  - GSE16561: **AUC=0.817**, p=1.1e-5, d=1.26
  - GSE58294: **AUC=0.888**, p=2.9e-8, d=1.65
- **Interpretation**: scRNAseq-derived FGS signature generalizes to bulk platforms, confirming biological validity
- Output: `results/external_validation/` (4 PNG + 4 CSV)

---

### Phase 11: CCI × DEG Cross-Reference ✅
- **CellChat pathways × DEG consensus**: 71 CCI genes, 16 (22.5%) overlap with DEG
- 7 unique DE genes in CCI pathways: TGFB1 (ligand), IL21R, IL2RA, CXCR4, CD44, FPR3, TNFRSF1B (receptors)
- 6 pathways with DE genes: IL2, MIF, ANNEXIN, GALECTIN, TGFb, TNF
- Output: `results/cci_deg_crossref/` (6 plots + 3 CSV + summary)

### Phase 12: cNMF (Gene Programs) ✅
- **16/16 cell types consensus completed**
- K values: 11 for most (12 types), k=8 for B_cell, Inflammatory Monocyte, MAIT, Plasma_cell
- Programs: 8-11 per cell type, 20,644 genes scored
- Total: 202,368 cells × 8-11 programs each
- **Fix**: Original script had k=10 hardcoded (not in factorized range [5,8,11,14,17,20]); resolved via k_selection plot inspection
- Output: `results/cnmf/{CELLTYPE}/` + `cnmf_consensus_summary.csv`

### Phase 13: MELD (Condition Likelihood) ✅
- **L1 (HC vs Stroke)**: Per-cell P(Stroke) computed (205K cells)
  - Highest P(Stroke): CD14+ Mono (0.889), Inflam Mono (0.869), ISG+ Myeloid (0.865)
  - Lowest: Platelet/PLA (0.393), ISG+ T_cell (0.518), pDC (0.597)
- **L2 (IS Good vs Bad)**: Per-cell P(Bad) computed (54K IS cells)
  - Highest P(Bad): CD14+ Mono (0.792), CD16+ Mono (0.760)
  - Lowest: pDC (0.595), Plasma_cell (0.610)
- **Concordant with MASC/MILO/scCODA**: Monocytes consistently most condition-associated
- Output: `results/meld/` (4 PNG + 6 CSV)

### Phase 14: Augur (Cell Type Prioritization) ✅
- **L1 (HC vs Stroke)**: Expression-based AUC per cell type
  - Top: cDC1 (0.859), CD16+ Mono (0.856), CD14+ Mono (0.855), CD4_S100A8_CSF3R (0.835)
  - Bottom: Platelet/PLA (0.609), CD4_S100A8_CD14 (0.617), Proliferating (0.622)
- **L2 (IS Good vs Bad)**: All AUCs near 0.5 (minimal expression separability)
  - Top: CD4_S100A8_CSF3R (0.624), CD14+ Mono (0.602)
  - Most types: 0.50-0.59 (near random)
- **Interpretation**: L1 has strong expression shifts (AUC 0.6-0.86); L2 shows subtle g3 effects
- Output: `results/augur/` (3 PNG + 3 CSV)

### Phase 15: propeller (Compositional DA) ✅
- **L1 + L2**: Logit-transformed proportions + limma
- Output: `frequency/propeller/` (2 PNG + 2 CSV)

### Phase 16: Subclustering ✅
- **3 compartments**: Monocyte (12 SC), CD4 T (17 SC), CD8 T (16 SC)
- **Subcluster biological annotation**: ✅ All subclusters named with biological identities
  - Monocyte: Classical (3 states), Inflammatory (2), Non-classical (3), S100A-high, contamination (2)
  - CD4: TEM, Naive, TCM, NK-like (2), CTL (2), Th1 (2), Treg, gdT-like, contamination (2)
  - CD8: TEM (4 states), Naive (2), TEMRA, NKT-like, NK-like, TCM, MAIT, gdT, contamination (2)
- **Subcluster frequency**: L1 30/45 significant, L2 0/45 (confirms g3 = expression not composition)
- Output: `subclustering/{Monocyte,CD4_T,CD8_T}/` (annotations, UMAPs, DotPlots, annotated .qs)

### Phase 17: cNMF GEP Naming ✅
- **155 GEPs across 16 cell types named** with biological pathway labels (GSEA-based)
- 92/155 significant in L1, 49/155 in L2
- Top Stroke-up: CD16+ Mono IFN response (logFC=3.25), ISG+ Myeloid TNFa (logFC=2.03)
- Top Stroke-down: CD14+ Mono TGFb/TNFa (logFC=-3.72)
- Output: `unbiased/cnmf/downstream/gep_named.csv`, `plots/gep_*.png`

### Phase 18: CCI Method Integration ✅
- **CellChat × MNN × LIANA cross-referenced**
- L1: 10 double-concordant interactions, 2 triple genes (ANXA1, FPR1)
- Top pathways (3-method): MIF (1.54), GALECTIN (1.39), RESISTIN (1.05), ANNEXIN (1.03)
- Pairwise: CellChat-LIANA r=0.578, MNN-LIANA r=0.545
- Output: `cci/plots/integrated_summary/` (9 PNG + 7 CSV + summary report)

### Phase 19: Folder Reorganization ✅
- **Entire stroke/ directory restructured** into 7 categories:
  - `frequency/` (MASC, milo, meld, augur, sccoda, propeller, scanvi)
  - `cci/` (cellchat, mnn, liana, cci_deg_crossref, integrated_summary)
  - `unbiased/` (cnmf, fgs, mofa, summary)
  - `deg/` (consensus, pathway, FindMarkers, external_validation)
  - `trajectory/` (v1, v2, v3)
  - `subclustering/` (Monocyte, CD4_T, CD8_T, frequency)
  - `pipeline/` (P0_PIPE, integration, annotation_qc)
- Scripts organized: `scripts/{cci,deg,descriptive,frequency,pipeline,trajectory,unbiased,claude}/`

### Phase 20: Paper Figure Strategy ✅
- **7 main figures + 10 supplementary + 11 tables** designed
- Figure strategy: `docs/FIGURE_STRATEGY.md`
- v3 paper figures: 🔄 Generation in progress (`figures/v3_paper/`)

---

## 보류 / 계획 (Pending)

| 분석 | 우선순위 | 비고 |
|------|---------|------|
| propeller / DirichletReg | Low | Additional compositional methods |
| ~~Subclustering~~ | ~~Medium~~ | ✅ 완료 — 3 compartments |
| ~~LIANA~~ | ~~Low~~ | ✅ 완료 — 5-method CCI consensus |
| ~~MOFA+~~ | ~~Low~~ | ✅ 완료 — 4 factors, Factor1=cohort |
| ~~scANVI~~ | ~~Low~~ | ✅ 완료 — 97.3% label agreement |
| ~~CCI × DEG cross-reference~~ | ~~High~~ | ✅ 완료 |
| ~~MELD~~ | ~~Low~~ | ✅ 완료 |
| ~~Augur~~ | ~~Low~~ | ✅ 완료 |
| ~~cNMF re-run~~ | ~~Medium~~ | ✅ 16/16 cell types consensus completed |
| ~~cNMF downstream~~ | ~~Medium~~ | ✅ 완료 — GSEA annotation + condition association |
| ~~External validation~~ | ~~Medium~~ | ✅ 완료 — 3 datasets |

---

## 시각화 결과물 위치

모든 결과물은 `results/` symlink → `/data/user3/sobj/stroke_hc_v8_2`

### CCI Plots
```
results/cci/plots/
├── cellchat_L1_cohort/          # 49 files (16 comparison + 33 pathway)
│   ├── 01-08_comparison_plots   # Circle, heatmap, boxplot (PNG + PDF)
│   ├── 11-12_bubble_*.png       # L-R bubble plots (merged)
│   └── pathways/                # pathway circle + contribution plots
├── cellchat_L2_g3/              # 53 files (16 comparison + 37 pathway)
│   └── (same structure)
├── mnn_L1_cohort_anno1/         # 6 PNG + 1 PDF + CSVs
│   ├── 01_de_summary_bar.png
│   ├── 02_volcano_grid.png
│   ├── 03_ligand_activity_grid.png
│   ├── 04_top_ligand_receiver.png
│   ├── 05_circos_ligand_receiver.{pdf,png}
│   └── 06_ligand_receiver_heatmap.png
├── mnn_L1_cohort_anno2/         # 11 files
├── mnn_L2_g3_anno1/             # 6 PNG + 1 PDF + CSVs
└── mnn_L2_g3_anno2/
```

### FGS Plots
```
results/fgs/plots/               # 8 PNG
├── 01_method_auc_comparison.png
├── 02_tml_model_comparison.png
├── 03_tml_all_metrics.png
├── 04_signature_importance.png
├── 05_cmgi_top_genes.png
├── 06_cmgi_gene_method_heatmap.png
├── 07_meta_score_distribution.png
└── 08_meta_score_by_celltype.png
```

### Trajectory Plots
```
results/trajectory_v2/           # 117 PNG (v2: GAM, no batch correction)
├── mono/
│   ├── trajectory_overview_scvi.png
│   ├── pseudotime_comparison.png
│   ├── pseudotime_by_celltype.png
│   └── gene_dynamics/           # 36 GAM plots
├── cd4/
└── cd8/

results/trajectory_v3/           # 216 PNG (v3: batch-corrected GAMM)
├── mono/
│   ├── gene_dynamics_cohort/    # 36 PNG + summary CSV
│   └── gene_dynamics_g3/        # 36 PNG + summary CSV
├── cd4/
│   ├── gene_dynamics_cohort/
│   └── gene_dynamics_g3/
└── cd8/
    ├── gene_dynamics_cohort/
    └── gene_dynamics_g3/
```

### scCODA
```
results/sccoda/
├── l1_composition.csv           # Patient × cell type count matrix (L1)
├── l2_composition.csv           # Patient × cell type count matrix (L2)
├── sccoda_L1_cohort_credible.csv  # Bayesian credible effects (14/21 TRUE)
├── sccoda_L1_cohort_fallback.csv  # Mann-Whitney fallback with FDR
├── sccoda_L2_g3_credible.csv      # Bayesian credible effects (1/21 TRUE)
└── sccoda_L2_g3_fallback.csv      # Mann-Whitney fallback with FDR
```

### MILO
```
results/milo/
├── L1_cohort/
│   ├── milo_object.rds           # 21GB milo object
│   ├── da_results.csv            # 11,739 nhood DA results
│   ├── 01_da_beeswarm.png
│   ├── 02_da_umap.png
│   ├── 03_da_volcano.png
│   └── 04_da_summary_barplot.png
├── L2_g3/                        # same structure (4,664 nhoods)
├── L1_concordance_3methods.csv   # MILO × MASC × scCODA cross-method
├── L1_milo_vs_masc_scatter.png
└── L1_concordance_heatmap.png
```

### Paper Figures v2
```
results/figures/v2_paper/
├── 01_frequency_da_concordance_heatmap.png  # 3-method DA concordance
├── 02_patient_fgs_score_cohort.png          # FGS score HC vs IS by compartment
├── 03_patient_fgs_score_g3.png              # FGS score Good vs Bad by compartment
├── 04_fgs_top30_patient_heatmap.png         # Top 30 genes pseudobulk heatmap
├── 05_cross_layer_upset_bar.png             # L1∩L2∩FGS overlap
├── 05_triple_overlap_genes.csv              # 12 genes in triple overlap
└── 06_milo_dual_layer_dotplot.png           # MILO L1+L2 dotplot
```

### External Validation
```
results/external_validation/
├── 01_validation_auc_comparison.png    # AUC bar chart (5 gene sets × 3 datasets)
├── 02_validation_score_boxplot.png     # FGS_TOP50 boxplot HC vs Stroke per dataset
├── 03_validation_roc_curves.png        # ROC curves (3 datasets overlaid)
├── 04_validation_forest_plot.png       # Effect size forest plot
├── validation_results_summary.csv      # Full results table (30 rows)
├── GSE16561_ssgsea_scores.csv          # Per-sample scores
├── GSE22255_ssgsea_scores.csv
└── GSE58294_ssgsea_scores.csv
```

### cNMF (16/16 cell types consensus) ✅
```
results/cnmf/
├── cnmf_consensus_summary.csv    # Summary: cell_type, k, n_cells, n_programs
├── CD4plus_T_Naive_Memory/       # k=11, 37,715 cells × 11 programs
│   ├── CD4plus_T_Naive_Memory.h5ad
│   └── CD4plus_T_Naive_Memory/
│       ├── *.k_selection.png
│       ├── *.usages.k_11.dt_0_1.consensus.txt
│       └── *.gene_spectra_score.k_11.dt_0_1.txt
├── Inflammatory_Monocyte/        # k=8, 29,226 cells × 8 programs
├── NK_cell/                      # k=11, 25,467 cells × 11 programs
└── ... (16 cell types total)
```

### Phase 15: cNMF Downstream (GEP Annotation + Condition Association) ✅
- **GSEA annotation**: fgsea on 4 gene set collections (HALLMARK, KEGG, REACTOME, GOBP)
  - 200,736 total GSEA results, **20,561 significant** enrichments (p.adj < 0.05)
  - 155/176 GEPs successfully annotated (≥1 significant pathway)
- **Condition association**: GEP usage ~ cohort (L1) and g3 (L2) via Wilcoxon test
  - L1 (cohort): 144 condition-associated GEPs
  - L2 (g3): 102 condition-associated GEPs
- **Key findings**:
  - Inflammatory Monocyte GEP5 = TNFa/NFkB signaling (HALLMARK NES=3.81, FDR=9.65e-50)
  - CD14+ Monocyte GEP4 = TNFa/NFkB signaling (HALLMARK NES=3.87, FDR=8.23e-56)
  - Multiple immune activation programs enriched in Stroke/Bad outcome
- Output: `results/cnmf/downstream/` (plots/, gsea/, condition_association/)
```
results/cnmf/downstream/
├── gsea/
│   ├── gsea_all_results.csv         # 200,736 total results
│   ├── gsea_significant.csv         # 20,561 significant enrichments
│   └── gsea_top_per_gep.csv         # Top pathway per GEP
├── condition_association/
│   ├── gep_condition_L1_cohort.csv  # 144 associated (Wilcoxon)
│   └── gep_condition_L2_g3.csv      # 102 associated
├── gep_annotations.csv              # 155 annotated GEPs
└── plots/
    ├── 01_gep_hallmark_heatmap.png  # NES heatmap (GEP × HALLMARK)
    ├── 02_gep_condition_L1.png      # GEP usage by cohort
    ├── 03_gep_condition_L2.png      # GEP usage by g3
    └── 04_top_genes_*.png           # Top spectra genes (5 cell types)
```

### Phase 16: MOFA+ (Patient-Level Factor Analysis) ✅
- **Multi-view**: 6 views (Bc, DC, Mono, NKc, Platelet-PLA, Tc) × 79 patients
- **4 active factors** (out of 15 requested; 11 removed as inactive)
- **Factor 1 = Cohort effect** (p=3.5e-11, Wilcoxon)
  - 26.5% variance in Mono, 17.3% in DC, 15.9% in NKc, 13.3% in Tc, 13.0% in Bc
  - Top Mono weights: LINC-PINT, MBP, CDKN1A (positive/HC); S100A8, S100A9, S100A12 (negative/Stroke)
  - Factor-age correlation: r=-0.345
- **Factor 2 = Platelet-PLA-specific** (19.8% variance; no cohort/g3 association)
  - NKG7, RPL3, RPS18 (positive); TUBB1, GP9 (negative) → platelet-leukocyte aggregate signature
- **No g3 association** for any factor (all p > 0.7)
- Output: `results/mofa/`
```
results/mofa/
├── mofa_model.hdf5                # Raw MOFA model
├── mofa_model.qs                  # R-ready MOFA model
├── mofa_L{1,2}_*.csv              # Factor values, weights, variance, associations
├── weights_{Bc,DC,Mono,NKc,Platelet-PLA,Tc}.csv  # Per-view gene weights
├── 01_variance_explained.png      # Variance explained per factor × view
├── 01b_total_variance.png         # Total variance per view
├── 02_factors_by_cohort.png       # Factor distributions HC vs Stroke
├── 03_factors_by_g3.png           # Factor distributions Good vs Bad
├── 04_top_weights_*.png           # Top gene weights per factor (6 views)
└── 05_factor_correlations.png     # Factor correlation matrix
```

### Phase 17: Subclustering (Deep Resolution) ✅
- **3 compartments**: Monocyte, CD4 T, CD8 T
- Uses original `umap.scvi` embeddings (no recomputation)
- Multiple resolutions tested (0.3, 0.5, 0.8, 1.0); working res=0.5

| Compartment | Cells | Subclusters | Markers | Top cohort subcluster |
|-------------|-------|-------------|---------|----------------------|
| **Monocyte** | 61,624 | 12 | 8,784 | SC1 (OR=956, HC-enriched), SC7 (OR=8.0, HC) |
| **CD4 T** | 51,214 | 17 | 9,362 | SC5 (OR=148, HC), SC6 (OR=44.7, HC) |
| **CD8 T** | 35,295 | 16 | 7,509 | SC1 (OR=177, HC), SC2 (OR=94.1, HC) |

- **g3 associations (IS only)**:
  - Mono: SC2 (OR=403, Bad), SC3 (OR=193, Bad), SC0 (OR=0.15, Good)
  - CD4 T: SC3 (OR=18.6, Bad), SC10 (OR=11.9, Bad), SC8 (OR=0.23, Good)
  - CD8 T: SC8 (OR=20.6, Bad), SC11 (OR=8.3, Bad), SC6 (OR=0.17, Good)
- Output: `results/subclustering/{Monocyte,CD4_T,CD8_T}/`
```
results/subclustering/{compartment}/
├── 01_umap_overview.png           # UMAP by subcluster + anno1 + cohort
├── 02_resolution_comparison.png   # res=0.3/0.5/0.8/1.0
├── 03_top_markers_heatmap.png     # Top 5 markers per subcluster
├── 04_frequency_barplot.png       # HC vs Stroke subcluster frequency
├── 05_frequency_g3_barplot.png    # Good vs Bad subcluster frequency
├── 06_dotplot_key_genes.png       # Compartment-specific marker DotPlot
├── 07_anno1_composition.png       # anno1 composition per subcluster
├── markers_all.csv                # FindAllMarkers results
└── {compartment}_subclustered.qs  # Seurat object with subclusters
```

### Phase 18: LIANA (Multi-Method CCI Consensus) ✅
- **5 methods**: natmi, connectome, sca, logfc, cellphonedb
- Rank aggregation via `liana_aggregate()`
- **L1 (HC vs Stroke)**:
  - HC: 6,755 interactions; Stroke: 3,825 interactions (Mast_cell removed, <5 cells)
  - Top HC: Tc→Platelet/PLA IL32-ITGB3, PF4-LRP1
  - Gained in Stroke: HLA-B→KLRD1 (DC→NKc), B2M→KIR2DL1 (Tc→NKc), S100A9→TLR4 (Mono→DC)
  - Lost in Stroke: Mast_cell interactions (removed), SPARC-ENG (Platelet→Mono)
- **L2 (Good vs Bad)**:
  - Good: 3,987 interactions; Bad: 3,845 interactions
  - Key differences: PF4-LDLR lost in Bad, SPARC-ENG gained in Bad, SNCA-LAG3 gained in Bad
- Output: `results/liana/`
```
results/liana/
├── L1_{HC,Stroke}_anno2_liana_agg.csv   # Aggregated ranks
├── L1_{HC,Stroke}_anno2_liana_raw.rds   # Full LIANA output
├── L1_{HC,Stroke}_anno2_top50.csv       # Top 50 interactions
├── L1_comparison_rank_diff.csv          # Rank difference analysis
├── L1_rank_comparison_scatter.png       # HC vs Stroke rank scatter
├── L2_{Good,Bad}_anno2_liana_agg.csv
├── L2_{Good,Bad}_anno2_liana_raw.rds
├── L2_{Good,Bad}_anno2_top50.csv
├── L2_comparison_rank_diff.csv
└── L2_rank_comparison_scatter.png
```

### Phase 19: scANVI (Label-Aware Integration) ✅
- **scVI base model**: 200 epochs, Final ELBO 4215.6 (n_latent=30, n_layers=2, NB likelihood)
- **scANVI**: 100 epochs from scVI base, unlabeled_category="Unknown"
- **Label agreement**: **97.3%** overall (scANVI predictions match anno1)
  - Highest: cDC1 (100%), cDC2 (98.9%), CD14+ Mono (98.7%), CD8+ T_Cytotoxic (98.2%)
  - Lowest: Mast_cell (86.4%, small N), CD4_S100A8_CSF3R (92.0%), CD4_S100A8_CD14 (94.0%)
- **Silhouette scores** (50K subsample):
  - Label: scVI=0.004, scANVI=**0.043** (10x better cell type separation)
  - Batch: scVI=-0.023, scANVI=-0.052 (both good, scANVI slightly better mixing)
- **Interpretation**: Annotation highly consistent — scANVI independently recovers 97.3% of manual labels
  - Low-agreement types (CD4_S100A8_*) confirm these are ambiguous/transitional populations
  - Validates anno1 quality for all downstream analyses
- Output: `results/scanvi/`
```
results/scanvi/
├── scanvi_model/                  # Trained scANVI model
├── adata_scanvi.h5ad              # 7.3 GB annotated h5ad
├── scanvi_latent.csv              # 205K × 30 latent dimensions
├── scanvi_predictions.csv         # Predicted vs actual labels
├── scanvi_per_type_agreement.csv  # Per cell type agreement
├── scanvi_metrics.csv             # Silhouette + agreement metrics
├── 01_scanvi_umap_overview.png    # UMAP: anno1, cohort, prediction
├── 02_scvi_vs_scanvi_umap.png    # scVI vs scANVI UMAP comparison
├── 03_confusion_matrix.png       # Prediction confusion matrix
└── 04_silhouette_comparison.png  # Label + batch silhouette bars
```

### scCODA Plots ✅
```
results/sccoda/plots/
├── 01_composition_stacked_L1.png    # Cell type proportions stacked bar
├── 02_credible_effects_dotplot.png  # L1 vs L2 credible effects
├── 03_volcano_L1.png               # Fold change × significance
├── 04_composition_heatmap_L1.png   # Per-patient composition heatmap
└── 05_composition_boxplot_L2.png   # Top 10 cell type boxplots (Good vs Bad)
```

### CCI × DEG Cross-Reference ✅
```
results/cci_deg_crossref/
├── cellchat_deg_crossref_L1.csv    # 6 crossref hits (pathway-gene-role)
├── pathway_deg_summary.csv         # 4 pathways with DE genes
├── gene_cci_involvement_summary.csv # 7 unique CCI-DE genes
├── crossref_summary.txt            # Full summary report
└── plots/
    ├── 01_pathway_deg_heatmap.png
    ├── 02_pathway_compartment_dotplot.png
    ├── 03_ligand_receptor_volcano.png
    ├── 04_top_de_ligands.png
    ├── 05_top_de_receptors.png
    └── 06_cci_deg_overlap.png
```

### MELD ✅
```
results/meld/
├── meld_L1_p_stroke.csv           # Per-cell P(Stroke) (205K cells)
├── meld_L1_summary_by_celltype.csv # Mean P(Stroke) per cell type
├── meld_L2_p_bad.csv              # Per-cell P(Bad) (54K IS cells)
├── meld_L2_summary_by_celltype.csv # Mean P(Bad) per cell type
└── plots/
    ├── 01_meld_L1_umap.png        # UMAP colored by P(Stroke)
    ├── 02_meld_L1_violin.png      # P(Stroke) violin per cell type
    ├── 03_meld_L2_umap.png        # UMAP colored by P(Bad)
    └── 04_meld_L2_violin.png      # P(Bad) violin per cell type
```

### Augur ✅
```
results/augur/
├── augur_L1_auc.csv               # L1 AUC per cell type (20 types)
├── augur_L2_auc.csv               # L2 AUC per cell type (17 types)
├── augur_combined_auc.csv         # Both layers combined
├── augur_L1_lollipop.png          # L1 AUC lollipop plot
├── augur_L2_lollipop.png          # L2 AUC lollipop plot
└── augur_L1_vs_L2_scatter.png     # L1 vs L2 AUC comparison scatter
```

---

## 핵심 발견 요약

1. **SIIS (Stroke-Induced Immunosuppression)**: HC 대비 Stroke에서 CCI -37%, DC hub 붕괴, cDC1 594 DE genes
2. **SIIS = Selective immunosuppression** (rankNet 신규 발견): 전체 CCI 감소하나 RESISTIN/MIF/CypA는 Stroke에서 증가 → 선택적 면역억제 + 염증 활성화 병존
3. **g3 Bad outcome**: 거의 모든 CCI pathway 증가 (IL1 Bad-only, RESISTIN Bad>>Good) → maladaptive broad immune activation
4. **Trajectory**: g3 effect는 Mono에서 최강 (p<2e-16), CD4 약함 (0.028), CD8 없음 (0.749)
5. **FGS meta-genes**: MT-CO1, HLA-DQA2 (DOWN in Bad), RPS26, IFI44L (UP in Bad)
6. **scCODA**: L1 14/21 cell types compositionally different; L2 only CD14+ Mono → g3 effect는 composition이 아닌 expression level
7. **Gene dynamics v3**: 3 compartments × 2 conditions, batch-corrected GAMM 전체 완료
8. **Trajectory effect sizes**: 10 genes confirmed by ALL 3 methods (ABC, Lamian, pseudobulk GAMM) in mono/cohort. Mono is ONLY compartment with g3 effect
9. **MILO 3-method concordance**: 14/19 cell types DA by ≥2 methods (MILO × MASC × scCODA) in L1. L2: only CD14+ Mono DA
10. **Cross-layer DEG overlap**: L1∩L2: 82 genes, L1∩FGS: 189 genes, triple overlap: 12 genes
11. **External validation**: FGS TOP50 signature (g3-derived) replicates in 3 independent bulk datasets (AUC 0.68-0.73); TOP25_DOWN achieves AUC up to 0.888
12. **MELD concordance**: P(Stroke) highest in CD14+ Mono (0.889), Inflam Mono (0.869) — consistent with MASC/MILO/scCODA. P(Bad outcome) also peaks in CD14+ Mono (0.792)
13. **Augur**: L1 strongly separable (AUC 0.61-0.86), L2 near random (AUC 0.50-0.62) — g3 effect is subtle, not global transcriptome shift
14. **cNMF downstream**: 155/176 GEPs annotated; Inflam Mono GEP5 = TNFa/NFkB (NES=3.81); 144 L1 + 102 L2 condition-associated programs
15. **CCI × DEG**: 22.5% of CCI pathway genes are DE; key DE signaling genes: TGFB1 (ligand), CXCR4/CD44/IL21R (receptors)
16. **MOFA+**: Factor1=Cohort (p=3.5e-11), 26.5% Mono variance; S100A8/A9/A12 drive Stroke signal; No g3 factor found
17. **Subclustering**: Bad-outcome–enriched subclusters found in all 3 compartments: Mono SC2/SC3 (OR>190), CD4 SC3 (OR=18.6), CD8 SC8/SC11 (OR>8)
18. **LIANA**: 5-method CCI consensus confirms HLA→KIR and S100A9→TLR4 as stroke-gained interactions; PF4-LDLR lost in bad outcome
19. **scANVI**: 97.3% label agreement validates annotation quality; scANVI latent gives 10x better cell type silhouette (0.043 vs 0.004) while maintaining batch mixing

---

## docs/claude/ 문서 목록

| # | 파일 | 내용 | 작성일 |
|---|------|------|--------|
| 1 | `1_STROKE_HC_V8_2_CONTEXT.md` | 전체 분석 컨텍스트, 데이터, 파이프라인 | 2026-02-13 |
| 2 | `2_gene_dynamics_approaches.md` | GAMM 모델 설계, two-stage pseudobulk 계획 | 2026-02-13 |
| 3 | `3_REFACTORING_PLAN.md` | 코드베이스 리팩토링 계획 | 2026-02-14 |
| 4 | `4_PROJECT_GOALS.md` | 논문 구성 분석 11개 항목, 의존성, 우선순위 | 2026-02-16 |
| 5 | `5_ANALYSIS_REPORT.md` | MASC, DEG, frequency 등 분석 결과 리포트 | 2026-02-17 |
| 6 | `6_ADDITIONAL_ANALYSES_PLAN.md` | 추가 7개 분석 실행 계획 | 2026-02-17 |
| 7 | `7_PROGRESS_STATUS.md` | **현재 문서** — 전체 진행 현황 | 2026-02-18 |
| — | `cci/CCI_ANALYSIS_REPORT.md` | CCI 전용: 결과, 해석, 문제점, 향후 계획 | 2026-02-18 |
