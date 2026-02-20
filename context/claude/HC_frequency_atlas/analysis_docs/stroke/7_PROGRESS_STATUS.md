# Progress Status — stroke_hc_v8_2

> Updated: 2026-02-18 | 전체 분석 진행 현황 및 시각화 결과물 정리

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

#### MNN (MultiNicheNet) — 근본적 한계 확인
- anno2 + anno1 수준 모두 실행 ✅, **permissive thresholds로도 재실행** ✅
- **`group_prioritization_tbl` = 0 rows** (모든 설정에서)
  - 근본 원인: `lr_target_prior_cor` step에서 cell type당 sample ≥5 불충족
  - multinichenetr 패키지가 이 데이터와 **근본적으로 비호환**
- Intermediate tables 활용 가능: DE genes (78K L1), ligand activities (1.2M L1)
- **MNN anno1 시각화** (intermediates 기반): L1 6 PNG + 1 PDF, L2 6 PNG + 1 PDF ✅
- **향후**: LIANA (multi-method CCI consensus) 대체 검토

### Phase 5: FGS (Feature Gene Signature)

#### Whole-dataset FGS (n=50)
- 10 methods: 9/10 succeeded (nmf_loadings FAIL) ✅
- **Ranger best**: AUC 0.713 (xgboost 0.591, elastic_net 0.552)
- **TML**: Ranger ROC 0.869, xgbTree 0.863, GLM 0.804 ✅
- **CMGI**: 241 genes, top: MT-CO1 (-9.21), RPS26 (+6.58), HLA-DQA2 (-6.52), IFI44L (+5.84) ✅
- **시각화**: 8 plots ✅

#### FGS sweep
- n=100: 완료 ✅
- n=200: 진행 중 🔄

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

### Phase 7: Descriptive Figures
- 63 files (UMAP, DotPlot, bars, frequency, heatmaps, pseudotime, FGS) ✅

### Phase 8: Compositional Analysis

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

---

## 진행 중 (In Progress)

| 분석 | 상태 | 세부 | 비고 |
|------|------|------|------|
| FGS n=200 sweep | 🔄 실행 중 | Background | 자동 진행 |
| Within-cell-type FGS/TML | 🔄 부분 완료 | CD4 T 완료 | 15 cell type subsets 준비됨, 나머지 실행 필요 |
| cNMF | 🔄 K selection 완료 | 3 cell types | CD4 T, Inflam Mono: k_selection done (311 files), NK: factorizing (230 files) |
| MASC anno1 × project_name | 🔄 실행 중 | Low priority | |

### cNMF Status Detail

| Cell Type | Prepare | Factorize | K Selection | Consensus |
|-----------|---------|-----------|-------------|-----------|
| CD4+ T_Naive/Memory | ✅ | ✅ (311 runs) | ✅ (plot generated) | Pending |
| Inflammatory Monocyte | ✅ | ✅ (311 runs) | ✅ (plot generated) | Pending |
| NK_cell | ✅ | 🔄 (230/311 runs) | Pending | Pending |

### Within-Cell-Type FGS Status

| Cell Type | Subset | FGS | TML | CMGI |
|-----------|--------|-----|-----|------|
| CD4+ T_Naive/Memory | ✅ | ✅ | ✅ | ✅ |
| CD14+ Monocyte | ✅ | Pending | | |
| CD16+ Monocyte | ✅ | Pending | | |
| Inflammatory Monocyte | ✅ | Pending | | |
| CD8+ T_Cytotoxic | ✅ | Pending | | |
| NK_cell | ✅ | Pending | | |
| B_cell | ✅ | Pending | | |
| ISG+ Myeloid | ✅ | Pending | | |
| CD4_S100A8_CD14 | ✅ | Pending | | |
| CD4_S100A8_CSF3R | ✅ | Pending | | |
| CD8+ Trm | ✅ | Pending | | |
| cDC2 | ✅ | Pending | | |
| MAIT | ✅ | Pending | | |
| Plasma_cell | ✅ | Pending | | |
| Treg | ✅ | Pending | | |

---

## 보류 / 계획 (Pending)

| 분석 | 우선순위 | 비고 |
|------|---------|------|
| MILO v2 (DA testing) | Medium | PCA fix 필요 (scVI→PCA mapping) |
| MNN permissive thresholds | Medium | group_prioritization empty → 더 관대한 filter 재실행 |
| CCI × DEG cross-reference | High | MNN DE genes와 DEG Consensus/NEBULA 교차 비교 |
| Subclustering (Mono, CD4, CD8) | Medium | anno1 세분화 |
| MELD (predicted likelihood) | Low | Per-cell condition density estimation |
| Augur (cell type prioritization) | Low | Expression-based discriminability |
| MOFA+ (multi-omics factor) | Low | Patient-level latent factors |
| scANVI (reference mapping) | Low | 추가 분석 |
| External validation | Low | 외부 데이터셋 검증 (GSE16561, GSE140275 등) |
| propeller / DirichletReg | Low | Additional compositional methods |

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

### cNMF
```
results/cnmf/
├── CD4plus_T_Naive_Memory/
│   ├── CD4plus_T_Naive_Memory.h5ad
│   └── CD4plus_T_Naive_Memory/
│       ├── CD4plus_T_Naive_Memory.k_selection.png
│       └── cnmf_tmp/             # 311 factorization runs
├── Inflammatory_Monocyte/        # same structure
└── NK_cell/                      # factorization in progress
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
