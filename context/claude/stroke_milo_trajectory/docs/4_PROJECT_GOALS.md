# PROJECT GOALS — stroke_hc_v8_2 논문 분석 계획

> 생성일: 2026-02-16 | 목적: 논문 구성에 필요한 분석 항목 명시 및 달성 추적

## 연구 질문 (Research Questions)

1. **Stroke가 PBMC 면역 landscape를 어떻게 바꾸는가?** (HC vs IS)
2. **같은 IS 환자 중 왜 어떤 환자는 예후가 나쁜가?** (g3==1 Good vs g3==2 Bad)
3. **HC vs IS 차이와 Good vs Bad 차이의 공통 면역 기전이 존재하는가?** (Cross-layer)
4. **핵심 면역 변화를 drive하는 gene signature는 무엇이며, 어떤 세포가 관여하는가?**

---

## 분석 항목 (Analysis Items)

### 1. Pre-processing & Integration

**목적**: scVI integration + annotation 완료, 깨끗한 분석 데이터 준비

| 항목 | 상태 | 비고 |
|------|------|------|
| CellBender + Souporcell + Solo | ✅ 완료 | Pipeline at mypy_re |
| scVI integration (n_latent=30, batch=GEM) | ✅ 완료 | GPU RTX4090 |
| Clustering (res=2, scvi_snn) | ✅ 완료 | 36 clusters |
| Annotation (anno1 26 types, anno2 8 compartments) | ✅ 완료 | |
| Annotation QC (doublet/LQ/RBC removal) | ✅ 완료 | 21,578 cells removed → 205,277 clean |
| Subset 생성 (5_1_hc_is.qs, 5_2_is_g3.qs) | ✅ 완료 | L1: 131K, L2: 54K |

**데이터 파일**: `/data/user3/sobj/stroke_hc_v8_2/5_strokev8_clean.qs`

---

### 2. Descriptive Figures (논문 Figure 1 후보)

**목적**: Integration 품질, annotation 근거, 전체 데이터 개요를 시각적으로 제시

#### 2-1. Integration & Annotation outcome
- [ ] **UMAP** (umap.scvi): anno1 colored, split by cohort (HC/Stroke)
- [ ] **DotPlot**: canonical markers × anno1 → annotation rationale 근거
- [ ] (Optional) UMAP colored by GEM/batch → integration 품질 확인

#### 2-2. Frequency analysis (stacked bar)
- [ ] **Sample-level stacked bar chart**: 모든 sample(~100명)에서 compartment(anno2) 비율
- [ ] **정렬 기준**: signature score 값 내림차순
- [ ] **시각적 설계**:
  - Monocyte → 바닥 (bottom), T cell → 위쪽 (top)
  - Score가 높을수록 해당 compartment 비율이 높다는 것을 시각적으로 암시
- [ ] HC/Stroke 또는 g3 그룹별 colour band 추가

#### 2-3. Signature score × condition (compartment별)
- [ ] **Compartment별로** signature score가 HC vs Stroke, Good vs Bad에서 다른 경우만 선별
- [ ] **Visualization**: boxplot + VlnPlot overlay
- [ ] **Two levels**: cell-level + pseudobulk(patient)-level 둘 다 제시
- [ ] 통계 검정: Wilcoxon (cell-level), patient-level t-test or mixed model

#### 2-4. Signature gene heatmap
- [ ] Signature score를 구성하는 **개별 gene**의 발현값 heatmap
- [ ] **Level options**:
  - Patient-level (pseudobulk mean per patient)
  - Compartment × patient level
  - Cell-level (sampled, 너무 크면)
- [ ] HC vs Stroke, Good vs Bad 패턴 시각적 확인
- [ ] Annotation bar: cohort, g3, anno2 등

#### 2-5. Group-level signature difference
- [ ] Patient-level / compartment-level에서 signature 발현량 boxplot + VlnPlot
- [ ] HC vs Stroke, Good vs Bad 각각
- [ ] 통계 검정 포함 (p-value annotation)

#### 2-6. Cross-layer overlap (HC-vs-IS ∩ Good-vs-Bad)
- [ ] HC vs Stroke DEG/signature와 Good vs Bad DEG/signature의 overlap
- [ ] **Venn diagram** 또는 **UpSet plot**
- [ ] Enrichment: shared genes의 pathway/function 분석
- [ ] 핵심 질문: "Stroke response와 bad outcome의 면역 기전이 겹치는가?"

**전제 조건**: FGS/TML 결과, DEG consensus 결과 필요

---

### 3. (번호 건너뜀 — 2번에 통합)

---

### 4. Frequency Analysis (Detail)

**목적**: HC vs Stroke, Good vs Bad에서 차이 나는 cell type을 robust하게 identify

| 방법 | 상태 | 비고 |
|------|------|------|
| **MASC** | ✅ 완료 | 6 comparisons (3 cluster levels × 2 comparisons) |
| **scCODA** | ✅ 완료 | L1: 14/21 credible, L2: CD14+ Mono only |
| **MILO** | ✅ 완료 | L1: 11,739 nhoods (63% sig), L2: 4,664 nhoods (19% sig) |
| **MELD** | 🔲 미시작 | Python, density estimation approach |
| **scANVI** | 🔲 미시작 | Cell type prediction + label transfer |

#### MILO Results Summary

**L1 (HC vs IS)**: 14/19 cell types DA by ≥2 methods (MILO×MASC×scCODA)
- HC-enriched (11): ISG+ T_cell, Platelet/PLA, CD16+ Mono, pDC, CD8+ Cyto, MAIT, NK, cDC2, Treg, CD8+ Trm, Infl. Mono
- IS-enriched (3): ISG+ Myeloid, B_cell, CD14+ Mono
- MILO–MASC direction correlation: Spearman rho=0.53, p=0.02

**L2 (g3 Good vs Bad)**: Only CD14+ Monocyte shows strong DA (MILO + scCODA concordant)
- CD14+ Mono: logFC=2.62, 57% of nhoods enriched in g3=2 (Bad)

- [x] 3가지 방법으로 frequency DA 수행 (MASC, scCODA, MILO) ✅
- [x] Concordance heatmap + scatter plot 생성 ✅
- [ ] MELD, scANVI (optional)

**스크립트 위치**:
- MASC: `Git_Repo/_wt/masc/scripts/masc/run_masc.R`
- MILO: `stroke_hc_v8_2/scripts/run_milo.R`
**출력**: `/data/user3/sobj/stroke_hc_v8_2/milo/`

---

### 5. CCI (Cell-Cell Interaction) Results

**목적**: HC vs Stroke, Good vs Bad에서 세포 간 소통 패턴 차이 규명

#### 현재 결과 요약

| 분석 | Layer | 결과 | 상태 |
|------|-------|------|------|
| CellChat L1 | HC vs IS | Stroke -37% interaction (SIIS) | ✅ 완료 |
| CellChat L2 | Good vs Bad | Bad +10% interaction, Mono→DC/Bc 증가 | ✅ 완료 |
| MNN L1 | HC vs IS | group_prioritization 빈 테이블 (anno2 too coarse) | ⚠️ 재실행 필요 |
| MNN L2 | Good vs Bad | 동일 문제 | ⚠️ 재실행 필요 |

#### 추가 필요 분석

- [ ] **MNN re-run at anno1** (26 types): DE가 제대로 나오도록 세분화
- [ ] **Circos plot** (메인): 주요 condition-specific LR pair 시각화
- [ ] **Incoming vs Outgoing** analysis: 어떤 cell type이 sender/receiver로 변화하는가
- [ ] **Ligand/Receptor level** analysis: 어떤 cell에서 어떤 signaling molecule을 많이 발현?
  - 실제 DEG에 해당하는가?
  - 그룹 간 차이가 실제로 나는가?
- [ ] **Cell count bias control**: 전체 세포 수 차이에 의한 편향 보정
  - CellChat의 `population.size = TRUE` 옵션
  - 또는 downsampling 후 재실행
- [ ] **Correlation scatter plot** (optional):
  - 예: 세포 내 TXNIP × CXCL8 발현 correlation
  - 또는 A세포 TXNIP → B세포 CXCL8 inter-cellular correlation
  - Pseudobulk level로 계산하면 patient 수준 상관관계

**스크립트 위치**:
- CellChat: `Git_Repo/_wt/cellchat/scripts/cellchat/run_cellchat_cli.R`
- MNN: `Git_Repo/_wt/cci/scripts/cci/mnn/run_multinichenet.R`
- CCI plots: `/data/user3/sobj/stroke_hc_v8_2/scripts/run_cci_plots.sh`

**출력 위치**: `/data/user3/sobj/stroke_hc_v8_2/cci/`

---

### 6. Pathway Analysis

**목적**: DEG 기반 functional enrichment로 생물학적 의미 부여

- [ ] **Gene set**: HALLMARK, KEGG, Reactome, GO-BP
- [ ] **Method**: fgsea / clusterProfiler (ORA + GSEA)
- [ ] **Input**: DEG consensus 결과 (multi-method), FGS signature genes
- [ ] **Comparison**: HC vs Stroke pathway vs Good vs Bad pathway overlap
- [ ] **Visualisation**: dotplot, barplot, network plot (enrichplot)
- [ ] Compartment별 pathway 차이 (Mono, CD4, CD8 등)

**전제 조건**: DEG Consensus 완료 필요

**스크립트**: `Git_Repo/_wt/deg-consensus/scripts/consensus/run_deg_consensus_cli.R` + 별도 enrichment script 필요

---

### 7. Pseudotime Distribution

**목적**: 분화 pseudotime 분포 차이로 면역 세포 maturation 변화 규명

#### 현재 결과

| Compartment | g3 p-value | Cohort p-value | 해석 |
|---|---|---|---|
| **Monocyte** | < 2e-16 | < 2e-16 | Stroke/Bad → 초기 pseudotime 집중 (분화 정체) |
| **CD4+ T** | 0.028 | < 2e-16 | 유사 패턴, effect size 작음 |
| **CD8+ T** | TBD | TBD | v2 실행 필요 |

#### 추가 필요 분석

- [ ] CD8 trajectory v2 완료
- [ ] **FeaturePlot + pseudotime overlay**: pseudotime gradient가 UMAP 어디에 위치?
  - 어느 subset에서 어디로 분화가 퍼져나가는지 추론
  - Canonical marker expression과 pseudotime 위치 대응
- [ ] HC vs Stroke의 pseudotime 분포 차이 → density plot / ridge plot
- [ ] **Branch analysis**: Slingshot lineage branch별 condition 분포 차이

**출력**: `/data/user3/sobj/stroke_hc_v8_2/trajectory_v2/`

---

### 8. Pseudotime-dependent Functional Change (Gene Dynamics)

**목적**: Pseudotime에 따른 유전자 발현 변화가 HC vs Stroke, Good vs Bad에서 다른지 검증

#### Gene Dynamics v3 (Batch-corrected GAMM) — ✅ ALL COMPLETE

**모델**: `expr ~ s(pt) + cond + s(pt, by=cond) + offset(log(nCount_RNA)) + percent.mt + s(GEM, bs="re")`
**Family**: NB(log), REML
**Test**: Interaction term `s(pt, by=cond)`의 유의성 (Chi-squared ANOVA)

| 항목 | 상태 |
|------|------|
| `analyze_gene_dynamics_v2()` 함수 | ✅ 구현 완료 |
| `run_gene_dynamics_v3.R` 실행 스크립트 | ✅ 완료 |
| Mono × g3 + cohort | ✅ 완료 |
| CD4 × g3 + cohort | ✅ 완료 |
| CD8 × g3 + cohort | ✅ 완료 |
| **Effect sizes (ABC/rABC/RMISE)** | ✅ Phase 1-2 완료 |
| **Lamian validation** | ✅ Phase 3 완료 |
| **Pseudobulk GAMM sensitivity** | ✅ Phase 4 완료 |

#### Key Results

| Analysis | GAMM Sig | Lamian Sig | PB Sig | 3-method Consensus |
|----------|----------|------------|--------|--------------------|
| mono/cohort | 21/24 | 14/24 | 16/32 | **10 genes** |
| mono/g3 | 17/23 | 8/23 | 16/30 | ~5 genes |
| cd4/cohort | 20/25 | 13/25 | 10/25 | ~8 genes |
| cd4/g3 | 0/25 | 0/25 | 1/28 | 0 |
| cd8/cohort | 18/23 | 7/23 | 3/25 | ~3 genes |
| cd8/g3 | 1/23 | 0/23 | 0/23 | 0 |

**핵심 발견**: Mono가 유일한 g3 effect compartment; CD8 cell-level 결과는 pseudoreplication artifact
**상세 문서**: `docs/claude/trajectory/4_TRAJECTORY_RESULTS_INTEGRATED.md`
**출력**: `trajectory_v3/{analysis,lamian,pseudobulk}/`

---

### 9. FGS/TML: Gene Signature Discovery

**목적**: HC vs IS, Good vs Bad를 구분하는 gene set을 다각도로 발굴하고, 그 생물학적 의미 탐색

#### 현재 상태
- n=50 sweep 실행 중 (method 5/10: pca_loadings), then n=100, n=200
- Target: g3 (IS only, 54K cells)
- Control vars: sex, age, GEM

#### 추가 분석 (FGS 완료 후)

- [ ] FGS 결과의 **gene overlap** across methods → consensus gene set
- [ ] **Gene expression pattern** 평가:
  - 어떤 cell type이 주로 발현하는가? (DotPlot / heatmap by anno1)
  - Pseudotime에 따른 발현 차이가 유의한가? (gene dynamics v3 결과와 교차)
  - 특정 subcelltype에서 더 높은 발현? → subclustering 후 확인
- [ ] **Functional role** 해석:
  - 염증 진행(pro-inflammatory)에 관여?
  - 면역 억제(immunosuppressive)에 관여?
  - 어떤 pathway에 enriched?
- [ ] HC vs IS signature와 Good vs Bad signature의 유사성 비교 (→ Item 10)

#### 핵심 질문
- 가장 강력한 signature signal은 무엇인가?
- 어떤 gene이 drive하는가?
- 그 gene은 거의 항상 top에 나오는가 (across methods)?
- 무슨 cell type에서 많이 발현되고, 어떤 일을 일으키나?
- 왜 그 세포가 늘어났을까? 증식인가? 어떻게?
  - Cell cycle scoring (Seurat `CellCycleScoring`)
  - Proliferating cell type의 비율 변화
  - Growth factor signaling (CCI 결과 참조)

**스크립트**: `Git_Repo/_wt/fgs/scripts/fgs/run_fgs_pipeline.R`
**출력**: `/data/user3/sobj/stroke_hc_v8_2/fgs/whole_IS_{50,100,200}/`

---

### 10. Cross-layer Concordance

**목적**: HC vs IS signature와 Good vs Bad signature의 유사성 검증

#### 핵심 질문
- HC vs IS signature ≈ Good vs Bad signature인가?
- 공통 gene이 있다면 → shared immune mechanism
- 차이가 있다면 → Stroke-specific vs outcome-specific pathway

#### 분석 계획
- [ ] **Gene-level overlap**: Venn/UpSet (FGS consensus genes)
- [ ] **Pathway-level overlap**: shared enriched pathways
- [ ] **Effect direction**: 동일 gene이 동일 방향으로 변화하는가?
  - Scatter plot: log2FC(HC-vs-IS) × log2FC(Good-vs-Bad)
  - Correlation coefficient + significance
- [ ] **Cell type specificity**: 어떤 cell type에서 공통 변화?
- [ ] **Signature transferability**: L1 signature로 L2 예측 가능한가? (ROC/AUC)

---

### 11. External Validation

**목적**: 주요 finding의 재현성 확인

#### 후보 전략
- [ ] **Public scRNAseq stroke data**: GEO/ArrayExpress에서 유사 cohort 검색
  - PBMC scRNAseq from stroke patients with outcome data
- [ ] **Bulk RNA-seq validation**: TCGA/GTEx 등에서 PBMC data 확보
- [ ] **Signature scoring**: FGS-derived signature를 external cohort에 적용
  - AddModuleScore 또는 AUCell
- [ ] **Clinical correlation**: external cohort에서 signature score ↔ outcome 연관성
- [ ] **Cross-species**: (optional) mouse stroke model transcriptome과 비교

**전제 조건**: 모든 주요 분석 완료, signature gene set 확정 후

---

## 분석 의존성 그래프 (Dependency)

```
1. Pre-processing (✅)
   ├── 2. Descriptive Figures (✅ 63 files, 추가 필요)
   ├── 4. Frequency Analysis
   │     ├── MASC (✅), scCODA (✅), MILO (✅)
   │     └── MELD, scANVI (🔲 optional)
   ├── 5. CCI
   │     ├── CellChat v2 (✅ L1+L2)
   │     └── MNN anno1 (✅ L1+L2 — 78K DEGs, 1.2M ligand activities)
   ├── 6. Pathway Analysis (✅ FGS genes HALLMARK/KEGG/GOBP)
   ├── 7. Pseudotime (✅ v2 all 3 compartments)
   ├── 8. Gene Dynamics (✅ v3 + effect sizes + Lamian + pseudobulk)
   ├── 9. FGS/TML (✅ n=50/100, n=200 🔄, within-celltype 🔄)
   ├── 10. Cross-layer (✅ DEG concordance 134 genes)
   └── 11. External Validation (🔲)
```

---

## 우선순위 (Priority Order)

| 순서 | 항목 | 상태 | 비고 |
|------|------|------|------|
| ~~1~~ | ~~Trajectory v2~~ | ✅ | 3 compartments 완료 |
| ~~2~~ | ~~Gene dynamics v3~~ | ✅ | + effect sizes + Lamian + pseudobulk |
| ~~3~~ | ~~FGS n=50/100~~ | ✅ | n=200 + within-celltype 진행중 |
| ~~4~~ | ~~MNN anno1~~ | ✅ | L1+L2 완료 |
| ~~5~~ | ~~DEG Consensus~~ | ✅ | L1: 3 methods, L2: NEBULA only |
| ~~6~~ | ~~scCODA~~ | ✅ | L1: 14/21 credible |
| ~~7~~ | ~~Frequency: MILO~~ | ✅ | L1+L2 완료, 3-method concordance |
| 8 | Within-celltype FGS | 🔄 | CD14 Mono 실행중 |
| 9 | cNMF | 🔄 | NK factorizing |
| 10 | Descriptive figures 추가 | 🔲 | signature-based plots |
| 11 | External validation | 🔲 | 모든 분석 후 |
