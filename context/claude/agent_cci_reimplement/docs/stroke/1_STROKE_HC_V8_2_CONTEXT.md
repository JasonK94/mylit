# Stroke_HC_v8_2 Analysis Context Document

> 생성일: 2026-02-12 | 목적: 분석 컨텍스트 요약 및 실행 계획

## 1. 프로젝트 개요

**Study**: Stroke vs Healthy Control scRNAseq PBMC analysis
- **Stroke**: 12 GEM (10X 5' v2 ST), ~80 patients (IS/ICH/IVH/SAH)
- **HC (Healthy Control)**: 5 GEM (10X 5' v2 HT from 4 HC_GEM), ~20 subjects
- **총 세포수**: 226,855 cells / 20,644 genes / 100 subjects

## 2. 데이터 출처 및 파이프라인

### Pipeline
- **Repository**: `/home/user3/data_user3/git_repo/mypy_re/`
- **Config**: `config/config_local_stroke_hc_v8_2.json`
- **Log**: `logs/pipeline_stroke_hc_v8_2.log` (Pipeline Complete ✓)

### Pipeline Steps (완료됨)
1. `demux_classify` - Souporcell (with known_genotypes VCF) → 모든 GEM1~8 rerun
2. `merger` - CellBender reuse + demux doublet filter
3. `qc_filter` - min_genes=200, min_cells=3, max_mt_pct=20, solo soft filter
4. `solo` - Doublet detection (100 epochs)
5. `integration` - scVI (GPU/RTX4090, n_hvgs=3000, n_latent=30, 200 epochs)

### Souporcell Demux 확인
- GEM1~8: VCF-based rerun, **Singlets only** (0% doublets in `_notroublet.csv`)
- GEM9~12: No demux path (pre-multiplexed/HT) → Direct barcode whitelist from CellRanger
- HC_GEM04~24: Souporcell with VCF, **doublets detected** (6.6~10.9%)
- **74_3 = GEM8의 sample_id '4-3'** → QC metrics에서 확인됨: raw 33,534 → solo_singlet 442 cells

## 3. 데이터 파일 위치

```
/data/user3/sobj/stroke_hc_v8_2/
├── 1_stroke_hc_integrated.h5ad     (7.6G) - 원본 scVI integration output
├── 2_strokev8_meta_add.h5ad        (7.7G) - metadata 추가
├── 2_strokev8_meta_add.qs          (533M) - Seurat 변환
├── 3_strokev8_clustered.qs         (9.3G) - Clustering 완료 (res=2, scvi_snn)
├── 4_strokev8_annotated.qs         (9.3G) - anno1/anno2 annotation
├── 5_strokev8_clean.qs             (~8G)  - LQ/RBC/doublet 제거 (205,277 cells)
├── 5_1_hc_is.qs                    (~5G)  - HC + IS only (131,552 cells) [Layer 1]
├── 5_2_is_g3.qs                    (~2G)  - IS only with g3 (54,002 cells) [Layer 2 + FGS]
├── results/
│   ├── FindMarkers/
│   │   ├── 1_fam.csv               - FindAllMarkers 결과
│   │   └── 1_fam_filtered.csv      - Filtered markers
│   └── annotation_qc/              - Annotation QC 결과 (2026-02-12)
│       ├── centroid_distances_scvi.csv
│       ├── cluster_nearest_neighbors.csv
│       ├── canonical_marker_dotplot_data.csv
│       ├── doublet_coexpression.csv
│       ├── suspicious_coexpression.csv
│       ├── cluster_cohort_distribution.csv
│       ├── cluster_g3_distribution.csv
│       └── patient_summary.csv
├── fgs/                             - FGS & TML 결과 (실행 중)
│   └── whole_IS_{50,100,200}/       - n_features sweep
├── cci/                             - CCI 결과 ✅
│   ├── cellchat/
│   │   ├── L1_cohort_anno2/         - HC vs IS (anno2 grouping)
│   │   │   └── merged/{Stroke,HC}/cellchat_merged.qs
│   │   └── L2_g3_anno2/            - g3==2 vs g3==1 (anno2 grouping)
│   │       └── merged/{X1,X2}/cellchat_merged.qs
│   └── mnn/
│       ├── L1_cohort_anno2/multinichenet_results.qs
│       └── L2_g3_anno2/multinichenet_results.qs
├── trajectory/                      - Trajectory 결과 ✅
│   ├── mono/                        - Monocyte (72K cells)
│   │   ├── subset.qs, slingshot_sce.rds, monocle3_cds.rds
│   │   └── gene_dynamics/*.png
│   ├── cd4/                         - CD4+ T (51K cells)
│   │   └── subset.qs, slingshot_sce.rds, monocle3_cds.rds
│   └── cd8/                         - CD8+ T (35K cells)
│       └── subset.qs, slingshot_sce.rds, monocle3_cds.rds
├── MASC/                            - MASC 결과 (6 comparisons done)
├── scripts/                         - 분석 스크립트
│   ├── annotation_qc.R
│   ├── create_subsets.R
│   ├── run_fgs_sweep.sh
│   ├── run_cci_all.sh
│   ├── run_trajectory.R
│   └── run_trajectory_part2.R
└── logs/                            - 실행 로그
```

## 4. Metadata 구조

### 핵심 컬럼
| Column | Description | Unique | NAs |
|--------|-------------|--------|-----|
| `patient_name` | 환자 ID (hos_no) | 100 | 0 |
| `GEM` | GEM batch ID | 17 | 0 |
| `SET` | SET1/2/3, HC_SET1/3 | 5 | 0 |
| `cohort` | Stroke / HC | 2 | 0 |
| `index_injury_hand` | IS/ICH/IVH/SAH/HC | 6 | 0 |
| `project_name` | Stroke / Healthy | 2 | 0 |
| `g3` | **Target variable** (1/2) | 3 | 167,373 ⚠️ |
| `sex` | M/F | 3 | 7,525 |
| `age` | 25~80+ | 48 | 7,525 |
| `day` | 1/2/3 (sampling day) | 4 | 80,823 |
| `anno1` | Fine annotation | 26 | 0 |
| `anno2` | Major compartments | 8 | 0 |
| `seurat_clusters` | Seurat clusters (res=2) | 36 | 0 |

### ⚠️ 중요: g3 NA 문제
- `g3`는 IS 환자 전용 → HC/ICH/IVH/SAH 세포는 NA
- FGS/TML에서 g3를 target_var로 사용 시, `index_injury_hand=="IS"` 필터 필요
- 또는 HC를 포함한 분석 설계 시 `cohort` 또는 `project_name` 사용

### anno1 categories (26개)
- T/NK: CD4+ T_Naive/Memory, CD8+ T_Cytotoxic, CD8+ Trm, Treg, MAIT, NK_cell, ISG+ T_cell, Proliferating
- Myeloid: CD14+ Monocyte, Inflammatory Monocyte, CD16+ Monocyte, ISG+ Myeloid, Macrophage-like
- DC: cDC1, cDC2, pDC
- B: B_cell, Plasma_cell
- Contam/LQ: Erythrocyte, Platelet/PLA, Mast_cell, T/NK_Low_quality_14, Mono_Low_quality_24, B_Low_quality_31
- Ambiguous: CD4_S100A8_CD14, CD4_S100A8_CSF3R

### anno2 categories (8개)
Tc, NKc, Mono, DC, Bc, Platelet/PLA, Erythrocyte, Mast_cell

## 5. 분석 설계 (Study Design)

### 비교 구조 (Dual-layer comparison)
```
Layer 1: HC vs IS (cohort-level)
  - Cell cluster frequency 차이
  - DEG 차이
  - CCI 차이
  → "Stroke가 면역 landscape를 어떻게 바꾸는가?"

Layer 2: IS 내 g3==1 (good) vs g3==2 (bad)
  - 같은 분석을 IS 내에서 반복
  - FGS/TML로 예후 인자(gene signature) 탐색
  → "같은 IS 환자 중 왜 어떤 환자는 예후가 나쁜가?"

Cross-layer insight:
  - Layer 1과 Layer 2에서 유사한 패턴이 발견되면 흥미로운 결과
  - 예: g3==2 (bad outcome)의 면역 프로필이 HC에서 더 벗어나 있는가?
  - "severe stroke" phenotype의 면역학적 특성 규명
```

### g3 정의
- **g3==1**: IS 환자 중 IAT 시행 후 예후 양호 (good outcome)
- **g3==2**: IS 환자 중 IAT 시행 후 예후 불량 (bad outcome)
- **NA**: HC, ICH, IVH, SAH 환자 → g3 해당 없음

### Annotation 전략 (확정)
- **PLT(Platelet) 유지**: PMA(Platelet-Monocyte Aggregate) 관찰 가능성 → PPBP+CD14 co-expression 확인
- **RBC(Erythrocyte) 제거**: 분석에 불필요
- **Low quality clusters 제거**: cluster 14(T/NK_LQ), 24(Mono_LQ), 31(B_LQ)
- **Ambiguous clusters 재평가 필요**:
  - CD4_S100A8_CD14 (cluster 19): T+Myeloid marker co-expression → doublet vs. 실제 세포
  - CD4_S100A8_CSF3R (cluster 25): 동일 이슈
- **Annotation QC 전략**:
  1. Canonical marker expression 확인 (절대적 발현 수준, 상대적 DEG가 아닌)
  2. Mutually exclusive marker co-expression 평가 (CD3 vs CD14 등) → doublet 판별
  3. Cluster centroid 간 거리 계산 → 이상 클러스터 탐지
  4. PERMANOVA / LISI 지표 평가 → batch/integration quality 확인

## 6. 이미 완료된 분석

### ✅ Pipeline (mypy_re)
- CellBender, Souporcell, Solo, scVI integration

### ✅ Clustering & Annotation
- FindAllMarkers (res=2) → `1_fam.csv`
- anno1 (26 types), anno2 (8 compartments) → `4_strokev8_annotated.qs`

### ✅ MASC (Frequency Analysis)
- 6 MASC runs completed:
  - cluster1(seurat_clusters) × comparison1(project_name) / comparison2(g3)
  - cluster2(anno1) × comparison1 / comparison2
  - cluster3(anno2) × comparison1 / comparison2
- Forest plots, p-value bar plots 생성됨

### 🔲 미완료 (실행 필요)
1. **Annotation QC & Refinement** - doublet/LQ 제거, ambiguous cluster 해결
2. **MILO** - 코드 비어있음 (plotting section만)
3. **FGS & TML** - 스크립트 존재, 새 데이터로 실행 필요
4. **CCI (CellChat + MultiNicheNet)** - 스크립트 존재, 새 데이터로 실행 필요
5. **Trajectory Analysis** - compartment별 설계 필요
6. **Subclustering** - 개별 compartment별 scVI re-integration 필요

## 7. 실행 계획

### Phase 1: Annotation QC & Data Preparation
1. **Annotation refinement**:
   - [ ] Cluster centroid 거리 계산 → 이상 클러스터 탐지
   - [ ] Canonical marker 절대 발현량 평가 (pct.exp + avg_exp)
   - [ ] Doublet co-expression 평가 (CD3 vs CD14, CD3 vs CD79A, etc.)
   - [ ] LQ cluster(14,24,31) + Erythrocyte(27,34) 제거
   - [ ] Ambiguous cluster(19,25) 판정: doublet이면 제거, 아니면 re-annotate
   - [ ] PLT 유지 (PMA 가능성)
2. **Subset 생성**:
   - [ ] `sobj_clean`: LQ/Erythrocyte 제거 후 전체 (HC + IS + 기타)
   - [ ] `sobj_hc_is`: HC + IS only (Layer 1 비교용)
   - [ ] `sobj_is`: IS only with valid g3 (Layer 2 비교용)
3. **기초 통계**:
   - [ ] g3 group 분포 (환자 수, 세포 수)
   - [ ] HC vs IS cell composition 비교
   - [ ] PERMANOVA / LISI batch quality 확인

### Phase 2: FGS & TML → 🔄 실행 중 (n=50: 5/10 methods done)
**Script**: `/home/user3/data_user3/git_repo/mylit/Git_Repo/_wt/fgs/scripts/fgs/run_fgs_pipeline.R`
**Sweep script**: `scripts/run_fgs_sweep.sh` (taskset -c 0-15, 16 cores)

#### Run 1: 전체 IS (g3 target) - 실행 중
- `--input`: `5_2_is_g3.qs` (54,002 cells)
- `--output`: `fgs/whole_IS_{50,100,200}/`
- `--target_var g3`, `--control_vars "sex,age,GEM"`
- `--n_features`: 50 (실행 중), 100, 200 (대기)

**진행 상황** (n_features=50):
| Method | Time | Status |
|---|---|---|
| random_forest_ranger | 128 min | ✅ |
| lasso | 91 min | ✅ |
| ridge | 385 min | ✅ |
| elastic_net | 93 min | ✅ |
| pca_loadings | ? | 🔄 실행 중 |
| nmf | - | 대기 |
| gam | - | 대기 |
| limma | - | 대기 |
| wilcoxon | - | 대기 |
| xgboost | - | 대기 |

**주의**: CPU limitation (taskset), NMF convergence, ridge가 매우 느림(6.5h)

### Phase 3: CCI (Dual-layer) → ✅ 완료 + 해석 (2026-02-13)

#### Layer 1: HC vs IS (CellChat)
- HC: mean 114.75 interactions, strength 7.35
- Stroke: mean 72.39 interactions, strength 5.00
- **Stroke에서 전반적 CCI 37% 감소** → Stroke-Induced Immunosuppression (SIIS)
- DC hub 중심 감소 (Mono→DC, DC→DC, NKc→DC 모두 감소)
- HC의 DC는 highest sender+receiver → Stroke에서 DC 기능 저하
- MNN 리간드 활성: IL-10, IL-4 (항염증), CCL11, TGFB2 → 면역억제 milieu
- NLGN1/3 (신경 리간드) → NKc 타겟: 뇌졸중 후 신경 손상 신호?

#### Layer 2: IS g3_Bad vs g3_Good (CellChat)
- g3_Bad: mean 76.24 interactions, strength 5.19
- g3_Good: mean 69.18 interactions, strength 4.71
- **Bad에서 ~10% 더 많은 interaction**
- **Mono→DC, Mono→Bc 증가** (myeloid→APC axis 강화)
- Bc autocrine 증가, NKc 관련 signaling 감소
- MNN 리간드 활성: IFN-I/II/III (IFNB1, IFNG, IL27), MMP1 → 인터페론 반응 signature
- 해석: Bad outcome = 비정상적 myeloid→APC 활성화 + NK 기능 저하

#### MNN 주의사항
- **anno2 레벨에서 group_prioritization_tbl = 0 rows** (both layers)
- DE 유의 유전자: L2에서 1개(Bc), L1에서 0개 → anno2가 너무 coarse
- **ligand_activities_target_de_tbl은 정상 작동** (L2: 105K rows, L1: 229K rows)
- **권장**: anno1 레벨 (26 types)으로 MNN 재실행

#### 출력 파일
```
cci/plots/
├── cellchat_L2_g3/         # 8 comparison plots + pathways/
├── cellchat_L1_cohort/     # 8 comparison plots + pathways/
├── mnn_L2_g3_anno2/        # DE volcanos, ligand activities (LR empty)
└── mnn_L1_cohort_anno2/    # DE volcanos, ligand activities (LR empty)
```

### Phase 4: Trajectory Analysis

#### v1 (trajectory/) → ✅ 완료 (2026-02-12) — **superseded by v2**
- subset 내 PCA/UMAP 재계산 → HC/Stroke batch effect 미보정 문제
- gene dynamics: raw counts, 기술적 보정 없음

#### v2 (trajectory_v2/) → ✅ 완료 (2026-02-13) — scVI UMAP 사용
- **핵심 수정**: 원본 `umap.scvi` embedding 사용 (batch-corrected)
- 모든 scatter plot → PNG (dpi=150)
- HC vs Stroke 비교 (Layer 1) 추가

| Compartment | Cells | Pseudotime | g3 Wilcoxon | Cohort Wilcoxon |
|---|---|---|---|---|
| **Monocyte** | 72,163 | 0 - 2.69 | **< 2e-16** | **< 2e-16** |
| **CD4+ T** | 50,996 | 0 - 2.06 | **0.028** | **< 2e-16** |
| **CD8+ T** | 35,295 | TBD | TBD | TBD |

**핵심 발견**:
- Mono: HC는 중간~후반 pseudotime에 분포 → 정상 분화. Stroke는 초기에 집중 → 분화 정체
- g3_Bad: pseudotime ~0에 집중 → 최초 분화 단계에서 멈춤
- CD4: 유사한 패턴이나 effect size 작음

#### Gene Dynamics v1/v2: `expr ~ s(pseudotime, by=cond)` — batch 보정 없음
- ⚠️ raw UMI counts 사용, GEM batch/library size/percent.mt 미보정
- outlier 때문에 y축 스케일 문제 (S100A8 등)

#### Gene Dynamics v3: Batch-corrected GAMM → 스크립트 준비 완료
- **함수**: `analyze_gene_dynamics_v2()` in `myR/R/analysis/pseudotime.R`
- **모델**: `expr ~ s(pt) + cond + s(pt, by=cond) + offset(log(nCount_RNA)) + percent.mt + s(GEM, bs="re")`
- family: NB(log), REML
- y축: 99th percentile clipping
- 예측 곡선: reference covariate values, batch RE marginalised
- **스크립트**: `scripts/run_gene_dynamics_v3.R`
- **출력**: `trajectory_v3/{comp}/gene_dynamics_{g3,cohort}/`
- **방법론 문서**: `docs/gene_dynamics_approaches.md`

#### Two-stage pseudobulk approach (향후 계획)
- Patient × pseudotime bin으로 aggregation → edgeR/limma-voom
- Patient RE 포함 가능 (pseudoreplication 해결)
- 상세: `docs/gene_dynamics_approaches.md`

#### 참고
- Slingshot: umap_scvi (2D) + approx_points=150
- Monocle3: 원본 umap_scvi 주입, use_partition=FALSE, root_cells 500개 sampling
- `dplyr::first()` masking 주의 (Bioconductor보다 dplyr 마지막 로드)

## 8. 답보상태 이슈

### Issue 1: Annotation 불명확성 → ✅ 해결 (2026-02-12)

**Annotation QC 수행** (centroid distance + canonical marker + co-expression 분석):

| Cluster | anno1 | Cells | 결정 | 근거 |
|---------|-------|-------|------|------|
| 13 | Macrophage-like | 7,435 | **제거 (doublet)** | CD3D+CD14=38%, CD3D+LYZ=48%, NKG7+CD14=46% |
| 22 | Plasma_cell | 3,020 | **제거 (doublet)** | MS4A1+CD14=29%, CD3D+MS4A1=24%, NKG7+CD14=23% |
| 14 | T/NK_LQ | 6,498 | **제거 (LQ)** | |
| 24 | Mono_LQ | 2,810 | **제거 (LQ)** | |
| 31 | B_LQ | 519 | **제거 (LQ)** | |
| 27,34 | Erythrocyte | 1,296 | **제거 (RBC)** | |
| 19 | CD4_S100A8_CD14 | 4,735 | **유지 → CD4+T_S100A8_stress** | CD3D+CD14=0.02%, CD3D=69%. S100A8 stress response |
| 25 | CD4_S100A8_CSF3R | 2,747 | **유지 → T/NK_activated** | CD3D+CD14=0%, NKG7=34%. CSF3R expressing T/NK |
| 28 | Platelet/PLA | 824 | **유지** | PMA 분석용 (PPBP+CD14 co-expression 관찰 가능) |

**제거 총계**: ~21,578 cells → 깨끗한 데이터 ~205,277 cells

### Issue 2: Souporcell/Demux 검증 → ✅ 확인됨
- GEM8의 74_3 (=4-3): 880 demux assigned → 877 QC → **442 solo singlets** (매우 적음)
- GEM9~12: No souporcell demux → CellRanger barcode whitelist 직접 사용
- HC GEMs: Doublet rate 6.6~10.9% → 정상 범위

### Issue 3: g3 NA 대응 → Dual-layer 설계로 해결
- Layer 1: HC vs IS (cohort 기준, g3 불필요)
- Layer 2: IS 내 g3==1 vs g3==2 (IS only subset)
- **g3 분포**: g3==1 (good) 11명, g3==2 (bad) 21명, 총 32 IS patients (59,482 cells)

## 9. 참고 파일

- **Master Rmd**: `/home/user3/GJC_KDW_250721/rmd/AfterPipe_260205.Rmd`
- **myR Package**: `/data/user3/git_repo/mylit/myR/R/`
- **FGS/TML Guide**: `docs/fgs/FGS_TML_INTEGRATED_GUIDE.md`
- **CCI Guide**: `docs/cci/CCI_INTEGRATED_GUIDE.md`
- **Pseudotime Guide**: `docs/pseudotime-dev/PSEUDOTIME_INTEGRATED_GUIDE.md`
- **Clinical Metadata**: `/home/user3/Clinical information/meta_data_prep_260121.csv`
- **QC Metrics**: `/home/user3/data_user3/git_repo/mypy_re/report/qc_sample_metrics.csv`
