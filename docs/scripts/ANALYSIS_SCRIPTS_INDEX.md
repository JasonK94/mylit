# Analysis Scripts Index — stroke_hc_v8_2

> 생성일: 2026-02-16 | 목적: 모든 분석 스크립트의 위치, 용도, 실행 방법 정리

## 1. 프로젝트 실행 스크립트

**위치**: `/data/user3/sobj/stroke_hc_v8_2/scripts/`

이 스크립트들은 stroke_hc_v8_2 데이터셋 전용으로, 파이프라인 함수를 호출하여 실제 분석을 실행함.

### Phase 1: Pre-processing & QC

| 스크립트 | 목적 | 상태 | 입력 | 출력 |
|----------|------|------|------|------|
| `annotation_qc.R` | Annotation QC: centroid distance, canonical marker, doublet co-expression | ✅ 완료 | `4_strokev8_annotated.qs` | `results/annotation_qc/` |
| `create_subsets.R` | LQ/RBC/doublet 제거 후 분석 subset 생성 | ✅ 완료 | `4_strokev8_annotated.qs` | `5_strokev8_clean.qs`, `5_1_hc_is.qs`, `5_2_is_g3.qs` |

### Phase 2: FGS & TML

| 스크립트 | 목적 | 상태 | 비고 |
|----------|------|------|------|
| `run_fgs_sweep.sh` | FGS n=50/100/200 순차 실행 (taskset 16 cores) | 🔄 실행중 | `5_2_is_g3.qs` → `fgs/whole_IS_{50,100,200}/` |
| `plot_fgs_results.R` | FGS 결과 시각화 | 🔲 대기 | n=50 결과 완료 후 |

```bash
# 실행 예시
cd /data/user3/sobj/stroke_hc_v8_2
nohup bash scripts/run_fgs_sweep.sh > logs/fgs_sweep.log 2>&1 &
```

### Phase 3: CCI

| 스크립트 | 목적 | 상태 | 비고 |
|----------|------|------|------|
| `run_cci_all.sh` | CellChat + MNN dual-layer 일괄 실행 | ✅ 완료 | L1(cohort) + L2(g3), anno2 grouping |
| `run_cci_plots.sh` | CellChat comparison + MNN 시각화 | ✅ 완료 | `cci/plots/` |
| `run_mnn_plots.R` | MNN-specific 시각화 (volcano, ligand activity) | ✅ 완료 | |

### Phase 4: Trajectory & Gene Dynamics

| 스크립트 | 목적 | 상태 | 비고 |
|----------|------|------|------|
| `run_trajectory.R` | v1 trajectory (deprecated: subset 내 PCA/UMAP) | ❌ superseded | `trajectory/` → v2로 교체 |
| `run_trajectory_v2.R` | v2 trajectory: scVI UMAP, HC vs Stroke 비교 | ✅/🔄 | Mono/CD4 완료, CD8 진행중 |
| `run_trajectory_part2.R` | CD8 completion + gene dynamics fix | 🔄 | v2의 후속 |
| `run_gene_dynamics_v3.R` | Batch-corrected GAMM (v3) | 🔲 대기 | → v2 CDS 완료 후 실행 |

```bash
# Gene dynamics v3 실행 명령 (v2 완료 후)
cd /data/user3/sobj/stroke_hc_v8_2
nohup taskset -c 16-31 Rscript scripts/run_gene_dynamics_v3.R > logs/gene_dynamics_v3.log 2>&1 &
```

### Phase 5: Frequency Analysis

| 스크립트 | 목적 | 상태 | 비고 |
|----------|------|------|------|
| `masc.sh` | MASC 6 comparisons (3 cluster × 2 contrast) | ✅ 완료 | `MASC/` |

---

## 2. 재사용 가능 파이프라인 스크립트

**위치**: `/data/user3/git_repo/mylit/Git_Repo/_wt/`

데이터셋에 무관하게 사용할 수 있는 범용 CLI 스크립트. 각 스크립트는 `optparse`로 CLI 인자를 받음.

### FGS (Feature Gene Signature)

**스크립트**: `_wt/fgs/scripts/fgs/run_fgs_pipeline.R`

```bash
Rscript run_fgs_pipeline.R \
  --input /path/to/seurat.qs \
  --output /path/to/output/ \
  --target_var g3 \
  --control_vars "sex,age,GEM" \
  --n_features 50 \
  --methods "all"
```

| 인자 | 설명 | 기본값 |
|------|------|--------|
| `--input` / `-i` | Seurat .qs 파일 | (필수) |
| `--output` / `-o` | 출력 디렉토리 | (필수) |
| `--target_var` | 예측 대상 변수 (e.g., g3) | (필수) |
| `--control_vars` | 보정 변수 (쉼표 구분) | NULL |
| `--n_features` | 선택할 feature 수 | 50 |
| `--methods` | 사용할 방법 목록 | "all" |

**Methods (10)**: random_forest_ranger, lasso, ridge, elastic_net, pca_loadings, nmf, gam, limma, wilcoxon, xgboost

### CellChat

**스크립트**: `_wt/cellchat/scripts/cellchat/run_cellchat_cli.R`

```bash
Rscript run_cellchat_cli.R \
  -i /path/to/seurat.qs \
  -a anno2 \
  --subset_aggregate sample_name \
  --group_var cohort \
  --output /path/to/output/
```

**비교 스크립트**: `run_cellchat_comparison_2a.R`
```bash
Rscript run_cellchat_comparison_2a.R \
  --control /path/to/HC/cellchat_merged.qs --control_label HC \
  --disease /path/to/Stroke/cellchat_merged.qs --disease_label Stroke \
  --output /path/to/plots/
```

### MultiNicheNet (MNN)

**스크립트**: `_wt/cci/scripts/cci/mnn/run_multinichenet.R`

```bash
Rscript run_multinichenet.R \
  -i /path/to/seurat.qs \
  -g cohort \
  -s patient_name \
  -c anno2 \
  -f "cohort + anno2" \
  --output /path/to/output/
```

| 인자 | 설명 |
|------|------|
| `-g` | Group variable (contrast column) |
| `-s` | Sample variable (patient ID) |
| `-c` | Celltype variable (annotation) |
| `-f` | Contrast formula |

### MASC

**스크립트**: `_wt/masc/scripts/masc/run_masc.R`

```bash
Rscript run_masc.R \
  --input /path/to/seurat.qs \
  --cluster_var anno1 \
  --contrast_var project_name \
  --random_effects patient_name \
  --fixed_effects "sex,age" \
  --prefix masc_anno1_cohort
```

### DEG Consensus

**스크립트**: `_wt/deg-consensus/scripts/consensus/run_deg_consensus_cli.R`

```bash
Rscript run_deg_consensus_cli.R \
  --input /path/to/seurat.qs \
  --output /path/to/output/ \
  --cluster anno1 \
  --group cohort \
  --contrast "Stroke-HC" \
  --covariates "sex,age,GEM" \
  --methods "MAST,NEBULA,muscat,edgeR,DESeq2" \
  --cores 16
```

---

## 3. myR 패키지 분석 함수

**위치**: `/data/user3/git_repo/mylit/myR/R/analysis/`

### 함수 인벤토리

| 파일 | 함수 | 줄 | Exported | 목적 |
|------|------|-----|----------|------|
| **pseudotime.R** | `run_slingshot_from_seurat()` | 54 | ❌ | Slingshot trajectory |
| | `analyze_gene_dynamics()` | 591 | ❌ | Gene dynamics v1 (raw counts GAM) |
| | `analyze_gene_dynamics_v2()` | 873 | ❌ | Gene dynamics v3 (batch-corrected GAMM) |
| | `process_gene_list_dynamics()` | 1290 | ❌ | Batch gene processor |
| | `run_monocle3_from_seurat()` | 1401 | ❌ | Monocle3 trajectory |
| **deg.R** | `runMAST()` | - | ✅ | MAST DEG |
| | `runMUSCAT()` | - | ✅ | muscat DEG |
| | `runNEBULA()` | - | ✅ | NEBULA DEG |
| **pseudobulk_deg.R** | `prepare_pseudobulk_edgeR()` | 820 | ✅ | Pseudobulk edgeR setup |
| | `run_edgeR_pseudobulk_advanced()` | 1062 | ✅ | edgeR pipeline |
| | `run_DESeq2_pseudobulk_advanced()` | 1247 | ✅ | DESeq2 pipeline |
| | `run_pseudobulk_deg()` | 1768 | ✅ | Meta-wrapper |
| **milo_pipeline.R** | `run_milo_pipeline()` | ~60 | ✅ | MILO DA pipeline |
| **nichenet_analysis.R** | `run_nichenet_analysis()` | - | ✅ | NicheNet workflow |
| **cellchat.R** | (various helpers) | - | - | CellChat wrappers |
| **pathway_enrichment.R** | (ORA/GSEA wrappers) | - | - | Pathway analysis |

### 주의사항

- **pseudotime.R 함수는 NAMESPACE에 export되지 않음** → `source()` 직접 로드 필요:
  ```r
  source("/data/user3/git_repo/mylit/myR/R/analysis/pseudotime.R")
  ```
- `.libPaths()` 설정을 `library()` 전에 반드시 실행:
  ```r
  .libPaths(c(
    "/home/user3/R/x86_64-pc-linux-gnu-library/4.3",
    "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu",
    .libPaths()
  ))
  ```
- dplyr은 Bioconductor 패키지들 이후에 로드 (masking 방지)

---

## 4. 실행 환경 & CPU 제한

```bash
# CPU 제한 (FGS, gene dynamics 등 heavy 작업)
taskset -c 16-31 Rscript script.R

# nohup + log 패턴
nohup taskset -c 16-31 Rscript scripts/run_xxx.R > logs/xxx.log 2>&1 &

# 로그 확인
tail -f /data/user3/sobj/stroke_hc_v8_2/logs/xxx.log

# R output buffering 주의: nohup 로그가 지연될 수 있음
# → 출력 파일의 mtime을 직접 확인하는 것이 더 정확
```

---

## 5. 출력 디렉토리 구조

```
/data/user3/sobj/stroke_hc_v8_2/
├── MASC/                          # Frequency analysis
│   └── {cluster}_{comparison}/    # Forest plots, p-value bars
├── fgs/                           # Gene signature discovery
│   └── whole_IS_{50,100,200}/     # n_features sweep results
├── cci/                           # Cell-cell interaction
│   ├── cellchat/
│   │   ├── L1_cohort_anno2/       # CellChat per-sample objects
│   │   └── L2_g3_anno2/
│   ├── mnn/
│   │   ├── L1_cohort_anno2/       # MNN results
│   │   └── L2_g3_anno2/
│   └── plots/
│       ├── cellchat_L1_cohort/    # CellChat comparison plots
│       ├── cellchat_L2_g3/
│       ├── mnn_L1_cohort_anno2/   # MNN plots
│       └── mnn_L2_g3_anno2/
├── trajectory/                    # v1 (deprecated)
├── trajectory_v2/                 # v2 (current, scVI UMAP)
│   ├── mono/                      # Monocyte trajectory
│   ├── cd4/                       # CD4+ T trajectory
│   └── cd8/                       # CD8+ T trajectory
├── trajectory_v3/                 # v3 gene dynamics output (pending)
│   └── {comp}/gene_dynamics_{g3,cohort}/
├── results/
│   ├── FindMarkers/               # FAM results
│   └── annotation_qc/            # QC matrices
├── scripts/                       # Execution scripts
└── logs/                          # Runtime logs
```
