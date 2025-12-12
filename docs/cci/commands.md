# MultiNicheNet Analysis Commands

실행 가능한 MultiNicheNet 분석 커맨드 모음입니다.

## NicheNet Analysis (Classic / Single Condition)

### 1. Run Analysis
To run a standard NicheNet analysis (Sender -> Receiver) for a specific condition or dataset:

```bash
Rscript scripts/cci/NicheNet/run_nichenet.R \
  -i /path/to/seurat_obj.qs \
  --senders "Macrophage,DC" \
  --receiver "CD8_T" \
  --output results/nichenet_run1 \
  --condition_col "group" \
  --condition_oi "Tumor" \
  --condition_ref "Normal" \
  --species human
```

### 2. Visualization
To regenerate plots (Heatmaps, Circos) from saved results (.qs):

```bash
Rscript scripts/cci/NicheNet/plot_nichenet.R \
  -f results/nichenet_run1/nichenet_results.qs \
  -o results/nichenet_run1/plots
```
**Note**: Plots are saved as both PDF and PNG.

---

## MultiNicheNet Analysis

### 1. Visualization (Circos & Diagnostics)

### Circos Plot (Comparison)
Generates comparative Circos plots with customizable filtering.

```bash
# Top 50 interactions, filtered by p-value < 0.05
Rscript scripts/cci/mnn/plot_mnn_circos.R \
  -f results/mnn/multinichenet_results.qs \
  -o results/mnn/plots \
  --top_n 50 \
  --sort_by p_val_adj \
  --p_val_threshold 0.05 \
  --logfc_threshold 0.10
```

### Diagnostic Plots (Bubble & Heatmaps)
Generates bubble plots of top interactions and ligand activity heatmaps.

```bash
Rscript scripts/cci/mnn/plot_mnn_diagnostics.R \
  -f results/mnn/multinichenet_results.qs \
  -o results/mnn/plots_diag
```

### Analysis Scripts (Run MultiNicheNet)
**Note**: We recommend using **1 core** (`--cores 1`) for standard datasets unless the object is extremely large, as parallelization overhead often outweighs the speed gain.

#### Permissive Mode (More Sensitivity)
Use this for detecting subtle interactions (p < 0.2, |logFC| > 0.05).
```bash
Rscript scripts/cci/mnn/run_multinichenet_IS2_IS3_permissive.R --cores 1
```

#### Standard Mode
Strict default thresholds.
```bash
Rscript scripts/cci/mnn/run_multinichenet_IS2_IS3.R --cores 1
```

### 범용 스크립트 (default parameters)
```bash
Rscript /home/user3/data_user3/git_repo/_wt/cci/scripts/cci/mnn/run_multinichenet.R
```

---

## 파라미터 조정

### 낮은 임계값 (더 많은 상호작용 발견)
```bash
Rscript /home/user3/data_user3/git_repo/_wt/cci/scripts/cci/mnn/run_multinichenet.R \
  -i /data/user3/sobj/is2_IS_3_1_plots.qs \
  -c anno3 \
  -f X2-X1 \
  --logfc_threshold 0.05 \
  --min_cells 5 \
  -o multinichenet_relaxed
```

### 특정 세포 타입만 분석
```bash
Rscript /home/user3/data_user3/git_repo/_wt/cci/scripts/cci/mnn/run_multinichenet.R \
  --senders "Mono_CD14_Classical,Mono_CD14_Inflammatory" \
  --receivers "CD4_Naive,CD8_Naive" \
  -o multinichenet_mono_to_tcells
```

### 병렬 처리 (코어 수 조정)
```bash
Rscript /home/user3/data_user3/git_repo/_wt/cci/scripts/cci/mnn/run_multinichenet.R \
  --cores 16 \
  -o multinichenet_fast
```

---

## 결과 확인

### Prioritization table 확인
```bash
Rscript scripts/cci/mnn/inspect_mnn_results.R
```

### Circos plot 생성 (Control vs Treatment) ⭐
*Best Practice Style: Perfect Inner Hint Track, Stable Interactions, Global Title*
```bash
Rscript scripts/cci/mnn/plot_comparison_circos.R
```

### Top interactions 요약
```bash 
Rscript scripts/cci/mnn/summarize_interactions.R
```

### 결과 검사 (Inspect Results)
```bash
Rscript scripts/cci/mnn/inspect_mnn_results.R
```

---

## Analysis Methodology
MultiNicheNet uses **Pseudobulk-based Differential Expression** via the `muscat` package.
- **Method**: `edgeR` or `limma-voom` on aggregated counts (sum of counts per sample).
- **Why**: This correctly accounts for biological variability across samples, avoiding false positives common in simple cell-level tests.
- **Covariates**: Can be adjusted using the `covariates` argument (e.g., `~ group + batch`).
- **Custom DEG**: You can supply your own DE results if needed.

---

## R Console에서 실행

### 기본 분석
```R
source('myR/R/cci_multinichenet_wrapper.R')
library(qs)

sobj <- qs::qread('/data/user3/sobj/is2_IS_3_1_plots.qs')

results <- run_multinichenet_analysis(
  sobj = sobj,
  sample_id = 'hos_no',
  group_id = 'g3',
  celltype_id = 'anno3',
  contrasts_oi = 'X2-X1',
  output_dir = 'my_analysis',
  verbose = TRUE
)
```

### 낮은 임계값으로 분석
```R
results <- run_multinichenet_analysis(
  sobj = sobj,
  sample_id = 'hos_no',
  group_id = 'g3',
  celltype_id = 'anno3',
  contrasts_oi = 'X2-X1',
  min_cells = 5,           # 10 → 5
  output_dir = 'analysis_relaxed'
)
```

### Circos plot 생성
```R
source('myR/R/cci_multinichenet_wrapper.R')

results <- readRDS('my_analysis/multinichenet_results.rds')
# 또는 (더 빠름)
results <- qs::qread('my_analysis/multinichenet_results.qs')

plot_multinichenet_circos(
  mnn_results = results,
  group_oi = 'X2',
  top_n = 50,
  output_file = 'circos_plot.pdf'
)
```

---

## 문제 해결

### Prioritization table이 비어있는 경우
```bash
# 임계값을 낮춰서 재실행
Rscript myR/scripts/run_multinichenet.R \
  --logfc_threshold 0.05 \
  --fraction_cutoff 0.03 \
  -o multinichenet_very_relaxed
```

### 특정 세포 타입 실패 시
```bash
# failed cell types 확인
Rscript myR/scripts/check_failed_celltypes.R

# min_cells를 낮춰서 재시도
Rscript myR/scripts/run_multinichenet.R \
  --min_cells 3 \
  -o multinichenet_lowcells
```

### 메모리 부족 시
```bash
# 세포 타입 제한
Rscript myR/scripts/run_multinichenet.R \
  --senders "Mono_CD14_Classical,NK_CD56dim" \
  --receivers "CD4_Naive" \
  --cores 4 \
  -o multinichenet_limited
```

---

## 고급 사용

### 배치 효과 포함
```R
results <- run_multinichenet_analysis(
  sobj = sobj,
  sample_id = 'hos_no',
  group_id = 'g3',
  celltype_id = 'anno3',
  batches = 'batch',      # 배치 컬럼
  contrasts_oi = 'X2-X1',
  output_dir = 'with_batch'
)
```

### Covariates 포함
```R
results <- run_multinichenet_analysis(
  sobj = sobj,
  sample_id = 'hos_no',
  group_id = 'g3',
  celltype_id = 'anno3',
  covariates = c('age', 'sex'),
  contrasts_oi = 'X2-X1',
  output_dir = 'with_covariates'
)
```

### 다중 대비 (Multiple contrasts)
```R
results <- run_multinichenet_analysis(
  sobj = sobj,
  sample_id = 'hos_no',
  group_id = 'g3',
  celltype_id = 'anno3',
  contrasts_oi = c('X2-X1', 'X3-X1', 'X3-X2'),
  output_dir = 'multiple_contrasts'
)
```

---

## 결과 파일 구조

```
output_dir/
├── multinichenet_results.qs        # 전체 결과 (qs 형식, 빠름)
├── multinichenet_results.rds       # 전체 결과 (rds 형식, fallback)
├── prioritization_table.csv        # Prioritization (CSV)
├── celltype_de.csv                 # DE results (CSV)
├── top_ligand_activities.csv       # Top ligands (CSV)
└── plots/
    ├── circos_X2.pdf               # Circos plot
    └── ...
```

---

## 파라미터 요약

| 파라미터 | 기본값 | 설명 | 영향 범위 |
|---------|--------|------|----------|
| `min_cells` | 10 | 최소 세포 수 | 모든 단계 |
| `logfc_threshold` | 0.10 | LogFC 임계값 | DE filtering, prioritization |
| `fraction_cutoff` | 0.05 | 발현 비율 임계값 | Expression filtering |
| `p_val_threshold` | 0.05 | P-value 임계값 | DE filtering |
| `p_val_adj` | FALSE | Adjusted p-value 사용 | DE filtering |
| `cores` | 8 | 병렬 처리 코어 수 | 실행 속도 |

---

## 팁

### 빠른 테스트
```bash
# 소수의 세포 타입으로 빠른 테스트
Rscript myR/scripts/run_multinichenet.R \
  --senders "Mono_CD14_Classical" \
  --receivers "CD4_Naive" \
  --cores 4 \
  -o quick_test
```

### 재현 가능한 분석
```R
# Seed 설정
set.seed(42)

results <- run_multinichenet_analysis(
  sobj = sobj,
  # ... parameters ...
  output_dir = 'reproducible_analysis'
)

# Git commit 기록
system('git rev-parse HEAD > reproducible_analysis/git_commit.txt')
```

### 결과 비교
```bash
# 두 가지 설정으로 분석
Rscript myR/scripts/run_multinichenet.R --logfc_threshold 0.10 -o strict
Rscript myR/scripts/run_multinichenet.R --logfc_threshold 0.05 -o relaxed

# R에서 비교
# results_strict <- readRDS('strict/multinichenet_results.rds')
# results_relaxed <- readRDS('relaxed/multinichenet_results.rds')
# compare results...
```

---

## 참고

- 개발 로그: `docs/cci/MULTINICHENET_DEV.md`
- 개발 로그 (한글): `docs/cci/MULTINICHENET_DEV_KR.md`
- README: `myR/scripts/README_MULTINICHENET.md`
