# CCI 분석 사용 예시

## 1. 기본 사용법

### 전체 파이프라인 실행
```r
# 함수 로드
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/prepare_cci_data.R")
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/utils_cci.R")
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/save_cci_results.R")
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/run_cci_analysis.R")
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/CCI.R")

# 데이터 로드
library(qs)
is5a <- qs::qread("/data/user3/sobj/IS6_sex_added_251110.qs")

# DEG 계산
deg_df <- FindMarkers(
  subset(is5a, anno3.scvi == "Platelets / Megakaryocytes"),
  group.by = "g3",
  ident.1 = "2"
) %>%
  marker_filter() %>%
  mutate(cluster = "Platelets / Megakaryocytes")

# CCI 분석 실행
results <- run_cci_analysis(
  sobj = is5a,
  cluster_col = "anno3.scvi",
  deg_df = deg_df,
  receiver_cluster = "Platelets / Megakaryocytes",
  condition_col = "g3",
  condition_oi = "2",
  condition_ref = "1",
  output_dir = "/data/user3/sobj/cci_output",
  verbose = TRUE
)

# 결과 저장
qs::qsave(results, "/data/user3/sobj/cci_results_full.qs")
```

## 2. 저장된 결과에서 Circos Plot 다시 그리기

### `nichenet_results.qs` 파일에서 plot 재생성

저장된 NicheNet 결과 파일에서 circos plot을 다시 그릴 수 있습니다:

```r
# replot_nichenet_circos 함수 로드
source("/home/user3/data_user3/git_repo/_wt/cci/myR/R/cci/utils_cci.R")

# 저장된 결과에서 plot 다시 그리기
replot_nichenet_circos(
  results_file = "/data/user3/sobj/run6/nichenet_results.qs",
  format = "png",
  output_file = "/data/user3/sobj/run6/circos_replot.png",
  circos_cex_text = 1.0,  # 텍스트 크기 조정 (기본값: 0.8)
  circos_show_legend = TRUE,  # 범례 표시
  circos_legend_position = "topright",  # 범례 위치
  circos_legend_size = 0.9,  # 범례 텍스트 크기
  circos_legend_inset = c(-0.15, 0),  # 범례 inset
  width = 10,
  height = 10,
  res = 100,
  verbose = TRUE
)

# 또는 화면에 직접 표시
replot_nichenet_circos(
  results_file = "/data/user3/sobj/run6/nichenet_results.qs",
  format = "screen",
  circos_cex_text = 1.0,
  circos_show_legend = TRUE
)
```

**참고**: `recordedplot` 객체는 파라미터를 변경할 수 없습니다. 다른 파라미터로 plot을 생성하려면 원본 분석을 다시 실행하거나, `run_cci_analysis()`에서 `circos_cex_text`, `circos_show_legend` 등의 파라미터를 조정해야 합니다.

## 3. 중간 체크포인트로부터 재개

### `cci_prepared_data_*.qs`로부터 plot 생성

`run_cci_analysis()`가 Step 5에서 저장한 `cci_prepared_data_*.qs` 파일로부터 NicheNet 분석을 재개할 수 있습니다:

```r
# 중간 데이터 로드 및 재개
results <- resume_cci_from_prepared(
  prepared_data_file = "/data/user3/sobj/cci_prepared_data_20251118_145434.qs",
  sobj = is5a,
  condition_col = "g3",
  condition_oi = "2",
  condition_ref = "1",
  output_dir = "/data/user3/sobj/cci_output",
  verbose = TRUE
  # top_n_targets_per_ligand는 자동 조정됨 (4637 DEGs → 50으로 설정)
)

# Plot 확인
plot(results$nichenet_results$plot_ligand_target_network)
plot(results$nichenet_results$plot_ligand_receptor_network)
plot(results$nichenet_results$plot_circos)
```

### 완료된 NicheNet 결과로부터 plot만 그리기

이미 완료된 분석 결과가 있다면, `.qs` 파일에서 직접 plot을 불러올 수 있습니다:

```r
library(qs)

# 최종 결과 로드
results <- qs::qread("/data/user3/sobj/cci_analysis_results_20251118_103840.qs")

# 또는 NicheNet 결과만 로드 (run_nichenet_analysis가 저장한 파일)
nichenet_results <- qs::qread("/data/user3/sobj/run1/nichenet_results.qs")

# Plot 그리기
plot(nichenet_results$plot_ligand_target_network)
plot(nichenet_results$plot_ligand_receptor_network)
plot(nichenet_results$plot_ligand_activity_hist)
plot(nichenet_results$plot_ligand_aupr_heatmap)
plot(nichenet_results$plot_circos)  # Circos plot
```

## 3. 대용량 데이터셋 최적화

### 자동 조정된 파라미터

`run_cci_analysis()`는 DEG 수에 따라 `top_n_targets_per_ligand`를 자동으로 조정합니다:

- **>3000 DEGs**: `top_n_targets_per_ligand = 50` (빠른 계산)
- **>1000 DEGs**: `top_n_targets_per_ligand = 100` (균형)
- **≤1000 DEGs**: `top_n_targets_per_ligand = 200` (기본값)

### 수동 조정

필요시 명시적으로 지정할 수 있습니다:

```r
results <- run_cci_analysis(
  sobj = is5a,
  cluster_col = "anno3.scvi",
  deg_df = deg_df,
  receiver_cluster = "Platelets / Megakaryocytes",
  condition_col = "g3",
  condition_oi = "2",
  condition_ref = "1",
  top_n_targets_per_ligand = 30,  # 더 빠르게 (정밀도는 낮아짐)
  verbose = TRUE
)
```

## 4. 진행 상황 모니터링

`verbose = TRUE`일 때 다음과 같은 로그가 출력됩니다:

```
Step 2: Extracting receiver DEGs...
  Found 4637 DEGs for receiver cluster after filtering
  Large DEG set detected (4637 genes). Auto-adjusting top_n_targets_per_ligand to 50 for faster computation.
Step 6: Running NicheNet analysis...
  Preparing NicheNet run with 23 sender cluster(s) and 4637 receiver DEGs.
  ...
Inferring active target genes for top ligands...
  Processing ligand 2/20: LFNG
  Processing ligand 4/20: FURIN
  ...
  Ligand-target inference completed (45.2 seconds, 850 links found)
```

## 5. 문제 해결

### Plot이 생성되지 않는 경우

1. **중간 체크포인트 확인**: `cci_prepared_data_*.qs`가 저장되었는지 확인
2. **NicheNet 결과 확인**: `output_dir/run*/nichenet_results.qs` 파일 존재 여부 확인
3. **재개 함수 사용**: `resume_cci_from_prepared()`로 중간부터 재개

### 실행 시간이 너무 긴 경우

1. **DEG 수 줄이기**: `p_val_adj_cutoff` 또는 `logfc_cutoff`를 더 엄격하게 설정
2. **Sender 수 줄이기**: `sender_clusters`를 명시적으로 지정
3. **Target 수 줄이기**: `top_n_targets_per_ligand`를 더 작게 설정 (예: 30)

### 메모리 부족

1. **다운샘플 데이터로 먼저 테스트**: `IS6_sex_added_0.1x_251110.qs` 사용
2. **Sender 수 줄이기**: 전체 클러스터 대신 관심 있는 몇 개만 지정
3. **중간 결과 저장 후 재시작**: `cci_prepared_data_*.qs`로부터 재개

