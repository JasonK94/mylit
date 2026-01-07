
# Frequency Analysis

## MASC

### Run (CLI)

```bash
logs# anno3 + g3 (권장 모델: age, sex, GEM만), random effect: hos_no
Rscript /home/user3/data_user3/git_repo/mylit/Git_Repo/_wt/masc/scripts/masc/run_masc.R \
  -i /data/user3/sobj/is2_IS_3_clustered.qs \
  -o /data/user3/sobj/masc/stroke_complex_cli2 \
  --cluster_var anno3 \
  --contrast_var g3 \
  --random_effects hos_no \
  --fixed_effects GEM,age,sex \
  --prefix masc_anno3_complex
```

### Plot only (from cached results)

분석 결과 파일(`*_results.qs`)이 이미 있는 경우, 플롯만 빠르게 생성:

```bash
Rscript /home/user3/data_user3/git_repo/mylit/Git_Repo/_wt/masc/scripts/masc/plot_masc.R \
  --results /data/user3/sobj/masc/stroke_complex_cli2/masc_anno3_complex_results.qs \
  -o /data/user3/sobj/masc/stroke_complex_cli2 \
  --prefix masc_anno3_complex \
  --cluster_var anno3 \
  --contrast_var g3 \
  --model_formula 'Frequency ~ g3 + age + sex + GEM + (1|hos_no)'
```

### Notes
- **Plot 저장**: PNG/PDF로 생성됨 (출력 디렉토리 확인).
- **변수 타입 자동 인식**: 실행 시 "Variable categories: numeric/categorical" 로그 출력됨.
- **중첩 구조 주의**: `hos_no` ⊂ `GEM` ⊂ `SET`. `GEM`과 `SET`을 동시에 포함하지 말 것. CLI는 자동으로 `SET`을 제거함.
- **bmi 제외 권장**: 결측값이 많아 기본적으로 제외하는 것을 권장.
- **renv**: 기본 `/home/user3/GJC_KDW_250721/renv` 사용. `--renv` 옵션으로 변경 가능.
- `hos_no`가 numeric으로 들어오면 내부에서 character로 변환됨.
- 너무 느리면 `--max_cells 20000` 같은 옵션으로 global downsample 가능.

## MILO


# DEG analysis

## Single Tools

## Consensus DEG

### run analysis
```bash
Rscript scripts/consensus/run_deg_consensus_cli.R \
  --input "/data/user3/sobj/is2_IS_3_clustered.qs" \
  --output "/data/user3/sobj/consensus/AG_run1_CLI" \
  --cluster "anno3big" \
  --group "g3" \
  --contrast "2 - 1" \
  --covariates "sex,age" \
  --methods "all" \
  --cores 8
```

### plotting
```bash
Rscript scripts/consensus/plot_consensus_cli.R \
  --input /data/user3/sobj/consensus/AG_run1/Broad_anno3big/results_consensus.qs \
  --output /data/user3/sobj/consensus/AG_run1/Broad_anno3big/plots_cli \
  --cluster "anno3big" \
  --group "g3" \
  --covariates "sex,age,GEM"
```

# Signature Analysis

## FGS

## TML

## CMGI

## AMSC

# CCI

## CellChat

### Run
```bash
# 1. Permissive 설정으로 데이터 다시 분석 (약 10-20분 소요 예상)
# -p 1: permissive p-value threshold (interaction cutoff)
# --thresh_p_expression 1: p-value threshold (expression cutoff)
# -m 5: min cells
# -c 1: core #
Rscript scripts/cellchat/run_cellchat_cli.R \
  -i /data/user3/sobj/is2_IS_3_1_plots.qs \
  -g anno3 \
  -s hos_no \
  -a g3 \
  --subset_aggregate "2,1" \
  -d "Secreted Signaling" \
  --thresh_p_expression 1 \
  -p 1 \
  -m 5 \
  -c 1 \
  -o /data/user3/sobj/cci/cellchat/permissive_run
```

### Plotting
```bash
Rscript scripts/cellchat/plot_cellchat_comparison.R \
  --file1 /data/user3/sobj/cci/cellchat/sample_wise1/merged/X1/cellchat_merged.qs \
  --name1 "g3_1" \
  --file2 /data/user3/sobj/cci/cellchat/sample_wise1/merged/X2/cellchat_merged.qs \
  --name2 "g3_2" \
  -o /data/user3/sobj/cci/cellchat/proper_full/comparison_plots
```
### Inspect data
```bash

Rscript scripts/cellchat/inspect_cellchat.R \
  -i /path/to/cellchat.qs \
  -o /path/to/interactions.csv
```

## NicheNet
###



## MNN (MultiNicheNet)

### Run
```bash
Rscript /home/user3/data_user3/git_repo/_wt/cci/scripts/cci/mnn/run_multinichenet.R \
  -i /data/user3/sobj/is2_IS_3_1_plots.qs \
  -c anno3 \
  -f X2-X1 \
  --logfc_threshold 0.10 \
  --min_cells 5 \
  -o multinichenet_relaxed\
  --cores 1
```
### Circos plot
```bash
Rscript scripts/cci/mnn/plot_comparison_circos.R   --top_n 100   --sort_by score   --max_per_sender 10   --max_per_receiver 10

Rscript scripts/cci/mnn/plot_comparison_circos.R   --top_n 100   --sort_by score   --max_per_sender 10   --max_per_receiver 10 -f /data/user3/sobj/cci/mnn/multinichenet_IS2_IS3_permissive/multinichenet_results.qs -o /data/user3/sobj/cci/mnn/multinichenet_IS2_IS3_permissive
```




# Trajectory Analysis
## Run analysis
```bash
Rscript scripts/pseudotime-dev/run_trajectory_subsets.R \
  -i /data/user3/sobj/is2_IS_3_1_plots.qs \
  -o /data/user3/sobj/pseudotime/pseudotime_results_run1 \
  -t "TXNIP,DDIT4,S100B,CCL4,HLA-B,UTY,XIST" \
  -s anno3big \
  -r umap \
  -c g3 \
  -f FALSE
```

## Plotting
### TXNIP
```bash
Rscript scripts/pseudotime-dev/plot_trajectory_results.R \
  -i /data/user3/sobj/pseudotime/pseudotime_results_run1/Tc_cds.qs \
  -o /data/user3/sobj/pseudotime/pseudotime_results_Tc2 \
  -t "TXNIP"
```
### flexible
```bash
Rscript scripts/pseudotime-dev/plot_trajectory_gene_list.R \
  -i /data/user3/sobj/pseudotime/pseudotime_results_run1/Tc_cds.qs \
  -o /data/user3/sobj/pseudotime/pseudotime_results_Tc_txnip_v1 \
  --target_list "txnip"

Rscript scripts/pseudotime-dev/plot_trajectory_gene_list.R \
  -i /data/user3/sobj/pseudotime/pseudotime_results_run1/Mono_cds.qs \
  -o /data/user3/sobj/pseudotime/pseudotime_results_Mono_txnip_v1 \
  --target_list "txnip"
```

## Differential Analysis
```bash
Rscript scripts/pseudotime-dev/analyze_trajectory_diff.R \
  -i /data/user3/sobj/pseudotime/pseudotime_results_run1/Tc_cds.qs \
  -o /data/user3/sobj/pseudotime/diff_analysis_Tc2 \
  -c g3 \
  -p 10
```

Rscript scripts/pseudotime-dev/analyze_trajectory_diff.R \
  -i /data/user3/sobj/pseudotime/pseudotime_results_run1/Mono_cds.qs \
  -o /data/user3/sobj/pseudotime/diff_analysis_Mono2 \
  -c g3 \
  -p 10

