
# Frequency Analysis

## MASC

### Run (CLI)

```bash
# anno3 + g3 (complex covariates), random effect: hos_no
Rscript /home/user3/data_user3/git_repo/mylit/Git_Repo/_wt/masc/scripts/masc/run_masc.R \
  -i /data/user3/sobj/is2_IS_3_clustered.qs \
  -o /data/user3/sobj/masc/stroke_complex_cli \
  --cluster_var anno3 \
  --contrast_var g3 \
  --random_effects hos_no \
  --fixed_effects GEM,SET,age,sex,bmi,hx_smok,hx_alcohol \
  --prefix masc_anno3_complex
```

### Notes
- Plot 저장은 PNG/PDF로 생성됨(출력 디렉토리 확인).
- `hos_no`가 numeric으로 들어오면 내부에서 character로 변환됨.
- 너무 느리면 `--max_cells 20000` 같은 옵션으로 global downsample 가능.

## MILO


# DEG analysis

## Single Tools

## Consensus DEG

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
Rscript scripts/cellchat/run_cellchat_cli.R \
  -i /data/user3/sobj/is2_IS_3_1_plots.qs \
  -g anno3 \
  -s hos_no \
  -a g3 \
  --subset_aggregate "2,1" \
  -d "Secreted Signaling" \
  --thresh_p_expression 0.1 \
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
```




# Trajectory Analysis

```bash

```


```bash

```


```bash

```


```bash

```

