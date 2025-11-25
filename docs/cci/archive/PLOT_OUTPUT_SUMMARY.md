# CCI 분석 플롯 출력 요약

## 생성된 플롯 파일

모든 플롯이 성공적으로 생성되었습니다.

### 출력 디렉토리
`/data/user3/sobj/cci_plots_output/run1/`

### 생성된 파일 목록

1. **NicheNet_Circos_LR.pdf** - Circos plot (PDF, 고해상도)
2. **NicheNet_Circos_LR.png** - Circos plot (PNG, 1000x1000)
3. **NicheNet_Circos_LR_legend.txt** - Circos plot 범례 정보
4. **NicheNet_Ligand_Target_Heatmap.png** - Ligand-Target 상호작용 히트맵
5. **NicheNet_Ligand_Receptor_Heatmap.png** - Ligand-Receptor 상호작용 히트맵
6. **NicheNet_Ligand_Activity_Histogram.png** - Ligand activity 히스토그램
7. **NicheNet_Ligand_AUPR_Heatmap.png** - Ligand AUPR 히트맵

## 분석 설정

### Cell Types
- **Receiver**: CD4+ T-cells (T cell)
- **Senders**: 
  - Monocytes / Macrophages
  - Activated Monocytes / Macrophages
  - IFN-Stimulated Monocytes
  - NK Cells
  - Effector T-cells / NKT

### 조건
- **Condition OI**: 2
- **Condition Ref**: 1
- **Group by**: g3

## Top Ligands

분석 결과 상위 리간드:
1. MIF
2. OSM
3. CXCL2
4. SPON2
5. FURIN
6. ADAM10
7. S100A9
8. TNFSF12
9. CTSD
10. NAMPT

## Circos Plot 설명

Circos plot은 ligand-receptor 상호작용을 시각화합니다:
- **Ligands** (왼쪽): Sender cell type별로 색상 구분
- **Receptors** (오른쪽): Receiver cell type의 수용체
- **화살표**: Ligand → Receptor 방향
- **두께**: 상호작용 강도 (weight)

### Sender Cell Type 색상
- 각 sender cell type은 서로 다른 색상으로 표시됩니다
- Receptor는 통일된 색상으로 표시됩니다
- 범례는 `NicheNet_Circos_LR_legend.txt` 파일에 저장되어 있습니다

## 사용 방법

### 플롯 확인
```bash
# 모든 플롯 파일 확인
ls -lh /data/user3/sobj/cci_plots_output/run1/

# Circos plot 확인
# PDF: 고해상도 인쇄용
# PNG: 웹/프레젠테이션용
```

### R에서 결과 로드
```r
library(qs)
results <- qs::qread("/data/user3/sobj/cci_analysis_results_20251114_172030.qs")

# Top ligands 확인
head(results$nichenet_results$ligand_activities, 10)

# Circos plot (recorded plot object)
results$nichenet_results$plot_circos
```

## 참고

- 모든 플롯은 `run_nichenet_analysis` 함수에 의해 자동 생성됩니다
- `output_dir` 파라미터를 지정하면 플롯이 자동으로 저장됩니다
- `run_circos = TRUE`로 설정하면 Circos plot이 생성됩니다

