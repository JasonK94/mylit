# 개발 요약 (항목 1-10) - 2025-12-16

이 문서는 DEG Consensus Analysis 파이프라인과 관련된 1번부터 10번까지의 요청 사항에 대한 변경 및 구현 내용을 요약합니다.

## 요청 및 구현 요약

1.  **상관관계 계산 메트릭 (Metric for Correlation Calculation)**
    *   **질문**: 상관관계 계산에 사용되는 메트릭은 무엇인가?
    *   **답변**: 각 방법론에서 추출한 **검정 통계량 (Test Statistics)** (예: t-value, F-value, Z-score, Wald statistic)을 기반으로 **스피어만 상관관계 (Spearman Correlation)**를 계산합니다. 이는 p-value나 adjusted p-value에서 발생하는 동률(tie) 문제를 방지합니다.

2.  **DESeq2 Wald 음의 상관관계 (DESeq2 Wald Negative Correlation)**
    *   **이슈**: DESeq2 Wald 결과가 다른 방법론들과 음의 상관관계를 보임.
    *   **원인**: 방법론 구현 시 대조군(contrast) 방향이 실수로 뒤집힘 (`2 vs 1` 대신 `1 vs 2`).
    *   **해결**: `runDESEQ2_Wald_v1` 함수의 contrast 인자를 `c("group", contrast_groups[1], contrast_groups[2])`로 수정하여 다른 방법론과 일치시킴.

3.  **PCA에서 Muscat-Limma-Trend의 이질성 (Muscat-Limma-Trend Heterogeneity in PCA)**
    *   **관찰**: `muscat-limma-trend`가 PCA 산점도에서 아웃라이어처럼 나타남.
    *   **배경**: 이는 `voom`이나 `edgeR` 방법론과 비교하여 `limma-trend`가 pseudobulked 데이터의 분산을 모델링하는 방식의 차이로 인한 것으로 보입니다. 별도의 코드 수정 요청은 없었으나, 이 관찰은 QC를 위한 Method PCA 플롯의 유용성을 입증합니다.

4.  **Meta P-value 아웃라이어 처리 및 시각화 (Meta P-value Outlier Handling & Visuals)**
    *   **이슈**: 극도로 작은 p-value로 인해 `-log10` 값이 무한대거나 너무 커져 플롯이 왜곡됨.
    *   **해결**:
        *   **Capping**: 유효한 `-log10(meta_p)` 값 분포의 **95분위수(95th percentile)**로 값을 제한(capping)하는 로직 구현.
        *   **시각화**: 산점도와 화산(Volcano) 플롯에서 이렇게 제한된 값들을 별도의 색상(예: 주황/빨강)으로 표시하여 구분함.

5.  **산점도 범례 (Scatter Plot Legend)**
    *   **이슈**: meta-p 산점도에 색상이 의미하는 바를 설명하는 범례가 없었음.
    *   **해결**: `plot_consensus_meta_distributions`에 명시적인 범례를 추가하여 색상 코드(예: "n_significant" 개수 vs "Capped Outlier")를 설명함.

6.  **UMAP 개선 (UMAP Improvements)**
    *   **요청**: UMAP에 대해 알아보기 쉬운 색상 구분 및 범례 추가.
    *   **해결**: `plot_gene_umap`을 업데이트하여 `color` (Consensus Score)와 `shape` (n_significant)에 대한 범례가 명확히 표시되도록 함.

7.  **히트맵 유전자 선택 (heatmap Gene Selection - Beta Matrix)**
    *   **이슈**: 히트맵에서 `muscat-limma-trend` 메서드만 색이 칠해져 보임.
    *   **원인**: 이전의 유전자 선택 방식이 특정 방법론에 치우치거나 겹치는 유전자가 적었을 가능성.
    *   **해결**: `plot_consensus_heatmap`을 업데이트하여 *각* 방법론에서 **|logFC| 기준 상위 10개 유전자**를 선택한 후 **합집합(Union)**을 취하도록 함. 이를 통해 모든 방법론의 유의한 유전자가 포함되도록 개선. 색상은 `logFC`로 매핑.

8.  **화산 플롯 개선 (Volcano Plot Enhancements)**
    *   **요청**: 아웃라이어 capping(5%), 아웃라이어 별도 색상 표시, 기준선 추가.
    *   **해결**:
        *   95% capping 로직 적용 (산점도와 동일).
        *   x = ±1 (logFC) 위치에 `geom_vline` 추가.
        *   y = -log10(0.05) 위치에 `geom_hline` 추가.
        *   아웃라이어는 별도의 색상으로 표시.

9.  **재현성 캡션 (Reproducibility Caption)**
    *   **요청**: 모든 플롯 하단에 분석 메타데이터 텍스트 추가.
    *   **해결**: 모든 플로팅 함수에 `caption` 인자를 추가함. 래퍼 스크립트는 다음과 같은 정보를 포함하는 캡션 문자열을 생성하여 전달:
        *   **Cluster**: (예: "Broad_anno3big")
        *   **Formula**: (예: "~ g3 + sex + age...")
        *   **Script**: (예: "scripts/consensus/run_AG_run1.R")
        *   **Time**: 실행 시간.
        *   **Commit**: 현재 Git 해시 (`d7e49ce...`).

10. **독립 플로팅 스크립트 (Independent Plotting Script)**
    *   **요청**: 결과 파일로부터 플롯을 재생성할 수 있는 스크립트 작성.
    *   **구현**: `scripts/consensus/plot_consensus.R` 작성.
    *   **사용법**: `Rscript scripts/consensus/plot_consensus.R [result_file.qs]`
    *   **기능**: `.qs` 객체를 로드하여 업데이트된 스타일과 캡션이 적용된 모든 컨센서스 플롯(Meta distribution, Volcano, Heatmap, PCA, UMAP, Similarity Heatmap)을 재생성.

## 수정/생성된 파일들
- `myR/R/deg_consensus/deg_methods_deseq2.R`: Contrast 방향 수정.
- `myR/R/deg_consensus/deg_consensus_pipeline.R`: 플로팅 함수 업데이트 (Capping, 범례, 캡션, 유전자 선택).
- `scripts/consensus/run_AG_run1.R`: 새로운 기능을 사용하도록 메인 실행 스크립트 업데이트.
- `scripts/consensus/plot_consensus.R`: 새로운 플로팅 유틸리티 생성.
- `docs/consensus/DEV_SUMMARY_1to10_251216.md`: 영문 요약 파일.
- `docs/consensus/DEV_SUMMARY_1to10_251216_KR.md`: 국문 요약 파일 (본 파일).
