# MultiNicheNet 분석 및 시각화 개발 보고서 (Final)

## 1. 프로젝트 개요
*   **목표**: IS2 vs IS3 데이터셋(Single-cell)을 사용한 세포 간 상호작용(Comparison Analysis) 분석 및 고품질 시각화(Circos Plot) 구현.
*   **핵심 성과**: MultiNicheNet 결과를 효과적으로 비교할 수 있는 **Best Practice Circos Plot (V10)** 개발 완료.

## 2. 파일 구조 (최종)
모든 관련 실행 스크립트는 `scripts/cci/mnn/` 하위로 통합되었습니다. 함수의 정의는 `myR/R`에 유지됩니다.

*   **메인 분석**: `scripts/cci/mnn/run_multinichenet_IS2_IS3.R`
    *   데이터 로드, 전처리, MultiNicheNet 분석 파이프라인.
    *   라이브러리 소스: `myR/R/cci_multinichenet_wrapper.R`
*   **시각화 (최종본)**: `scripts/cci/mnn/plot_comparison_circos.R`
    *   V10 Perfection Logic 적용.
    *   중복 제거(`distinct`), 결정론적 정렬(`arrange`), Inner Hint Track 구현 포함.
*   **결과 검사**: `scripts/cci/mnn/inspect_mnn_results.R`
*   **아카이브**: `scripts/cci/mnn/archive/` (V2~V9 등 개발 중간 버전 보관)

## 3. 개발 여정 및 주요 의사결정

### A. Prioritization Table 공백 문제 (`prioritization_tables` is Empty)
*   **문제**: `generate_prioritization_tables` 함수 실행 후 `group_prioritization_tbl`이 비어있는 현상 지속.
*   **해결**: 표준 테이블을 우회(Bypass)하고, `ligand_activities_target_de_tbl`, `celltype_de`, `lr_network`를 직접 조인하여 상호작용 테이블을 재구성함.

### B. 시각화 고도화 (Circos Plot Evolution)
*   **V1-V4**: 기본 기능 구현. 라벨 겹침 및 레이아웃 불안정 해결 시도.
*   **V5-V7**: **Fixed Layout** 도입. Control(X1)과 Treatment(X2)가 동일한 순서를 가지도록 하여 1:1 비교 가능하게 함.
*   **V8-V9**: **Inner Hint Track** 아이디어 구현.
*   **V10 (Final)**: 중복 제거 및 재현성 확보 (Final Polish).

## 4. 해결된 주요 이슈 (Solved Issues)

### 1. Interaction Count Instability (상호작용 수 변동)
*   **증상**: 코드 실행 시마다 상호작용 수가 미세하게 변함.
*   **원인**: 동점자 정렬 기준 부재(Non-deterministic sort).
*   **해결**: `set.seed(1234)` 및 `arrange(..., ligand, receptor, sender, receiver)`로 정렬 순서 고정.

### 2. Duplicate Links (중복 링크)
*   **증상**: 동일한 Ligand-Receptor Pair에 대해 링크가 겹쳐 그려짐.
*   **원인**: 데이터 조인 과정에서의 중복 발생.
*   **해결**: `distinct(ligand, receptor, sender, receiver)` 함수를 사용하되, **반드시 `arrange(desc(score))`를 먼저 수행**하여 점수가 가장 높은 Target 정보를 선별(Pick Best)한 후 중복을 제거함. 순서가 바뀌면 랜덤한 Target이 선택되어 결과가 달라짐.

### 3. Permissive Mode for Rare Celltypes (pDC 등)
*   **증상**: 기본 설정(Strict) 사용 시, pDC와 같이 세포 수가 적거나 변화가 미세한 세포 타입에서 유의한 상호작용이 거의 검출되지 않음.
*   **해결**: `run_multinichenet_IS2_IS3_permissive.R` 생성 및 Wrapper 수정.
    *   **Threshold 완화**: `p_val_thresh` 0.05 -> 0.20 (Raw), `logFC` 0.10 -> 0.05, `min_cells` 10 -> 5.
    *   이로 인해 잠재적인(suggestive) 상호작용을 더 많이 포착하게 됨.

---

## 5. 최신 업데이트 (2025-12-12)

### A. Circos Plot 유연성 및 균형 확보
1.  **Command Line Arguments (optparse) 도입**:
    *   하드코딩된 경로와 설정을 제거하고 `-n`, `-s` 등 파라미터로 제어 가능하게 변경.
2.  **Greedy Selection for Balancing (`--max_per_sender/receiver`)**:
    *   **문제**: 단순 필터링 시 Top N개가 특정 우점 세포(예: pDC)에 의해 독점되거나, 교차 필터링으로 인해 결과 개수가 급감하는 문제 발생.
    *   **해결**: Greedy 알고리즘 도입. 전체 랭크 리스트를 순회하며 Sender/Receiver별 할당량(Quota)이 남은 경우에만 선택, Top N을 채울 때까지 반복.
    *   이로써 pDC와 같은 dominant source를 억제하면서도 전체 interaction 수를 보장함.

### B. 중요 버그 수정 (Troubleshooting)
1.  **sce Object Not Found**: Wrapper 함수 수정 중 Seurat->SCE 변환 로직 누락 -> 복구 완료.
2.  **Invalid Color Name**: 색상 코드 생성 시 `paste0(..., "CC")` 사용으로 이름 있는 색상(예: "blue") 처리 불가 -> `adjustcolor()`로 교체.
3.  **Duplicate Factor Levels**: Circos 초기화 시 Sector 순서 중복 -> `unique()` 처리.
4.  **Object 'direction' Not Found**: Greedy Selection 로직 변경 시 `direction` 컬럼 생성 누락 -> 로직 후단으로 이동하여 해결.

---

## 6. 개발 회고 및 Best Practice 가이드 (Lessons Learned)
이번 프로젝트에서 시각화 품질을 높이기 위해 수많은 시행착오를 겪었습니다. 다음 개발자를 위해 핵심 노하우를 정리합니다.

### A. Link 가시성: Gradient 대신 "Inner Hint Track"
*   **실패 사례**: 링크 선 자체에 그라데이션을 넣어 방향성을 표현하려 했으나, 렌더링 시 색이 깨지거나(White-out), 시각적으로 너무 복잡해 알아보기 힘들었습니다.
*   **성공 사례 (V10)**: **"Inner Hint Track"** 도입.
    *   링크 시작점(Sender Block) 바로 안쪽에 **도착지(Receiver)의 색상**을 가진 얇은 띠를 두릅니다.
    *   사용자는 선을 따라갈 필요 없이, 출발점의 색상 띠만 보고도 "아, 이건 B세포로 가는 신호구나"라고 직관적으로 알 수 있습니다.
    *   구현 팁: `circos.link` 외에 별도의 `circos.track`을 생성하고, `panel.fun` 내부에서 링크 좌표와 동일한 위치에 `rect`를 그립니다.

### B. 레이아웃: "Clustered Grouping"으로 정돈하기
*   **문제**: 세포 타입(Sender/Receiver)과 상관없이 유전자들이 나열되면 매우 산만해 보입니다.
*   **해결**: 3단계 Gap 전략을 사용합니다.
    1.  **Tiny Gap (0.5)**: 같은 세포 타입 내 유전자 사이. (거의 붙어 보임)
    2.  **Medium Gap (5)**: 서로 다른 세포 타입 사이.
    3.  **Large Gap (30)**: Sender 그룹과 Receiver 그룹 사이.
    *   이렇게 하면 시각적으로 '덩어리감'이 생겨 정보 파악이 쉬워집니다.

### C. 데이터 정확성: 재현성과 중복 제거 보장
*   **Deterministic Sorting**: `top_n`은 동점자가 있을 때 무작위로 자를 수 있습니다. 이를 방지하기 위해 정렬 기준에 식별자(ID)를 포함시켜 순서를 고정해야 합니다.
*   **Deduplication**: 시각화 전 데이터 정제(Cleaning) 단계에서 `distinct()`를 사용하여 중복된 Link가 겹쳐 그려지는 것을 방지해야 합니다. 이는 시각적 품질뿐 아니라 분석 결과의 신뢰도와도 직결됩니다.

이 가이드가 향후 Agent나 개발자가 CCI 시각화 도구를 개발할 때 유용한 참고자료가 되기를 바랍니다.

**최종 수정일**: 2025-12-12
**작성자**: Antigravity (AI Assistant)
