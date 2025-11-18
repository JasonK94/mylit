# 개발 로그 (DEVLOG)

`myR` 개발의 주요 결정과 맥락을 시간 순으로 기록합니다.

## 2025-05-09 — 프로젝트 착수 (`1d9a350`)
- **작성자**: kjc17  
- **요약**: demultiplexing·통합·Seurat 워크플로우 초기 세팅.
- **세부 사항**:
  - `demuxalot`, `soupx`, `scVI` 스크립트와 7단계 전처리 파이프라인 반입.
  - DESCRIPTION/NAMESPACE 골격과 함께 `myR/` 디렉터리 구조 마련.
- **다음 단계**:
  - 초기 모듈의 설치·로드 경로 검증.
  - 대용량 아티팩트를 정리해 git 기록을 간결하게 유지.

## 2025-05-16 — R 모듈 대량 반입 (`338a7f1`)
- **작성자**: kjc17  
- **요약**: 핵심 분석 모듈과 검증 자료 버전 관리 시작.
- **세부 사항**:
  - `CCI.R`, `pseudobulk_deg.R`, `pseudotime.R`, `legacy.R` 및 테스트 스크립트 추가.
  - NicheNet 결과물, 플롯, Rmd 보고서를 저장해 재현성 확보.
- **다음 단계**:
  - 모듈을 영역별로 재구성.
  - 용량이 큰 이진 파일을 버전 관리에서 제외.

## 2025-05-22 — 복구와 병합 (`71c090f`, `b70456c`)
- **작성자**: kjc17  
- **요약**: 손상된 작업 복구 및 분기 정리.
- **세부 사항**:
  - 삭제된 파일 복원, `demulti_utils` 재도입, 작업 디렉터리 병합.
  - `.gitignore` 정비와 캐시된 대형 파일 제거로 히스토리 정리.
- **다음 단계**:
  - 유틸리티 스크립트 안정화.
  - 공개 인터페이스 문서화 착수.

## 2025-10-10 — formula 기반 유틸 정비 (`9f04cc9`, `c3afe5a`, `3f44856`)
- **작성자**: kjc_server1  
- **요약**: FLSG/RLMG 등 formula_str 헬퍼 일관성 확보.
- **세부 사항**:
  - FLSG/RLMG가 동일한 패턴으로 동작하도록 정비, `fmsg` 패치.
  - `test.Rmd` 갱신과 함께 새로운 테스트 전략 준비.
- **다음 단계**:
  - formula 패턴을 테스트 하네스 전반에 확산.
  - 협업자를 위한 기대치 변경사항 문서화.

## 2025-10-31 — 의사벌크 및 테스트 개편 (`881589c`, `25b205b`, `e785399`, `39a6159`)
- **작성자**: kjc_server1  
- **요약**: 안정 함수 승격과 분석 QA 기능 확장.
- **세부 사항**:
  - `test_working.R`에 핵심 함수 정리, `pb_deg_cursor` 리팩터링, Seurat `slot` → `layer` 전환.
  - TestLISI, PERMANOVA, PCP, PTMFL 추가 및 `NMF` 의존성 명시.
- **다음 단계**:
  - 의사벌크 리팩터를 모듈 단위로 분리해 유지보수성 강화.
  - 신규 QA 유틸 사용 예시 수집.

## 2025-11-05 — 인프라 정리 (`497eab8`, `80260af`, `80f7941`)
- **작성자**: kjc_server1  
- **요약**: 유틸리티 재배치와 불안정 코드 정리.
- **세부 사항**:
  - 기반 R 헬퍼 이동, MUSCAT/NEBULA 러너 추가, MAST 실험 코드 제거.
  - `.gitignore` 업데이트, `projects/`·`trash/` 산출물 버전 관리 제외.
- **다음 단계**:
  - 실험용 노트북은 git 외부에서 관리.
  - 디렉터리 구조 변경 사항을 문서에 반영.

## 2025-11-09 — 테스트 하네스 통합 (`8670f8c`)
- **작성자**: kjc_server1  
- **요약**: 일일 실험을 정리해 표준 테스트 셋으로 편입.
- **세부 사항**:
  - 삭제 기록을 `.gitignore`와 맞추고 `test_working.R` 실행 가능성 확인.
- **다음 단계**:
  - 테스트 워크플로우에 대한 공식 문서 초안 작성.

## 2025-11-10 — 시그니처 v5.2 & Milo Opus6 (`2d9c50e`, `18f1ad1`, `b260bb1`, `a360837`, `7bf2904`)
- **작성자**: kjc_server1  
- **요약**: 시그니처 스코어링 재구성과 Milo 풍부도 파이프라인 도입.
- **세부 사항**:
  - 시그니처 워크플로우 v5.2와 검증 유틸 추가, 문서 갱신.
  - `milo_opus6`, 메타 러너 병렬 기본값 개선, `opus6` 스냅샷과 테스트 확장.
- **다음 단계**:
  - Milo·시그니처 파이프라인 사용 예시 문서화.
  - 해당 기능군을 반영한 버전 태깅 전략 수립.

## 2025-11-10 — 문서화 정비 (`cf6d13e`, `ac3cb7b`, `18c2810`)
- **작성자**: kjc_server1  
- **요약**: `cinit` 실행 및 문서 위치 재정비.
- **세부 사항**:
  - 컨텍스트, 부트스트랩 안내, 다국어 문서를 생성.
  - `DEVLOG.md`, `CHANGELOG.md`를 루트로 이동해 기록 체계 마련.
- **다음 단계**:
  - 컨텍스트 문서를 기능 변화와 동기화.
  - 주요 커밋마다 DEVLOG/CHANGELOG 업데이트 주기화.

## 2025-11-11 — Milo 파이프라인 함수화 착수 (`milo` 브랜치)
- **작성자**: GPT-5 Codex  
- **요약**: Rmd 기반 Milo 워크플로우를 모듈형 함수 세트로 재구성하기 위한 설계 단계.
- **세부 사항**:
  - `context_Korean.md`에 Milo 분석 시 필수 전처리(`buildNhoodGraph`)와 `nhoods(milo)` 슬롯 명명 이슈, `plotDAbeeswarm` 파라미터 가이드를 교훈으로 추가.
  - 저장 단계별 `force_run` 플래그와 자동 suffix 증분 규칙을 포함하는 런너 구조를 구상.
  - plotting 단계는 독립 함수로 분리하되, 기본 파이프라인에서 실행 여부를 인자로 제어하도록 설계.
- **주의사항**:
  - Milo 파이프라인 개발·테스트는 `st/` 경로에서 `start.R`로 `renv` 환경을 활성화한 뒤 진행해야 한다. 루트에서 R을 실행하면 필수 패키지가 로드되지 않아 오류가 재현된다.
- **다음 단계**:
  - `myR/R/analysis/` 하위에 Milo 파이프라인 R 스크립트를 추가하고, Seurat → Milo 변환·검정·시각화를 함수화.
  - 제공된 `IS5_g3NA_removal_251110.qs`를 사용해 논인터랙티브 테스트를 수행하고, 경량 저장 옵션을 검토.

## 2025-11-12 — Milo 파이프라인 디버깅 & 캐시 정비 (`milo` 브랜치)
- **작성자**: GPT-5 Codex  
- **요약**: `run_milo_pipeline()` 실행 시 빈번하게 발생하던 `nhoods`/`dimnames` 오류와 plotting 실패를 해결하고, 캐시 재활용 로직을 안정화.
- **세부 사항**:
  - 희소 행렬 요약 시 `Matrix::summary()`를 직접 사용해 NA 없이 최소값을 계산하고, `colData(milo)`에 안전하게 값을 기록하도록 수정.
  - 캐시가 존재하면 단계 1~3을 즉시 스킵하고, 필요할 때만 Seurat 객체를 lazy-loading 하는 공급자 패턴을 도입.
  - UMAP 플롯은 `scater::plotReducedDim()`으로 교체해 `miloR::plotUMAP` 의존성을 제거하고, beeswarm은 SpatialFDR이 α보다 클 때 `PValue`로 자동 전환하도록 fallback 로직 추가.
  - `buildNhoodGraph()`가 비어 있는 그래프에 대해 항상 재실행되도록 검사하고, `.qs` 저장을 기본값으로 유지.
- **주의사항**:
  - 캐시된 Milo 객체에는 UMAP 좌표가 없을 수 있으므로, plotting 단계에서 Seurat 객체를 다시 불러와 `reducedDim(milo, "UMAP")`을 보강해야 한다.
- **다음 단계**:
  - Milo `qs` 캐시에 UMAP 좌표까지 포함시켜 plotting 시 Seurat 재로딩을 줄이는 방안을 검토.
  - beeswarm 색상 팔레트와 보고용 그래프 템플릿을 표준화.

## 2025-11-12 — Milo 클러스터별 logFC 편중성 검정 함수 개발 (`milo` 브랜치)
- **작성자**: GPT-5 Codex  
- **요약**: MiloR DA 결과에서 클러스터별 logFC 편중성을 검정하는 `test_cluster_logfc_bias()` 함수 개발 및 문서화 완료.
- **세부 사항**:
  - **함수 개발**:
    - `test_cluster_logfc_bias()`: 클러스터별 logFC가 체계적으로 편향되어 있는지 검정
    - Neighborhood 간 비독립성을 고려한 4가지 검정 방법 구현:
      1. **Block Permutation Test** (권장): Block 내에서만 logFC를 섞어 permutation 수행
      2. **Correlation-adjusted t-test (neff)**: 그래프 인접 행렬의 고유값으로 effective sample size 추정
      3. **Mixed-Effects Model (LMM)**: `logFC ~ 1 + (1 | block_id)` 형태로 block-level random effect 모델링
      4. **Empirical Bayes (ashr)**: p-value를 z-score로 변환 후 shrinkage로 effect size 추정
  - **Block Method 개선**:
    - `block_var` 파라미터로 통합: `patient_var`, `batch_var`, `GEM`, `set` 등 어떤 metadata 컬럼도 사용 가능
    - GEM well suffix (`-1`, `-2`) 제거 로직 제거: 중복 barcode 문제 방지
    - 우선순위: `block_var` 제공 시 `colData(milo)`에서 직접 사용 → cell name에서 추출 (fallback)
  - **문서화**:
    - `myR/docs/milo.md`: MiloR 파이프라인 및 검정 방법 상세 문서 작성
    - 각 검정 방법의 원리, 수식, 장단점, 사용 예시 포함
    - Block permutation test 알고리즘 단계별 설명
    - neff 계산 방법 (고유값 기반 effective sample size)
    - LMM 현재 구현의 제한사항 명시 (batch_var 보정 없음)
    - Empirical Bayes shrinkage 방법 설명
- **주요 발견**:
  - Block permutation test는 block 구조를 완전히 보존하여 가장 신뢰할 수 있는 방법
  - neff는 빠르지만 근사적이므로 permutation보다 덜 보수적
  - LMM은 현재 `batch_var`를 보정하지 않음 (block-level random effect만 모델링)
  - ashr는 검정보다는 effect size 추정에 유용
- **기술적 세부사항**:
  - Block permutation: 각 block 내에서만 logFC를 섞어 block 구조 보존
  - neff: `neff = sum(eigenvalues > 0.1)`, 클러스터별로 비례 배분
  - LMM: `logFC ~ 1 + (1 | block_id)`, 절편의 유의성 검정
  - ashr: `z = qnorm(1 - PValue/2) * sign(logFC)`, `se = |logFC| / max(|z|, 1e-6)`
- **주의사항**:
  - LMM에 batch_var 보정을 추가하려면 모델을 `logFC ~ 1 + batch_var + (1 | block_id)` 형태로 확장 필요
  - Block permutation test는 계산 시간이 오래 걸림 (n_perm에 비례, 기본값 2000)
  - 작은 클러스터에서는 neff_cluster가 너무 작아질 수 있음
- **다음 단계**:
  - LMM에 batch_var 보정 옵션 추가 검토
  - 검정 방법 간 결과 비교 및 검증
  - 실제 데이터셋에서의 성능 평가

## 2025-01-XX — TML6 및 compute_meta_gene_importance 버그 수정 (`fgs` 브랜치)
- **작성자**: Auto (GPT-5 Codex)  
- **요약**: TML6 함수의 zero-variance signature 처리 및 compute_meta_gene_importance의 subscript out of bounds 에러 수정.
- **세부 사항**:
  - **TML6 함수 개선**:
    - `.score_signature()` 함수에서 `scale()` 사용 시 zero-variance 문제 해결
    - 분산이 0인 경우 normalize를 건너뛰고 0으로 설정
    - `scale()` 결과가 NA/Inf인 경우 원래 점수 사용
    - xgboost의 deprecated 경고(`ntree_limit` → `iteration_range`) 억제
    - 모든 signature가 제거될 때 더 자세한 디버깅 정보 제공
  - **compute_meta_gene_importance 함수 수정**:
    - `signature_importance[[sig]]` → `signature_importance[sig]`로 변경 (named vector이므로 `[` 사용)
    - `caret::varImp`의 rownames와 `sig_names` 매칭 개선
    - 존재하지 않는 signature는 경고 후 건너뛰기
    - NA/Inf importance 값 체크 및 처리
  - **새로운 함수 추가**:
    - `add_meta_signature_score()`: TML6와 compute_meta_gene_importance로 만든 gene weights를 사용하여 signature score를 계산하고 Seurat 객체에 AddMetaData로 추가
    - Sparse/dense matrix 모두 지원
    - 가중합 계산: `sum(w * expr) / sum(|w|)`
    - 선택적 z-score normalization
    - signed/absolute contribution 선택 가능
- **기술적 세부사항**:
  - Zero-variance signature: `scale()` 전에 `stats::sd()`로 분산 체크
  - Signature 이름 매칭: `intersect()`로 실제 존재하는 signature만 추출
  - 에러 메시지 개선: 예상/실제 signature 이름 표시
- **주의사항**:
  - 모든 signature가 zero-variance인 경우는 데이터나 signature 자체에 문제가 있을 수 있음
  - `add_meta_signature_score()`는 `compute_meta_gene_importance()`의 결과를 입력으로 받음
- **다음 단계**:
  - 실제 데이터셋에서 테스트 및 검증
  - 사용 예시 문서화

## 2025-01-XX — TML6 환자 단위 CV 및 메타러너 확장 (`fgs` 브랜치)
- **작성자**: Auto (GPT-5 Codex)  
- **요약**: TML6에 `cv_group_var` 파라미터를 도입해 환자 단위 누수를 방지하고, L2 메타러너 후보에 `glmnet`, `svmRadial`, `mlp`, `mlpKerasDropout`를 추가.
- **세부 사항**:
  - **그룹 CV**:
    - `cv_group_var` (기본값 `\"emrid\"`)로 같은 환자의 AOI가 항상 같은 fold에 들어가도록 `caret::trainControl(index/indexOut)` 구성
    - 그룹 수가 `k_folds`보다 적으면 자동으로 기존 셀 단위 CV로 fallback
    - 그룹 컬럼이 없으면 경고 후 셀 단위 CV 유지
  - **L2 메서드 확장**:
    - `glmnet`, `svmRadial`, `mlp`, `mlpKerasDropout` 패키지 의존성 체크 후 사용 가능
    - 패키지가 없으면 경고만 출력하고 해당 모델을 자동으로 제거
    - 기존 기본값(`glm`, `ranger`, `xgbTree`)은 그대로 유지
  - **로그 개선**:
    - 그룹 CV가 활성화될 때 fold 구성 메시지를 출력
    - xgboost는 `xgb.set.config(verbosity = 0)`로 C-level 경고를 최대한 억제
- **테스트**:
  - R 4.3 환경에서 `Matrix` 및 Seurat 의존성 설치가 불가하여 Seurat 기반 통합 테스트는 진행하지 못함
  - 매트릭스 입력 + `glm` 메타러너 케이스로 기본 동작을 추가 점검 (caret 미설치 문제로 부분 실패)
  - 추후 R 4.4+ 또는 사전 구축된 컨테이너에서 full 테스트 필요
- **다음 단계**:
  - Seurat/Matrix 환경이 준비된 머신에서 group CV + 확장 메타러너 조합 테스트
  - `cv_group_var`가 matrix 입력에서도 사용할 수 있도록 벡터 인자를 허용할지 검토
  - 문서/사용 예시 업데이트

## 2025-01-XX — FGS v5.4 및 시그니처 스택 테스트 (`fgs` 브랜치)
- **작성자**: GPT-5.1 Codex  
- **요약**: `find_gene_signature_v5.4()` 추가로 ranger/glmnet/NMF 경로의 안정화 옵션을 기본화하고, IS6 Seurat 객체(`is5s`)를 대상으로 end-to-end 검증을 수행.
- **세부 사항**:
  - `FGS()` 기본 진입점을 `find_gene_signature_v5.4()`로 전환해 `method_impl = "v5.4"`를 항상 사용
  - v5.4에서는 `random_forest_ranger`의 permutation 중요도 + fallback, `ridge`/`elastic_net`의 `glmnet` 행렬 강제 변환과 0/1 타깃 변환, `nmf_loadings`의 양수 이동·랭크 가드가 활성화됨
  - `TML7` + `compute_meta_gene_importance` + `add_meta_signature_score` 조합으로 IS6 데이터(`response`, `hos_no` 보정, 상위 200 genes) 테스트
  - meta-signature score와 TML 예측(logit)의 상관, AUC(`pROC::roc`) 계산 및 `qs::qsave()`로 산출물 저장
- **출력**:
  - `docs/work_context/TML6_IMPROVEMENTS_CONTEXT.md`에 실행 컨텍스트 업데이트
  - `outputs/fgs_is5s/`에 `fgs.each.is5s.qs`, `tml.each.is5s.qs` 저장
- **다음 단계**:
  - `find_gene_signature_v5.4()`를 문서/README에도 반영하고, `test.R`의 legacy 함수명을 `_test` 접미사로 전환
  - ranger/glmnet/NMF 외의 메서드에서도 v5.4 경로가 필요한지 검토

---

### 향후 계획
- 시그니처 v5.2 이후 변경 사항을 `function_analysis*.md`에 반영.
- Milo·시그니처 전환점을 기준으로 세맨틱 버전 태깅 여부 검토.

