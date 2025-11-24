# FGS/TML 코드 정리 계획서 (Code Refactoring Plan)

**작성일**: 2025-11-24  
**버전**: 1.0  
**목적**: signature.R 파일 모듈화 및 코드베이스 정리

---

## 📊 현재 상태 분석

### 문제점
1. **signature.R 비대화**
   - 현재: 3,567줄, 149개 함수/아이템
   - 너무 길어서 유지보수 및 협업 어려움
   - FGS, TML7, 유틸리티 함수들이 모두 섞여 있음

2. **중복 파일 존재**
   - `signature.R.bak`: 백업 파일 (132,228 bytes)
   - `signature_dev.R`: 개발/실험 함수 (325줄)
   - 역할이 불명확하고 관리되지 않음

3. **버전 관리 문제**
   - FGS 함수들이 여러 버전으로 나뉘어 있을 가능성
   - 어떤 함수가 실제로 사용되는지 불명확

### 목표
- **사용자 인터페이스 함수**만 `signature.R`에 남기기
- **백엔드 구현**을 별도 모듈로 분리
- **개발 함수** 및 **레거시 코드** 정리
- **명확한 디렉토리 구조** 확립

---

## 🎯 정리 작업 계획

### Phase 1: 함수 분류 및 인벤토리 작성

#### 1.1 현재 signature.R 함수 분석
```bash
# 모든 함수 정의 추출
grep -n "^[A-Za-z_][A-Za-z0-9_.]* <- function" myR/R/signature.R > /tmp/signature_functions.txt

# 분류:
# - 핵심 사용자 함수 (FGS, TML7, etc.)
# - FGS 내부 함수 (find_gene_signature_*, etc.)
# - TML 내부 함수 (train_tml_*, etc.)
# - 유틸리티 함수 (AddModuleScore, plot_*, etc.)
# - 레거시/사용하지 않는 함수
```

**예상 분류**:
```
📦 핵심 API 함수 (사용자가 직접 호출)
├── FGS()                          # Gene signature 탐지
├── TML7()                         # Transfer meta-learning
├── compute_meta_gene_importance() # 결과 해석
└── AddMultipleModuleScores()      # Seurat 통합

📦 백엔드 구현 함수 (내부 로직)
├── find_gene_signature_*()        # FGS 메서드 구현
├── train_tml_model_*()            # TML 학습 로직
├── create_cv_folds_*()            # CV 생성 로직
└── validate_*()                   # 검증 함수들

📦 유틸리티 함수
├── plot_*()                       # 시각화
├── analyze_*()                    # 분석 도구
└── helper functions               # 내부 헬퍼

📦 레거시/미사용 함수
├── 이전 버전 FGS 함수들
├── 실험적 함수들
└── 더 이상 호출되지 않는 함수
```

#### 1.2 실제 사용 현황 조사
```bash
# 프로젝트 전체에서 함수 호출 확인
for func in FGS TML7 compute_meta_gene_importance; do
  echo "=== $func ==="
  grep -r "$func(" --include="*.R" --include="*.Rmd" scripts/ docs/ myR/
done
```

---

### Phase 2: 디렉토리 구조 설계

#### 2.1 제안하는 새 구조

```
myR/R/
├── signature.R                    # 사용자 API만 포함 (≤500줄 목표)
│   ├── FGS()
│   ├── TML7()
│   ├── compute_meta_gene_importance()
│   └── AddMultipleModuleScores()
│
├── fgs/                           # FGS 백엔드 구현
│   ├── find_gene_signature.R      # FGS 메인 로직
│   ├── methods_supervised.R       # Lasso, Ridge, RF, etc.
│   ├── methods_unsupervised.R     # NMF, PCA, etc.
│   └── utils_fgs.R                # FGS 전용 유틸리티
│
├── tml/                           # TML 백엔드 구현
│   ├── train_tml.R                # TML7 메인 로직
│   ├── cv_methods.R               # CV 생성 (cv, LOGO, repeatedcv)
│   ├── l2_models.R                # L2 모델 학습
│   └── utils_tml.R                # TML 전용 유틸리티 (기존 파일)
│
├── visualization/                 # 시각화 함수들
│   ├── plot_signatures.R          # FGS 시각화
│   ├── plot_tml_results.R         # TML 시각화
│   └── plot_gene_importance.R     # 유전자 importance 플롯
│
├── legacy/                        # 레거시 코드 보관
│   ├── signature_v1.R             # 이전 버전 FGS
│   ├── signature_v2.R             # 중간 버전
│   └── deprecated_functions.R     # 더 이상 사용 안 함
│
└── dev/                           # 개발/실험 함수
    └── experimental.R             # signature_dev.R 내용 이동
```

#### 2.2 NAMESPACE 및 패키지 로딩 전략

**문제**: 
- 현재 모든 함수가 한 파일에 있어 자동으로 export됨
- 모듈화하면 internal 함수와 exported 함수 구분 필요

**해결책**:
```R
# myR/R/signature.R
#' @export
FGS <- function(...) {
  # source 대신 패키지 내부 함수 호출
  .find_gene_signature_core(...)
}

# myR/R/fgs/find_gene_signature.R
# (No @export - internal function)
.find_gene_signature_core <- function(...) {
  # 실제 구현
}
```

또는 기존 방식 유지 (source 방식):
```R
# myR/R/signature.R
source(file.path(.pkgroot, "R", "fgs", "find_gene_signature.R"))
source(file.path(.pkgroot, "R", "tml", "train_tml.R"))

#' @export
FGS <- function(...) {
  find_gene_signature_core(...)
}
```

---

### Phase 3: 단계적 마이그레이션

#### Step 1: 백업 및 테스트 환경 구축 (완료 예정: Day 1)

1. **현재 signature.R 백업**
   ```bash
   cp myR/R/signature.R myR/R/legacy/signature_complete_v5.5.R
   git add myR/R/legacy/signature_complete_v5.5.R
   git commit -m "Backup: signature.R before refactoring"
   ```

2. **테스트 스크립트 준비**
   - `test_real_data.R` 실행 결과 저장
   - 리팩토링 후 동일한 테스트 재실행하여 결과 비교

#### Step 2: FGS 모듈 분리 (완료 예정: Day 2)

1. **디렉토리 생성**
   ```bash
   mkdir -p myR/R/fgs
   ```

2. **FGS 관련 함수 이동**
   - `find_gene_signature_*()` 함수들 → `myR/R/fgs/find_gene_signature.R`
   - Method 구현 함수들 분류:
     - Supervised: `myR/R/fgs/methods_supervised.R`
     - Unsupervised: `myR/R/fgs/methods_unsupervised.R`

3. **signature.R 수정**
   - FGS() wrapper 함수만 남기고, 실제 구현은 source로 불러오기
   ```R
   # Source FGS implementation
   source(file.path(.pkgroot, "R", "fgs", "find_gene_signature.R"))
   source(file.path(.pkgroot, "R", "fgs", "methods_supervised.R"))
   source(file.path(.pkgroot, "R", "fgs", "methods_unsupervised.R"))
   ```

4. **테스트**
   ```bash
   Rscript -e "source('myR/R/signature.R'); FGS(sobj, 'g3', methods='nmf_loadings')"
   ```

#### Step 3: TML 모듈 분리 (완료 예정: Day 3)

1. **디렉토리 생성**
   ```bash
   mkdir -p myR/R/tml
   ```

2. **TML 관련 함수 이동**
   - `TML7()` 내부 로직 → `myR/R/tml/train_tml.R`
   - CV 생성 함수들 → `myR/R/tml/cv_methods.R`
   - L2 모델 학습 → `myR/R/tml/l2_models.R`

3. **기존 utils_tml.R 통합**
   - 현재 `myR/R/utils_tml.R`을 `myR/R/tml/` 아래로 이동

4. **테스트**
   ```bash
   Rscript scripts/fgs/test_nmf_l2.R
   ```

#### Step 4: 유틸리티 및 시각화 함수 정리 (완료 예정: Day 4)

1. **시각화 함수 분리**
   ```bash
   mkdir -p myR/R/visualization
   ```
   - `plot_*()` 함수들 → `myR/R/visualization/`

2. **유틸리티 함수 정리**
   - `AddMultipleModuleScores()` 등 Seurat 관련 유틸리티는 `signature.R`에 유지
   - 또는 별도 `myR/R/seurat_utils.R` 생성

#### Step 5: 레거시 파일 정리 (완료 예정: Day 5)

1. **signature_dev.R 처리**
   ```bash
   # 실험적 함수들 확인
   grep "^[A-Za-z_].*<- function" myR/R/signature_dev.R
   
   # compute_meta_gene_importance_v2가 실제로 사용되는지 확인
   grep -r "compute_meta_gene_importance_v2" scripts/ docs/
   
   # 사용되지 않으면 dev/로 이동, 사용되면 signature.R에 통합
   ```

2. **signature.R.bak 제거**
   ```bash
   # Git에서 완전히 제거 (히스토리는 남음)
   git rm myR/R/signature.R.bak
   git commit -m "Remove: signature.R.bak (redundant backup)"
   ```

3. **legacy/ 디렉토리 정리**
   - 이전 버전 FGS 함수들 식별 및 문서화
   - 실제로 사용되지 않으면 제거 또는 별도 보관

---

### Phase 4: 문서화 및 검증

#### 4.1 함수 문서 업데이트

각 파일에 명확한 헤더 추가:
```R
# ==============================================================================
# FGS Core Implementation
# ==============================================================================
# File: myR/R/fgs/find_gene_signature.R
# Purpose: Core logic for finding gene signatures (FGS framework)
# Dependencies: methods_supervised.R, methods_unsupervised.R
# Public API: Called by FGS() in signature.R
# ==============================================================================
```

#### 4.2 API 문서 생성

```markdown
# FGS/TML API Reference

## User Functions (Public API)
- `FGS()`: Find Gene Signatures
- `TML7()`: Transfer Meta Learning
- `compute_meta_gene_importance()`: Gene importance from TML results

## Internal Modules (Do not call directly)
- `myR/R/fgs/`: FGS implementation
- `myR/R/tml/`: TML implementation
- `myR/R/visualization/`: Plot functions
```

#### 4.3 전체 테스트 재실행

```bash
# 1. 실제 데이터 테스트
Rscript scripts/fgs/test_real_data.R

# 2. 더미 데이터 테스트
Rscript scripts/fgs/test_nmf_l2.R

# 3. 결과 비교
Rscript -e "
  old <- readRDS('logs/fgs/test_real_data_results_before_refactor.rds')
  new <- readRDS('logs/fgs/test_real_data_results.rds')
  
  # Compare metrics
  all.equal(old$tmla_full$results[[1]]$models$glm$performance$ROC,
            new$tmla_full$results[[1]]$models$glm$performance$ROC)
"
```

---

## 📋 체크리스트

### Phase 1: 분석 및 준비
- [ ] signature.R 함수 목록 추출
- [ ] 각 함수 사용처 조사
- [ ] 레거시 함수 식별
- [ ] 테스트 결과 사전 저장

### Phase 2: FGS 모듈화
- [ ] `myR/R/fgs/` 디렉토리 생성
- [ ] FGS 함수들 이동
- [ ] signature.R에서 source 연결
- [ ] FGS 테스트 통과 확인

### Phase 3: TML 모듈화
- [ ] `myR/R/tml/` 디렉토리 생성
- [ ] TML 함수들 이동
- [ ] utils_tml.R 통합
- [ ] TML 테스트 통과 확인

### Phase 4: 정리
- [ ] signature_dev.R 처리 (통합 or 이동)
- [ ] signature.R.bak 제거
- [ ] 시각화 함수 분리
- [ ] 문서화 완료

### Phase 5: 검증
- [ ] 전체 테스트 재실행
- [ ] 결과 비교 (리팩토링 전/후)
- [ ] 코드 리뷰
- [ ] Git commit 및 태깅

---

## ⚠️ 위험 요소 및 대응 방안

### 위험 1: 패키지 로딩 순서 문제
- **증상**: Source된 파일에서 함수를 찾지 못함
- **대응**: 
  - 모든 source 문을 signature.R 최상단에 배치
  - 또는 .onLoad() 함수에서 일괄 로딩

### 위험 2: 네임스페이스 충돌
- **증상**: 동일한 함수명이 여러 파일에 존재
- **대응**:
  - Internal 함수는 `.function_name()` 형식 사용
  - 또는 명확한 prefix 사용 (e.g., `fgs_internal_*()`)

### 위험 3: 의존성 체인 끊김
- **증상**: A 함수가 B 함수를 호출하는데 B가 다른 파일로 이동됨
- **대응**:
  - 의존성 그래프 작성
  - 관련 함수들은 같은 파일에 유지

### 위험 4: 테스트 결과 불일치
- **증상**: 리팩토링 후 성능 메트릭이 달라짐
- **대응**:
  - `set.seed()` 확인
  - 함수 로딩 순서 확인
  - 코드 로직 변경 없이 구조만 변경했는지 검증

---

## 🎯 성공 기준

1. **코드 가독성**
   - signature.R이 500줄 이하로 줄어듦
   - 각 모듈이 명확한 책임을 가짐

2. **유지보수성**
   - 함수 위치를 쉽게 찾을 수 있음
   - 새 기능 추가 시 어디에 넣을지 명확함

3. **테스트 통과**
   - 모든 기존 테스트가 동일한 결과 반환
   - 실제 데이터 테스트 통과

4. **문서화**
   - 각 모듈의 역할이 문서화됨
   - API vs Internal 구분이 명확함

---

## 📅 예상 일정

| 단계 | 작업 | 예상 소요 시간 |
|------|------|----------------|
| Phase 1 | 분석 및 인벤토리 | 2-3시간 |
| Phase 2 | FGS 모듈 분리 | 3-4시간 |
| Phase 3 | TML 모듈 분리 | 3-4시간 |
| Phase 4 | 정리 및 레거시 처리 | 2-3시간 |
| Phase 5 | 문서화 및 검증 | 2-3시간 |
| **총계** | | **12-17시간** (2-3일) |

---

## 📝 참고사항

- 이 계획서는 실제 코드 분석 후 조정될 수 있습니다
- 각 Phase 완료 후 Git commit을 권장합니다
- 문제 발생 시 이전 commit으로 롤백 가능하도록 단계를 작게 나눕니다
- 실제 데이터 테스트 결과를 기준으로 검증합니다
