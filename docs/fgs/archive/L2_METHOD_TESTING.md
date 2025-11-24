# L2 Method 개별 테스트 가이드

TML7의 L2 methods를 개별적으로 테스트하고 디버깅하는 방법입니다.

**✅ 테스트 완료**: 스크립트는 다양한 데이터셋과 target variable로 테스트되었으며 정상 작동합니다.

## 개별 Method 테스트

### 방법 1: R 스크립트 직접 실행 (권장)

**기본 사용법** (data_seurat_251104.qs, target_var="response"):
```bash
# taskset을 사용하여 CPU 코어 제한 (xgboost 문제 해결)
taskset -c 0-7 Rscript /home/user3/data_user3/git_repo/_wt/fgs/scripts/fgs/test_l2_method_individual.R glm

# xgbTree 테스트 (taskset 필수)
taskset -c 0-7 Rscript /home/user3/data_user3/git_repo/_wt/fgs/scripts/fgs/test_l2_method_individual.R xgbTree
```

**is5 데이터 사용** (target_var="g3"):
```bash
# g3를 target으로 사용
taskset -c 0-7 Rscript /home/user3/data_user3/git_repo/_wt/fgs/scripts/fgs/test_l2_method_individual.R glm g3

# 다른 데이터 파일과 cv_group_var 지정
taskset -c 0-7 Rscript /home/user3/data_user3/git_repo/_wt/fgs/scripts/fgs/test_l2_method_individual.R glm g3 /data/user3/sobj/is5.qs emrid
```

### 방법 2: 환경 변수 사용

```bash
export METHOD=glm
taskset -c 0-7 Rscript /home/user3/data_user3/git_repo/_wt/fgs/scripts/fgs/test_l2_method_individual.R
```

### 방법 3: 스크립트 내에서 수정

`test_l2_method_individual.R` 파일을 열어서 다음 줄의 주석을 해제하고 method 이름을 지정:

```r
METHOD <- "glm"  # 이 줄의 주석을 해제하고 method 이름 지정
```

## 모든 Methods 일괄 테스트

```bash
# 기본: 8코어 사용 (0-7)
cd /home/user3/data_user3/git_repo/_wt/fgs
./scripts/fgs/test_all_l2_methods.sh

# 다른 코어 수 지정 (예: 4코어)
./scripts/fgs/test_all_l2_methods.sh 0-3
```

## 지원되는 Methods

- **기본 methods**: `glm`, `ranger`, `xgbTree`
- **확장 methods**: `glmnet`, `svmRadial`, `mlp`, `mlpKerasDropout`, `nnet`, `earth`

## 출력 파일

각 method 테스트 결과는 다음 위치에 저장됩니다:

```
/data/user3/sobj/tml2_v7_<method_name>.qs
```

예:
- `/data/user3/sobj/tml2_v7_glm.qs`
- `/data/user3/sobj/tml2_v7_xgbTree.qs`

전체 테스트 요약은 다음 파일에 저장됩니다:

```
/data/user3/sobj/tml2_v7_test_summary.txt
```

## 중요 사항

### xgboost (xgbTree) CPU 제한 필수

**⚠️ xgbTree는 반드시 taskset을 사용해야 합니다.**

```bash
# 올바른 방법 (taskset 사용)
taskset -c 0-7 Rscript test_l2_method_individual.R xgbTree

# 잘못된 방법 (CPU 완전 점유 발생)
Rscript test_l2_method_individual.R xgbTree
```

### 다른 Methods

- `glm`, `ranger`, `glmnet` 등은 taskset 없이도 작동할 수 있지만, 안정성을 위해 taskset 사용 권장
- 특히 큰 데이터셋(is5s 수준)에서는 모든 method에 taskset 사용 권장

## 테스트 데이터

**기본 설정**:
- **FGS 결과**: `/data/user3/sobj/fgs/fgs2.qs`
- **Seurat 데이터**: `/data/user3/sobj/data_seurat_251104.qs` (기본값)
- **target_var**: `response` (기본값, data_seurat_251104.qs용)
- **cv_group_var**: `hos_no` (기본값)
- **k_folds**: 5
- **metric**: AUC

**is5 데이터 사용 시**:
- **Seurat 데이터**: `/data/user3/sobj/is5.qs` (또는 다른 경로)
- **target_var**: `g3` (is5 데이터용)
- **cv_group_var**: `emrid` (또는 다른 그룹 변수)

스크립트는 `target_var`가 데이터에 없으면 자동으로 `response` 또는 `g3`를 찾아서 사용합니다.

## 결과 확인

각 method 테스트 후 다음 정보가 출력됩니다:

- ✓/✗ 성공/실패 여부
- 실행 시간 (초, 분)
- Best model 이름
- CV ROC (best)
- Training AUC
- 저장된 파일 경로

## 디버깅 팁

1. **개별 method로 테스트**: 문제가 있는 method만 개별적으로 테스트
2. **작은 데이터로 먼저 테스트**: 큰 데이터셋 대신 작은 샘플로 먼저 확인
3. **에러 메시지 확인**: 실패한 경우 에러 메시지를 자세히 확인
4. **CPU 사용량 모니터링**: `htop`으로 CPU 사용량 확인 (taskset이 제대로 작동하는지)

## 예시

```bash
# 1. glm 테스트
taskset -c 0-7 Rscript test_l2_method_individual.R glm

# 2. ranger 테스트
taskset -c 0-7 Rscript test_l2_method_individual.R ranger

# 3. xgbTree 테스트 (taskset 필수)
taskset -c 0-7 Rscript test_l2_method_individual.R xgbTree

# 4. 모든 methods 일괄 테스트
./scripts/fgs/test_all_l2_methods.sh
```

