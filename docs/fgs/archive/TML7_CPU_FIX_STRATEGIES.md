# TML7 CPU 점유 문제 해결 전략

## 문제 분석

TML7에서 `caret::train`이 여전히 전체 CPU를 점유하는 원인:

1. **caret 내부 병렬 처리**: `allowParallel=FALSE`로 설정했지만, caret이 내부적으로 다른 병렬 백엔드를 사용할 수 있음
2. **모델 패키지 내부 병렬**: xgboost, ranger 등이 C/C++ 레벨에서 자체 병렬 처리
3. **BLAS/LAPACK 멀티스레딩**: 환경 변수가 제대로 전달되지 않거나 런타임에 재설정됨
4. **doParallel/foreach 백엔드**: caret이 등록된 백엔드를 사용할 수 있음

## 해결 전략

### 전략 1: caret 호출 전 강제 백엔드 정리 (권장)

```r
# caret::train 호출 직전에 모든 병렬 백엔드 강제 정리
if (requireNamespace("doParallel", quietly = TRUE)) {
  tryCatch({
    clusters <- getDoParWorkers()
    if (clusters > 0) {
      stopCluster(getDoParWorkers())
    }
  }, error = function(e) {})
  registerDoSEQ()
}
if (requireNamespace("foreach", quietly = TRUE)) {
  registerDoSEQ()
}
```

**장점**: 간단하고 즉시 적용 가능  
**단점**: caret이 내부적으로 다시 등록할 수 있음

### 전략 2: 모델별 스레드 제한 강화

각 모델 패키지의 스레드 수를 명시적으로 제한:

```r
# xgboost
options(xgboost.nthread = 1)

# ranger  
options(ranger.num.threads = 1)

# glmnet (OpenMP)
Sys.setenv(OMP_NUM_THREADS = "1")

# nnet
options(nnet.nthreads = 1)
```

**장점**: 모델 레벨에서 제어 가능  
**단점**: 모든 모델에 대해 설정 필요

### 전략 3: caret 래퍼 함수로 격리

caret::train을 래핑하여 환경을 완전히 격리:

```r
train_with_limits <- function(...) {
  # 환경 변수 재설정
  old_env <- Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS"))
  on.exit(do.call(Sys.setenv, old_env))
  
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1"
  )
  
  # 병렬 백엔드 정리
  if (exists("registerDoSEQ", envir = asNamespace("foreach"))) {
    foreach::registerDoSEQ()
  }
  
  # caret 호출
  caret::train(...)
}
```

**장점**: 완전한 격리  
**단점**: 코드 복잡도 증가

### 전략 4: tuneLength 감소

CV fold 수는 유지하되, 하이퍼파라미터 튜닝 범위를 줄임:

```r
# tuneLength를 5에서 3으로 감소
caret::train(..., tuneLength = 3)
```

**장점**: 실행 시간 단축으로 CPU 사용 시간 감소  
**단점**: 모델 성능에 영향 가능

### 전략 5: 모델별 순차 실행 + 강제 GC

각 모델 사이에 가비지 컬렉션과 백엔드 정리:

```r
for (m in l2_methods) {
  # 이전 모델 정리
  gc()
  if (requireNamespace("doParallel", quietly = TRUE)) {
    tryCatch(doParallel::stopImplicitCluster(), silent = TRUE)
  }
  foreach::registerDoSEQ()
  
  # 모델 학습
  fit <- caret::train(...)
}
```

**장점**: 모델 간 격리  
**단점**: 약간의 오버헤드

### 전략 6: cgroups 또는 nice를 사용한 시스템 레벨 제한

R 스크립트 실행 시 CPU 사용량을 시스템 레벨에서 제한:

```bash
# cgroups 사용 (Linux)
cgcreate -g cpu:/fgs_limit
echo 50000 > /sys/fs/cgroup/cpu/fgs_limit/cpu.cfs_quota_us  # 50% CPU
cgexec -g cpu:/fgs_limit Rscript script.R

# nice 사용
nice -n 19 Rscript script.R  # 낮은 우선순위
```

**장점**: 시스템 레벨 제어  
**단점**: 시스템 권한 필요

### 전략 7: 모델 학습을 별도 R 프로세스로 격리

각 모델을 별도 Rscript 프로세스로 실행:

```r
train_model_isolated <- function(method, data, ...) {
  # 임시 스크립트 생성
  script_file <- tempfile(fileext = ".R")
  writeLines(sprintf("
    Sys.setenv(OMP_NUM_THREADS='1')
    library(caret)
    fit <- train(...)
    saveRDS(fit, '%s')
  ", output_file), script_file)
  
  # 별도 프로세스로 실행
  system2("Rscript", script_file, env = c("OMP_NUM_THREADS=1"))
  
  # 결과 로드
  readRDS(output_file)
}
```

**장점**: 완전한 프로세스 격리  
**단점**: 오버헤드 큼, 구현 복잡

## 구현된 조합 (TML7에 적용됨)

현재 TML7 함수에 다음 전략들이 조합되어 적용되었습니다:

### 적용된 전략들:

1. **전략 1 (백엔드 정리)**: 각 모델 학습 전후에 모든 병렬 백엔드 강제 정리
   - `doParallel::stopImplicitCluster()` 호출
   - `doParallel::registerDoSEQ()` 강제 등록
   - `foreach::registerDoSEQ()` 강제 등록
   - `doMC::registerDoMC(cores = 1)` 설정

2. **전략 2 (환경 변수 재설정)**: 각 모델 학습 전에 BLAS/LAPACK 스레드 제한
   - `OMP_NUM_THREADS=1`
   - `OPENBLAS_NUM_THREADS=1`
   - `MKL_NUM_THREADS=1`
   - `VECLIB_MAXIMUM_THREADS=1`
   - `NUMEXPR_NUM_THREADS=1`

3. **전략 3 (모델별 스레드 제한)**: 모델 패키지별 옵션 설정
   - `xgboost`: `options(xgboost.nthread = 1)` + `nthread = 1` 파라미터
   - `ranger`: `options(ranger.num.threads = 1)` + `num.threads = 1` 파라미터
   - `glmnet`, `nnet`, `earth`: `OMP_NUM_THREADS=1` 재설정

4. **전략 5 (GC 및 정리)**: 각 모델 학습 전후에 가비지 컬렉션 수행
   - 학습 전: `gc(verbose = FALSE)`
   - 학습 후: 백엔드 정리 + `gc(verbose = FALSE)`

### 추가 권장사항:

만약 여전히 CPU 점유 문제가 발생한다면 (특히 모델 패키지가 C/C++ 레벨에서 병렬 처리하는 경우):

1. **시스템 레벨 제한 (전략 6) - 강력 권장**: 
   ```bash
   # 방법 1: nice를 사용한 우선순위 낮추기
   nice -n 19 Rscript script.R
   
   # 방법 2: cgroups를 사용한 CPU 제한 (Linux)
   # 50% CPU로 제한 (예: 8코어 시스템에서 4코어만 사용)
   cgcreate -g cpu:/fgs_limit
   echo 50000 > /sys/fs/cgroup/cpu/fgs_limit/cpu.cfs_quota_us
   cgexec -g cpu:/fgs_limit Rscript script.R
   
   # 방법 3: taskset을 사용한 CPU 코어 제한
   taskset -c 0-7 Rscript script.R  # 코어 0-7만 사용
   ```

2. **tuneLength 감소 (전략 4)**: `tuneLength = 3`으로 변경 (성능 영향 있음)

3. **프로세스 격리 (전략 7)**: 각 모델을 별도 R 프로세스로 실행 (오버헤드 큼)

### 중요: xgboost의 C/C++ 레벨 병렬 처리 (특별 주의)

**⚠️ 확인된 문제**: **xgboost (xgbTree)만** 다른 모델들과 달리 R 레벨 제어만으로는 CPU 점유를 완전히 막을 수 없습니다.

**증상**:
- 다른 모델들(glm, ranger, glmnet 등)은 R 레벨 제어로 충분
- **xgboost만** 즉시적, 동시적으로 모든 CPU를 점유
- `nthread=1` 파라미터, `xgb.set.config(nthread=1)`, 환경 변수 모두 설정해도 부족

**해결책**:
1. **전역 설정**: `xgb.set.config(nthread=1)` 전역 호출 (함수 시작 시) - 부분적 효과
2. **환경 변수**: `XGBOOST_NTHREAD=1`, `OMP_NUM_THREADS=1` 설정 - 부분적 효과
3. **시스템 레벨 제한 (필수)**: **taskset 사용이 유일한 확실한 해결책**
   ```bash
   taskset -c 0-7 Rscript script.R
   ```

**즉시적, 동시적 CPU 점유 현상**:
- 이는 xgboost가 C/C++ 레벨에서 병렬 처리를 시작했다는 의미입니다
- R 레벨 제어만으로는 해결이 불가능합니다
- **taskset을 사용한 시스템 레벨 제한이 필수입니다**

**권장사항**:
- xgboost를 포함한 TML7 실행 시 **반드시 taskset 사용**
- 다른 모델들만 사용하는 경우는 R 레벨 제어로 충분할 수 있음

