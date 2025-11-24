# CPU 설정 가이드

FGS 프로젝트에서 CPU 사용량을 유연하게 제어하는 방법을 설명합니다.

## 환경 변수 설정

`init_fgs_env.R`을 source하기 전에 환경 변수를 설정하여 CPU 사용량을 제어할 수 있습니다.

### 기본 사용법 (병렬 처리 비활성화)

```r
# 기본: 병렬 처리 완전 비활성화 (sequential only)
source('/home/user3/data_user3/git_repo/_wt/fgs/scripts/fgs/init_fgs_env.R')
```

### 병렬 처리 활성화 (제한된 CPU 사용)

```r
# 최대 8개 코어 사용, BLAS는 4 스레드
Sys.setenv(
  FGS_MAX_CPU_CORES = "8",
  FGS_BLAS_THREADS = "4",
  FGS_DISABLE_PARALLEL = "FALSE"
)
source('/home/user3/data_user3/git_repo/_wt/fgs/scripts/fgs/init_fgs_env.R')
```

### 메모리 한도 조정

```r
# 메모리 한도를 100GB로 설정
Sys.setenv(FGS_MEMORY_GB = "100")
source('/home/user3/data_user3/git_repo/_wt/fgs/scripts/fgs/init_fgs_env.R')
```

## 환경 변수 설명

| 변수명 | 기본값 | 설명 |
|--------|--------|------|
| `FGS_MAX_CPU_CORES` | `16` | 최대 사용할 CPU 코어 수 |
| `FGS_BLAS_THREADS` | `1` | BLAS/LAPACK 라이브러리 스레드 수 (1=sequential) |
| `FGS_DISABLE_PARALLEL` | `TRUE` | 병렬 처리 비활성화 여부 |
| `FGS_MEMORY_GB` | `200` | Future 메모리 한도 (GB) |

## 사용 예시

### 예시 1: 완전 순차 실행 (기본)

```r
source('/home/user3/data_user3/git_repo/_wt/fgs/scripts/fgs/init_fgs_env.R')
# 병렬 워커: 0, BLAS 스레드: 1
```

### 예시 2: 제한된 병렬 처리

```r
Sys.setenv(
  FGS_MAX_CPU_CORES = "4",
  FGS_BLAS_THREADS = "2",
  FGS_DISABLE_PARALLEL = "FALSE"
)
source('/home/user3/data_user3/git_repo/_wt/fgs/scripts/fgs/init_fgs_env.R')
# 최대 4개 코어, BLAS는 2 스레드
```

### 예시 3: 스크립트 실행 시 환경 변수 설정

```bash
# bash에서 실행
export FGS_MAX_CPU_CORES=8
export FGS_BLAS_THREADS=2
export FGS_DISABLE_PARALLEL=FALSE
cd /home/user3/GJC_KDW_250721
Rscript /home/user3/data_user3/git_repo/_wt/fgs/scripts/fgs/benchmark_l2_methods.R
```

## 주의사항

1. **자식 프로세스**: 환경 변수는 자식 프로세스에 상속되므로, `MYLIT_DISABLE_PARALLEL=TRUE`로 설정하면 모든 워커가 sequential로 실행됩니다.

2. **BLAS 스레드**: `FGS_BLAS_THREADS=1`로 설정하면 각 프로세스가 단일 스레드로 실행되어 연쇄 병렬화를 방지합니다.

3. **메모리**: 큰 데이터셋의 경우 메모리 한도를 충분히 설정하세요 (기본 200GB).

## start.R과의 관계

`init_fgs_env.R`은 내부적으로 `start.R`을 호출하지만, 환경 변수를 통해 병렬 설정을 제어합니다:
- `MYLIT_DISABLE_PARALLEL`: start.R의 future plan 설정 제어
- `MYLIT_FUTURE_MEMORY_GB`: 메모리 한도 설정

