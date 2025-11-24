# CPU 사용량 제어 가이드

## 빠른 참조

### 기본 사용 (병렬 처리 완전 비활성화)

```r
source('/home/user3/data_user3/git_repo/_wt/fgs/scripts/fgs/init_fgs_env.R')
```

### 병렬 처리 활성화 (제한된 CPU 사용)

```r
# 최대 8개 코어, BLAS는 4 스레드
Sys.setenv(
  FGS_MAX_CPU_CORES = "8",
  FGS_BLAS_THREADS = "4",
  FGS_DISABLE_PARALLEL = "FALSE"
)
source('/home/user3/data_user3/git_repo/_wt/fgs/scripts/fgs/init_fgs_env.R')
```

### 스크립트 실행 시 환경 변수 설정

```bash
# bash에서
export FGS_MAX_CPU_CORES=8
export FGS_BLAS_THREADS=2
export FGS_DISABLE_PARALLEL=FALSE
cd /home/user3/GJC_KDW_250721
Rscript /path/to/script.R
```

## 환경 변수 설명

| 변수명 | 기본값 | 설명 | 권장값 |
|--------|--------|------|--------|
| `FGS_MAX_CPU_CORES` | `16` | 최대 CPU 코어 수 | 4-16 |
| `FGS_BLAS_THREADS` | `1` | BLAS/LAPACK 스레드 (1=sequential) | 1-4 |
| `FGS_DISABLE_PARALLEL` | `TRUE` | 병렬 처리 비활성화 | TRUE/FALSE |
| `FGS_MEMORY_GB` | `200` | 메모리 한도 (GB) | 100-200 |

## 사용 시나리오

### 시나리오 1: 서버 공유 환경 (기본)

```r
# 다른 사용자와 서버를 공유하는 경우
# 병렬 처리 완전 비활성화 (기본 설정)
source('/home/user3/data_user3/git_repo/_wt/fgs/scripts/fgs/init_fgs_env.R')
```

### 시나리오 2: 전용 서버 (제한된 병렬 처리)

```r
# 서버를 혼자 사용하지만 CPU 사용량을 제한하고 싶은 경우
Sys.setenv(
  FGS_MAX_CPU_CORES = "8",      # 최대 8코어
  FGS_BLAS_THREADS = "2",       # BLAS는 2스레드
  FGS_DISABLE_PARALLEL = "FALSE"
)
source('/home/user3/data_user3/git_repo/_wt/fgs/scripts/fgs/init_fgs_env.R')
```

### 시나리오 3: 빠른 실행이 필요한 경우

```r
# 더 많은 CPU를 사용하여 빠르게 실행 (주의: 다른 프로세스에 영향)
Sys.setenv(
  FGS_MAX_CPU_CORES = "16",
  FGS_BLAS_THREADS = "4",
  FGS_DISABLE_PARALLEL = "FALSE"
)
source('/home/user3/data_user3/git_repo/_wt/fgs/scripts/fgs/init_fgs_env.R')
```

## 동작 원리

1. **환경 변수 설정**: `init_fgs_env.R`이 환경 변수를 읽어 설정
2. **start.R 호출**: 환경 변수를 전달하여 start.R이 병렬 설정 제어
3. **자식 프로세스 상속**: 모든 자식 프로세스가 동일한 환경 변수 상속
4. **연쇄 병렬화 방지**: 각 프로세스가 단일 스레드로 실행되거나 제한된 코어만 사용

## 확인 방법

실행 후 다음 메시지를 확인하세요:

```
=== FGS Environment Ready ===
  - start.R: sourced with MYLIT_DISABLE_PARALLEL=TRUE
  - Memory limit: 200GB
  - Parallel workers: 0 (sequential only)
```

또는 병렬 처리 활성화 시:

```
=== FGS Environment Ready ===
  - start.R: sourced with MYLIT_DISABLE_PARALLEL=FALSE
  - Memory limit: 200GB
  - Max CPU cores: 8
  - BLAS/LAPACK threads: 4
```

