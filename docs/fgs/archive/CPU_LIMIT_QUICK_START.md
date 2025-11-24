# CPU 제한 빠른 시작 가이드

TML7에서 CPU 완전 점유 문제가 발생할 때 사용하는 방법입니다.

## 문제 증상

- `FGS_DISABLE_PARALLEL=FALSE`로 설정해도 CPU가 완전히 점유됨
- htop에서 CPU 활용 개수가 즉시적, 동시적으로 증가
- **특히 xgboost (xgbTree)가 주요 원인**: 다른 모델들은 괜찮지만 xgboost만 문제 발생
- 모델 패키지(xgboost, ranger 등)가 C/C++ 레벨에서 병렬 처리 시작

## 해결 방법

### 방법 1: taskset 사용 (권장, 가장 간단)

특정 CPU 코어만 사용하도록 제한:

```bash
# 8개 코어만 사용 (코어 0-7)
taskset -c 0-7 Rscript /path/to/script.R

# 또는 헬퍼 스크립트 사용
cd /home/user3/data_user3/git_repo/_wt/fgs
./scripts/fgs/run_with_cpu_limit.sh scripts/fgs/benchmark_l2_methods.R taskset 8
```

### 방법 2: nice 사용

프로세스 우선순위를 낮춰서 다른 프로세스에 영향 최소화:

```bash
nice -n 19 Rscript /path/to/script.R

# 또는 헬퍼 스크립트 사용
./scripts/fgs/run_with_cpu_limit.sh scripts/fgs/benchmark_l2_methods.R nice
```

### 방법 3: cgroups 사용 (root 권한 필요)

CPU 사용량을 백분율로 제한:

```bash
# 50% CPU로 제한 (예: 8코어 시스템에서 4코어만 사용)
sudo cgcreate -g cpu:/fgs_limit
echo 50000 > /sys/fs/cgroup/cpu/fgs_limit/cpu.cfs_quota_us
echo 100000 > /sys/fs/cgroup/cpu/fgs_limit/cpu.cfs_period_us
sudo cgexec -g cpu:/fgs_limit Rscript /path/to/script.R

# 또는 헬퍼 스크립트 사용 (root 권한 필요)
sudo ./scripts/fgs/run_with_cpu_limit.sh scripts/fgs/benchmark_l2_methods.R cgroups 4
```

## 헬퍼 스크립트 사용법

```bash
cd /home/user3/data_user3/git_repo/_wt/fgs

# 기본 사용 (taskset, 8코어)
./scripts/fgs/run_with_cpu_limit.sh scripts/fgs/benchmark_l2_methods.R

# nice 방법 사용
./scripts/fgs/run_with_cpu_limit.sh scripts/fgs/benchmark_l2_methods.R nice

# taskset으로 4코어만 사용
./scripts/fgs/run_with_cpu_limit.sh scripts/fgs/benchmark_l2_methods.R taskset 4

# cgroups로 50% CPU 제한 (root 필요)
sudo ./scripts/fgs/run_with_cpu_limit.sh scripts/fgs/benchmark_l2_methods.R cgroups 4
```

## 권장 설정

### 개발/테스트 환경
```bash
# 4-8 코어만 사용
taskset -c 0-7 Rscript script.R
```

### 프로덕션 환경
```bash
# nice로 우선순위 낮추기 + taskset으로 코어 제한
nice -n 19 taskset -c 0-7 Rscript script.R
```

## xgboost 특별 주의사항

**⚠️ 중요**: xgboost (xgbTree)는 다른 모델들과 달리 R 레벨 제어만으로는 CPU 점유를 완전히 막을 수 없습니다.

- `nthread=1` 파라미터 설정
- `xgb.set.config(nthread=1)` 전역 설정
- `XGBOOST_NTHREAD=1` 환경 변수 설정

위의 모든 설정을 해도 xgboost는 C/C++ 레벨에서 병렬 처리를 시작할 수 있습니다.

**해결책**: **taskset을 사용한 시스템 레벨 제한이 필수입니다.**

```bash
# xgboost를 포함한 TML7 실행 시 반드시 taskset 사용
taskset -c 0-7 Rscript script.R
```

## 왜 R 레벨 제어만으로는 부족한가?

1. **xgboost의 C/C++ 레벨 병렬 처리 (특히 문제)**
   - xgboost는 C/C++로 작성되어 내부적으로 스레드를 생성
   - `nthread=1` 파라미터, 전역 설정, 환경 변수 모두 설정해도 부족할 수 있음
   - **시스템 레벨 제한(taskset)이 유일한 확실한 해결책**

2. **다른 모델 패키지들**
   - `ranger`, `glmnet` 등은 R 레벨 제어로 대부분 제어 가능
   - xgboost만 특별히 문제가 되는 경우가 많음

3. **caret의 CV fold 처리**
   - caret이 CV fold를 처리할 때 각 fold마다 모델 학습
   - 모델 패키지가 각 fold마다 병렬 처리를 시작할 수 있음

4. **시스템 레벨 제한의 장점**
   - 운영체제 레벨에서 강제로 CPU 사용량 제한
   - R 코드나 패키지 설정과 무관하게 작동
   - 가장 확실한 방법

## 추가 참고

- 상세한 전략은 `docs/fgs/TML7_CPU_FIX_STRATEGIES.md` 참고
- CPU 설정 가이드는 `docs/fgs/CPU_CONFIGURATION.md` 참고

