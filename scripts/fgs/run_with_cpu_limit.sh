#!/bin/bash
# CPU 제한을 적용하여 R 스크립트 실행
# 사용법: ./run_with_cpu_limit.sh <R_script> [method] [cores]
# method: nice, taskset, cgroups (기본값: taskset)
# cores: 사용할 CPU 코어 수 (기본값: 8)

SCRIPT="$1"
METHOD="${2:-taskset}"
CORES="${3:-8}"

if [ -z "$SCRIPT" ]; then
    echo "사용법: $0 <R_script> [method] [cores]"
    echo "  method: nice, taskset, cgroups (기본값: taskset)"
    echo "  cores: 사용할 CPU 코어 수 (기본값: 8)"
    exit 1
fi

if [ ! -f "$SCRIPT" ]; then
    echo "오류: 스크립트 파일을 찾을 수 없습니다: $SCRIPT"
    exit 1
fi

# CPU 코어 수 확인
TOTAL_CORES=$(nproc)
if [ "$CORES" -gt "$TOTAL_CORES" ]; then
    echo "경고: 요청한 코어 수($CORES)가 시스템 코어 수($TOTAL_CORES)보다 큽니다."
    CORES=$TOTAL_CORES
fi

case "$METHOD" in
    nice)
        echo "=== nice를 사용하여 낮은 우선순위로 실행 ==="
        nice -n 19 Rscript "$SCRIPT"
        ;;
    taskset)
        echo "=== taskset을 사용하여 CPU 코어 $CORES개로 제한 ==="
        CORES_LIST=$(seq -s, 0 $((CORES-1)))
        taskset -c "$CORES_LIST" Rscript "$SCRIPT"
        ;;
    cgroups)
        echo "=== cgroups를 사용하여 CPU 제한 ==="
        # cgroups 설정 (root 권한 필요)
        if [ "$EUID" -ne 0 ]; then
            echo "오류: cgroups 사용은 root 권한이 필요합니다."
            echo "대신 'taskset' 방법을 사용하세요: $0 $SCRIPT taskset $CORES"
            exit 1
        fi
        
        # cgroup 생성 (이미 존재하면 무시)
        cgcreate -g cpu:/fgs_limit 2>/dev/null || true
        
        # CPU 할당량 설정 (예: 50% = 50000/100000)
        QUOTA=$((CORES * 100000 / TOTAL_CORES * 100000))
        echo "$QUOTA" > /sys/fs/cgroup/cpu/fgs_limit/cpu.cfs_quota_us
        echo "100000" > /sys/fs/cgroup/cpu/fgs_limit/cpu.cfs_period_us
        
        # 스크립트 실행
        cgexec -g cpu:/fgs_limit Rscript "$SCRIPT"
        
        # 정리 (선택사항)
        # cgdelete -g cpu:/fgs_limit
        ;;
    *)
        echo "오류: 알 수 없는 방법: $METHOD"
        echo "사용 가능한 방법: nice, taskset, cgroups"
        exit 1
        ;;
esac

