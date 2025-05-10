#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n'

#------------------------------------------------------------------------------
# Usage 함수: 잘못된 사용법 전달 시 안내
usage() {
  cat <<EOF
Usage: $0 <input_file>

Supported formats: .csv, .tsv, .xlsx

Columns must be:
  env,kernel,script,option1,option2,...

Requires:
  - conda (for 'conda run')
  - xlsx2csv or in2csv (for .xlsx)
EOF
  exit 1
}

# 인자 체크
[ $# -eq 1 ] || usage
input="$1"
[[ -f "$input" ]] || { echo "[ERROR] File not found: $input"; exit 1; }

# 로그 함수: 타임스탬프 포함
log() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

# 임시 파일 cleanup 설정
tmpcsv=""
trap '[[ -n "$tmpcsv" && -f "$tmpcsv" ]] && rm -f "$tmpcsv"' EXIT

# 확장자 추출 (소문자로)
ext="${input##*.}"
ext="${ext,,}"
case "$ext" in
  csv)   delimiter=','; processor_cat="cat" ;;  
  tsv)   delimiter=$'\t'; processor_cat="cat" ;;  
  xls|xlsx)
    # 중간 CSV 생성
    tmpcsv="$(mktemp --suffix=.csv)"
    if command -v xlsx2csv &>/dev/null; then
      xlsx2csv "$input" "$tmpcsv"
    elif command -v in2csv &>/dev/null; then
      in2csv "$input" > "$tmpcsv"
    else
      log "[ERROR] xlsx2csv or in2csv required to process .$ext"
      exit 1
    fi
    input="$tmpcsv"
    delimiter=','; processor_cat="cat" ;;  
  *)
    log "[ERROR] Unsupported extension: .$ext"
    exit 1 ;;  
esac

# CSV/TSV 처리: 헤더 건너뛰고, 빈줄·주석 제거
$processor_cat "$input" \
  | tail -n +2 \
  | sed -e '/^\s*$/d' -e '/^\s*#/d' \
  | while IFS="$delimiter" read -ra cols; do

    env_name="${cols[0]}"
    kernel="${cols[1]}"
    script="${cols[2]}"
    options=("${cols[@]:3}")

    # 빈 옵션 필터링
    cmd_opts=()
    for o in "${options[@]}"; do
      [[ -n "$o" ]] && cmd_opts+=("$o")
    done

    log "[INFO] env=$env_name, kernel=$kernel, script=$script"
    log "[INFO] opts=${cmd_opts[*]}"

    # conda 실행 prefix
    if command -v conda &>/dev/null; then
      conda_prefix=(conda run -n "$env_name")
    else
      log "[WARN] conda not in PATH; fallback activation"
      source ~/.bashrc
      conda_prefix=(conda activate "$env_name" '&&')
    fi

    # R 스크립트는 Rscript로
    if [[ "${kernel,,}" == "r" ]]; then
      interpreter=(Rscript)
    else
      interpreter=("$kernel")
    fi

    # 커맨드 조립 및 실행
    cmd=( "${conda_prefix[@]}" "${interpreter[@]}" "$script" "${cmd_opts[@]}" )
    log "[INFO] Running: ${cmd[*]}"
    "${cmd[@]}"
  done
