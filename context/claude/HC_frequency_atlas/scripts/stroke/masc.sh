#!/bin/bash

# 공통 설정
INPUT="/data/user3/sobj/stroke_hc_v8_2/4_strokev8_annotated.qs"
OUTPUT_DIR="/data/user3/sobj/stroke_hc_v8_2/MASC/"
RSCRIPT_PATH="/home/user3/data_user3/git_repo/mylit/Git_Repo/_wt/masc/scripts/masc/run_masc.R"

# 반복할 리스트 정의
CLUSTERS=("scvi_snn_res.2" "anno1" "anno2")
COMPARISONS=("project_name" "g3")

# 출력 디렉토리가 없으면 생성
mkdir -p $OUTPUT_DIR

# 중첩 루프 시작
for CLUSTER in "${CLUSTERS[@]}"; do
    for COMP in "${COMPARISONS[@]}"; do
        
        # 파일명 접두사 동적 생성 (예: MASC_anno1_g3)
        PREFIX="MASC_${CLUSTER}_${COMP}"
        
        echo ">>> Running: $PREFIX"
        
        # Rscript 실행
        Rscript $RSCRIPT_PATH \
          -i "$INPUT" \
          -o "$OUTPUT_DIR" \
          --cluster_var "$CLUSTER" \
          --contrast_var "$COMP" \
          --random_effects patient_name \
          --fixed_effects GEM,age,sex \
          --prefix "$PREFIX"
          
    done
done

echo ">>> All tasks completed!"