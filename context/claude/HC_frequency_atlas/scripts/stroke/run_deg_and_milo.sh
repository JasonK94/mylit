#!/bin/bash
# ==============================================================================
# DEG Consensus + MILO Pipeline
# Runs both Layer 1 (HC vs IS) and Layer 2 (g3==1 vs g3==2)
# ==============================================================================

set +e  # Don't exit on error, allow subsequent steps to continue

BASE="/data/user3/sobj/stroke_hc_v8_2"
DEG_SCRIPT="/home/user3/data_user3/git_repo/mylit/Git_Repo/_wt/deg-consensus/scripts/consensus/run_deg_consensus_cli.R"
DEG_SRC="/home/user3/data_user3/git_repo/mylit/Git_Repo/_wt/deg-consensus/myR/R/deg_consensus"
MASC_SCRIPT="/home/user3/data_user3/git_repo/mylit/Git_Repo/_wt/masc/scripts/masc/run_masc.R"

# ---- DEG Consensus: Layer 1 (HC vs IS at anno1) ----
# NOTE: GEM removed from covariates — GEM is completely confounded with cohort
# (HC GEMs and Stroke GEMs are entirely separate). Including GEM causes collinearity.
echo ">>> [1/5] DEG Consensus: Layer 1 (cohort, anno1)"
mkdir -p "$BASE/deg_consensus/L1_cohort_anno1"
taskset -c 16-31 Rscript "$DEG_SCRIPT" \
  --input "$BASE/5_1_hc_is.qs" \
  --output "$BASE/deg_consensus/L1_cohort_anno1" \
  --cluster anno1 \
  --group cohort \
  --contrast "Stroke - HC" \
  --covariates "sex,age" \
  --sampleid patient_name \
  --methods "muscat-edgeR,muscat-DESeq2,nebula" \
  --cores 16 \
  --src "$DEG_SRC" \
  2>&1 | tee "$BASE/logs/deg_L1_cohort_anno1.log"

# ---- DEG Consensus: Layer 2 (g3 at anno1) ----
echo ">>> [2/5] DEG Consensus: Layer 2 (g3, anno1)"
mkdir -p "$BASE/deg_consensus/L2_g3_anno1"
taskset -c 16-31 Rscript "$DEG_SCRIPT" \
  --input "$BASE/5_2_is_g3.qs" \
  --output "$BASE/deg_consensus/L2_g3_anno1" \
  --cluster anno1 \
  --group g3 \
  --contrast "2 - 1" \
  --covariates "sex,age,GEM" \
  --sampleid patient_name \
  --methods "muscat-edgeR,muscat-DESeq2,nebula" \
  --cores 16 \
  --src "$DEG_SRC" \
  2>&1 | tee "$BASE/logs/deg_L2_g3_anno1.log"

# ---- MASC: Missing comparisons (anno1/seurat_clusters × project_name) ----
echo ">>> [3/5] MASC: anno1 × project_name"
Rscript "$MASC_SCRIPT" \
  -i "$BASE/5_strokev8_clean.qs" \
  -o "$BASE/MASC/" \
  --cluster_var anno1 \
  --contrast_var project_name \
  --random_effects patient_name \
  --fixed_effects "GEM,age,sex" \
  --prefix "MASC_anno1_project_name" \
  2>&1 | tee "$BASE/logs/masc_anno1_project.log"

echo ">>> [4/5] MASC: seurat_clusters × project_name"
Rscript "$MASC_SCRIPT" \
  -i "$BASE/5_strokev8_clean.qs" \
  -o "$BASE/MASC/" \
  --cluster_var "scvi_snn_res.2" \
  --contrast_var project_name \
  --random_effects patient_name \
  --fixed_effects "GEM,age,sex" \
  --prefix "MASC_seurat_project_name" \
  2>&1 | tee "$BASE/logs/masc_seurat_project.log"

echo ">>> [5/5] All DEG + MASC tasks completed!"
echo "Output: $BASE/deg_consensus/ and $BASE/MASC/"
