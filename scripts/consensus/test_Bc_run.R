#!/usr/bin/env Rscript
# scripts/consensus/test_Bc_run.R
# Test DEG Consensus CLI on Bc dataset

renv_lib <- "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu"
if (dir.exists(renv_lib)) .libPaths(c(renv_lib, .libPaths()))

# Input/Output
input_file <- "/data/user3/sobj/local_stroke_snp1/1_obj/4_5_3_Bc_processed.qs"
output_dir <- "/data/user3/sobj/consensus/test_Bc"

# Command
cmd <- sprintf(
    "Rscript scripts/consensus/run_deg_consensus_cli.R \\
    --input '%s' \\
    --output '%s' \\
    --cluster 'anno3' \\
    --group 'g3' \\
    --covariates 'sex,age,GEM' \\
    --methods 'edgeR-LRT,DESeq2-Wald,limma-voom' \\
    --contrast '2 - 1' \\
    --sampleid 'sample_name' \\
    --cores 4",
    input_file, output_dir
)

message("Running Test Command:")
message(cmd)
system(cmd)
