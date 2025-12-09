#!/bin/bash
set -e

CONFIG="config/manifest_stroke.csv"
RUN_ID="run_stroke_test_v1"
DOWNSAMPLE="0.05"

echo "Starting pipeline test with Run ID: $RUN_ID and Downsample: $DOWNSAMPLE"

echo "Running Step 1: Read and Demultiplex"
Rscript scripts/pipe/pipe1_read_demulti.R --config $CONFIG --run_id $RUN_ID --input_step 0 --output_step 1 --downsample $DOWNSAMPLE

echo "Running Step 2: Normalize and Cluster"
Rscript scripts/pipe/pipe2_nmz_clustering.R --config $CONFIG --run_id $RUN_ID --input_step 1 --output_step 2

echo "Running Step 3: Ambient RNA Removal"
Rscript scripts/pipe/pipe3_ambient_removal.R --config $CONFIG --run_id $RUN_ID --input_step 2 --output_step 3

echo "Running Step 4: SCTransform"
Rscript scripts/pipe/pipe4_sctransform.R --config $CONFIG --run_id $RUN_ID --input_step 3 --output_step 4

echo "Running Step 5: DoubletFinder"
Rscript scripts/pipe/pipe5_doubletfinder.R --config $CONFIG --run_id $RUN_ID --input_step 4 --output_step 5

echo "Running Step 6: Integration"
Rscript scripts/pipe/pipe6_integration.R --config $CONFIG --run_id $RUN_ID --input_step 5 --output_step 6

echo "Pipeline test completed successfully!"
