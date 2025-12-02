#!/bin/bash
# Pipeline wrapper script
# Usage: ./pipe_wrapper.sh --config <config_path> [--run_id <run_id>] [--skip_steps <step1,step2,...>]

set -e  # Exit on error

# Default values
CONFIG=""
RUN_ID="run1"
SKIP_STEPS=""
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --config|-c)
      CONFIG="$2"
      shift 2
      ;;
    --run_id|-r)
      RUN_ID="$2"
      shift 2
      ;;
    --skip_steps|-s)
      SKIP_STEPS="$2"
      shift 2
      ;;
    --help|-h)
      echo "Usage: $0 --config <config_path> [--run_id <run_id>] [--skip_steps <step1,step2,...>]"
      echo ""
      echo "Options:"
      echo "  --config, -c       Path to config.csv file (required)"
      echo "  --run_id, -r       Run ID (default: run1)"
      echo "  --skip_steps, -s   Comma-separated list of steps to skip (e.g., '1,3')"
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

if [ -z "$CONFIG" ]; then
  echo "Error: --config argument is required"
  exit 1
fi

if [ ! -f "$CONFIG" ]; then
  echo "Error: Config file not found: $CONFIG"
  exit 1
fi

# Function to check if step should be skipped
should_skip_step() {
  local step=$1
  if [ -z "$SKIP_STEPS" ]; then
    return 1  # Don't skip
  fi
  IFS=',' read -ra STEPS <<< "$SKIP_STEPS"
  for s in "${STEPS[@]}"; do
    if [ "$s" == "$step" ]; then
      return 0  # Skip
    fi
  done
  return 1  # Don't skip
}

# Function to run step
run_step() {
  local step_num=$1
  local script_name=$2
  shift 2
  local extra_args="$@"
  
  if should_skip_step "$step_num"; then
    echo "Skipping step $step_num: $script_name"
    return 0
  fi
  
  echo "=========================================="
  echo "Running step $step_num: $script_name"
  echo "=========================================="
  
  local script_path="$SCRIPT_DIR/$script_name"
  if [ ! -f "$script_path" ]; then
    echo "Error: Script not found: $script_path"
    exit 1
  fi
  
  Rscript "$script_path" \
    --config "$CONFIG" \
    --run_id "$RUN_ID" \
    $extra_args
  
  if [ $? -ne 0 ]; then
    echo "Error: Step $step_num failed"
    exit 1
  fi
  
  echo "Step $step_num completed successfully"
  echo ""
}

# Step 0: Validation
echo "=========================================="
echo "Step 0: Validation"
echo "=========================================="
Rscript "$SCRIPT_DIR/pipe_validate.R" \
  --config "$CONFIG" \
  --run_id "$RUN_ID"

if [ $? -ne 0 ]; then
  echo "Error: Validation failed"
  exit 1
fi
echo "Validation passed"
echo ""

# Step 1: Read and demultiplex
run_step 1 "pipe1_read_demulti.R" \
  --input_step 0 \
  --output_step 1

# Step 2: LogNormalize normalization and clustering (for SoupX)
run_step 2 "pipe2_nmz_clustering.R" \
  --input_step 1 \
  --output_step 2

# Step 3: SoupX ambient RNA removal
run_step 3 "pipe3_ambient_removal.R" \
  --input_step 2 \
  --output_step 3

# Step 4: SCTransform normalization (after SoupX)
run_step 4 "pipe4_sctransform.R" \
  --input_step 3 \
  --output_step 4

# Step 5: Doublet detection
run_step 5 "pipe5_doubletfinder.R" \
  --input_step 4 \
  --output_step 5

# Step 6: Integration (RPCA)
run_step 6 "pipe6_integration.R" \
  --input_step 5 \
  --output_step 6 \
  --method RPCA

# Step 6: Integration (scVI) - optional, can be skipped if not needed
# Uncomment to run scVI integration as well
# run_step 6 "pipe6_integration.R" \
#   --input_step 5 \
#   --output_step 6 \
#   --method scVI

echo "=========================================="
echo "Pipeline completed successfully!"
echo "Run ID: $RUN_ID"
echo "Config: $CONFIG"
echo "=========================================="

