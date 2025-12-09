# Pipeline Archive

This directory contains historical files that are no longer actively used but preserved for reference.

## Contents

### Old Scripts (`old_scripts/`)
Helper and utility scripts that have been superseded or are no longer needed for regular pipeline execution:
- `prep_metadata.R` - Metadata preprocessing (now integrated into main pipeline)
- `fix_manifest_hto.R` - One-time manifest correction script
- `update_manifest.R` - Legacy manifest update utility
- `update_config_columns.R` - Configuration synchronization (deprecated)
- `verify_meta.R` - Metadata verification (for debugging)
- `create_downsampled_data.R` - Testing data generator
- `vars_config.R` - Variable configuration utility
- `load_renv.R` - renv loading script
- `start_pipe.R` - Old startup script (replaced by `scripts/pipe/start_pipe.R`)

### Run Logs and Test Scripts
Historical execution logs and test run scripts:
- `run_pipeline_*.log` - Pipeline execution logs from various test runs
- `run_pipeline_*.sh` - Pipeline execution wrapper scripts
- `install_log*.txt` - Installation and setup logs

## Note
These files are kept for historical reference and debugging purposes. For current pipeline execution, refer to the main documentation in `docs/pipe/README.md`.
