# Pipeline Context & Status

## Overview
This pipeline processes single-cell RNA-seq data, including demultiplexing, normalization, ambient RNA removal, doublet detection, and integration.

## Current Status & Issues

### Step 1: Read & Demultiplexing (`pipe1_read_demulti.R`)
- **Status**: Previous debugging consumed significant effort due to incomplete demultiplexing data processing.
- **Action**: Verify if the current logic handles demultiplexing correctly.

### Step 3: SoupX Ambient RNA Removal (`pipe3_ambient_removal.R`)
- **Status**: The limit for parameter adjustment is harsh.
- **Current Logic**: Decrements parameter by 0.1 until a plot is generated.
- **Issue**: If the parameter gets too low, it should pass `adjustCount`.
- **Action**: Check and modify the parameter adjustment logic to be more robust.

### Step 4 & 5: SCTransform & Doublet Detection
- **Status**: Believed to be working correctly.
- **Potential Improvements**: 
    - Consider alternative ambient RNA removal methods besides SoupX.
    - Consider alternative doublet finder modules.
    - Consider converting `qs` to `h5ad` to use Python modules for these steps.

### Step 6: Integration (`pipe6_integration.R`)
- **Status**: Suspicious results with scVI.
- **Issue**: Batch correction is poor. GEM1~4 and GEM5~8 (different sets) are separated in UMAP.
- **Hypothesis**: Batch factor "GEM" might be too granular; might need a larger "SET" factor.
- **Observation**: Within-batch clustering is good. Previously, removing doublets improved results significantly.

## Goals
1. Run step-by-step checks with `run_id=run_ag1`.
2. Debug and fix Step 1 if necessary.
3. Refine Step 3 parameter adjustment logic.
4. Investigate Step 6 batch correction issues (scVI parameters, batch factors).
5. Explore `qs` -> `h5ad` conversion for Python module integration.

## References
- [COMMANDS.md](COMMANDS.md) - Common commands and usage
- [DATA_FLOW.md](DATA_FLOW.md) - Data flow and directory structure
- [PIPE_INTEGRATED_GUIDE_KR.md](PIPE_INTEGRATED_GUIDE_KR.md) - 통합 가이드 (국문)
- [PIPELINE_REVIEW.md](PIPELINE_REVIEW.md) - Code review notes
- [README.md](README.md) - Documentation index

