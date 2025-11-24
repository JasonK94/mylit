# 2025-11-21 Debugging & Development Plan

## Objective
Debug FGS/TML7 issues (ranger, nmf, earth, cmgi) and verify robustness across all clusters.

## Workflow

### 1. CMGI Worktree (Branch: `cmgi`)
- [ ] **Fix `earth` package missing issue**:
    - Use `renv::activate("/home/user3/GJC_KDW_250721")`.
    - Verify `earth` package availability.
    - Install if missing.
- [ ] **Verify CMGI with `earth`**:
    - Ensure `compute_meta_gene_importance` works correctly for `earth` models.
- [ ] **Merge to `fgs`**:
    - Once verified, merge `cmgi` branch into `fgs`.

### 2. FGS Worktree (Branch: `fgs`)
- [ ] **Resolve Failed Clusters**:
    - Create a script to iterate FGS over all clusters using full data (`is5a`).
    - Log errors for failed clusters.
    - Debug and fix issues.

### 3. L1 Ranger/NMF Worktree (Branch: `l1_ranger/nmf`)
- [ ] **Verify NMF & Ranger**:
    - Run FGS with `nmf_loadings` and `random_forest_ranger` for all clusters.
    - Confirm valid signature generation.
- [ ] **Compare Ranger vs Random Forest**:
    - Check consistency between `random_forest_ranger` (ranger) and `random_forest` (randomForest).
- [ ] **Merge to `fgs`**:
    - If successful, merge into `fgs`.

### 4. TML Worktree (Branch: `tml`)
- [ ] **Review Documentation**:
    - Analyze `docs/tml` content.
- [ ] **Develop TML Improvements**:
    - Implement planned features.
- [ ] **Verify CMGI for L2 Ranger**:
    - Ensure CMGI works for L2 ranger models (with `importance` parameter).

## Execution Log
- 2025-11-21: Plan created.
