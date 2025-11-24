# Overview
Total package for scRNAseq / Spatial Transcriptomics analysis tool development.

## Details

### 1. Test Environment
When executing scripts with `cd /home/user3/GJC_KDW_250721 && Rscript`, there are no package dependency issues. (All installed in renv, start.R automatically loads them)
Cursor should directly test & debug loops and save analysis results to /data/user3/sobj. Overriding should be avoided as much as possible.

### 2. Data Used
- **stroke PBMC full data**: is5, location: /data/user3/sobj/IS6_sex_added_251110.qs
- **stroke PBMC downsampled data**: is5s, location: /data/user3/sobj/IS6_sex_added_0.1x_251110.qs
- **mIBD GeoMx data**: data_seurat, location: /data/user3/sobj/data_seurat_251104.qs

### 3. Variable Descriptions

#### 1) stroke
`g3` is the target variable ("2", "1", "NA"), patients are in `$hos_no` (8-digit number), cluster information is in `$anno3.scvi` (approximately 22 types). `$GEM` and `$set` are random effects other than patients, but since there are 23 patients, `$GEM` has a maximum of 12, and `$set` has a maximum of 3, they can be treated as fixed effects.

However, since `$hos_no < $GEM < $set` is completely nested, collinearity issues occur when all are included simultaneously, so be careful and select only some to use.

Examples:
* `gene_expression ~ g3 + anno3.scvi + (1|hos_no)`
* (After subsetting, e.g., `subset(anno3.scvi == "CD4+ T-cells")`) `gene_expression ~ g3 + (1|hos_no)`

#### 2) mIBD GeoMx
`$emrid` is the patient, with 114 samples (much fewer). With `$drug` (3 drugs), `$ck` (tissue type; PanCK positive or negative; TRUE, FALSE), `$treatment` ("pre", "post"), `$response` ("R", "NR"), there are many analysis considerations and the data is small, making it more challenging.

However, since the amount is small, testing may be faster. The number of unique AOIs (cells) per `$emrid` is about 2-8, so it may not be necessary to consider it as a random effect.

Examples:
* simple: (After subsetting by drug, tissue type) `gene_expression ~ treatment*response`
* full model: `gene_expression ~ treatment*ck*drug*response`

### 4. Notes
1) Each worktree should only commit its work content in its own branch, and should not affect other worktrees, especially the main or merger branch. They exist in their own folders in `_wt/`.
2) When merging each branch, docs and scripts should be managed separately in subfolders under `mylit/docs/` and `mylit/scripts/` to avoid conflicts.

### 5. Objectives
1) **Statistically rigorous yet powerful analysis (DEG analysis) tool development**: Mainly in myR/analysis.R, being developed in the analysis branch, `_wt/analysis` worktree. Previously the main2 branch/worktree.
2) **Tool development for DEG analysis consensus**: Source code mainly in `/data/user3/git_repo/_wt/deg-consensus/myR/R/deg_consensus`, test scripts in `/data/user3/git_repo/_wt/deg-consensus/scripts/`, being developed in the deg-consensus-dev branch, `_wt/deg-consensus` worktree.
3) **Signature and meta-learner tool development for group discrimination**: Source code managed and developed in `/data/user3/git_repo/_wt/fgs/myR/R/signature.R`, some unwanted duplicate functions exist in `/data/user3/git_repo/_wt/fgs/myR/R/test.R`.
4) **lds branch/worktree**: Measure and visualize variance explained by covariates through SVA after limma-dream. As of 251118-12:55, it is thought that other deg lists besides limma-dream should be passable. Patch not yet started.
   - `/data/user3/git_repo/_wt/lds/myR/R/lds_corrplot.R`, `/data/user3/git_repo/_wt/lds/myR/R/lds_08_heatmaps.R`, `/data/user3/git_repo/_wt/lds/myR/R/lds.R`, `/data/user3/git_repo/_wt/lds/scripts/lds`
5) **pt.umap branch/worktree**: Aggregate cell-level scRNAseq to patient-level data and analyze similarities and differences between patient groups through dimension reduction.
6) **Package usage convenience for CCI (Nichenetr) / MILO / trajectory analysis**: milo, cci, pseudotime-dev branches/worktrees
7) **Plotting function convenience**: plots-dev branch/worktree

### 6. Additional References
Refer to this for documentation rules:

@DOCS_ORGANIZATION_RULE.md
