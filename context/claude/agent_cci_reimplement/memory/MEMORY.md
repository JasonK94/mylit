# Project Memory - mylit / myR

## Project Structure
- **myR**: R package for Seurat-based scRNAseq analysis at `/data/user3/git_repo/mylit/myR/R/`
- **Docs**: Comprehensive bilingual (EN/KR) integrated guides at `docs/`
- **Master Rmd**: `/home/user3/GJC_KDW_250721/rmd/AfterPipe_260205.Rmd`
- **Context doc**: `docs/claude/1_STROKE_HC_V8_2_CONTEXT.md` - comprehensive analysis state
- **Project goals**: `docs/claude/4_PROJECT_GOALS.md` - 11-item analysis plan
- **Progress**: `docs/claude/7_PROGRESS_STATUS.md` - master progress tracker
- **CCI report**: `docs/claude/cci/CCI_ANALYSIS_REPORT.md` - CCI analysis details
- **Scripts index**: `docs/scripts/ANALYSIS_SCRIPTS_INDEX.md`

## Dataset: hc_only_v1 (HC-only independent analysis)
- 549,340 cells / 96 HC patients / 24 GEMs
- Seurat object: `/data/user3/sobj/hc_only_v1/5_hc_clean.qs`
- Integration: scVI (n_latent=30, batch=GEM) + MrVI
- 32 anno1 types, 10 anno2 compartments
- Demographics: age 19-66, sex M(39)/F(57)
- Docs: `docs/claude/HC/HC_ANALYSIS_CONTEXT.md` (master), `FGS_AGE_RESULTS.md`, `TREG_EXHAUSTION_ANALYSIS.md`
- Completed: MrVI, MASC, MILO, DEG Consensus, FGS continuous (6 methods), CCI (CellChat+MNN for age/sex), Treg exhaustion markers
- **Multi-method frequency analysis COMPLETE**: MELD, MILO, MASC(downsampled), scCODA, CellChat, MNN, Literature Validation
  - 5-method concordance: CD8 Naive↓, CD8 Effector↑, NK CD16hi↑, MAIT↓, CD8 TEM↑ — ALL 5/5 methods agree
  - MASC age FDR-sig: 8 cell types (CD8 Naive OR=0.47 top); sex: 1 nominal (NK CD16hi)
  - scCODA age credible: 6 cell types; sex credible: 1 (CD4 Th17/Th22)
  - Novel: Inflammatory Mono decline (anti-inflammaging), MAIT top-2 decline, ISG-hi age×sex interaction, male-specific naive decline, MIF sex dimorphism (FDR=0.006)
- Key FGS finding: PLAC8, EPSTI1, JAML, ISG15 decline with age; MSC, CDKN2A increase
- Exhaustion markers: No age association in Treg (too sparse, wrong cell context)
- FGS by compartment: T cell strongest (898 sig), then Myeloid (117), NK (111)
- CCI outputs: `/data/user3/sobj/hc_only_v1/cellchat/{sex_v1,age_v1}/`, `mnn/{sex_v1,age_v1}/`
- Results doc: `docs/claude/HC/HC_RESULTS_INTERIM.md` (comprehensive, 7 methods)

## Dataset: stroke_hc_v8_2
- 226,855 cells / 100 patients (Stroke 12GEM + HC 5GEM)
- Clean: 205,277 cells (`5_strokev8_clean.qs`)
- Layer 1 subset: `5_1_hc_is.qs` (131K: HC 73K + IS 58K)
- Layer 2 subset: `5_2_is_g3.qs` (54K: g3=1 15K + g3=2 39K)
- Integration: scVI (n_latent=30, batch=GEM)
- **g3** target variable has 167K NAs (IS-only); filter by `index_injury_hand=="IS"`
- anno1 (21 types after QC), anno2 (8 major compartments)

## Key Script Locations
- FGS: `Git_Repo/_wt/fgs/scripts/fgs/run_fgs_pipeline.R`
- CellChat: `Git_Repo/_wt/cellchat/scripts/cellchat/run_cellchat_cli.R`
- MultiNicheNet: `Git_Repo/_wt/cci/scripts/cci/mnn/run_multinichenet.R`
- MASC: `Git_Repo/_wt/masc/scripts/masc/run_masc.R`
- DEG Consensus: `Git_Repo/_wt/deg-consensus/scripts/consensus/run_deg_consensus_cli.R`
- Pseudotime functions: `myR/R/analysis/pseudotime.R` (NOT exported in NAMESPACE)
- Analysis scripts: `/data/user3/sobj/stroke_hc_v8_2/scripts/`

## User Preferences
- **Plots**: Heavy scatter plots (DimPlot etc.) → PNG (not PDF). PDF only for simple/vector plots
- **Reductions**: ALWAYS use `integrated.scvi` / `umap.scvi` from the original object. Do NOT recompute PCA/UMAP within subsets
- **Trajectory design**: Include HC vs Stroke comparison (Layer 1), not just g3 within IS
- **Stacked bars**: Sort by FGS signature score, Monocyte at bottom, T cell at top

## Known Issues & Fixes
- View() bug fixed: `import(stats)` in NAMESPACE caused BioGenerics override
- S4Vectors moved to Suggests to prevent as.data.frame override
- FGS: CPU limitation needed (taskset), NMF convergence issues
- R lib paths: `/home/user3/R/x86_64-pc-linux-gnu-library/4.3` + renv
- **dplyr::first() masking**: When sourcing pseudotime.R with Bioconductor packages, load dplyr LAST
- **Slingshot on large datasets**: Use umap (2D) not PCA (30D), add approx_points=150
- **myR pseudotime functions not exported**: Must `source()` directly, not `library(myR)`
- **R output buffering**: nohup logs may be delayed; check file mtime and output files directly
- **DEG Consensus muscat filter**: Changed filter="edgeR"→"none" (filter="genes" too aggressive for sparse pseudobulk)
- **NEBULA "cells must be grouped"**: Fixed by sorting count matrix/design by subject ID before nebula() call
- **GEM × cohort confounding**: GEM completely separated from cohort in L1. Remove GEM from covariates for L1 DEG. L2 (within IS) is OK.
- **Gene dynamics anova column name**: `anova()` returns `Pr(>Chi)` not `P(>|Chi|)` — fixed in pseudotime.R
- **MILO testNhoods needs PCA**: Add `reducedDim(milo, "PCA") <- reducedDim(milo, "SCVI")` before testNhoods
- **MILO plotNhoodGraphDA needs buildNhoodGraph()**: Must call `buildNhoodGraph(milo)` before plotting
- **MASC cluster names use dots**: `CD14..Monocyte` not `CD14+ Monocyte` — need explicit name mapping
- **L2 muscat too sparse**: muscat-edgeR/DESeq2 fail on L2 (g3) data; use NEBULA only for L2
- **CellChat v2 mergeCellChat drops @data**: Build CellChat from Seurat per condition, not from per-sample merge. Script: `scripts/reimplement_cellchat_v2.R`
- **CellChat ComplexHeatmap + png()**: `netVisual_heatmap`/`netAnalysis_signalingRole_heatmap` return HeatmapList; must explicitly `draw()` after capture
- **CellChat centrality lost on merge**: Compute centrality on individual objects, then lift→merge for scatter. Don't expect centrality on merged object
- **MNN fundamentally incompatible**: lr_target_prior_cor needs ≥5 samples/celltype/condition. Never met → 0 prioritized interactions
- **GEM gene vs GEM metadata**: Seurat DotPlot confuses "GEM" gene with GEM metadata column; exclude "GEM" from features
- **FGS compartment patient alignment**: After pseudobulk filtering, always align `patient_meta <- patient_meta[colnames(pb_counts), ]`
- **model.matrix single-level factor**: Use manual dummy construction when a factor may have only 1 level after filtering
- **sprintf %d for age**: Age is numeric/double, use `%.0f` not `%d`
- **MASC glmer convergence**: 549K cells → all NA. Must downsample to ~100-120K. Even then ~25/32 cell types may fail
- **NB-GAMM sparse genes**: CD4/CD8 T cells fitting S100A8/IL1B etc → convergence hang. Filter: skip genes with <1% expressing cells per condition. Also add setTimeLimit(300s)
- **Lamian patient vs GEM**: Lamian `cellanno$sample` must be patient_name (biological replicate), NOT GEM (batch). GEM contains 4-8 patients each
- **Lamian chisq vs permutation**: `test.method="chisq"` is ~100x faster than `"permutation"` with similar results. Use for 36-gene targeted tests
- **Lamian CDS has raw counts**: Monocle3 CDS `exprs()` returns raw counts (integers). Must log-normalize before Lamian: `log1p(count/lib_size * 10000)`
- **%||% not in base R**: Must define `\`%||%\` <- function(a,b) if(is.null(a)) b else a` in standalone scripts
- **slice_sample(n=n())**: `n()` can't be used in `slice_sample(n=)` argument; use `split/lapply/sample.int` instead
- **CellChat merged objects @data empty**: Per-sample merged objects have empty @data slot; aggregate from individual sample-level objects instead

## Completed Analyses (stroke_hc_v8_2)
- Pipeline (CellBender/Souporcell/Solo/scVI) ✓
- Clustering (res=2) + Annotation (anno1/anno2) ✓
- Annotation QC (doublet removal: clusters 13,22; LQ: 14,24,31; RBC: 27,34) ✓
- MASC (6 comparisons) ✓
- FindAllMarkers ✓
- CCI: CellChat v2 re-implementation ✓ — L1: 59 PNG, L2: 59 PNG (rankNet, bubble, scatter all working)
- CCI: MNN abandoned (0 prioritized interactions) — intermediate DE/ligand tables still usable
- Trajectory v2: Slingshot + Monocle3, scVI UMAP, 3 compartments ✓
- DEG Consensus L1: muscat-edgeR + muscat-DESeq2 + NEBULA (all 3 succeeded) ✓
  - 15 clusters with consensus, 5 skipped (insufficient methods)
- DEG L2: NEBULA only (203K results) ✓ — muscat too sparse for g3 data
- Cross-layer concordance: 134 genes significant in both L1 and L2 ✓
- Descriptive figures: 63 files (UMAP, DotPlot, bars, frequency, heatmaps, pseudotime, FGS) ✓
- Pathway enrichment (FGS genes): HALLMARK, KEGG, GOBP ✓
- FGS × DEG overlap: 218/241 FGS genes significant in L1 DEG ✓

## Running / Pending Analyses
- Gene dynamics v3: ALL COMPLETE (3 compartments × 2 conditions) ✓
- **Trajectory effect sizes**: ALL PHASES COMPLETE ✓ (1-4)
  - Phase 1-2: ABC/rABC/RMISE + comparison figures → `trajectory_v3/analysis/`
  - Phase 3: Lamian validation (chisq method) → `trajectory_v3/lamian/`
  - Phase 4: Pseudobulk GAMM sensitivity → `trajectory_v3/pseudobulk/`
  - Scripts: `stroke/scripts/{compute_trajectory_effect_sizes,run_lamian_validation,run_pseudobulk_gamm}.R`
  - Docs: `docs/claude/trajectory/` (4 documents, #4 = integrated results)
  - Key finding: 10 genes confirmed by all 3 methods (mono/cohort); mono is ONLY compartment with g3 effect
- scCODA: COMPLETE — L1 14/21 credible, L2 only CD14+ Mono ✓
- FGS n=100: COMPLETE ✓, n=200: running
- Within-cell-type FGS: CD4 T done, 14 more pending
- cNMF: K selection done for CD4 T + Inflam Mono, NK factorizing
- MILO: COMPLETE ✓ — L1: 11,739 nhoods (63% sig), L2: 4,664 nhoods (19% sig)
  - 3-method concordance (MILO×MASC×scCODA): 14/19 cell types DA by ≥2 methods (L1)
  - L2: CD14+ Mono only DA (MILO + scCODA concordant)
  - Output: `/data/user3/sobj/stroke_hc_v8_2/milo/{L1_cohort,L2_g3}/`
- MASC: anno1×project_name (running)
- Subclustering, MELD, Augur, MOFA+, scANVI, external validation still pending

## Environment
- renv path: `/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu`
- monocle3 is in renv, not default R library
- Set .libPaths() BEFORE library() calls
