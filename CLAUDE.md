# CLAUDE.md — mylit Repository

## Repository Overview

**mylit** is a monorepo for single-cell RNA-seq (scRNAseq) analysis, centred on
the `myR` R package and associated analysis workflows. The primary dataset is
`stroke_hc_v8_2` (226K cells, 100 patients: ischemic stroke vs healthy controls).

## Package: myR

An R package (`myR/`) providing Seurat-based scRNAseq utilities.

### Directory Structure
```
myR/
├── R/
│   ├── analysis/
│   │   ├── pseudotime.R       # Slingshot, Monocle3, gene dynamics (v1 + v2 GAMM)
│   │   ├── deg.R              # DEG methods (MAST, muscat, NEBULA, pseudobulk)
│   │   ├── cellchat.R         # CellChat wrappers
│   │   ├── milo.R             # MILO neighbourhood DA
│   │   └── nichenet.R         # NicheNet analysis
│   ├── fgs/                   # Feature/Gene Signature discovery
│   ├── plots.R                # Plotting utilities (DimPlot wrappers, heatmaps)
│   ├── seurat_to_scanpy.R     # h5ad ↔ Seurat conversion
│   ├── demultiplex.R          # Demultiplexing utilities
│   ├── validation.R           # Input validation helpers
│   └── zzz.R                  # Package load hooks
├── DESCRIPTION
├── NAMESPACE
└── man/                       # roxygen2 docs
```

### Key Functions

| Function | File | Purpose |
|----------|------|---------|
| `run_slingshot_from_seurat()` | pseudotime.R:54 | Slingshot trajectory on Seurat object |
| `run_monocle3_from_seurat()` | pseudotime.R:1401 | Monocle3 trajectory with custom reduction |
| `analyze_gene_dynamics()` | pseudotime.R:591 | Gene dynamics GAM (v1, no batch correction) |
| `analyze_gene_dynamics_v2()` | pseudotime.R:873 | **Batch-corrected GAMM** (NB + offset + GEM RE) |
| `runMAST()` / `runMUSCAT()` | deg.R | DEG analysis wrappers |
| `load_h5ad_to_seurat_qs()` | seurat_to_scanpy.R | h5ad → Seurat conversion |

### Known Issues
- Pseudotime functions (`run_slingshot_from_seurat`, `run_monocle3_from_seurat`,
  `analyze_gene_dynamics*`) are NOT exported in NAMESPACE. Must `source()` directly.
- `S4Vectors` in Suggests (not Imports) to prevent `as.data.frame` override.
- `dplyr::first()` masked by Bioconductor packages — load dplyr LAST.
- `myR/R/cci_multinichenet_wrapper.R` has duplicate `p_val_thresh` param (line 53) — known bug, use CLI script instead.

## Analysis Workflows (Git_Repo/_wt/)

Standalone CLI scripts for each analysis type:

| Analysis | Script | Key Args |
|----------|--------|----------|
| FGS | `_wt/fgs/scripts/fgs/run_fgs_pipeline.R` | `--input`, `--target_var`, `--n_features` |
| CellChat | `_wt/cellchat/scripts/cellchat/run_cellchat_cli.R` | `-i`, `-a`, `--subset_aggregate` |
| MultiNicheNet | `_wt/cci/scripts/cci/mnn/run_multinichenet.R` | `-i`, `-f` (contrast) |
| MASC | `_wt/masc/scripts/masc/run_masc.R` | `--input`, `--cluster_col`, `--comparison_col` |
| DEG Consensus | `_wt/deg-consensus/scripts/consensus/run_deg_consensus_cli.R` | multi-method DEG |

## Environment

- R 4.3.2 at `/opt/R/4.3.2/`
- User libraries: `/home/user3/R/x86_64-pc-linux-gnu-library/4.3`
- renv libraries: `/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu`
- monocle3 is only in renv — must set `.libPaths()` BEFORE `library()` calls
- GPU: RTX 4090 (used by scVI in Python pipeline)
- Python: conda env `scvi-tools` for Python-based analyses (scCODA, cNMF, MELD)

### R Script Boilerplate
```r
.libPaths(c(
  "/home/user3/GJC_KDW_250721/renv/library/R-4.3/x86_64-pc-linux-gnu",
  "/home/user3/R/x86_64-pc-linux-gnu-library/4.3",
  .libPaths()
))
```

### Python Script Boilerplate
```bash
conda activate scvi-tools
# pertpy (for scCODA), cnmf, meld, scanpy, anndata available
```

## Conventions

- **Plots**: Heavy scatter plots (DimPlot, trajectory) → PNG. Simple/vector → PDF.
- **Reductions**: Always use `integrated.scvi` / `umap.scvi` from the original object.
  Never recompute PCA/UMAP within subsets — HC/Stroke batch effects corrupt de novo embeddings.
- **Dual-layer design**: Layer 1 = HC vs Stroke (cohort), Layer 2 = IS g3==1 vs g3==2.
- **Commit style**: Short imperative subject + optional body with details.
- **Language**: User communicates in Korean, documentation is bilingual (EN primary).
- **CPU control**: Use `taskset -c 0-15` or `16-31` to limit CPU affinity for long-running jobs.
- **Background jobs**: `nohup Rscript script.R > log.log 2>&1 &` — check output files (not just log) for completion.

## Current Analysis State

**Master progress tracker**: `docs/claude/7_PROGRESS_STATUS.md`

### Dataset: stroke_hc_v8_2 (primary)

See `docs/claude/1_STROKE_HC_V8_2_CONTEXT.md` for full state.

**Data locations**:
- Clean object: `/data/user3/sobj/stroke_hc_v8_2/5_strokev8_clean.qs` (205K cells)
- L1 subset: `/data/user3/sobj/stroke_hc_v8_2/5_1_hc_is.qs` (131K: HC 73K + IS 58K)
- L2 subset: `/data/user3/sobj/stroke_hc_v8_2/5_2_is_g3.qs` (54K: g3=1 15K + g3=2 39K)
- Results symlink: `results/` → `/data/user3/sobj/stroke_hc_v8_2`

**Completed**: Pipeline, Clustering, Annotation, QC, MASC (6 runs), FindAllMarkers,
DEG Consensus (L1 3 methods, L2 NEBULA), Trajectory v2+v3, Gene dynamics v3 (all),
CCI CellChat v2 (re-implemented, rankNet+bubble working), FGS (n=50, n=100), scCODA,
Descriptive figures (63 files)

**In Progress**: FGS n=200, Within-cell-type FGS (CD4 done, 14 pending), cNMF (3 cell types)

**Pending**: MNN v3 visualization, MILO v2, Subclustering,
MELD, Augur, MOFA+

**Recently Fixed**: MNN (multinichenetr) — contrast string mismatch bug fixed via apostrophe-wrapping.
v3 results: L1 anno1 148,958 / L1 anno2 8,244 / L2 anno1 91,111 / L2 anno2 7,600 prioritized interactions.
Output: `cci/mnn/{L1_cohort_anno1_v3, L1_cohort_anno2_v3, L2_g3_anno1_v3, L2_g3_anno2_v3}/`

### Dataset: hc_only_v1 (HC-only independent atlas)
See `docs/claude/HC/HC_ANALYSIS_CONTEXT.md` for full state.
- **549K cells / 96 patients / 24 GEMs** (Korean HC, age 19-66, sex M39/F57)
- **32 anno1 types, 10 anno2 compartments**
- Seurat object: `/data/user3/sobj/hc_only_v1/2_hc_annotated.qs`
- Integration: scVI + MrVI

---

## Analysis Scripts & Results (stroke_hc_v8_2)

### Executed Analysis Scripts
All scripts in `/data/user3/sobj/stroke_hc_v8_2/scripts/`:

| Script | Purpose | Status |
|--------|---------|--------|
| `run_cci_all.sh` | CellChat + MNN dual-layer (anno2) | ✅ Done (v1) |
| `run_mnn_anno1.sh` | MNN re-run at anno1 (19 types) | ✅ Done |
| `reimplement_cellchat_v2.R` | **CellChat v2 from Seurat** (rankNet+bubble) | ✅ Done |
| `rerun_mnn_all_v3.sh` | MNN v3 with apostrophe fix (all 4 configs) | ✅ Done (256K interactions) |
| `plot_cci_comprehensive.R` | CellChat native plots (pathway circles, bubbles) | ✅ Done |
| `plot_mnn_anno1.R` | MNN DE volcanos, circos, heatmaps | ✅ Done |
| `plot_fgs_results.R` | FGS 8-panel visualization | ✅ Done |

### CCI Status (UPDATED)

**CellChat v2** (re-implemented, USE THIS):
- **Condition-level CellChat built from Seurat directly** → `@data` populated, all comparison functions work
- `rankNet()`, `netVisual_diffInteraction()`, `netVisual_bubble()` all **WORKING** ✅
- Output: `results/cci/plots/cellchat_{L1_cohort,L2_g3}_v2/` (7 main + 36 pathway per layer)
- Saved comparison objects: `cellchat_comparison.qs` + per-condition `cellchat_{HC,Stroke}.qs`
- Still missing: `netAnalysis_signalingRole_heatmap` (centrality lost in merge)
- See `docs/claude/cci/CCI_AUDIT_AND_REIMPLEMENTATION.md` for full audit

**MultiNicheNet** (FIXED ✅, v3):
- Contrast string bug fixed: apostrophe-wrapping `"'Stroke-HC'"` preserves column name in makeContrasts
- **v3 results**: L1 anno1: 148,958 / L1 anno2: 8,244 / L2 anno1: 91,111 / L2 anno2: 7,600 interactions
- Output: `cci/mnn/{L1_cohort_anno1_v3, L1_cohort_anno2_v3, L2_g3_anno1_v3, L2_g3_anno2_v3}/`
- Top L1: Mono→Mono S100A8-CD36 (score=0.978), TYROBP-KLRD1 (Mono→NKc)
- See `docs/claude/cci/MNN_FAILURE_ANALYSIS.md` for full diagnosis and fix

### Per-Sample CellChat Objects
```
results/cci/cellchat/L1_cohort_anno2/samples/
├── h19457220/cellchat.qs   # HC sample
├── h19877398/cellchat.qs   # HC sample
├── ...                     # 20 HC + 36 Stroke = 56 total
results/cci/cellchat/L2_g3_anno2/samples/
├── ...                     # 32 IS samples
```

---

## Additional Analyses: Implementation Guide

These analyses are planned or in progress. Reference implementations for other agents.

### scCODA (Bayesian Compositional) — COMPLETED ✅
- Output: `results/sccoda/`
- L1: 14/21 cell types credible; L2: only CD14+ Mono → g3 effect is expression, not composition
- **Implementation**: Python (`pertpy`)
```python
import pertpy as pt
import anndata as ad
# Load composition table (exported from R: patient × cell_type counts)
# → sccoda.load() → run_nuts() → credible_effects()
```

### cNMF (Gene Programs) — IN PROGRESS 🔄
- Output: `results/cnmf/{CD4plus_T_Naive_Memory,Inflammatory_Monocyte,NK_cell}/`
- K selection done for CD4 T + Inflam Mono; NK factorizing
- **Implementation**: Python (`cnmf`)
```python
from cnmf import cNMF
# Step 1: Convert Seurat subset → h5ad
# Step 2: cnmf.prepare(components=np.arange(5,31,5), n_iter=200, num_highvar_genes=2000)
# Step 3: cnmf.factorize()
# Step 4: cnmf.combine()
# Step 5: cnmf.k_selection_plot()  # visual inspect
# Step 6: cnmf.consensus(k=K, density_threshold=0.1)
```

### MELD (Condition Likelihood) — PLANNED
- **Purpose**: Per-cell condition density estimation using KNN graph
- **Implementation**: Python (`meld`)
```python
import meld
import scanpy as sc
adata = sc.read_h5ad("stroke_hc_v8_2_scvi.h5ad")
sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=15)
meld_op = meld.MELD()
meld_op.fit(adata.obsp['connectivities'])
# L1: meld_densities = meld_op.transform(sample_labels)  # → P(Stroke) per cell
# L2: filter IS only, transform(g3_labels) → P(Bad) per cell
```

### Augur (Cell Type Prioritization) — PLANNED
- **Purpose**: Per-cell-type AUC for condition separability (expression-based)
- **Implementation**: R (`Augur`)
```r
# install: remotes::install_github("neurorestore/Augur")
library(Augur)
sobj <- qread("5_strokev8_clean.qs")
# L1:
augur_l1 <- calculate_auc(sobj, cell_type_col="anno1", label_col="cohort",
                           n_subsamples=50, subsample_size=20, min_cells=50)
# L2 (IS only):
sobj_is <- subset(sobj, index_injury_hand == "IS" & !is.na(g3))
augur_l2 <- calculate_auc(sobj_is, cell_type_col="anno1", label_col="g3",
                           n_subsamples=50, subsample_size=20, min_cells=50)
plot_lollipop(augur_l1)
```

### MOFA+ (Patient-level Latent Factors) — PLANNED
- **Purpose**: Multi-view factor analysis from pseudobulk per cell type
- **Implementation**: R (`MOFA2`) + Python backend (`mofapy2`)
```r
library(MOFA2)
# Build pseudobulk: patient × gene per cell type → create_mofa()
# → prepare_mofa(num_factors=15) → run_mofa()
# → plot_variance_explained(), plot_factor(color_by="cohort")
```

### propeller (Compositional DA) — PLANNED
- **Purpose**: Logit-transformed proportions + limma for DA testing
- **Implementation**: R (`speckle`)
```r
library(speckle)
results <- propeller(clusters=sobj$anno1, sample=sobj$patient_name, group=sobj$cohort)
```

### LIANA (Multi-method CCI Consensus) — PLANNED (CCI alternative)
- **Purpose**: Replace/complement CellChat+MNN with multi-method consensus CCI
- **Implementation**: R (`liana`)
```r
library(liana)
liana_res <- liana_wrap(sce, method=c("natmi","connectome","sca","cellchat"))
liana_agg <- liana_aggregate(liana_res)
# Condition-specific with LIANA+
```

---

## Documentation

### Project-specific (Claude context)
| Doc | Path | Content |
|-----|------|---------|
| **Progress status** | `docs/claude/7_PROGRESS_STATUS.md` | 전체 진행 현황, 시각화 결과물 위치, 핵심 발견 |
| **Project goals** | `docs/claude/4_PROJECT_GOALS.md` | 논문 구성 분석 11개 항목, 의존성, 우선순위 |
| Analysis context | `docs/claude/1_STROKE_HC_V8_2_CONTEXT.md` | Full analysis state, data locations, results |
| Analysis report | `docs/claude/5_ANALYSIS_REPORT.md` | MASC, DEG, frequency 분석 결과 |
| Additional plans | `docs/claude/6_ADDITIONAL_ANALYSES_PLAN.md` | 추가 분석 실행 계획 |
| Gene dynamics methods | `docs/claude/2_gene_dynamics_approaches.md` | GAMM model design, two-stage pseudobulk plan |

### CCI-specific documentation
| Doc | Path | Content |
|-----|------|---------|
| **CCI results** | `docs/claude/cci/CCI_ANALYSIS_REPORT.md` | CellChat/MNN 결과 해석, 핵심 발견 |
| **CCI audit** | `docs/claude/cci/CCI_AUDIT_AND_REIMPLEMENTATION.md` | 감사 결과, 문제점, 재구현 계획 |

### HC-only analysis context
| Doc | Path | Content |
|-----|------|---------|
| **HC master context** | `docs/claude/HC/HC_ANALYSIS_CONTEXT.md` | 완료된 분석, 데이터 위치, 핵심 결과 |
| **HC novel plan** | `docs/claude/HC/HC_NOVEL_ANALYSIS_PLAN.md` | Literature-guided 신규 분석 계획 |
| FGS age results | `docs/claude/HC/FGS_AGE_RESULTS.md` | FGS continuous 결과 (6 methods, compartment별) |
| Treg exhaustion | `docs/claude/HC/TREG_EXHAUSTION_ANALYSIS.md` | Treg exhaustion marker 분석 (null result) |

### Shared scripts (HC analysis)
| Script | Location | Purpose |
|--------|----------|---------|
| FGS continuous | `scripts/hc/run_fgs_continuous.R` | 6-method FGS for continuous target |
| FGS by compartment | `scripts/hc/run_fgs_by_compartment.R` | Per-compartment FGS pipeline |
| Treg deep analysis | `scripts/hc/run_treg_exhaustion_deep.R` | Subclustering + age/sex stratification |
| MELD | `scripts/hc/run_meld.py` | Continuous perturbation scoring |
| scCODA | `scripts/hc/run_sccoda.py` | Compositional DA |
| MASC | `scripts/hc/run_masc_hc.R` | Mixed-effects abundance testing |
| MILO age | `scripts/hc/run_milo_age.R` | Neighbourhood DA for age |

### Module guides (workflow documentation)
| Doc | Path | Content |
|-----|------|---------|
| CCI guide | `docs/cci/CCI_INTEGRATED_GUIDE.md` | CellChat + MNN workflow |
| FGS guide | `docs/fgs/FGS_TML_INTEGRATED_GUIDE.md` | Feature selection pipeline |
| Pseudotime guide | `docs/pseudotime-dev/PSEUDOTIME_INTEGRATED_GUIDE.md` | Trajectory analysis |
| Docs rules | `docs/DOCS_ORGANIZATION_RULE.md` | Documentation structure conventions |

### Scripts & execution
| Doc | Path | Content |
|-----|------|---------|
| **Scripts index** | `docs/scripts/ANALYSIS_SCRIPTS_INDEX.md` | 모든 스크립트 위치, CLI 인자, 실행 방법, 출력 구조 |
