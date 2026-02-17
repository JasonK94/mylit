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

## Conventions

- **Plots**: Heavy scatter plots (DimPlot, trajectory) → PNG. Simple/vector → PDF.
- **Reductions**: Always use `integrated.scvi` / `umap.scvi` from the original object.
  Never recompute PCA/UMAP within subsets — HC/Stroke batch effects corrupt de novo embeddings.
- **Dual-layer design**: Layer 1 = HC vs Stroke (cohort), Layer 2 = IS g3==1 vs g3==2.
- **Commit style**: Short imperative subject + optional body with details.
- **Language**: User communicates in Korean, documentation is bilingual (EN primary).

## Current Analysis State

See `docs/claude/STROKE_HC_V8_2_CONTEXT.md` for full analysis state.

### Completed
- Pipeline (CellBender/Souporcell/Solo/scVI) + Clustering + Annotation + QC
- MASC (6 comparisons), FindAllMarkers
- CCI: CellChat + MNN (dual-layer) + interpretation
- Trajectory v2: Slingshot + Monocle3, scVI UMAP, 3 compartments, HC vs Stroke

### In Progress / Pending
- FGS/TML: n=50 sweep (method 5/10), then n=100, n=200
- Gene dynamics v3 (batch-corrected GAMM): script ready, pending v2 completion
- MNN at anno1 level (26 types) for proper LR prioritisation
- Subclustering, MILO, DEG Consensus

## Documentation

### Project-specific (Claude context)
| Doc | Path | Content |
|-----|------|---------|
| **Project goals** | `docs/claude/PROJECT_GOALS.md` | 논문 구성 분석 11개 항목, 의존성, 우선순위 |
| Analysis context | `docs/claude/STROKE_HC_V8_2_CONTEXT.md` | Full analysis state, data locations, results |
| Gene dynamics methods | `docs/claude/gene_dynamics_approaches.md` | GAMM model design, two-stage pseudobulk plan |

### Scripts & execution
| Doc | Path | Content |
|-----|------|---------|
| **Scripts index** | `docs/scripts/ANALYSIS_SCRIPTS_INDEX.md` | 모든 스크립트 위치, CLI 인자, 실행 방법, 출력 구조 |

### Module guides (workflow documentation)
| Doc | Path | Content |
|-----|------|---------|
| CCI guide | `docs/cci/CCI_INTEGRATED_GUIDE.md` | CellChat + MNN workflow |
| FGS guide | `docs/fgs/FGS_TML_INTEGRATED_GUIDE.md` | Feature selection pipeline |
| Pseudotime guide | `docs/pseudotime-dev/PSEUDOTIME_INTEGRATED_GUIDE.md` | Trajectory analysis |
| Docs rules | `docs/DOCS_ORGANIZATION_RULE.md` | Documentation structure conventions |
