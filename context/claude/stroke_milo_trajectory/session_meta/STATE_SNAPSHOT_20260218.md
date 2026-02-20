# State Snapshot — 2026-02-18 14:17 KST

Quick-start document for a new agent resuming this project.
Read this + MEMORY.md + CLAUDE.md to reconstruct full context.

## Datasets

### stroke_hc_v8_2 (primary)
- **226,855 cells → 205,277 clean** / 100 patients (HC 5GEM + Stroke 12GEM)
- **Seurat objects**:
  - Full: `/data/user3/sobj/stroke_hc_v8_2/5_strokev8_clean.qs`
  - L1 (HC vs IS): `5_1_hc_is.qs` (131K: HC 73K + IS 58K)
  - L2 (g3): `5_2_is_g3.qs` (54K: g3=1 15K + g3=2 39K)
- **Results**: `/data/user3/sobj/stroke_hc_v8_2/` (symlinked as `results/`)
- **21 anno1 types, 8 anno2 compartments**

### hc_only_v1 (HC-only atlas)
- **549,340 cells** / 96 HC patients / 24 GEMs
- Seurat: `/data/user3/sobj/hc_only_v1/5_hc_clean.qs`
- 32 anno1, 10 anno2

## Completed Analyses Summary

### Frequency DA (3-method concordance)
- **L1**: 14/19 cell types DA by ≥2 methods (MASC × scCODA × MILO)
  - HC-enriched: ISG+ T_cell, Platelet/PLA, CD16+ Mono, pDC, CD8+ Cyto, MAIT, NK, cDC2, Treg, CD8+ Trm, Infl. Mono
  - IS-enriched: ISG+ Myeloid, B_cell, CD14+ Mono
- **L2**: Only CD14+ Monocyte (MILO + scCODA concordant)
- Output: `milo/{L1_cohort,L2_g3}/`, `MASC/`, `sccoda/`

### DEG Consensus
- L1: muscat-edgeR + muscat-DESeq2 + NEBULA (15 clusters)
- L2: NEBULA only (203K results)
- Cross-layer: 134 genes sig in both L1 and L2
- Output: `deg_consensus/`

### CCI
- CellChat v2: L1+L2 working (rankNet, bubble, scatter)
- MNN: ABANDONED (0 prioritized interactions)
- Output: `cci/`

### Trajectory + Gene Dynamics
- Slingshot + Monocle3, 3 compartments (mono, cd4, cd8) × 2 conditions
- Gene dynamics v3 (NB-GAMM) + 3 validation methods
- **Key finding**: 10 genes confirmed 3/3 methods (mono/cohort)
- Mono = ONLY compartment with g3 effect; CD8 cell-level = pseudoreplication
- Output: `trajectory_v2/`, `trajectory_v3/`

### FGS
- n=50, n=100 complete; n=200 running
- Pathway enrichment: HALLMARK, KEGG, GOBP done
- FGS × DEG overlap: 218/241 genes significant in L1 DEG
- Output: `fgs/`

## Running Background Processes (at backup time)
1. FGS n=200 (whole IS) — PID 2711001
2. FGS within-celltype CD14 Mono n=50 — PID 3634494
3. cNMF v2 — PID 3544661

## Key Scripts
| Script | Location |
|--------|----------|
| MILO | `/data/user3/sobj/stroke_hc_v8_2/scripts/run_milo.R` |
| FGS pipeline | `Git_Repo/_wt/fgs/scripts/fgs/run_fgs_pipeline.R` |
| CellChat v2 | `/data/user3/sobj/stroke_hc_v8_2/scripts/reimplement_cellchat_v2.R` |
| DEG Consensus | `Git_Repo/_wt/deg-consensus/scripts/consensus/run_deg_consensus_cli.R` |
| Trajectory effect sizes | `stroke/scripts/compute_trajectory_effect_sizes.R` |
| Lamian validation | `stroke/scripts/run_lamian_validation.R` |
| Pseudobulk GAMM | `stroke/scripts/run_pseudobulk_gamm.R` |
| Gene dynamics v3 | `stroke/scripts/run_gene_dynamics_v3.R` |

## Remaining Work (priority order)
1. Within-celltype FGS (CD14 Mono running, 14 more)
2. cNMF (NK factorizing)
3. Descriptive figures (signature-based plots for paper)
4. External validation
5. Optional: MELD, Augur, MOFA+, scANVI

## Critical Known Issues
- `buildNhoodGraph(milo)` required before `plotNhoodGraphDA()`
- MASC cluster names use dots: `CD14..Monocyte` → need name mapping
- `reducedDim(milo, "PCA") <- reducedDim(milo, "SCVI")` before testNhoods
- MNN fundamentally incompatible (lr_target_prior_cor needs ≥5 samples)
- GEM × cohort confounded in L1 → remove GEM from L1 DEG covariates
- Lamian: use patient_name (not GEM) as sample; chisq >> permutation speed
- CDS `exprs()` = raw counts → must log-normalize before Lamian
