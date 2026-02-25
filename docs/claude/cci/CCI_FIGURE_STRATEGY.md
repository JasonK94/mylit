# CCI Figure Strategy — CellChat + MultiNicheNet

## Overview

Two complementary CCI methods are used to analyze cell-cell interactions:

| Method | Resolution | Directionality | Key Strength |
|--------|-----------|----------------|--------------|
| **CellChat** | Compartment (anno2, 8 types) | Bidirectional (HC and Stroke computed independently) | Pathway-level decomposition, curated L-R database |
| **MultiNicheNet** | Fine cell-type (anno1, 21 types) | Condition-stratified via logFC direction | Integrated prioritization (DE + ligand activity + L-R) |

### MNN Bidirectionality — How It Works

MNN's `group_prioritization_tbl` is **uni-directional** (only Stroke-enriched interactions for L1).
However, **condition-stratified** visualization is achieved by:

1. Using `ligand_activities_target_de_tbl` (NOT `group_prioritization_tbl`)
2. Joining with `celltype_de` to get sender info via differential expression logFC
3. Splitting by `direction = ifelse(logFC > 0, "Stroke", "HC")`
4. Filtering for **concordant** sender-receiver pairs (both ligand DE and activity in same direction)

Scripts: `plot_comparison_circos.R` (basic), `plot_mnn_circos_balanced.R` (balanced, preferred)

---

## Data Locations

### CellChat
- Objects: `cci/plots/cellchat_{L1_cohort,L2_g3}_v2/{cellchat_HC,Stroke,Good,Bad}.qs`
- Comparison: `cellchat_comparison.qs` per layer
- Per-sample: `cci/cellchat/{L1_cohort,L2_g3}_anno2/samples/*/cellchat.qs` (88 total)
- Plots: `cci/plots/cellchat_{L1_cohort,L2_g3}_v2/` (62 files per layer)

### MultiNicheNet (v3)
- Results: `cci/mnn/{L1_cohort,L2_g3}_{anno1,anno2}_v3/multinichenet_results.qs`
- Interaction counts: L1a1=148,958 / L1a2=8,244 / L2a1=91,111 / L2a2=7,600
- Plots: `cci/plots/mnn_v3_{L1,L2}_{anno1,anno2}/` (13 files each)

### Cross-Method
- Concordance: `cci/plots/concordance/` (5 figures + 4 CSV)
- Compartment summary: `cci/plots/compartment_summary/` (5 figures + 3 CSV)
- Pathway deep-dive: `cci/plots/pathway_deepdive/` (5+ figures + 3 CSV)
- Inflammation scoring: `cci/plots/inflammation_score/` (CellChat) + `mnn_inflammation/` (MNN)

---

## Paper Figure Composition

### Main Figures

| Figure | Title | Panels | Source |
|--------|-------|--------|--------|
| **A** | CCI Overview: Who talks to whom | A1: Diff circle (L1), A2: Diff heatmap (L1), A3: Same for L2 | Existing cellchat_v2 |
| **B** | Pathway Information Flow | B1: RankNet L1, B2: RankNet L2, B3: Inflammation dotplot | Existing cellchat_v2 + inflammation_score |
| **C** | Inflammatory Balance Shift | C1: Pro/Anti ratio bars, C2: Pathway waterfall | Existing inflammation_score |
| **D** | MNN Condition-Stratified Interactions | D1: Balanced circos L1, D2: Balanced circos L2 | Existing mnn_balanced |
| **E** | MNN Top Molecular Drivers | E1: Top 30 L-R dotplot (L1), E2: SR network | Existing mnn_v3 |
| **F** | MNN Functional Annotation | F1: Category composition, F2: SR heatmap (inflammatory) | Existing mnn_inflammation |
| **G** | Cross-Method Validation | G1: Ligand overlap, G2: CC vs MNN scatter, G3: SR agreement | NEW: concordance |
| **H** | Compartment-Level Integration | H1: SR heatmap, H2: Difference heatmap, H3: Pathway × compartment | NEW: compartment_summary |
| **I** | MIF Pathway Deep-Dive | I1: CC strength per SR pair, I2: MNN top L-R, I3: Expression dotplot | NEW: pathway_deepdive |

### Supplementary

| Figure | Content | Status |
|--------|---------|--------|
| S1 | All pathway chord diagrams (L1+L2) | ✅ |
| S2 | All MNN top 30 dotplots (4 configs) | ✅ |
| S3 | Signaling role scatter plots | ✅ |
| S4 | MNN score component heatmaps | ✅ |
| S5 | Full inflammation scoring (5+5 figures) | ✅ |
| S6 | Per-sender bubble plots | ✅ |
| S7 | CCI × DEG overlap (16 genes) | ✅ |

---

## Key Findings

### CellChat
- **Total CCI count**: HC > Stroke (diverse homeostatic signaling maintained)
- **Pro/Anti ratio**: HC 2.10 → Stroke 3.50 (inflammatory bias ↑ despite fewer total interactions)
- **L2**: Bad > Good in absolute weight (intense inflammatory interactions)
- **Dominant pathways**: MIF, GALECTIN, CypA (pro); TGFb, ANNEXIN (anti)

### MNN
- **Top interaction**: S100A8→CD36 (Mono→Mono), score=0.978
- **L1 top 500**: 73% pro-inflammatory (S100A8/A9, TNF, ANXA1)
- **L2 top 500**: 24% pro-inflammatory, 15% anti-inflammatory (TGFB1, TIMP1 emerge)
- **38/40 pro-inflammatory ligands shared** between L1 and L2
- Uni-directional group_prioritization_tbl; bidirectional via logFC stratification

### Concordance
- Cross-method: ligand overlap Jaccard ~0.3-0.4
- Sender-receiver agreement stronger at compartment (anno2) than cell-type (anno1) level
- MIF pathway consistently dominant in both methods

---

## Scripts Index

| Script | Purpose | Output Directory |
|--------|---------|-----------------|
| `reimplement_cellchat_v2.R` | CellChat from Seurat (primary) | `cellchat_{L1,L2}_v2/` |
| `plot_cci_comprehensive.R` | CellChat native + MNN combined | `cci/plots/` |
| `plot_pathway_inflammation_score.R` | CellChat inflammation scoring | `inflammation_score/` |
| `plot_mnn_inflammation_score.R` | MNN inflammation scoring | `mnn_inflammation/` |
| `plot_mnn_circos_balanced.R` | MNN balanced circos | `mnn_balanced/` |
| `plot_mnn_v3.R` | MNN standard visualization (4 configs) | `mnn_v3_*/` |
| `plot_cci_concordance.R` | CellChat × MNN cross-validation | `concordance/` |
| `plot_cci_compartment_summary.R` | Integrated compartment view | `compartment_summary/` |
| `plot_cci_pathway_deepdive.R` | MIF/TNF/TGFb mechanistic detail | `pathway_deepdive/` |
