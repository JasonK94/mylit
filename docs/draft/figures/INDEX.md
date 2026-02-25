# Figure Index — stroke_hc_v8_2 Paper Draft

> Created: 2026-02-25
> Notation: C1 = coarse annotation (anno2: Tc, Mono, etc.), C2 = fine annotation (anno1: CD4+Trm, etc.)
>           L1 = HC vs Stroke (cohort), L2 = Good vs Bad (g3, within IS)

---

## Overview

| Folder | Topic | Files | Narrative Role |
|--------|-------|:-----:|----------------|
| `fig1_integration_annotation/` | Integration QC, UMAP, DotPlot | 15 | Data quality & cell type landscape |
| `fig2_celltype_frequency/` | Composition bar, sample-level distribution | 12 | Frequency overview (sample & group) |
| `fig3_differential_abundance/` | MASC, scCODA, MILO, MELD, Augur | 24 | Statistical evidence for DA (C1/C2, L1/L2) |
| `fig4_deg_signature/` | Volcano, FGS, cNMF, MOFA+, validation | 30 | Expression differences & unbiased discovery |
| `fig5_cci/` | CellChat, MNN, LIANA, integrated | 25 | Cell-cell communication (C1/C2, L1/L2) |
| `fig6_trajectory/` | Pseudotime UMAP, density, gene dynamics | 19 | Trajectory analysis & prognosis genes |
| `fig7_subclustering_programs/` | Mono/CD4/CD8 SC, cNMF, frequency | 30 | Deep-dive: pathogenic subpopulations |
| **Total** | | **155** | |

---

## Fig 1: Integration & Annotation

**Purpose**: Data quality, batch mixing, cell type annotation의 무결성 제시

| File | Content | Key Message |
|------|---------|-------------|
| `01_umap_anno1_C2.png` | UMAP colored by anno1 (21 types) | C2-level cell type landscape |
| `02_umap_anno2_C1.png` | UMAP colored by anno2 (8 compartments) | C1-level overview |
| `03_umap_anno1_split_cohort.png` | UMAP split HC / Stroke | Batch mixing across conditions |
| `04_umap_batch_GEM.png` | UMAP colored by GEM (batch) | GEM mixing after scVI integration |
| `05_dotplot_canonical_markers.png` | DotPlot: mutually exclusive markers per C2 type | Annotation validity |
| `06-09_scanvi_*` | scANVI: UMAP, scVI comparison, confusion matrix, silhouette | Independent label validation (97.3%) |
| `10_v3_umap_anno1.png` | High-res UMAP (v3 paper version) | Publication-ready |
| `11_v3_dotplot_markers.png` | v3 DotPlot | Publication-ready |

**Discussion points**:
- C1 (anno2) vs C2 (anno1): 논문에서 먼저 C1으로 overview를 보여주고, C2로 detail을 보여줄지, 아니면 C2만 보여줄지 결정 필요
- DotPlot: `features = list(Tc = c("CD3D","CD3E"), Mono = c("CD14","LYZ"), ...)` 형태의 mutually exclusive marker 확인
- GEM mixing: scVI integration 후 GEM 간 잘 섞여 있음을 보여주는 것이 필수

---

## Fig 2: Cell Type Frequency

**Purpose**: Sample-level, group-level composition 차이를 시각적으로 보여줌

| File | Content | Key Message |
|------|---------|-------------|
| `01_stacked_bar_C1_by_patient.png` | Cumulative bar (C1, all patients) | Patient-level composition overview |
| `02_stacked_bar_C1_by_patient_cohort.png` | Same, split by cohort | HC vs Stroke pattern |
| `03_stacked_bar_C2_by_patient_cohort.png` | Same at C2 level | Fine-grained composition |
| `04_stacked_bar_fgs_score_sorted.png` | Sorted by FGS signature score | FGS score ↔ composition relationship |
| `05_stacked_bar_C1_fgs_sorted.png` | C1-level, FGS-sorted | Mono ↑ = FGS score ↑ trend |
| `06-07_freq_boxplot_C1_L1/L2.png` | Boxplots C1, L1 and L2 | Group-level frequency comparison |
| `08-09_freq_boxplot_C2_L1/L2.png` | Boxplots C2, L1 and L2 | Fine-grained frequency |
| `10_v3_composition_bar.png` | v3 paper composition bar | Publication-ready |
| `11-12_sccoda_*` | scCODA composition & heatmap | Bayesian compositional analysis |

**Key narrative**:
- FGS score 순 정렬 시 Monocyte ↑ / T cell ↓ trend가 보임
- Sample-level variation이 있지만 overall trend 일관

---

## Fig 3: Differential Abundance

**Purpose**: 5-method concordance로 어떤 cell type이 DA인지 통계적으로 보여줌

| File | Content | Key Message |
|------|---------|-------------|
| `01_3method_concordance_heatmap.png` | MASC×scCODA×MILO concordance | 14/19 types DA by ≥2 methods (L1) |
| `02-03_v3_concordance_L1/L2.png` | v3 concordance tiles | L1 robust, L2 mono-only |
| `04_masc_forest.png` | MASC forest plot (OR) | Effect size visualization |
| `05-09_milo_*` | MILO beeswarm, UMAP, summary (L1/L2) | Neighbourhood-level DA |
| `10-11_milo_concordance/scatter` | MILO×MASC×scCODA cross-method | Method agreement |
| `12-14_sccoda_*` | scCODA credible effects, volcano, boxplot | Bayesian DA |
| `15-18_meld_*` | MELD UMAP & violin (L1/L2) | Continuous per-cell probability |
| `19-21_augur_*` | Augur lollipop (L1/L2), L1 vs L2 scatter | Expression-based separability |
| `22-23_v3_multimethod/augur` | v3 paper panels | Publication-ready |
| `24_milo_dual_layer_dotplot.png` | MILO L1+L2 dual-layer dotplot | L1 vs L2 comparison |

**Key narrative**:
- L1: 14/19 cell types DA (robust, multi-method)
- L2: CD14+ Mono = sole DA type → prognosis = monocyte-specific
- Augur L2 near random (AUC ~0.5) → g3 effect is subtle

---

## Fig 4: DEG & Signature

**Purpose**: Expression-level 차이, unbiased signature discovery, external validation

| File | Content | Key Message |
|------|---------|-------------|
| `01_volcano_4celltypes.png` | Volcano (CD14+ Mono, Inflam Mono, CD4 T, NK) | Key cell type DEG patterns |
| `02-03_deg_count*` | DEG count per cluster | Massive asymmetry (99% down) |
| `04-05_deg_sharing/concordance` | Method sharing, edgeR-NEBULA concordance | Multi-method robustness |
| `06-08_cross_layer_*` | L1∩L2 overlap, upset, scatter | 12 triple-overlap genes |
| `09-12_fgs_*` | FGS method AUC, top genes, heatmap, celltype | Unbiased feature selection |
| `13_v3_fgs_top30.png` | Top 30 FGS genes | Key signature genes |
| `14-17_fgs_score_*` | FGS score distribution (C1, L1/L2, patient-level) | Score differences by condition |
| `18_fgs_top30_heatmap.png` | Patient × gene heatmap | Individual-level patterns |
| `19_fgs_genes_in_deg.png` | FGS∩DEG overlap heatmap | 218/241 FGS genes in DEG |
| `20-23_external_*` | AUC, boxplot, ROC, forest | AUC 0.82-0.89 validation |
| `24-27_cnmf_*` | cNMF hallmark, condition, overview | Gene program discovery |
| `28-30_mofa_*` | MOFA+ variance, cohort, weights | Patient-level factor analysis |

**Key narrative**:
- DEG: massive transcriptional suppression (99% downregulated)
- FGS: g3-derived signature validated externally (AUC 0.89)
- cNMF: TNFa/NFkB = dominant monocyte program
- MOFA+: Factor 1 = cohort (26.5% Mono variance, S100A8/A9 drivers)

---

## Fig 5: CCI

**Purpose**: Cell-cell communication changes (L1: SIIS pattern, L2: maladaptive activation)

| File | Content | Key Message |
|------|---------|-------------|
| `01-06_cellchat_L1_*` | CellChat L1: interaction count, diff, rankNet, bubble, scatter | CCI -37%, selective pathway activation |
| `07-09_cellchat_L2_*` | CellChat L2: diff, rankNet, bubble | All pathways ↑ in Bad; IL1 Bad-only |
| `10-13_mnn_L1_*` | MNN: DE summary, volcano, circos, heatmap | Ligand-receptor level analysis |
| `14-15_liana_L1/L2_scatter` | LIANA rank comparison | 5-method consensus interactions |
| `16-20_cci_integrated_*` | 3-method concordance: UpSet, pathway, sender-receiver | MIF (1.54), GALECTIN (1.39) top |
| `21-22_cci_deg_*` | CCI×DEG crossref | 22.5% overlap (TGFB1, CD44, CXCR4) |
| `23-25_v3_*` | v3 paper panels | Publication-ready |

**Key narrative**:
- L1: Global CCI collapse + selective RESISTIN/MIF/CypA activation = SIIS
- L2: Bad = maladaptive broad hyperactivation (ALL pathways ↑)
- 3-method concordance: MIF, GALECTIN, RESISTIN, ANNEXIN

---

## Fig 6: Trajectory

**Purpose**: Pseudotime dynamics — monocyte = only g3-responsive compartment

| File | Content | Key Message |
|------|---------|-------------|
| `01-03_pseudotime_umap_*` | Pseudotime UMAP (Mono, CD4, CD8) | Trajectory structure |
| `04-06_pseudotime_*_cohort` | Split by cohort | HC vs Stroke distribution |
| `07-09_density_*_L1` | Density plots, cohort comparison | Pseudotime distribution shift |
| `10-12_density_*_L2` | Density plots, g3 comparison | Mono p<2e-16, CD8 NS |
| `13-19_v3_*` | v3 panels: UMAP, density, gene dynamics, effect sizes | Publication-ready |
| `20_genedyn_mono_cohort_*.png` | Individual gene dynamics (S100A8/A9, VCAN, IL1B, HIF1A, TNF) | Key prognosis genes |
| `21_genedyn_mono_g3_*.png` | Same, g3 comparison | Outcome-associated dynamics |

**Key narrative**:
- Monocyte = ONLY compartment with g3 trajectory effect (p<2e-16)
- 10 genes confirmed by 3 methods (ABC, Lamian, pseudobulk GAMM)
- CD8 T: no g3 effect at all (p=0.749)

---

## Fig 7: Subclustering & Gene Programs

**Purpose**: Deep-dive into pathogenic subpopulations (C2 → subclusters)

| File | Content | Key Message |
|------|---------|-------------|
| `01-07_mono_*` | Mono SC: UMAP, identity, markers, frequency L1/L2 | SC2 (OR=403), SC3 (OR=193) Bad-enriched |
| `08-12_cd4_*` | CD4 T SC: UMAP, identity, frequency | SC3 (OR=18.6) Bad-enriched |
| `13-17_cd8_*` | CD8 T SC: UMAP, identity, frequency | SC8 (OR=20.6) Bad-enriched |
| `18_combined_volcano.png` | 3-compartment frequency volcano | Overview of all SC changes |
| `19-23_*_boxplot/stacked_bar` | Per-compartment frequency details | Individual SC patterns |
| `24-26_cnmf_*` | cNMF top genes, named programs | TNFa/NFkB = pathogenic program |
| `27-30_v3_*` | v3 paper panels | Publication-ready |

**Key narrative**:
- Mono SC2/SC3: extreme Bad-enrichment (OR > 100)
- cNMF GEP4/5: TNFa/NFkB signaling (NES = 3.81-3.87)
- L2 subclustering frequency: 0/45 significant → composition은 안 다르지만, specific SC는 OR > 100

---

## Narrative Flow Summary

```
Fig 1: "Here is our data" — 226K cells, 21 types, well-integrated, well-annotated
    ↓
Fig 2: "Here is what the composition looks like" — sample-level bars, Mono↑/Tc↓ trend
    ↓
Fig 3: "Which cell types are statistically different?" — 5-method DA concordance
    → L1: 14/19 DA, L2: CD14+ Mono only
    ↓
Fig 4: "What genes/programs are different?" — DEG (99% down), FGS (AUC 0.89), cNMF (TNFa/NFkB)
    ↓
Fig 5: "How does cell communication change?" — SIIS (CCI -37% + selective activation)
    → L2: Bad = all ↑ (maladaptive hyperactivation)
    ↓
Fig 6: "How does differentiation trajectory change?" — Mono = only g3-responsive
    → 10 genes × 3 methods
    ↓
Fig 7: "Which specific subpopulations drive this?" — Mono SC2/SC3 (OR > 100)
    → TNFa/NFkB program = pathogenic mechanism
```
