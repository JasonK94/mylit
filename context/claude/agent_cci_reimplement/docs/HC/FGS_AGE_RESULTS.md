# FGS Continuous (Age) — Results Summary

> Created: 2026-02-18 | Dataset: hc_only_v1 (96 HC patients, 549K cells)

## 1. Overview

Feature Gene Signature discovery for **continuous age** in healthy PBMCs. Unlike the binary FGS pipeline used in stroke_hc_v8_2, this pipeline uses regression-based methods appropriate for a continuous target variable.

### Why Continuous FGS?

The standard FGS pipeline (`run_fgs_pipeline.R`) uses 10 methods designed for binary classification (Wilcoxon, DESeq2, logistic regression, etc.). For age as a continuous variable, we adapted to:
- **Correlation-based**: Spearman rank correlation
- **Penalized regression**: LASSO (alpha=1), Elastic Net (alpha=0.5)
- **Ensemble**: Random Forest (permutation importance)
- **Linear modelling**: limma-voom (pseudobulk, `~ age + sex`)
- **Nonlinear**: GAM (mgcv, `expr ~ s(age) + sex`)

Other candidate methods not yet implemented: mutual information, distance correlation, Boruta (RF wrapper), PLS regression.

---

## 2. Treg × Age (Primary Analysis)

### 2.1 v1: 4 Methods, Top-50

| Method | Key Metric | Top Gene | Genes Selected |
|---|---|---|---|
| Spearman | 0 sig (padj<0.05) | PLAC8 (rho=−0.443) | 50 |
| LASSO | R²=0.815*, 40 non-zero | JAML (coef=−2.79) | 50 |
| Random Forest | OOB R²=0.033 | PSMD6 (imp=0.554) | 50 |
| limma-voom | 2 sig (adj.P<0.05) | PLAC8 (logFC=−0.025) | 50 |

*LASSO R² inflated by overfitting (40 coefs for 96 samples)

**Consensus**: 45 genes in top-50 of ≥2 methods (expected by chance: 0.82 → 55x enrichment, p=5.5e-61)

### 2.2 v2: 6 Methods, Top-200

| Method | Sig (padj<0.05) | Top Gene | Non-zero / R² |
|---|---|---|---|
| Spearman | 0 | PLAC8 (rho=−0.443) | — |
| LASSO | — | PLAC8 | 10 / R²=0.539 |
| RF | — | PSMD6 | OOB R²=0.033 |
| limma-voom | 2 | PLAC8 | — |
| **GAM** | 2 | PLAC8 (dev.expl=0.233) | — |
| **Elastic Net** | — | PLAC8 | 4 / R²=0.437 |

#### v2 Consensus

| Threshold | Genes | Fold Enrichment | p-value |
|---|---|---|---|
| ≥2 methods | 362 | 11.3x | 3.1e-240 |
| ≥3 methods | **84** | 177.9x | 7.9e-155 |
| ≥4 methods | **32** | 8,157x | 3.8e-113 |
| ≥5 methods | **7** | 402,127x | 9.6e-38 |
| All 6 | **4** | 1.24×10⁸x | 4.5e-32 |

#### Top Consensus Genes (v2, 4+ methods)

| Rank | Gene | #Methods | Mean Rank | Direction | Notes |
|---|---|---|---|---|---|
| 1 | **PLAC8** | 6 | 1 | −age | ISG, immune regulation |
| 2 | **EPSTI1** | 6 | 2 | −age | IFN-stimulated |
| 3 | **ENSG00000274460** | 6 | 3 | +age | Unannotated locus |
| 4 | **JAML** | 6 | 4 | −age | Junctional adhesion molecule |
| 5 | ENSG00000284624 | 5 | 5 | +age | Unannotated |
| 6 | ISG15 | 5 | 6 | −age | Ubiquitin-like ISG |
| 7 | TRAPPC8 | 5 | 7 | −age | Vesicle trafficking |
| 8 | AGMAT | 4 | 8 | mixed | Agmatinase |
| 9 | GPA33 | 4 | 9 | −age | Glycoprotein |
| 10 | C2orf42 | 4 | 10 | −age | Uncharacterized |
| 11 | TMEM115 | 4 | 11 | −age | Transmembrane |
| 12 | ABCF2 | 4 | 12 | −age | ABC transporter |
| 13 | TIGD2 | 4 | 13 | −age | Transposase-derived |
| 14 | IL15 | 4 | 14 | +age | T/NK homeostatic cytokine |
| 15 | IRF2 | 4 | 15 | −age | IFN regulatory factor |
| 16 | LIN7C | 4 | 16 | −age | PDZ scaffold |
| 17 | WDR55 | 4 | 17 | −age | WD-repeat protein |
| 18 | MSC | 4 | 18 | **+age** | Musculin (Treg/NK senescence?) |
| 19 | ENSG00000288043 | 4 | 19 | +age | Unannotated |
| 20 | DGKA | 4 | 20 | −age | DAG kinase → T cell anergy |

#### Pairwise Method Overlap (Top-200)

| | Spearman | LASSO | RF | limma | GAM | ENet |
|---|---|---|---|---|---|---|
| Spearman | 200 | 11 | 46 | 74 | **107** | 7 |
| LASSO | 11 | 200 | 9 | 12 | 13 | **194** |
| RF | 46 | 9 | 200 | 43 | 41 | 6 |
| limma | 74 | 12 | 43 | 200 | **97** | 7 |
| GAM | **107** | 13 | 41 | **97** | 200 | 7 |
| ENet | 7 | **194** | 6 | 7 | 7 | 200 |

- **Spearman-GAM**: highest overlap (107) — both test marginal associations
- **LASSO-ENet**: nearly identical (194) — same glmnet framework
- **RF is most orthogonal**: only 6–46 overlaps with other methods

---

## 3. Compartment × Age

### Signal Strength by Compartment

| Compartment | Spearman sig | limma sig | Consensus (2+) | Top Gene |
|---|---|---|---|---|
| **T cell** | **898** | **935** | 42 | SOX4 (−0.68) |
| **Myeloid** | 117 | 60 | 43 | SPI1 (+0.51) |
| **NK** | 111 | 70 | 36 | CDKN2A (+0.61) |
| DC | 59 | 11 | 21 | SEMA5A (−0.49) |
| B cell | 9 | 20 | 31 | RHEX (−0.46) |
| Platelet | 0 | 0 | 40 | RHOC (+0.46) |

### Biological Themes per Compartment

**T cell** — Naive/developmental decline:
- SOX4 (−0.68): developmental TF, critical for T cell differentiation
- BLK (−0.62): B lymphoid kinase (ectopic in T?)
- ROBO1 (−0.62): Slit-Robo pathway, guidance/migration
- CASC15 (−0.63): lncRNA, developmental regulator

**NK** — Senescence/maturation:
- CDKN2A (+0.61): p16INK4a, canonical senescence marker
- IGFBP3 (+0.51): growth factor regulation, senescence-associated
- LAG3 (+0.46): exhaustion/senescence checkpoint
- MSC (+0.55): musculin, also +age in T cells

**Myeloid** — Activation/inflammation:
- SPI1 (+0.51): PU.1, master myeloid TF → more "committed" myeloid with age
- EMP3 (+0.49): epithelial membrane protein, activation marker
- APH1B (+0.48): γ-secretase subunit, Notch/APP pathway

**B cell** — Signaling decline:
- RHEX (−0.46): erythroid expansion factor
- VAV1 (−0.44): Rho GEF, BCR signaling
- TNFRSF13B (−0.43): TACI, B cell survival receptor

### Pan-Compartment Gene

**ENSG00000237471**: consensus in NK, Platelet, T cell (mean rho=+0.47, increases with age). Unannotated lincRNA — warrants further investigation.

---

## 4. Output Locations

```
fgs_continuous_treg/                    # v1 (4 methods, top-50)
├── pseudobulk_data.rds                 # Reusable for downstream
├── consensus_genes_2plus_methods.csv   # 45 genes
├── consensus_all_genes.csv             # Full ranking (18,182 genes)
├── fgs_continuous_results.rds          # Per-method detailed results
├── heatmap_consensus_genes_vs_age.png
└── scatter_top9_genes_vs_age.png

fgs_continuous_treg_v2/                 # v2 (6 methods, top-200)
├── consensus_genes_2plus_methods.csv   # 362 genes
├── consensus_genes_3plus_methods.csv   # 84 genes
├── consensus_genes_4plus_methods.csv   # 32 genes
├── pairwise_overlap_top200.csv
├── heatmap_consensus_top60.png
└── scatter_top12_consensus.png

fgs_compartment/                        # By compartment
├── fgs_compartment_results.rds         # All compartment results
├── fgs_top_genes_across_compartments.csv
├── heatmap_fgs_top_genes_across_compartments.png
└── {B_cell,DC,HSPC,...,T_cell}/
    ├── consensus_genes.csv
    ├── spearman_full.csv
    └── limma_full.csv
```
