# HC-Only PBMC Atlas: Novel Analysis Plan

> Created: 2026-02-18 | Dataset: hc_only_v1 | 549K cells, 96 Korean HC, age 19-66

## 1. Literature Context & Novelty Claims

### 1.1 Key Prior Art
| Study | Cohort | Cells | Key Contribution |
|-------|--------|-------|------------------|
| Luo et al. 2022 (Nat Aging) | 166 donors, 25-85y | 2M | 12 age-changing subpops, NKG2C+GZMB- CD8 Tmem decline, frailty monocytes |
| Yazar/OneK1K 2022 (Science) | 982 donors | 1.27M | Cell-type eQTLs; Connolly 2024 re-analyzed for aging |
| Sopena-Rios 2024 (bioRxiv) | OneK1K re-analysis | 1.27M | Sex-specific immunosenescence: female CD8 TEM expansion, male CD5+ B |
| AIDA 2025 (Cell) | 619 Asian donors (7 pops) | 1.27M | **Korean-specific: lower Treg proportions**; MAIT~ancestry |
| Filippov 2024 (npj Aging) | 7 datasets integrated | 1M+ | Consensus: CD8 naive decline, MAIT decline, non-classical mono expansion |
| Connolly 2024 (Aging Cell) | 678+982 donors | 4M | Loss of cell identity in aging effector/cytotoxic cells |

### 1.2 What Makes This Cohort Novel
1. **First Korean-specific aging PBMC atlas** — AIDA has Korean donors but cross-sectional, no deep per-donor analysis
2. **First systematic CCI aging atlas in PBMC** from any population
3. **First East Asian age×sex interaction analysis** — Sopena-Rios 2024 only European
4. **Technical**: scVI/MrVI integration (state-of-art), CellBender + Souporcell + Solo QC pipeline
5. **Korean Treg depletion mechanism** — AIDA flagged low Tregs, we can characterize with subclustering + age stratification
6. **Single-hospital Korean cohort**: homogeneous population, controlled for geography

### 1.3 Confirmatory vs Novel Findings
**Confirm** (expected from literature):
- CD8 naive decline, CD8 TEM/effector expansion with age
- NK CD56dim expansion, pDC decline
- Non-classical monocyte (CD16+) increase
- MAIT decline with age
- ISG/IFN signaling decline with age

**Novel** (where we can add value):
- Korean-specific cell type frequencies vs published references
- Age×Sex interaction effects (especially Treg, NK, B cells)
- CCI remodeling with aging (completely unexplored in PBMC)
- MrVI sample-level heterogeneity in aging
- Korean Treg depletion: deeper subclustering, functional characterization

---

## 2. Analysis Plan

### Phase A: Frequency / Differential Abundance (Multi-Method)

#### A1. MASC (Mixed-effects Association of Single Cells) — R
- **Status**: NOT YET RUN for HC
- **Comparisons** (6 total):
  - anno1 × sex, anno1 × age_group, anno1 × age (continuous via quartiles)
  - anno2 × sex, anno2 × age_group
  - anno1 × age×sex interaction
- **Script**: `scripts/hc/run_masc_hc.R`
- **Command**: Uses existing `_wt/masc/scripts/masc/run_masc.R` CLI
- **Covariates**: age + sex (for sex comparison: age only; for age: sex only)
- **Expected runtime**: ~10-15 min per comparison

#### A2. MILO (Neighbourhoods) — R
- **Status**: Sex DONE, Age NOT DONE
- **New runs needed**:
  - MILO age: `~ age` design (continuous), reuse existing nhood object
  - MILO age_group: `~ age_group` (categorical Young/Middle/Old)
  - MILO age×sex interaction: `~ age * sex`
- **Script**: `scripts/hc/run_milo_age.R`
- **Note**: Must add `reducedDim(milo, "PCA") <- reducedDim(milo, "SCVI")` before testNhoods

#### A3. MELD (Manifold-Enhanced Likelihood of Differential abundance) — Python
- **Status**: NOT YET RUN
- **What it does**: Assigns per-cell continuous score for condition enrichment
- **Application**: MELD score for age (binned) or sex, visualize on UMAP
- **Script**: `scripts/hc/run_meld.py`
- **Package**: `meld` (Python, pip install meld)
- **Input**: scVI latent space + age/sex metadata from h5ad
- **Key advantage**: Shows gradients on UMAP, not just cluster-level DA

#### A4. scCODA (Compositional DA) — Python
- **Status**: NOT YET RUN
- **What it does**: Bayesian compositional analysis accounting for compositionality constraint
- **Application**: age_group (Young/Old) and sex as conditions
- **Script**: `scripts/hc/run_sccoda.py`
- **Package**: `sccoda` (part of scverse/pertpy)
- **Reference cell type**: Largest or most stable type (CD14 Classical Mono or CD4 Naive)
- **Key advantage**: Proper compositional analysis (proportions sum to 1)

#### A5. MrVI Differential Expression — Python
- **Status**: MrVI trained, DE NOT YET EXTRACTED
- **What it does**: Sample-aware DE using MrVI's learned sample representations
- **Application**: `mrvi.differential_expression()` for age/sex effects per cell type
- **Script**: `scripts/hc/run_mrvi_de.py`
- **Key advantage**: Sample-level variation modeling → more accurate DE than pseudobulk

### Phase B: CCI Interpretation & Novel Analysis

#### B1. CellChat Age Comparison — R
- **Status**: Objects READY (Young/Old merged), interpretation PENDING
- **Tasks**:
  - Load merged CellChat objects (Young, Old)
  - `netAnalysis_computeCentrality()`, `compareInteractions()`
  - Identify pathways gained/lost with aging
  - Focus on: Treg-myeloid crosstalk, NK-T cell interactions
  - Circos plots, pathway-level heatmaps
- **Script**: `scripts/hc/interpret_cellchat_age.R`

#### B2. CellChat Sex Comparison — R
- **Status**: Objects READY (M/F merged), interpretation PENDING
- **Tasks**: Same as B1 but for sex contrast
- **Focus**: BAFF/APRIL (expected female-enriched), IFN pathways

#### B3. MNN Age/Sex Interpretation — R
- **Status**: Results READY (~690K sex + ~890K age LR activities), interpretation PENDING
- **Tasks**:
  - Extract top differential LR pairs by cell type
  - Focus on Treg sender/receiver interactions
  - Cross-reference with CellChat findings
  - Pathway enrichment of differential LR target genes
- **Script**: `scripts/hc/interpret_mnn_results.R`

### Phase C: Korean-Specific & Novel Discoveries

#### C1. Korean Treg Characterization
- **Status**: Subclustered (7 clusters), basic age analysis done
- **New analyses**:
  - Compare Treg proportion vs AIDA published values
  - FOXP3/IL2RA/CTLA4/TIGIT expression per donor → age/sex correlation
  - Treg:Tconv ratio analysis
  - Treg CCI (specifically Treg → DC, Treg → T effector suppression)

#### C2. Age×Sex Interaction Analysis
- **Rationale**: Sopena-Rios 2024 found female-specific CD8 TEM expansion and male-specific CD5+ B expansion in Europeans. Replicate in Korean cohort.
- **Methods**:
  - Interaction term in MASC: `~ age_group * sex`
  - Per-cell-type: `lm(proportion ~ age * sex)` interaction p-value
  - Visualization: 2D heatmap (age bins × sex) for each cell type
- **Key cell types to watch**: CD8 TEM, NK, B cells, Treg, inflammatory monocytes

#### C3. Aging Signature Validation
- **Cross-reference our FGS top genes with published signatures**:
  - Luo 2022: 12 age-changing subpopulations gene lists
  - Connolly 2024: "loss of identity" gene program
  - Filippov 2024: 330 "global aging genes"
  - Tabula Muris Senis: ribosomal/EIF aging signature
- **Method**: GSEA / Fisher's overlap test
- **Korean-specific deviation**: genes that are age-associated in our data but NOT in published European/Chinese cohorts

#### C4. MAIT Cell Population Analysis
- **Rationale**: AIDA found MAIT~genetic ancestry; our data shows MAIT decline with age (rho=−0.40, FDR=5.6e-5)
- **Tasks**:
  - MAIT proportion × age × sex interaction
  - MAIT marker genes (TRAV1-2, SLC4A10, KLRB1) expression with age
  - Compare MAIT decline rate with published values

#### C5. Ribosomal Gene Aging Clock
- **Rationale**: Camara-Lemarroy 2023 (Sci Adv) showed ribosome-to-inflammation ratio as aging clock
- **Method**: Score RPL/RPS genes vs S100A8/A9/A12 per cell → ratio per patient → age correlation
- **Advantage**: Novel metric for Korean cohort

### Phase D: Advanced Methods

#### D1. scANVI (semi-supervised integration)
- **Rationale**: Train scANVI with known cell types → better rare cell type identification
- **Application**: Confirm annotation, discover novel subtypes

#### D2. Subclustering Priority Cell Types
- **Priority order**: CD8 TEM, NK, CD14 Mono, B cells
- **Per subcluster**: age/sex frequency + DEG + FGS
- **Expected yield**: age-related subpopulation shifts within major types

---

## 3. Execution Priority

| Priority | Analysis | Est. Time | Script |
|----------|----------|-----------|--------|
| **1** | MASC (6 comparisons) | 1-2 hrs | `run_masc_hc.R` |
| **2** | MILO age (continuous) | 30 min | `run_milo_age.R` |
| **3** | MELD (Python, GPU) | 30 min | `run_meld.py` |
| **4** | scCODA (Python) | 30 min | `run_sccoda.py` |
| **5** | CellChat interpretation | 1 hr | `interpret_cellchat_age.R` |
| **6** | MNN interpretation | 1 hr | `interpret_mnn_results.R` |
| **7** | Age×Sex interaction | 1 hr | `run_age_sex_interaction.R` |
| **8** | Aging signature validation | 1 hr | `validate_aging_signatures.R` |
| **9** | MrVI DE | 1 hr | `run_mrvi_de.py` |
| **10** | scCODA, subclustering, scANVI | 2+ hrs | various |

---

## 4. Output Structure

```
/data/user3/sobj/hc_only_v1/
├── masc/                        # MASC results (6 comparisons)
├── milo_age/                    # MILO age DA
├── meld/                        # MELD scores
├── sccoda/                      # scCODA results
├── mrvi_de/                     # MrVI DE results
├── cellchat_interpretation/     # CellChat parsed results
├── mnn_interpretation/          # MNN parsed results
├── age_sex_interaction/         # Age×Sex analysis
├── aging_signatures/            # Literature validation
└── figures/                     # Publication-quality figures
```

## 5. Key Genes to Track Across Analyses

### Age-declining (from our FGS + literature)
PLAC8, EPSTI1, ISG15, JAML, SOX4, BLK, ROBO1, CCR7, LEF1, TCF7, SELL

### Age-increasing (from our FGS + literature)
MSC, CDKN2A, IGFBP3, SPI1, EMP3, LAG3, GZMA, GZMB, S100A8, S100A9

### Korean-specific candidates (from AIDA)
FOXP3 (Treg depletion), TRAV1-2/SLC4A10 (MAIT ancestry), MAIT frequency

### X-escape / sex-dimorphic
XIST, RPS4X, DDX3X, EIF2S3, JPX, CD99, TLR7, KDM6A
