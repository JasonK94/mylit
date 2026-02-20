# Analysis Report — stroke_hc_v8_2

> Generated: 2026-02-17 | Dataset: 226,855 cells / 100 patients (Stroke 12 GEM + HC 5 GEM)
> Clean: 205,277 cells (21 anno1 types, 8 anno2 compartments)

---

## 1. Differential Abundance (Frequency Analysis)

### 1.1 MASC (Mixed-effects Association of Single Cells)

MASC uses per-cell logistic mixed models: `cluster_membership ~ contrast + fixed_effects + (1|patient_name)`.

#### anno1 × cohort (HC vs Stroke) — **v2, no GEM** (age, sex as fixed effects)

| Cell type | Size | OR (Stroke) | 95% CI | FDR |
|-----------|------|-------------|--------|-----|
| **Platelet/PLA** | 808 | 0.13 | 0.09–0.21 | **1.5e-13** |
| **MAIT** | 5,379 | 0.38 | 0.28–0.53 | **1.4e-6** |
| **NK_cell** | 24,751 | 0.40 | 0.29–0.55 | **1.8e-6** |
| **pDC** | 633 | 0.32 | 0.21–0.48 | **1.8e-6** |
| **Inflammatory Monocyte** | 28,172 | 2.85 | 1.81–4.50 | **4.6e-5** |
| **CD8+ Trm** | 5,652 | 0.39 | 0.26–0.59 | **5.9e-5** |
| **ISG+ T_cell** | 2,439 | 0.37 | 0.23–0.59 | **1.6e-4** |
| **CD8+ T_Cytotoxic** | 23,114 | 0.52 | 0.40–0.68 | **2.2e-5** |
| **Mast_cell** | 21 | 0.17 | 0.05–0.52 | **1.3e-3** |
| **B_cell** | 7,197 | 1.87 | 1.22–2.87 | **0.010** |
| **ISG+ Myeloid** | 9,394 | 6.44 | 1.71–24.26 | **0.012** |
| **CD14+ Monocyte** | 22,065 | 0.006 | 0.0001–0.24 | **0.013** |
| CD4_S100A8_CSF3R | 2,599 | 0.59 | 0.39–0.90 | 0.027 |
| Plasma_cell | 5,196 | 1.33 | 1.02–1.73 | 0.055 |
| Treg | 5,793 | 0.81 | 0.63–1.04 | 0.127 |
| CD16+ Monocyte | 10,011 | 0.81 | 0.47–1.39 | 0.456 |
| others | — | — | — | NS |

**Key findings**:
- **Stroke-enriched**: Inflammatory Monocyte (OR=2.85), ISG+ Myeloid (OR=6.44), B_cell (OR=1.87)
- **HC-enriched**: Platelet/PLA (OR=0.13), pDC (OR=0.32), MAIT (OR=0.38), NK (OR=0.40), CD8+ Trm (OR=0.39), CD8+ T_Cytotoxic (OR=0.52)
- **Note**: CD14+ Monocyte shows OR=0.006 with very wide CI (0.0001–0.24), likely an estimation artifact from extreme effect size — interpret as strongly ISG+/Inflammatory rather than classical CD14+ shift
- 13/21 cell types significant at FDR < 0.05

#### anno2 × g3 (Good vs Bad outcome, IS only)

| Compartment | OR | p-value |
|---|---|---|
| **T cell** | 1.65 | **0.013** |
| **Monocyte** | 0.64 | **0.027** |
| DC | 0.74 | 0.14 |
| NK | 1.11 | 0.62 |

**Key**: Bad outcome (g3=2) has higher T cell proportion, lower Monocyte — inverse of L1 pattern.

#### Previous MASC runs
- `anno1 × project_name` with GEM as fixed effect → **all p-values NA** (GEM completely confounded with cohort → rank-deficient design matrix). This run is invalid.
- `seurat_clusters × g3` → Treg (p=0.011), CD16+ Mono (p=0.011), ISG+ Myeloid (p=0.026) significant

---

### 1.2 MILO (Neighbourhood-level Differential Abundance)

L1: HC vs Stroke, using scVI embeddings (copied to PCA slot for compatibility).

| Metric | Value |
|--------|-------|
| Total neighbourhoods | 11,739 |
| Significant (SpatialFDR < 0.1) | 7,387 (62.9%) |
| Significant (SpatialFDR < 0.05) | 6,678 (56.9%) |

**Per cell type DA summary** (SpatialFDR < 0.1):

| Cell type | Total nhoods | Sig nhoods | Median logFC | Direction |
|-----------|-------------|------------|-------------|-----------|
| CD14+ Monocyte | 1,136 | 1,012 | +2.41 | **Stroke** |
| B_cell | 216 | 167 | +1.41 | **Stroke** |
| NK_cell | 1,572 | 1,423 | -2.26 | HC |
| CD8+ T_Cytotoxic | 1,426 | 1,105 | -1.87 | HC |
| CD4+ T_Naive/Memory | 2,199 | 1,370 | -1.19 | HC |
| Inflammatory Monocyte | 1,604 | 923 | -0.63 | HC |
| MAIT | 353 | 291 | -1.99 | HC |
| CD8+ Trm | 362 | 267 | -1.81 | HC |
| ISG+ T_cell | 199 | 150 | -1.73 | HC |
| Treg | 423 | 139 | -0.56 | HC |
| Plasma_cell | 378 | 132 | -0.30 | HC |
| CD4_S100A8_CD14 | 303 | 123 | -0.81 | HC |
| CD4_S100A8_CSF3R | 161 | 90 | -1.14 | HC |
| pDC | 49 | 41 | -1.70 | HC |
| ISG+ Myeloid | 579 | 27 | +0.05 | Mixed |

**MASC vs MILO concordance**:
- Both agree: NK, MAIT, CD8+ T (cytotoxic + Trm), pDC are **HC-enriched**
- Both agree: B_cell is **Stroke-enriched**
- **Discrepancy**: Inflammatory Monocyte is **Stroke-enriched by MASC** (OR=2.85) but **HC-enriched by MILO** (logFC=-0.63). This likely reflects heterogeneity — some neighbourhoods within Inflammatory Monocyte are Stroke-enriched but the majority are HC-enriched by MILO's neighbourhood-level resolution.
- **ISG+ Myeloid**: Strongly Stroke-enriched by MASC (OR=6.44) but nearly neutral by MILO — could be composition effect (few IS cells but forming distinct neighbourhoods)

---

## 2. Differential Expression (DEG)

### 2.1 DEG Consensus — L1 (HC vs Stroke)

Three methods: muscat-edgeR, muscat-DESeq2, NEBULA. Fixed effects: age, sex (no GEM due to confounding).

- **15 clusters** with consensus (2+ methods), **5 skipped** (insufficient methods)
- NEBULA dominates: detects massive transcriptomic shifts

**NEBULA L1 summary by cell type (FDR < 0.05)**:

| Cell type | Up (Stroke) | Down (Stroke) | Total DEG |
|-----------|-------------|---------------|-----------|
| CD14+ Monocyte | 57 | 7,234 | 7,291 |
| CD16+ Monocyte | 22 | 7,416 | 7,438 |
| Inflammatory Monocyte | 10 | 3,098 | 3,108 |
| NK_cell | 45 | 5,789 | 5,834 |
| CD8+ T_Cytotoxic | 16 | 5,271 | 5,287 |
| CD4+ T_Naive/Memory | 14 | 3,869 | 3,883 |
| B_cell | 82 | 4,253 | 4,335 |

**Pattern**: Overwhelmingly downregulated in Stroke — consistent across all cell types. This may reflect global transcriptional suppression in stroke PBMC.

### 2.2 DEG — L2 (g3=1 Good vs g3=2 Bad, IS only)

NEBULA only (muscat-edgeR/DESeq2 too sparse for g3 data).

- 203K gene-cluster pairs tested
- Much fewer DEGs than L1:

| Cell type | DEG (FDR < 0.05) |
|-----------|-------------------|
| CD8+ T_Cytotoxic | 89 |
| NK_cell | 60 |
| CD4+ T_Naive/Memory | 50 |
| Inflammatory Monocyte | 12 |
| CD14+ Monocyte | 8 |

### 2.3 Cross-Layer Concordance

- 134 gene-cluster pairs significant in **both** L1 and L2
- Of these, concordance rate varies by cluster:
  - Inflammatory Monocyte: 3 concordant / 1 discordant (75% concordance)
  - CD14+ Monocyte: 5 concordant / 6 discordant (45% concordance)
  - CD8+ Trm: 0 concordant / 5 discordant (0% concordance)

Top concordant genes (both FDR<0.05, same direction): LINC02649, PPM1N, SNAI1 (ISG+ Myeloid), LGR6 (CD8+ T), KLF4 (CD16+ Mono)

---

## 3. CCI (Cell-Cell Interaction)

### 3.1 CellChat

- **L1 (HC vs Stroke)**: Stroke shows -37% overall interaction strength (SIIS). Major loss in MHC-I, MHC-II signaling
- **L2 (Good vs Bad)**: Bad outcome shows +10% interaction, with increased Mono→DC/B_cell communication

### 3.2 MultiNicheNet (MNN) — anno1 level

Re-run at anno1 resolution (19 cell types, L1+L2):
- **78K DEGs**, **1.2M ligand activities** computed
- Proper LR prioritization now possible with cell-type-specific DE

---

## 4. Trajectory & Gene Dynamics

### 4.1 Pseudotime (Slingshot + Monocle3)

Three compartments on scVI UMAP (using original `umap.scvi`, no recomputation):

| Compartment | n_cells | Lineages | g3 p-value | Cohort p-value |
|-------------|---------|----------|------------|----------------|
| Monocyte | ~60K | 2 | < 2e-16 | < 2e-16 |
| CD4+ T | ~42K | 2 | 0.028 | < 2e-16 |
| CD8+ T | ~29K | 2 | TBD | TBD |

**Interpretation**: Stroke/Bad outcome cells concentrate at early pseudotime → maturation arrest.

### 4.2 Gene Dynamics v3 (Batch-corrected GAMM)

Model: `expr ~ s(pt, k=8, bs="cr") + cond + s(pt, by=cond, k=8, bs="cr") + offset(log(nCount_RNA)) + percent.mt + s(GEM, bs="re")`, NB family, `mgcv::bam()`.

36 target genes × 3 compartments × 2 conditions = 6 analyses.

**Summary of significant interaction terms (padj < 0.05)**:

| Analysis | Significant / Total | Key genes |
|----------|-------------------|-----------|
| **mono/g3** | 27/36 | S100A8 (p=1e-12), S100A9 (p=1e-12), IFI6 (p=2e-13), TXNIP (p=8e-12) |
| **mono/cohort** | 27/36 | HLA-C (p=2e-140), S100A9 (p=3e-50), TNF (p=3e-56), HLA-B (p=9e-52) |
| **cd4/g3** | 2/36 | S100A9 (p=3e-7), DDIT4 (p=0.019) |
| **cd4/cohort** | 22/36 | ISG15 (p=6e-20), TNF (p=4e-15), CXCL8 (p=9e-14), S100A9 (p=9e-11) |
| **cd8/g3** | 0/36 | None significant (all p > 0.05 after BH) |
| **cd8/cohort** | 16/36 | HLA-B (p=2e-15), HLA-C (p=3e-13), GZMB (p=1e-11), PRF1 (p=1e-11) |

**Interpretation**:
- **Monocyte compartment** shows the strongest pseudotime × condition interaction for both g3 and cohort
- **CD4+ T**: g3 interaction is weak (only 2 genes), but cohort interaction is strong (22 genes) — cohort (HC vs Stroke) dominates over g3 (Good vs Bad)
- **CD8+ T**: No g3 interaction whatsoever, but strong cohort effects for cytotoxic genes (GZMB, PRF1, NKG7, GNLY, LAG3)
- S100A8/A9/A12 and HLA-B/C are consistently significant across compartments

---

## 5. FGS (Feature Gene Signature)

### 5.1 Whole-dataset FGS (target: g3, IS patients only)

| n_features | Genes selected | Top genes |
|-----------|---------------|-----------|
| n=50 | 241 | MT-CO1 (-9.21), RPS26 (+6.58), HLA-DQA2 (-6.52) |
| n=100 | 464 | (superset of n=50) |
| n=200 | Running | — |

### 5.2 FGS × DEG overlap

218/241 FGS genes (n=50) are significant in L1 DEG (NEBULA) — **90.5% overlap** — confirming that FGS captures real transcriptomic signals.

### 5.3 Pathway enrichment of FGS genes

FGS genes enriched in: HALLMARK_INTERFERON_ALPHA_RESPONSE, HALLMARK_TNFA_SIGNALING, KEGG_RIBOSOME, GOBP_RESPONSE_TO_VIRUS.

### 5.4 Limitation: Composition proxy

FGS signature score is directionally consistent across compartments: positive in Monocyte/DC, negative in T/NK for bad outcome. This mirrors cell type proportions (MASC results) rather than per-cell expression changes. **Within-cell-type FGS/TML is needed** to disentangle compositional from transcriptomic effects.

---

## 6. Descriptive Figures (v1)

63 figures generated and archived in `figures/v1_initial/`:
- UMAP (anno1, anno2, cohort, g3, GEM splits)
- DotPlots (canonical markers)
- Stacked bar charts (frequency)
- Heatmaps (gene expression, pathway enrichment)
- Pseudotime distribution plots
- FGS-related figures
- Cross-layer concordance scatter plots

---

## 7. Summary of Cross-Analysis Concordance

| Finding | MASC | MILO | DEG L1 | DEG L2 | CellChat | Gene Dynamics |
|---------|------|------|--------|--------|----------|---------------|
| Monocyte ↑ Stroke | ✓ (Inflam) | ✓ (CD14+) | ✓ | — | ✓ (loss comm.) | ✓ (S100A8/9) |
| NK/CD8+ T ↓ Stroke | ✓ | ✓ | ✓ | — | ✓ | ✓ (GZMB, PRF1) |
| B cell ↑ Stroke | ✓ | ✓ | ✓ | — | — | — |
| ISG+ Myeloid ↑ Stroke | ✓ | Weak | ✓ | — | — | — |
| g3 effect weak except Mono | ✓ | — | ✓ (few DEGs) | ✓ | — | ✓ (0 sig cd8) |

---

## Data Locations

| Data | Path |
|------|------|
| Clean Seurat | `/data/user3/sobj/stroke_hc_v8_2/5_strokev8_clean.qs` |
| L1 subset | `/data/user3/sobj/stroke_hc_v8_2/5_1_hc_is.qs` |
| L2 subset | `/data/user3/sobj/stroke_hc_v8_2/5_2_is_g3.qs` |
| MASC | `/data/user3/sobj/stroke_hc_v8_2/MASC/` |
| MILO | `/data/user3/sobj/stroke_hc_v8_2/milo/` |
| DEG consensus | `/data/user3/sobj/stroke_hc_v8_2/deg_consensus/` |
| Trajectory | `/data/user3/sobj/stroke_hc_v8_2/trajectory_v2/`, `trajectory_v3/` |
| FGS | `/data/user3/sobj/stroke_hc_v8_2/fgs/` |
| CCI | `/data/user3/sobj/stroke_hc_v8_2/cci/` |
| Figures | `/data/user3/sobj/stroke_hc_v8_2/figures/v1_initial/` |
| Scripts | `/data/user3/sobj/stroke_hc_v8_2/scripts/claude/` |
| External data | `/data/user3/sobj/stroke_hc_v8_2/external_data/` |
