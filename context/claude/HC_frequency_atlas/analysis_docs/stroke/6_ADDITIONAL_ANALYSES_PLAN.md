# Additional Analyses Plan — stroke_hc_v8_2

> Generated: 2026-02-17 | Status: Execution plans for 7 additional analyses

---

## Priority 1: Within-Cell-Type FGS/TML

### Rationale
Whole-dataset FGS captures **cell composition** as its primary signal (monocyte↑ = stroke signature positive, T cell↑ = negative). This is biologically meaningful (showing composition shifts predict outcome) but does not reveal **per-cell transcriptomic predictors**. Running FGS within individual cell types isolates expression-level changes from composition effects.

### Execution Plan

**Cell types to run** (from L2 subset `5_2_is_g3.qs`, IS patients only):

| Cell type | Approx cells (L2) | Patients | Feasible |
|-----------|-------------------|----------|----------|
| CD14+ Monocyte | ~8K | ~50 | Yes |
| CD16+ Monocyte | ~4K | ~50 | Yes |
| Inflammatory Monocyte | ~10K | ~50 | Yes |
| CD4+ T_Naive/Memory | ~13K | ~50 | Yes |
| CD8+ T_Cytotoxic | ~8K | ~50 | Yes |
| NK_cell | ~9K | ~50 | Yes |
| B_cell | ~3K | ~50 | Marginal |
| ISG+ Myeloid | ~3K | ~50 | Marginal |

**Script**: `scripts/claude/12_fgs_within_celltype.sh`

```bash
#!/bin/bash
# Within-cell-type FGS/TML sweep
# Step 1: Subset Seurat object by anno1, save as individual .qs files
# Step 2: Run FGS pipeline on each subset

INPUT="/data/user3/sobj/stroke_hc_v8_2/5_2_is_g3.qs"
FGS_SCRIPT="/home/user3/data_user3/git_repo/mylit/Git_Repo/_wt/fgs/scripts/fgs/run_fgs_pipeline.R"
OUTPUT_BASE="/data/user3/sobj/stroke_hc_v8_2/fgs/within_celltype"

CELL_TYPES=("CD14+ Monocyte" "Inflammatory Monocyte" "CD4+ T_Naive/Memory"
            "CD8+ T_Cytotoxic" "NK_cell" "CD16+ Monocyte" "B_cell" "ISG+ Myeloid")

# For each cell type: subset → FGS with n_features=50
for CT in "${CELL_TYPES[@]}"; do
    SAFE_NAME=$(echo "$CT" | tr ' /+' '___')
    taskset -c 0-15 Rscript "$FGS_SCRIPT" \
        --input "${OUTPUT_BASE}/subsets/${SAFE_NAME}.qs" \
        --output "${OUTPUT_BASE}/${SAFE_NAME}_50" \
        --target_var g3 \
        --control_vars "sex,age,GEM" \
        --n_features 50 \
        --cores 16
done
```

**Pre-step**: Create subset .qs files via R:
```r
sobj <- qread("5_2_is_g3.qs")
cell_types <- c("CD14+ Monocyte", "Inflammatory Monocyte", ...)
for (ct in cell_types) {
    sub <- subset(sobj, anno1 == ct)
    qsave(sub, paste0("fgs/within_celltype/subsets/", sanitize(ct), ".qs"))
}
```

**Expected output per cell type**: `{ct}_50_fgs.qs`, `{ct}_50_tml.qs`, `{ct}_50_cmgi.qs`, `{ct}_50_cmgi_genes.csv`, `{ct}_50_sobj.qs`

**Analysis after completion**:
1. Compare top genes across cell types → shared vs cell-type-specific predictors
2. Compare within-cell-type signatures with whole-dataset FGS → composition vs expression
3. Pathway enrichment of cell-type-specific predictor genes
4. Cross-reference with DEG L2 (NEBULA) results

**Estimated runtime**: ~2-4 hours per cell type × 8 = 16-32 hours total

---

## Priority 2: scCODA (Bayesian Compositional Analysis)

### Rationale
MASC uses per-cell logistic regression which ignores the **compositional constraint** (cell type proportions must sum to 1). scCODA uses a Bayesian model (logistic-normal multinomial) that properly handles compositional data, providing more statistically rigorous abundance testing.

### Execution Plan

**Environment**: Python (scvi-tools conda env)

**Installation**:
```bash
conda activate scvi-tools
pip install pertpy  # scCODA is now part of pertpy
# or: pip install sccoda
```

**Script**: `scripts/claude/13_sccoda.py`

```python
import pertpy as pt
import scanpy as sc
import pandas as pd
import numpy as np

# Option 1: Create compositional table from metadata
# (no need to load full AnnData — just need counts per sample × cell type)

# Load metadata exported from R
meta = pd.read_csv("sccoda_input_counts.csv")
# Columns: patient_name, anno1, cohort/g3, n_cells

# Pivot to wide format: patients × cell types
comp_table = meta.pivot_table(index='patient_name', columns='anno1',
                               values='n_cells', fill_value=0)

# Create AnnData for scCODA
import anndata as ad
adata = ad.AnnData(X=comp_table.values,
                    obs=patient_meta,  # patient-level covariates
                    var=pd.DataFrame(index=comp_table.columns))

# Run scCODA
sccoda = pt.tl.Sccoda()
sccoda_data = sccoda.load(adata, type="cell_level",
                           cell_type_identifier="anno1",
                           sample_identifier="patient_name",
                           covariate_obs=["cohort"])  # or "g3"

# Set reference cell type (choose stable type)
sccoda.set_reference(sccoda_data, reference_cell_type="Proliferating")

# Run model
sccoda.run_nuts(sccoda_data, num_warmup=500, num_samples=10000)
sccoda.summary(sccoda_data)
sccoda.credible_effects(sccoda_data)
```

**Two runs**:
1. L1: `cohort` (HC vs Stroke) — all patients
2. L2: `g3` (Good vs Bad) — IS patients only

**Pre-step (R)**: Export cell count matrix per patient × cell type:
```r
sobj <- qread("5_strokev8_clean.qs")
counts <- sobj@meta.data %>%
    group_by(patient_name, anno1) %>%
    summarise(n = n()) %>%
    pivot_wider(names_from = anno1, values_from = n, values_fill = 0)
write.csv(counts, "sccoda_input.csv", row.names = FALSE)
```

**Expected output**: Credible intervals for each cell type's log-fold change, posterior inclusion probabilities, credible effects table.

**Estimated runtime**: 10-30 minutes

---

## Priority 3: cNMF (Consensus NMF)

### Rationale
cNMF discovers **transcriptional programs** (gene expression programs, GEPs) in an unsupervised manner. Unlike FGS (which finds genes discriminating conditions), cNMF finds co-regulated gene modules that may correspond to biological processes, cell states, or pathways. Running cNMF within cell types reveals condition-associated programs beyond simple DE.

### Execution Plan

**Environment**: Python (scvi-tools conda env)

**Installation**:
```bash
conda activate scvi-tools
pip install cnmf
```

**Script**: `scripts/claude/14_cnmf.py`

**Workflow per cell type**:
```python
from cnmf import cNMF
import scanpy as sc

# Step 1: Prepare AnnData (from Seurat via h5ad or direct conversion)
adata = sc.read_h5ad("cnmf_input_{celltype}.h5ad")

# Step 2: Initialize cNMF
cnmf_obj = cNMF(output_dir="cnmf_output/{celltype}", name="{celltype}")

# Step 3: Prepare (normalize, select HVGs)
cnmf_obj.prepare(counts=adata, components=np.arange(5, 31, 5),  # K = 5,10,15,20,25,30
                  n_iter=200, seed=42, num_highvar_genes=2000)

# Step 4: Factorize (embarrassingly parallel)
cnmf_obj.factorize(worker_i=0, total_workers=1)  # or parallelize

# Step 5: Combine results
cnmf_obj.combine()

# Step 6: Select K (using stability plot)
cnmf_obj.k_selection_plot()  # visual inspection

# Step 7: Consensus (at selected K)
cnmf_obj.consensus(k=K_selected, density_threshold=0.1)

# Step 8: Get results
usage, spectra_score, spectra_tpm, top_genes = cnmf_obj.load_results(K=K_selected)
```

**Cell types to run**: Same as within-cell-type FGS (CD14+ Mono, Inflam Mono, CD4+ T, CD8+ T, NK)

**Pre-step**: Convert Seurat subsets to h5ad format using `myR::seurat_to_h5ad()` or `SeuratDisk`

**Expected output per cell type**:
- Gene programs (spectra): genes × programs matrix with weights
- Usage matrix: cells × programs, showing program activity per cell
- Top genes per program
- Program-condition association: test if program usage differs by g3/cohort

**Post-analysis**:
1. Correlate program usage with g3/cohort (Wilcoxon test per program)
2. Pathway enrichment of top genes per program
3. Compare condition-associated programs across cell types
4. Identify programs that overlap with FGS genes or DEG results

**Estimated runtime**: 1-3 hours per cell type (K selection + consensus)

---

## Priority 4: MELD (Manifold Enhancement of Latent Dimensions)

### Rationale
MELD computes per-cell **condition likelihood** using the KNN graph, producing continuous density estimates. Unlike binary condition assignment, MELD provides a nuanced view of how cells position relative to conditions on the manifold. Useful for identifying transitional/mixed-state cells.

### Execution Plan

**Installation**:
```bash
conda activate scvi-tools
pip install meld
```

**Script**: `scripts/claude/15_meld.py`

```python
import meld
import scanpy as sc
import numpy as np

# Load AnnData with scVI embeddings
adata = sc.read_h5ad("stroke_hc_v8_2_scvi.h5ad")

# Use scVI latent space for KNN graph
sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=15)

# Run MELD
meld_op = meld.MELD()
meld_op.fit(adata.obsp['connectivities'])

# For L1: HC vs Stroke
sample_labels = adata.obs['cohort'].values  # or use patient_name for sample-level
meld_densities = meld_op.transform(sample_labels)
adata.obs['meld_stroke_likelihood'] = meld_densities[:, 1]  # P(Stroke)

# For L2: g3 (IS only)
is_mask = adata.obs['index_injury_hand'] == 'IS'
adata_is = adata[is_mask].copy()
meld_op2 = meld.MELD()
meld_op2.fit(adata_is.obsp['connectivities'])
g3_densities = meld_op2.transform(adata_is.obs['g3'].values)
adata_is.obs['meld_bad_likelihood'] = g3_densities[:, 1]  # P(g3=2, Bad)
```

**Expected output**:
- Per-cell condition likelihood scores
- UMAP colored by MELD likelihood → identify regions of condition mixing
- Correlation of MELD scores with FGS signature scores
- Cell-type-level aggregation of MELD scores → compare with MASC/MILO

**Estimated runtime**: 10-30 minutes (graph-based, fast)

---

## Priority 5: Augur (Cell Type Prioritization)

### Rationale
Augur quantifies how well each cell type separates conditions using cross-validated classification (random forests). The AUC per cell type reflects **discriminability** — higher AUC = more transcriptionally different between conditions. Complementary to MASC (abundance) because Augur captures **expression** differences.

### Execution Plan

**Installation**:
```r
# R package
remotes::install_github("neurorestore/Augur")
```

**Script**: `scripts/claude/16_augur.R`

```r
library(Augur)
library(Seurat)
library(qs)

sobj <- qread("5_strokev8_clean.qs")

# L1: HC vs Stroke
augur_l1 <- calculate_auc(sobj,
                           cell_type_col = "anno1",
                           label_col = "cohort",
                           n_subsamples = 50,
                           subsample_size = 20,
                           min_cells = 50)

# L2: g3 (IS only)
sobj_is <- subset(sobj, index_injury_hand == "IS" & !is.na(g3))
augur_l2 <- calculate_auc(sobj_is,
                           cell_type_col = "anno1",
                           label_col = "g3",
                           n_subsamples = 50,
                           subsample_size = 20,
                           min_cells = 50)

# Plot results
plot_lollipop(augur_l1)
plot_lollipop(augur_l2)
```

**Expected output**: AUC per cell type, lollipop plots, rankings

**Estimated runtime**: 30-60 minutes

---

## Priority 6: MOFA+ (Multi-Omics Factor Analysis)

### Rationale
MOFA+ discovers **patient-level latent factors** from multi-view data. By creating pseudobulk expression profiles per patient × cell type, MOFA+ identifies factors that explain variation across cell types simultaneously. This bridges cell-type-specific effects to patient-level phenotypes.

### Execution Plan

**Installation**:
```r
# R interface
BiocManager::install("MOFA2")
# Python backend
pip install mofapy2
```

**Script**: `scripts/claude/17_mofa.R`

```r
library(MOFA2)
library(Seurat)
library(qs)

sobj <- qread("5_strokev8_clean.qs")

# Create pseudobulk: patient × gene per cell type (views)
# Select top variable genes per cell type (e.g., top 2000 HVG)
cell_types <- c("CD14+ Monocyte", "Inflammatory Monocyte",
                "CD4+ T_Naive/Memory", "CD8+ T_Cytotoxic", "NK_cell")

# Build MOFA input: list of matrices (views × patients × genes)
mofa_data <- list()
for (ct in cell_types) {
    sub <- subset(sobj, anno1 == ct)
    pb <- AggregateExpression(sub, group.by = "patient_name",
                               slot = "counts", return.seurat = FALSE)
    mofa_data[[ct]] <- log1p(pb$RNA)  # or use normalized
}

# Create MOFA object
mofa <- create_mofa(mofa_data)

# Define options
data_opts <- get_default_data_options(mofa)
model_opts <- get_default_model_options(mofa)
model_opts$num_factors <- 15

train_opts <- get_default_training_options(mofa)
train_opts$convergence_mode <- "medium"
train_opts$seed <- 42

mofa <- prepare_mofa(mofa, data_options = data_opts,
                      model_options = model_opts,
                      training_options = train_opts)

# Train
mofa <- run_mofa(mofa)

# Analyze factors
plot_variance_explained(mofa)
plot_factor(mofa, factor = 1, color_by = "cohort")
plot_weights(mofa, view = 1, factor = 1, nfeatures = 20)
```

**Expected output**:
- Latent factors explaining patient-level variation
- Factor-condition association (which factors correlate with cohort/g3)
- Gene weights per factor per cell type → shared/specific genes
- Variance explained per view per factor

**Estimated runtime**: 30-60 minutes

---

## Priority 7: MixedPools (Compositional Mixed Models)

### Rationale
MixedPools / mixed-effects compositional models address the same question as MASC and scCODA but use Dirichlet-multinomial mixed models that handle repeated measures (multiple samples per patient) and compositional constraints simultaneously.

### Execution Plan

**Options**:
1. **DirichletReg** (R): Dirichlet regression on patient-level proportions
2. **corncob** (R): Beta-binomial regression for differential abundance
3. **propeller** (R, from limma authors): Logit-transformed proportions + limma

**Script**: `scripts/claude/18_mixedpools.R`

```r
# Using propeller (speckle package)
library(speckle)
library(limma)

# Create proportion table
props <- getTransformedProps(sobj$anno1, sobj$patient_name)

# Run propeller
results <- propeller(clusters = sobj$anno1,
                      sample = sobj$patient_name,
                      group = sobj$cohort)

# Or using corncob
library(corncob)
# ... beta-binomial models per cell type
```

**Estimated runtime**: 10-30 minutes

---

## Execution Priority & Dependencies

```
Phase 1 (Immediate — can run in parallel):
├── Within-cell-type FGS/TML (16-32 hrs, CPU)
├── scCODA (30 min, Python)
├── Augur (30-60 min, R)
└── MELD (10-30 min, Python)

Phase 2 (After Phase 1):
├── cNMF per cell type (8-15 hrs, Python)
├── MOFA+ (30-60 min, R/Python)
└── MixedPools / propeller (10-30 min, R)

Phase 3 (After all above):
├── Integration of results
├── Consensus frequency analysis (MASC + MILO + scCODA + Augur + MELD)
├── Comparison: within-cell-type FGS vs cNMF programs vs DEG
└── External validation with top findings
```

---

## External Validation Datasets

| Dataset | Type | Comparison | Status |
|---------|------|------------|--------|
| GSE16561 | Illumina array, 39 IS + 24 HC (PBMC) | HC vs IS validation | Downloaded |
| GSE140275 | RNA-seq FPKM, stroke blood | HC vs IS validation | Downloaded |
| GSE22255 | Microarray, stroke blood | HC vs IS validation | Downloading |
| GSE58294 | Microarray, stroke blood (time series) | HC vs IS temporal | Downloading |
| BACTRAC | Blood clot SCFA data | Supplementary | Downloaded |
| UKB-PPP | Proteomics (4,907 proteins × 54K) | Protein-level validation | Access via Synapse |

**Validation strategy**: Apply FGS signatures (whole + within-cell-type) to external bulk datasets using ssGSEA / GSVA scoring, then test association with stroke status/outcome.
