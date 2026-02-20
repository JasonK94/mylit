# Gene Dynamics Along Pseudotime: Statistical Approaches

## Context
- Dataset: stroke_hc_v8_2 (205K cells, 112 patients, 17 GEMs)
- Dual-layer: L1 (HC vs Stroke), L2 (IS: g3==1 Good vs g3==2 Bad)
- Pseudotime: Slingshot/Monocle3 on scVI UMAP (batch-corrected embedding)
- Expression: Raw UMI counts from Seurat RNA assay

## Approach 1: Single-stage GAMM (Implemented — `analyze_gene_dynamics_v2`)

### Model
```
expr ~ s(pseudotime, k=6, bs="cr") +     # shared baseline trajectory
       disease +                           # disease main effect
       s(pseudotime, by=disease, k=6) +    # disease-specific trajectory (TEST TARGET)
       offset(log(nCount_RNA)) +           # library size normalisation
       percent.mt +                        # mitochondrial fraction (technical)
       s(GEM, bs="re")                     # batch random effect (17 levels)
family = nb(link="log"), method = "REML"
```

### Key design decisions
1. **Raw counts + NB family**: Statistically appropriate for count data. The NB
   distribution naturally handles overdispersion in scRNAseq.
2. **offset(log(nCount_RNA))**: Standard exposure offset in NB-GLM. Converts
   the model to effectively model rates (counts/library_size) without destroying
   the count nature. Does NOT absorb patient-level variance like a patient RE would.
3. **s(GEM, bs="re")**: Random intercept for batch (17 GEMs, 4-8 patients each).
   Controls for capture efficiency, sequencing depth differences across GEMs.
   Since patients are nested within GEMs but not 1:1, the GEM RE adjusts batch
   while leaving patient-level disease variance intact.
4. **No patient random effect**: disease is a patient-level variable, so
   `(1|patient)` would be collinear with disease, absorbing the effect of interest.
5. **percent.mt as fixed effect**: Controls for cell quality/stress without
   confounding with disease.
6. **Interaction test**: Full model vs null (without `s(pseudotime, by=disease)`)
   via Chi-squared ANOVA. Tests whether trajectory shape differs between conditions.

### Visualisation
- Y-axis clipped at 99th percentile to handle outlier counts
- Prediction curves at reference covariate values (median lib size, mean percent.mt)
- Batch RE excluded from predictions (marginalised)

### Limitations
- Cell-level model with 50-70K observations per compartment → computationally heavy
- Pseudoreplication: multiple cells from same patient treated independently
  (partially mitigated by GEM RE, but not fully by patient RE)
- NB dispersion estimated globally, not per-patient

---

## Approach 2: Two-stage Pseudobulk (Planned — NOT YET IMPLEMENTED)

### Rationale
Addresses pseudoreplication properly by aggregating to patient level first.
This is the Nature/Cell standard for multi-sample scRNAseq DE analysis
(Squair et al. 2021, Zimmerman et al. 2021).

### Stage 1: Pseudobulk aggregation along pseudotime
```r
# 1. Bin pseudotime into quantile-based bins (e.g., 10 bins per patient)
# 2. For each patient × bin: sum raw counts across cells
# 3. Result: patient × bin × gene matrix (pseudobulk)
# 4. Apply edgeR TMM normalisation across all patient-bin samples
```

### Stage 2: Mixed-effect model on pseudobulk
```r
# Option A: NB-GLMM (edgeR/glmmTMB)
logCPM ~ s(pseudotime_bin, by=disease) + disease +
         (1|patient) + offset(log(lib_size))

# Option B: limma-voom on pseudobulk
# Standard limma-voom pipeline treating patient-bins as samples
# Design: ~ 0 + disease:pseudotime_bin (or spline basis)
# Random effect: duplicateCorrelation(block=patient)
```

### Advantages over Approach 1
- Properly handles pseudoreplication (patient is the statistical unit)
- Can include patient random effect without disease confounding
  (because aggregation creates multiple observations per patient across bins)
- More conservative, less prone to inflated significance
- Familiar to reviewers (standard pseudobulk framework)

### Disadvantages
- Loses continuous pseudotime (binned)
- Requires sufficient cells per patient-bin (sparse patients may drop out)
- More complex implementation
- May miss fine-grained temporal patterns visible at cell level

### Implementation plan
1. Determine optimal bin count (5, 10, or 20 bins by pseudotime quantile)
2. Filter: minimum 10 cells per patient-bin
3. Pseudobulk sum per patient × bin
4. edgeR TMM normalisation
5. Fit glmmTMB or limma-voom with patient RE
6. Extract bin-wise coefficients → reconstruct trajectory curves
7. Visualise: patient-level trajectory lines (thin) + condition mean (thick)

### Status: PLANNED for future implementation

---

## Comparison

| Feature | Approach 1 (GAMM) | Approach 2 (Pseudobulk) |
|---------|-------------------|------------------------|
| Statistical unit | Cell | Patient-bin |
| Pseudoreplication | Partial control (GEM RE) | Fully addressed |
| Pseudotime | Continuous smooth | Binned |
| Patient RE | Not included (collinear) | Can include |
| Batch correction | GEM RE + offset | TMM + patient RE |
| Sensitivity | Higher (more power) | More conservative |
| Reviewer acceptance | Moderate | High |
| Implementation | Done | Planned |

## Recommendation
- Use Approach 1 for exploratory analysis and hypothesis generation
- Use Approach 2 for publication-grade results and sensitivity analysis
- Report both in supplementary if discordant results arise
