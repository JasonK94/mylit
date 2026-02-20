# Differential Gene Dynamics Methods — Detailed Comparison

> Created: 2026-02-18

## Overview

This document compares established methods for testing and quantifying differential gene expression dynamics along pseudotime between conditions in multi-sample scRNA-seq.

---

## 1. tradeSeq (Van den Berge et al., Nature Communications 2020)

### Model
Negative binomial GAM per gene:
```
g(E[Ygi]) = Σ_l Σ_j s_lj(T_li) × Z_li + U_i × α_g
```
- Smooth splines per lineage, NB family
- Uses `mgcv::gam()` backend
- `fitGAM()` fits per-gene models

### Tests provided
| Test | Null hypothesis | Use case |
|------|-----------------|----------|
| `associationTest` | Expression constant along pseudotime | Gene associated with trajectory |
| `startVsEndTest` | Start expression = end expression | Activation/repression |
| `diffEndTest` | End expression equal across lineages | Differential fate |
| `patternTest` | Expression pattern identical across lineages | Different dynamics shape |
| `conditionTest` | Expression pattern identical across conditions | **Differential dynamics between conditions** |
| `earlyDETest` | Early expression identical across lineages | Early divergence |

### conditionTest (key for our use case)
- Wald test on GAM spline coefficients between conditions
- Returns: Wald statistic, df, p-value per gene
- `l2fc` argument: minimum log2 fold-change threshold
- Can test global (all conditions) or pairwise

### Effect size
- **Wald statistic**: test statistic, not directly interpretable as effect size
- **meanLogFC**: average log2 fold-change across pseudotime (available for some tests)
- **No built-in area between curves** — must compute post-hoc

### Pseudoreplication handling
- **Not handled** in base tradeSeq — treats cells as independent
- Can be partially addressed by fitting on pseudobulk or using `condiments` wrapper
- The original paper acknowledges this limitation

### Strengths
- Well-established, widely cited (>500 citations)
- Flexible GAM framework
- Multiple targeted tests for different biological questions
- Active maintenance (Bioconductor)

### Weaknesses
- No built-in sample/patient random effects
- Pseudoreplication inflates significance
- Effect size reporting limited to Wald stat and meanLogFC
- Computationally heavy for large datasets (70K+ cells)

### References
- Van den Berge et al. (2020). Trajectory-based differential expression analysis for single-cell sequencing data. *Nature Communications*, 11, 1201.

---

## 2. Lamian (Hou et al., Nature Communications 2023)

### Model
Mixed-effects B-spline model:
```
Y_sc = φ(t_sc)ᵀ(B·x_s + u_s) + ε_sc
```
- φ(t): B-spline basis functions (k knots)
- B: fixed effect matrix (covariate → spline coefficient mapping)
- x_s: sample-level covariates (e.g., disease status)
- u_s: sample-level random effects (captures inter-sample variability)
- ε_sc: cell-level noise

### Tests provided
| Test | What it detects | Null hypothesis |
|------|----------------|-----------------|
| **TDE** (Trajectory DE) | Genes varying along pseudotime | Expression flat along trajectory |
| **XDE overall** | Genes with any covariate effect | No difference between conditions |
| **XDE mean** | Uniform shift (intercept difference) | Same mean expression across conditions |
| **XDE trend** | Shape difference (slope/curvature change) | Same trajectory shape across conditions |

### XDE gene classification
Lamian classifies XDE genes into three categories:
1. **Mean shift only**: Parallel curves at different levels
2. **Trend difference only**: Same mean but different shapes
3. **Both mean and trend**: Both shifted and differently shaped

This is exactly what we need — distinguishes S100A8 (probably mean + trend) from a gene with only a level shift.

### Effect size
- **B-spline coefficients** (β): how covariate modifies each basis function weight
- **Fitted curve differences**: Lamian produces group-specific fitted curves
- **No built-in ABC/RMISE metric** — but trivially computable from fitted curves
- Genes clustered by temporal pattern (hierarchical clustering on fitted differences)

### Pseudoreplication handling
- **Properly handled** via sample-level random effects u_s
- EM algorithm jointly estimates fixed effects, random effects, and variances
- Permutation-based p-values (100 permutations by default) at sample level
- This is the gold standard approach

### Strengths
- Properly handles multi-sample design (the only method that does this natively)
- Distinguishes mean shift from trend difference
- FDR well-calibrated (demonstrated in simulations)
- Comprehensive framework (topology + expression + density)
- Published in Nature Communications 2023

### Weaknesses
- Requires aligned pseudotime across samples (assumes shared trajectory)
- Permutation-based testing is slow (100× model fitting)
- B-spline basis may not capture sharp expression changes as well as penalized splines
- Less widely adopted than tradeSeq (newer)
- Documentation less mature

### References
- Hou et al. (2023). A statistical framework for differential pseudotime analysis with multiple single-cell RNA-seq samples. *Nature Communications*, 14, 7286.

---

## 3. condiments (Roux de Bézieux et al., Nature Communications 2024)

### Focus
Differential **topology** (branching), **progression** (pseudotime distribution), and **differentiation** (fate selection) across conditions. Gene expression is secondary.

### Three-step workflow
1. **Differential topology**: Are trajectory structures (branches) different between conditions?
2. **Differential progression**: Are cells distributed differently along pseudotime? (KS test on pseudotime distributions)
3. **Differential differentiation**: Do cells in different conditions choose different fates?

### Gene expression analysis
- Wraps tradeSeq's `conditionTest` for differential expression
- Adds topology-aware corrections
- Does not add its own expression model

### Pseudoreplication handling
- Sample weights for topology/progression tests
- Gene expression still uses tradeSeq (no sample RE)

### Relevance to our analysis
- **Differential progression** (step 2) is directly relevant — we already showed pseudotime distribution differences (KS test p < 2e-16 for monocytes)
- Gene expression analysis is just tradeSeq under the hood
- Not a replacement for proper differential dynamics testing

### References
- Roux de Bézieux et al. (2024). Trajectory inference across multiple conditions with condiments. *Nature Communications*, 15, 833.

---

## 4. Our GAMM v3 (Current Implementation)

### Model
```
expr ~ s(pseudotime, k=8, bs="cr") +    # shared baseline smooth
       cond +                              # condition main effect
       s(pseudotime, by=cond, k=8, bs="cr") +  # condition-specific smooth (TEST)
       offset(log(nCount_RNA)) +           # library size
       percent.mt +                        # QC covariate
       s(GEM, bs="re")                     # batch RE (17 levels)
```
Family: NB(log), method: REML, engine: `mgcv::bam()`

### Test
Chi-squared ANOVA: full model vs model without `s(pseudotime, by=cond)` interaction term.

### Current output
- interaction_p, interaction_padj (BH correction)
- dev_explained (deviance explained by full model)
- TV (total variation of fitted smooth per condition)
- DR (derivative range per condition)

### Pseudoreplication handling
- GEM random effect (17 levels) — partially mitigates
- **No patient random effect** (collinear with condition for single disease variable)
- Cell-level model → inflated significance (known limitation)

### What's missing
- **No ABC/RMISE effect size** — only p-values and ad-hoc TV/DR
- **No systematic gene ranking** by effect magnitude
- **No cross-layer comparison framework**
- **No density weighting**

---

## 5. Pseudobulk GAMM (Planned Approach 2)

### Model
```
logCPM ~ s(pseudotime_bin, by=condition, k=5, bs="cr") +
         condition +
         (1|patient) +
         offset(log(lib_size))
```
Engine: `glmmTMB` or `limma-voom + duplicateCorrelation`

### Key design
1. Bin pseudotime into 10 quantile-based bins per patient
2. Sum raw counts per patient × bin (pseudobulk)
3. TMM normalize across all patient-bin samples
4. Fit mixed model with **patient random effect**

### Advantages
- Statistical unit = patient-bin (not cell) → no pseudoreplication
- Patient RE properly handles within-patient correlation
- More conservative, reviewer-friendly
- Effect size from fitted curves directly meaningful at patient level

### Disadvantages
- Loses continuous pseudotime resolution
- Sparse patients may drop bins → unbalanced design
- Implementation complexity higher

---

## Method Selection Rationale

### For this project, we recommend a two-pronged approach:

**Primary**: Enhance current GAMM v3 with post-hoc effect sizes (ABC, rABC, RMISE)
- Already computed for all 6 analyses
- Just need to extract prediction curves and compute metrics
- Fast, no re-fitting needed

**Validation**: Implement Lamian as the formal statistical framework
- Properly handles pseudoreplication
- Provides mean shift vs trend distinction
- Publication-ready (Nature Communications 2023)
- If results concordant with GAMM v3 → strong evidence
- If discordant → Lamian results take precedence

**Sensitivity**: Pseudobulk GAMM (Approach 2)
- Third line of evidence
- Patient-level aggregation is most conservative
- Addresses reviewer concern about pseudoreplication directly

---

## Curve Difference Metrics — Detailed Discussion

### Why p-values alone are insufficient

Consider two genes:
- Gene A: p = 1e-50, but curves differ by 0.1 expression units everywhere
- Gene B: p = 0.01, but curves differ by 5 expression units in a specific region

Gene B is biologically more interesting, but Gene A has a smaller p-value (due to high sample size → extreme power). Effect sizes are essential for biological interpretation.

### L1 Norm (ABC) vs L2 Norm (RMISE)

| Property | L1 (ABC = ∫|Δ|dt) | L2 (RMISE = √∫Δ²dt) |
|----------|--------------------|-----------------------|
| Sensitivity to outlier regions | Robust | Sensitive |
| Interpretation | "Total area of disagreement" | "RMS disagreement" |
| Metric space properties | Proper metric | Proper metric |
| Common in FDA | Yes (robust) | Yes (standard) |
| Favours | Broad moderate differences | Sharp large differences |

### Practical example

If f₁(t) - f₂(t) = 1 for all t ∈ [0,1] (constant shift):
- ABC = 1.0, RMISE = 1.0 (same)

If f₁(t) - f₂(t) = 10 for t ∈ [0, 0.1] and 0 elsewhere (spike):
- ABC = 1.0, RMISE = √10 ≈ 3.16

L2 amplifies the spike scenario → RMISE prefers sharp localized differences.
For biological interpretation, ABC (L1) is generally preferred because a gene that's consistently slightly different throughout pseudotime is more robust/reproducible than one with a single spike.

### Normalization for cross-gene comparison

Raw ABC is scale-dependent: a gene expressed at count=100 will have larger ABC than one at count=1. Normalize by average expression level:

```
rABC = ABC / mean_expr
     = ∫|f₁(t) - f₂(t)|dt / (0.5 × ∫(f₁(t) + f₂(t))dt)
```

This gives a "relative difference" interpretation — rABC = 0.5 means the curves differ by 50% of the mean expression level on average.

### Conventional precedent

- **Functional data analysis (FDA)** literature uses L2 norm most commonly (Ramsay & Silverman, 2005)
- **scRNA-seq trajectory** papers don't have a standardized effect size — most report only p-values
- **Lamian** implicitly uses coefficient-level tests (B-spline coefficients), which corresponds to a weighted L2-like metric in function space
- **tradeSeq** uses Wald statistics on spline coefficients — similar to L2
- **Our proposal** to use ABC (L1) + RMISE (L2) as complementary metrics is **reasonable and defensible**, though not conventional (because convention = no effect size at all!)

---

## References

1. Van den Berge et al. (2020). Trajectory-based differential expression analysis for single-cell sequencing data. *Nature Communications*, 11, 1201.
2. Hou et al. (2023). A statistical framework for differential pseudotime analysis with multiple single-cell RNA-seq samples. *Nature Communications*, 14, 7286.
3. Roux de Bézieux et al. (2024). Trajectory inference across multiple conditions with condiments. *Nature Communications*, 15, 833.
4. Song & Li (2021). PseudotimeDE: inference of differential gene expression along cell pseudotime with well-calibrated p-values. *Genome Biology*, 22, 124.
5. Ramsay, J. O., & Silverman, B. W. (2005). *Functional Data Analysis*. Springer.
6. Squair et al. (2021). Confronting false discoveries in single-cell differential expression. *Nature Communications*, 12, 5692.
