# Trajectory Analysis Strategy — stroke_hc_v8_2

> Created: 2026-02-18 | Status: Active development

## Problem Statement

Trajectory plotting (Slingshot + Monocle3, scVI UMAP, 3 compartments) and gene dynamics GAMM fitting are complete. However, **systematic quantitative analysis** is lacking:

1. **Curve difference quantification**: How to measure "how different" two fitted expression curves are along pseudotime?
2. **Effect size ranking**: Which genes show the *largest* condition-dependent trajectory differences?
3. **Cross-comparison**: HC-vs-Stroke differences vs g3 Good-vs-Bad differences — are they concordant?
4. **Compartment comparison**: Which compartment shows the strongest differential dynamics?
5. **Biological interpretation**: What do differential dynamics *mean* — maturation arrest? Activation? Exhaustion?

---

## Current State (v3 Gene Dynamics)

### What we have
- 36 genes × 3 compartments × 2 conditions = 6 GAMM analyses (all complete)
- Per-gene: interaction p-value (Chi-sq ANOVA), deviance explained, TV/DR metrics
- Fitted GAM prediction curves (visualized in 216 PNG plots)

### What's missing
1. **No effect size metric** — p-values alone don't tell us *how much* the curves differ
2. **No systematic ranking** — can't compare S100A8 "difference magnitude" vs TXNIP
3. **No cross-layer comparison** — can't say "mono/cohort S100A8 difference > mono/g3 S100A8 difference"
4. **No compartment-level summary** — can't say "monocyte trajectory is more perturbed than CD8"
5. **TV/DR metrics** are ad-hoc (total variation / derivative range) — not standardized, not interpretable

---

## Research Question Hierarchy

### Level 1: Per-gene differential dynamics
> "Does gene X follow a different expression trajectory in condition A vs B?"

**Current answer**: interaction p-value from GAMM ANOVA — this works, but only tells us significance, not magnitude.

### Level 2: Gene ranking by effect size
> "Which genes show the *biggest* condition-dependent trajectory changes?"

**Needed**: A scalar effect size metric that:
- Captures both mean shifts and shape changes
- Is comparable across genes (normalized)
- Works with our GAMM framework

### Level 3: Cross-condition comparison
> "Is the HC→Stroke trajectory shift concordant with the Good→Bad shift?"

**Needed**: Paired effect sizes for each gene across L1 and L2, then correlation/concordance analysis.

### Level 4: Compartment-level summary
> "Which compartment is most perturbed by stroke / outcome?"

**Needed**: Aggregate metric (e.g., median effect size across genes) per compartment × condition.

### Level 5: Biological interpretation
> "What biological process does the trajectory perturbation reflect?"

**Needed**: Link differential dynamics genes to pathways, cross-reference with DEG/CCI/FGS results.

---

## Proposed Effect Size Metrics

### Option A: Area Between Curves (ABC) — ∫|f₁(t) - f₂(t)|dt

**Definition**: Integrate the absolute difference between two fitted prediction curves over the pseudotime range [0, 1].

```
ABC = ∫₀¹ |f̂_cond1(t) - f̂_cond2(t)| dt
```

**Properties**:
- Intuitive: "total area of disagreement"
- Scale-dependent: larger for highly-expressed genes
- L1 norm of the difference function
- Easy to compute from discrete predictions (trapezoidal rule)
- Sensitive to both mean shifts and shape differences

**Normalization**: Divide by mean expression to get relative ABC:
```
rABC = ABC / (0.5 × ∫₀¹ (f̂_cond1(t) + f̂_cond2(t)) dt)
```

### Option B: Root Mean Integrated Squared Error (RMISE) — √∫(f₁(t) - f₂(t))²dt

**Definition**: L2 norm of the difference function.

```
RMISE = √(∫₀¹ (f̂_cond1(t) - f̂_cond2(t))² dt)
```

**Properties**:
- Penalizes large local differences more than ABC
- Standard in functional data analysis
- More sensitive to sharp peaks of divergence
- Less robust to outlier regions

### Option C: Maximum Pointwise Difference (MPD)

```
MPD = max_t |f̂_cond1(t) - f̂_cond2(t)|
```

**Properties**:
- Simple, interpretable: "biggest gap between curves"
- Ignores the breadth of difference (narrow spike = same as broad difference)
- Useful as complementary metric

### Option D: Weighted ABC (wABC)

Weight by cell density at each pseudotime position — regions with more cells contribute more to the effect size.

```
wABC = ∫₀¹ |f̂_cond1(t) - f̂_cond2(t)| × d(t) dt
```

where d(t) is the combined cell density at pseudotime t.

**Rationale**: A large difference at pseudotime=0.9 where few cells exist is less biologically meaningful than a moderate difference at pseudotime=0.3 where many cells reside.

### Recommendation

**Use ABC (L1 norm) as primary metric**, with:
- **rABC** (relative) for cross-gene comparison
- **RMISE** (L2 norm) as sensitivity analysis
- **wABC** (density-weighted) for biological interpretation

**Rationale**:
- ABC is more robust than RMISE to isolated spikes
- ABC has a direct interpretation as "total shift"
- L1 norm is conventional in functional data analysis for comparing curves
- Neither Lamian, tradeSeq, nor condiments provide a built-in effect size — this is commonly computed post-hoc

---

## Conventional Methods Comparison

See [2_DIFFERENTIAL_DYNAMICS_METHODS.md](2_DIFFERENTIAL_DYNAMICS_METHODS.md) for detailed comparison.

### Summary

| Method | Statistical unit | Pseudoreplication | Effect size | Our applicability |
|--------|-----------------|-------------------|-------------|-------------------|
| **tradeSeq** | Cell | Not handled | Wald stat + l2fc | ⚠️ No sample RE |
| **Lamian** | Sample (mixed) | Properly handled | B-spline coeff, permutation p | ✅ Best fit |
| **condiments** | Cell (topology focus) | Partial | KS statistic | ⚠️ Not expression-focused |
| **Our GAMM v3** | Cell | GEM RE (partial) | interaction p (no effect size) | ✅ Already done, needs effect size |
| **Pseudobulk GAMM** | Patient-bin | Fully handled | ABC / RMISE on fitted curves | ✅ Planned Approach 2 |

---

## Proposed Analytical Framework

### Stage 1: Extract effect sizes from existing v3 results (immediate)

From the already-fitted GAMM models, extract prediction curves and compute:
1. **ABC** (area between curves) per gene × compartment × condition
2. **rABC** (normalized)
3. **Max pointwise difference** and its pseudotime location
4. **Direction**: which condition has higher expression at which pseudotime regions

**Implementation**: Predict from saved models at 200 evenly-spaced pseudotime points, compute metrics.

### Stage 2: Systematic comparison tables (immediate)

1. **Gene ranking table**: Rank genes by rABC within each compartment × condition
2. **Cross-layer scatter**: rABC(cohort) vs rABC(g3) per gene — are the same genes perturbed in both comparisons?
3. **Compartment summary**: Median rABC per compartment — which compartment is most affected?
4. **Heatmap**: genes × (compartment × condition), colored by signed ABC

### Stage 3: Pseudobulk validation (Approach 2, future)

Implement two-stage pseudobulk GAMM (from `2_gene_dynamics_approaches.md`):
1. Bin pseudotime into quantiles per patient
2. Pseudobulk sum per patient × bin
3. Fit mixed model with patient RE
4. Compare ABC from pseudobulk vs cell-level — do rankings agree?

### Stage 4: Formal method (Lamian or tradeSeq + condiments)

If reviewer requires established method:
- **Lamian**: Best for multi-sample, properly handles pseudoreplication, provides XDE test
- **tradeSeq + condiments**: More established ecosystem, but weaker on pseudoreplication
- Run as sensitivity analysis to confirm our GAMM findings

---

## Integration with Other Analyses

| Analysis | How trajectory integrates |
|----------|--------------------------|
| **DEG consensus** | Genes with large ABC should overlap with L1/L2 DEGs |
| **FGS** | FGS signature genes — do they show differential dynamics? |
| **CCI** | Ligand/receptor genes with differential dynamics → mechanism |
| **MASC/MILO** | Cell types with frequency shifts → trajectory maturation arrest |
| **Pathway** | ABC-ranked genes → pathway enrichment (do pro-inflammatory pathways dominate?) |

---

## File Index

| File | Content |
|------|---------|
| `1_TRAJECTORY_ANALYSIS_STRATEGY.md` | This document — overall strategy and framework |
| `2_DIFFERENTIAL_DYNAMICS_METHODS.md` | Detailed methods comparison (tradeSeq, Lamian, condiments, GAMM) |
| `3_CURRENT_STATE_AND_PLAN.md` | Current results assessment, gaps, and concrete action items |
