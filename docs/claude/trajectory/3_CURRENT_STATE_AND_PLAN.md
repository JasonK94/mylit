# Current State Assessment & Action Plan — Trajectory Analysis

> Created: 2026-02-18

## Current Results Summary

### Gene Dynamics v3 Results (6 analyses complete)

| Analysis | Sig genes (padj<0.05) / 36 | Top gene (p-value) | Interpretation |
|----------|----------------------------|---------------------|----------------|
| **mono/g3** | 27/36 (75%) | IFI6 (2e-13) | Strong g3 effect in monocyte trajectory |
| **mono/cohort** | 27/36 (75%) | HLA-C (2e-140) | Very strong cohort effect |
| **cd4/g3** | 2/36 (6%) | S100A9 (3e-7) | Weak g3 effect in CD4 |
| **cd4/cohort** | 22/36 (61%) | ISG15 (6e-20) | Strong cohort effect |
| **cd8/g3** | 0/36 (0%) | None | No g3 effect in CD8 |
| **cd8/cohort** | 16/36 (44%) | HLA-B (2e-15) | Moderate cohort effect |

### Key patterns
1. **Monocyte >> CD4 > CD8** for differential dynamics strength
2. **Cohort >> g3** for all compartments (HC-vs-Stroke is a larger signal than Good-vs-Bad)
3. **g3 effect concentrated in monocytes** — consistent with MASC (Inflammatory Mono OR=2.85 in Stroke)
4. **CD8 has zero g3 signal** — outcome prediction is not mediated by CD8 trajectory
5. **HLA class I genes (HLA-B, HLA-C) consistently top** in cohort comparisons across all compartments

---

## Identified Gaps

### Gap 1: No effect size / magnitude metric
**Problem**: We have 27 significant genes in mono/g3, but can't rank them by "how much the curves differ."
**Impact**: Cannot identify the biologically most important genes vs marginally significant ones.
**Solution**: Compute ABC (area between curves) from existing GAMM predictions.

### Gap 2: No cross-layer comparison
**Problem**: S100A8 is significant in both mono/g3 (p=1e-12) and mono/cohort (p=1e-40), but we can't compare the *magnitudes*.
**Impact**: Cannot answer "is the HC→Stroke shift in S100A8 dynamics larger than the Good→Bad shift?"
**Solution**: ABC for same gene across L1 and L2, then scatter plot.

### Gap 3: No compartment-level aggregate
**Problem**: Monocyte has 27 sig genes, CD8 has 0 — but p-value count is crude.
**Impact**: Need a quantitative "trajectory perturbation score" per compartment.
**Solution**: Median rABC across all 36 genes per compartment × condition.

### Gap 4: No direction analysis
**Problem**: We know curves differ, but not *how* — which condition has higher expression at which pseudotime?
**Impact**: Cannot say "g3=2 monocytes fail to upregulate S100A8 at late pseudotime."
**Solution**: Signed difference curves + identify regions of divergence.

### Gap 5: No pseudotime-resolved difference
**Problem**: ABC gives a single number per gene — loses temporal information.
**Impact**: Cannot identify "early vs late" differential dynamics.
**Solution**: Compute pointwise difference curves, identify pseudotime regions of peak divergence.

### Gap 6: TV/DR metrics are not interpretable
**Problem**: Current "TV" (total variation) and "DR" (derivative range) are ad-hoc metrics without clear biological meaning.
**Impact**: Confusing, not publishable.
**Solution**: Replace with ABC/rABC/RMISE in the summary table.

### Gap 7: Pseudoreplication concern
**Problem**: Cell-level GAMM with GEM RE only partially addresses pseudoreplication.
**Impact**: P-values may be inflated (especially for cohort comparisons where p ~ 1e-140).
**Solution**: (a) Lamian as validation, (b) Pseudobulk GAMM as sensitivity analysis.

---

## Action Plan

> **Status (2026-02-18)**: Phases 1–4 COMPLETE. See `4_TRAJECTORY_RESULTS_INTEGRATED.md` for full results.

### Phase 1: Post-hoc Effect Sizes from v3 — COMPLETE

**Goal**: Extract effect sizes from existing GAMM results without re-fitting models.

**Steps**:
1. Load saved GAMM model objects (or re-predict from saved models)
2. Generate predictions at 200 evenly-spaced pseudotime points for each condition
3. Compute per-gene metrics:
   - ABC = ∫|f₁(t) - f₂(t)|dt (trapezoidal rule)
   - rABC = ABC / mean_expression
   - RMISE = √∫(f₁(t) - f₂(t))²dt
   - MPD = max|f₁(t) - f₂(t)| and its pseudotime location
   - Direction: proportion of pseudotime where cond1 > cond2
4. Generate updated summary CSV with effect size columns
5. Create ranking plots

**Output**: `gene_dynamics_v3_effect_sizes.csv` per compartment × condition

**Script**: `scripts/compute_trajectory_effect_sizes.R`

### Phase 2: Systematic Comparison Tables & Figures (Priority: HIGH, Effort: 2-3 hours)

**Goal**: Cross-layer, cross-compartment comparison figures.

**Steps**:
1. **Gene ranking heatmap**: 36 genes × 6 analyses, colored by rABC (signed)
2. **Cross-layer scatter**: rABC(cohort) vs rABC(g3) for each gene (3 panels = 3 compartments)
   - Correlation coefficient + significance
   - Quadrant analysis: concordant vs discordant
3. **Compartment summary barplot**: Median rABC per compartment × condition (6 bars)
4. **Pointwise difference curves**: For top 10 genes, plot f₁(t) - f₂(t) with CI
5. **Integration with DEG**: Volcano overlay — genes with high rABC highlighted on DEG volcano

**Output**: `/stroke/trajectory_v3/analysis/` directory with figures + tables

### Phase 3: Lamian Validation (Priority: MEDIUM, Effort: 4-6 hours)

**Goal**: Run Lamian as formal multi-sample method to validate GAMM v3 findings.

**Steps**:
1. Install Lamian (available on GitHub: `Winnie09/Lamian`)
2. Prepare input: expression matrix + pseudotime + sample labels + condition
3. Run XDE test for each compartment × condition (6 runs)
4. Compare Lamian XDE results with GAMM v3:
   - Gene-level concordance (Jaccard of significant sets)
   - Effect direction concordance
   - Mean-shift vs trend classification
5. Report discrepancies

**Decision point**: If Lamian and GAMM v3 are concordant → report both. If discordant → investigate why, prioritize Lamian.

### Phase 4: Pseudobulk GAMM (Priority: MEDIUM, Effort: 3-4 hours)

**Goal**: Address pseudoreplication by aggregating to patient level.

**Steps**:
1. Bin pseudotime into 10 quantile-based bins per patient
2. Pseudobulk sum per patient × bin (filter: ≥10 cells per bin)
3. Fit per-gene: `edgeR` or `glmmTMB` with patient RE
4. Extract fitted curves, compute ABC
5. Compare rankings with cell-level GAMM v3

### Phase 5: Biological Interpretation (Priority: HIGH, but depends on Phase 1-2)

**Goal**: Translate quantitative trajectory differences into biological narrative.

**Steps**:
1. **Pathway enrichment of top-rABC genes** (by compartment × condition)
2. **Cross-reference with DEG/FGS**:
   - How many top-rABC genes are also DEGs?
   - Are FGS signature genes enriched among high-rABC genes?
3. **Maturation arrest narrative**:
   - Monocytes: Are classical→non-classical transition genes (FCN1→FCGR3A) differentially dynamic?
   - CD4: Naive→Memory transition genes?
   - CD8: Cytotoxic function genes (GZMB, PRF1)?
4. **Direction interpretation**:
   - "Stroke monocytes fail to upregulate S100A8/A9 at late pseudotime" = maturation block
   - "Bad outcome monocytes show premature ISG activation" = aberrant activation
5. **Figure for paper**: Composite figure with trajectory UMAP + density + top gene dynamics + ABC barplot

---

## Timeline & Dependencies

```
Week 1:
  Phase 1 (effect sizes) ──→ Phase 2 (comparisons) ──→ Phase 5 partial (narrative)

Week 2:
  Phase 3 (Lamian) ──→ Phase 5 (integration with Lamian)
  Phase 4 (Pseudobulk) parallel with Phase 3

Week 3:
  Final figures + paper section draft
```

---

## Key Questions to Resolve

1. **Saved models**: Are GAMM model objects saved, or just predictions? If only summaries, we need to re-fit for prediction curves.
2. **Pseudotime scale**: Is pseudotime [0,1] normalized or raw? Need consistent scale for ABC.
3. **Expression scale**: Predictions on log(count) or count scale? ABC interpretation depends on this.
4. **Gene set expansion**: Current 36 genes are curated targets. Should we expand to genome-wide for Lamian? (Computational cost: ~10x higher, but more comprehensive.)
5. **Presentation**: Should we present ABC as the primary metric or keep p-values primary with ABC as supplement?

---

## Expected Deliverables

| Deliverable | Format | Location |
|-------------|--------|----------|
| Effect size summary table (per analysis) | CSV | `trajectory_v3/analysis/effect_sizes_*.csv` |
| Gene ranking heatmap (6 analyses) | PNG | `trajectory_v3/analysis/heatmap_rABC.png` |
| Cross-layer scatter (3 compartments) | PNG | `trajectory_v3/analysis/crosslayer_scatter_*.png` |
| Compartment summary barplot | PNG | `trajectory_v3/analysis/compartment_summary.png` |
| Pointwise difference curves (top genes) | PNG | `trajectory_v3/analysis/diff_curves_*.png` |
| Lamian validation results | CSV + PNG | `trajectory_v3/lamian/` |
| Pseudobulk validation results | CSV + PNG | `trajectory_v3/pseudobulk/` |
| Paper figure composite | PNG/PDF | `figures/trajectory_composite.png` |
