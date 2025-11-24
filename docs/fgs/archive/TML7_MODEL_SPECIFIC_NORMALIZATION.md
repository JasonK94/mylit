# Phase 1: Model-Specific Normalization

## Overview
This document details the implementation of model-specific normalization for `compute_meta_gene_importance` in TML7.

## Problem
Different L2 models produce feature importances with different scales and distributions:
- **GLM**: Log-odds coefficients (can be negative, unbounded).
- **Ranger**: Permutation importance (usually positive, but can be negative for noisy features).
- **XGBoost**: Gain (positive, unbounded).

Simply summing these (AMSC approach) without normalization can lead to one model dominating the ensemble or incorrect gene ranking.

## Solution
We implemented a flexible normalization strategy with the following methods:

### 1. `max_abs` (Default)
- **Formula**: $x' = x / \max(|x|)$
- **Range**: $[-1, 1]$
- **Use case**: Preserves sign and relative magnitude. Good for GLM and general use.

### 2. `min_max`
- **Formula**: $x' = (x - \min(x)) / (\max(x) - \min(x))$
- **Range**: $[0, 1]$
- **Use case**: When only ranking/magnitude matters, not sign. Good for tree-based models if negative importance is noise.

### 3. `softmax`
- **Formula**: $x' = \frac{e^{x - \max(x)}}{\sum e^{x - \max(x)}}$
- **Range**: $[0, 1]$, sums to 1.
- **Use case**: Interpreting importance as probabilities. Highlights top features strongly.

### 4. `rank`
- **Formula**: $x' = \text{rank}(x) / N$
- **Range**: $(0, 1]$
- **Use case**: Robust against outliers.

### 5. `z_score`
- **Formula**: $x' = (x - \mu) / \sigma$
- **Range**: Unbounded (typically $[-3, 3]$)
- **Use case**: Statistical standardization.

## Implementation
The `compute_meta_gene_importance` function was updated to accept a `normalization_method` parameter.

```r
compute_meta_gene_importance(meta_result, normalization_method = "max_abs", ...)
```

## Validation
Tested on `CD4+ T-cells` dataset with TML7 results.
- `max_abs` and `min_max` yielded similar top genes.
- `softmax` highlighted top genes more aggressively.
- `z_score` significantly altered rankings, suggesting it should be used with caution for importance scores that are not normally distributed.
