
# DEG Consensus Methods Usage

This document provides concise examples for running each implemented differential expression method.

## Setup
```{r}
sobj <- qs::qread("/data/user3/sobj/is2_IS_3_1_plots.qs")

# Parameters
group_col <- "g3"
contrast <- "2 - 1"
cluster_col <- "anno3big"
fixed_effects <- c("age", "sex", "GEM") # For fixed-effect models
random_effect_var <- "hos_no"           # For mixed models
mixed_effects <- c("age", "sex")        # Use "hos_no" as random intercept
```

## 1. Pseudobulk Methods (Fixed Effects)

### edgeR (LRT)
```{r}
# Uses fixed effects (GEM included as batch proxy)
res_edger <- runEDGER_LRT_v1(
  sobj, cluster_id=cluster_col, group_id=group_col,
  covar_effects=fixed_effects, contrast=contrast
)
```

### DESeq2 (LRT)
```{r}
# Uses fixed effects
res_deseq2 <- runDESEQ2_LRT_v1(
  sobj, cluster_id=cluster_col, group_id=group_col,
  covar_effects=fixed_effects, contrast=contrast
)
```

### limma-voom
```{r}
# Uses fixed effects
res_voom <- runLIMMA_voom_v1(
  sobj, cluster_id=cluster_col, group_id=group_col,
  covar_effects=fixed_effects, contrast=contrast
)
```

### limma-trend (muscat-style)
```{r}
# Uses fixed effects
res_trend <- runLIMMA_trend_v1(
  sobj, cluster_id=cluster_col, group_id=group_col,
  covar_effects=fixed_effects, contrast=contrast
)
```

## 2. Mixed Model Methods

### dream (Pseudobulk LMM)
```{r}
# Uses (1|hos_no) random effect. Removing GEM from fixed effects to avoid over-parameterization.
res_dream <- runDREAM_v1(
  sobj, cluster_id=cluster_col, group_id=group_col,
  sample_id=random_effect_var,  # Random Effect ID
  covar_effects=mixed_effects,  # age, sex
  contrast=contrast, use_batch_as_random=TRUE 
)
```

### Nebula (Single-Cell NBMM)
```{r}
# Fast single-cell mixed model
res_nebula <- runNebula_v1(
  sobj, cluster_id=cluster_col, group_id=group_col,
  sample_id=random_effect_var,  # Random Effect
  covar_effects=mixed_effects,  # age, sex
  contrast=contrast
)
```

### MAST (Single-Cell Hurdle LMM)
```{r}
# Single-cell model. Can use random effects (method='glmer') or fixed (method='bayesglm').
# Warning: glmer is slow.
res_mast <- runMAST_v1(
  sobj, cluster_id=cluster_col, group_id=group_col,
  sample_id=random_effect_var,  # Used if random_effect is set
  covar_effects=mixed_effects,
  contrast=contrast,
  random_effect=random_effect_var # Triggers glmer
)
```