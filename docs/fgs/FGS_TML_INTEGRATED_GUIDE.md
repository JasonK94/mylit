# FGS & TML Integrated Guide

## 1. Introduction
This document provides a comprehensive guide to the **Find Gene Signature (FGS)** and **Targeted Meta Learner (TML)** frameworks, including their methodology, usage, and recent improvements.

- **FGS (L1 Layer)**: Identifies gene signatures that distinguish sample groups (e.g., cell clusters) using various feature selection methods.
- **TML (L2 Layer)**: An ensemble meta-learner that combines L1 signatures to predict a target variable, improving robustness and accuracy.
- **CMGI (Compute Meta Gene Importance)**: A method to derive gene-level importance from the TML model, explaining *which* genes drive the prediction.
- **AMSC (Aggregated Meta Signature Contribution)**: The specific metric used in CMGI.

## 2. Workflow Visualization

```mermaid
flowchart TD
    subgraph FGS["FGS (L1 Layer)"]
        direction TB
        Input[(Count Matrix<br/>Genes x Cells)] --> Preprocess[Preprocessing<br/>(Filter, Norm, Scale)]
        Preprocess --> Methods{Select Methods}
        Methods -->|Tree| RF[Random Forest / Ranger]
        Methods -->|Reg| Lasso[Lasso / Ridge]
        Methods -->|DimRed| NMF[NMF / PCA]
        RF --> Sig1[Signature 1]
        Lasso --> Sig2[Signature 2]
        NMF --> Sig3[Signature 3]
        Sig1 & Sig2 & Sig3 --> L1Sigs[(L1 Signatures)]
    end

    subgraph TML["TML (L2 Layer)"]
        direction TB
        L1Sigs --> Score[Score Signatures<br/>on Holdout Data]
        Score --> L2Input[(L2 Input Matrix)]
        L2Input --> MetaLearner{Meta Learner<br/>(GLM, Ranger, Earth)}
        MetaLearner --> TrainedModel[Trained L2 Model]
        TrainedModel --> CMGI[CMGI<br/>(Feature Importance)]
        CMGI --> GeneImp[(Gene Importance<br/>AMSC Score)]
    end

    FGS --> TML
```

---

## 3. Development Log & Recent Improvements

### Version 5.4 (Current)
*   **L1 Method Fixes**:
    *   **`random_forest_ranger`**: Fixed name mapping issue where importance scores were assigned to wrong genes. Verified to match `random_forest` results (Jaccard Index 1.0).
    *   **`nmf_loadings`**: Fixed `do.call` scope issue by explicitly passing the `NMF::nmf` function.
*   **TML & CMGI Improvements**:
    *   **L2 Model Support**: Fixed `compute_meta_gene_importance` for `ranger` (added importance extraction) and `earth` (added package dependency).
    *   **Normalization**: Implemented model-specific normalization (`max_abs`, `min_max`, `softmax`, etc.) in TML7 to handle different scales of L1 weights.
    *   **Robustness**: Improved name matching between L1 signatures and L2 importance names (handling `make.names` and backticks).

### Key Commit History
*   `Fix compute_meta_gene_importance function definition syntax` (Recent)
*   `Integrate TML7 Phase 1 (Normalization) into signature.R`
*   `Update progress report and scripts (Ranger vs RF verified)`

---

## 3. User Guide & Best Practices

### ⚠️ Critical Warnings

1.  **CPU Resource Management (Taskset)**
    *   **Issue**: Many L1/L2 methods (e.g., `xgboost`, `ranger`) use C++ backends that may ignore R-level CPU limits or OpenMP settings. This can lead to CPU saturation and system unresponsiveness.
    *   **Solution**: **ALWAYS** use `taskset` when running scripts to strictly limit CPU core usage at the OS level.
    *   **Example**:
        ```bash
        # Limit to cores 0-15 (16 cores)
        taskset -c 0-15 Rscript run_analysis.R
        ```

2.  **Environment Management (renv)**
    *   **Issue**: Missing or conflicting package versions can break the pipeline (e.g., `earth` package missing).
    *   **Solution**: Always activate the project environment before running analysis.
    *   **Example**:
        ```r
        renv::activate("/home/user3/GJC_KDW_250721")
        ```

### Usage Guide

#### FGS (Find Gene Signature)
Input: Count matrix (Genes x Cells), Sample ID, Target Variable (e.g., 'g3'), Control Vars.
Output: A list of gene signatures (named vectors of weights) for each method.

```r
fgsa <- find_gene_signature_v5.4(
  sobj,
  target_var = "g3",
  control_vars = "hos_no",
  n_features = 200,
  method = c("random_forest_ranger", "lasso", "nmf_loadings", ...)
)
```

#### TML (Targeted Meta Learner)
Input: FGS result (`l1_signatures`), Holdout Data.
Output: Trained meta-model, L2 importance, and performance metrics.

```r
tmla <- TML7(
  l1_signatures = fgsa,
  holdout_data = sobj,
  target_var = "g3"
)
```

#### CMGI (Compute Meta Gene Importance)
Input: TML result.
Output: Gene-level contribution scores.

```r
# Default (max_abs normalization)
cmgi_res <- compute_meta_gene_importance(tmla)

# With specific normalization (e.g., for probability interpretation)
cmgi_res_soft <- compute_meta_gene_importance(tmla, normalization_method = "softmax")
```

---

## 4. Methodology Details

### CMGI & AMSC
**CMGI** calculates the contribution of each gene to the final prediction by propagating importance back from the L2 model to the L1 signatures and then to the genes.

**AMSC (Aggregated Meta Signature Contribution)** is the metric:
$$ \text{AMSC}_g = \sum_{s \in S} (\text{Imp}_s \times w_{g,s}) $$
Where:
*   $\text{Imp}_s$: Importance of signature $s$ in the L2 model (normalized).
*   $w_{g,s}$: Weight of gene $g$ in signature $s$.

**Reliability**:
*   Provides a transparent explanation of the "Black Box" ensemble.
*   High AMSC indicates a gene is consistently important across multiple high-performing L1 methods.

---

## 5. Appendix: Methods Reference

### FGS Methods (L1 Layer)

| Method | Description | Pros | Cons |
| :--- | :--- | :--- | :--- |
| **random_forest** | Traditional Random Forest (randomForest pkg). | Robust, non-linear. | Slow on large data. |
| **random_forest_ranger** | Fast Random Forest implementation (ranger pkg). | **Very Fast** (C++), parallel. | - |
| **xgboost** | Gradient Boosting Trees. | High performance. | Tuning required. |
| **lasso** | L1-regularized Logistic Regression. | Sparse (feature selection). | Linear only. |
| **ridge** | L2-regularized Logistic Regression. | Handles multicollinearity. | Keeps all features. |
| **elastic_net** | Combination of Lasso and Ridge. | Balanced. | - |
| **pca_loadings** | Principal Component Analysis loadings. | Unsupervised, captures variance. | Not target-specific. |
| **nmf_loadings** | Non-negative Matrix Factorization. | Parts-based representation. | Slow. |
| **gam** | Generalized Additive Models. | Non-linear relationships. | Computationally expensive. |
| **limma** | Linear Models for Microarray/RNA-seq. | Standard for DE analysis. | Linear. |
| **wilcoxon** | Wilcoxon Rank Sum Test. | Simple, non-parametric. | Univariate only. |

### TML Methods (L2 Layer)

| Method | Description | Key Features |
| :--- | :--- | :--- |
| **glm** | Logistic Regression. | Interpretable (Log-odds). Baseline. |
| **ranger** | Random Forest (ranger). | Non-linear, captures interactions. |
| **svmRadial** | Support Vector Machine (RBF kernel). | Complex boundaries. |
| **xgbTree** | Gradient Boosting. | High accuracy. |
| **earth** | Multivariate Adaptive Regression Splines (MARS). | Captures non-linearities and interactions explicitly. |
| **nnet** | Neural Network. | Non-linear. |

---

## 6. Helper Functions Assessment
*   **Current Status**: Helper functions (e.g., `standardise_signature`) are embedded within `signature.R`.
*   **Assessment**: While functional, extracting them to a dedicated `utils_signature.R` would improve code readability and testability, especially as the logic for normalization and name mapping grows more complex.
