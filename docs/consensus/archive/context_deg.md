# Multi-Model DEG Consensus Engine — Specification Document

## Overview
Apply multiple DEG models (limma/voom/dream/edgeR/DESeq2/muscat) to the same dataset, extract gene-level statistics, integrate outputs using unsupervised learning, and produce a final consensus DEG signature.

## parameters
- test data = /data/user3/sobj/IS_scvi_251107_ds2500.qs
- data = /data/user3/sobj/IS_scvi_251107.qs
- format: Seurat object
- meta.data specifications
* sample: $hos_no
* target variable: $g3 (1, 2, "NA"; character, but numeric after removal of NA)
* batch: $GEM
* confounding: $GEM perfectly include $hos_no, and some $g3 contains whole $GEM==n
* biological factor: sex

## Input Specification
- counts: raw count matrix (genes × samples)
- meta: sample metadata
- contrast: contrast string
- models: list of DEG engines

## Supported DEG Engines
- limma
- voom
- dream
- limma-wt
- limma-trend
- edgeR-LRT
- edgeR-QLF
- edgeR-robust
- DESeq2-Wald
- DESeq2-LRT
- muscat
- nebula
- mast (if possible)

## Per-Model Statistical Extraction
Extract per gene and model:
- Effect size (beta)
- Standard error
- Test statistic
- p-value and FDR
- Optional weights/dispersion

## Matrix Construction
- Effect size matrix B[g, m]
- P-value matrix P[g, m]
- -log10(p) matrix
- Significance indicator S[g, m]

## Consensus Analysis
Agreement score:
A_g = mean(S[g, m] across models)

## Unsupervised Integration
Using beta, -log10(p), and S:
- PCA
- NMF
- UMAP / Diffusion map
- Clustering (k-means, hierarchical, HDBSCAN)

## Outputs
- Per-model DEG tables
- Gene × model matrices
- Consensus DEG list
- PCA/NMF/UMAP plots
- Stability metrics

## Constraints
- All models must use same genes and contrast
- Pseudobulk required for scRNA
- dream requires random effect variables

## Goal
Produce a robust consensus DEG signature from multiple DEG engines.
