#!/usr/bin/env python3
"""
myeloid_v1_TNBC.h5ad 파일 구조를 상세히 확인하는 스크립트
"""

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import scipy.sparse

# h5ad 파일 로드
h5ad_path = "/data/ARPAH/250305_TNBC_26/Xenium_analysis/transfer_to_python/TNBC_only/myeloid_v1_TNBC.h5ad"
print(f"Loading h5ad file: {h5ad_path}")
adata = ad.read_h5ad(h5ad_path)

print("=" * 80)
print("BASIC INFO")
print("=" * 80)
print(f"Shape: {adata.shape} (n_obs × n_vars)")
print(f"n_obs (cells): {adata.n_obs}")
print(f"n_vars (genes): {adata.n_vars}")

print("\n" + "=" * 80)
print("OBS (Cell Metadata)")
print("=" * 80)
print(f"Shape: {adata.obs.shape}")
print(f"Columns: {list(adata.obs.columns)}")
print(f"\nData types:")
print(adata.obs.dtypes)

# 필수 컬럼 확인
required_cols = ['batch', 'Core', 'intra_T', 'TIL_2', 'Annot3']
print(f"\n=== Required Columns Check ===")
for col in required_cols:
    if col in adata.obs.columns:
        print(f"✓ {col}: Found")
        print(f"  Unique values: {sorted(adata.obs[col].unique())}")
        print(f"  Value counts:\n{adata.obs[col].value_counts()}")
    else:
        print(f"✗ {col}: MISSING")

print(f"\n=== Sample values for key columns ===")
for col in required_cols + ['leiden']:
    if col in adata.obs.columns:
        print(f"\n{col}:")
        print(f"  Unique values: {sorted(adata.obs[col].unique())}")
        counts = adata.obs[col].value_counts()
        print(f"  Value counts:\n{counts}")

# intra_T와 TIL_2의 low/high 분포 확인
print("\n" + "=" * 80)
print("COMPARISON GROUP DISTRIBUTION")
print("=" * 80)
if 'intra_T' in adata.obs.columns:
    print("\nintra_T distribution:")
    print(adata.obs['intra_T'].value_counts())
    print("\nintra_T by Core (sample):")
    print(pd.crosstab(adata.obs['Core'], adata.obs['intra_T']))

if 'TIL_2' in adata.obs.columns:
    print("\nTIL_2 distribution:")
    print(adata.obs['TIL_2'].value_counts())
    print("\nTIL_2 by Core (sample):")
    print(pd.crosstab(adata.obs['Core'], adata.obs['TIL_2']))

print("\n" + "=" * 80)
print("OBSM (Cell Embeddings/Reductions)")
print("=" * 80)
if adata.obsm:
    print(f"Available keys: {list(adata.obsm.keys())}")
    for key in adata.obsm.keys():
        emb = adata.obsm[key]
        print(f"\n{key}:")
        print(f"  Shape: {emb.shape}")
        print(f"  Data type: {emb.dtype}")
        print(f"  Contains NaN: {np.isnan(emb).any()}")
        print(f"  Contains Inf: {np.isinf(emb).any()}")
else:
    print("No obsm found")

print("\n" + "=" * 80)
print("LAYERS")
print("=" * 80)
if adata.layers:
    print(f"Available layers: {list(adata.layers.keys())}")
else:
    print("No layers found")

print("\n" + "=" * 80)
print("SUMMARY FOR MILO ANALYSIS")
print("=" * 80)
print(f"Total cells: {adata.n_obs}")
print(f"Total genes: {adata.n_vars}")
if 'Core' in adata.obs.columns:
    print(f"Number of samples (Core): {adata.obs['Core'].nunique()}")
    print(f"Samples: {sorted(adata.obs['Core'].unique())}")
if 'batch' in adata.obs.columns:
    print(f"Number of batches: {adata.obs['batch'].nunique()}")
    print(f"Batches: {sorted(adata.obs['batch'].unique())}")
if 'Annot3' in adata.obs.columns:
    print(f"Number of clusters (Annot3): {adata.obs['Annot3'].nunique()}")
    print(f"Clusters: {sorted(adata.obs['Annot3'].unique())}")

# Reduction 확인
print("\nAvailable reductions (obsm keys):")
if adata.obsm:
    print(f"  {list(adata.obsm.keys())}")
    if 'X_integrated_scvi' in adata.obsm.keys():
        print("  ✓ X_integrated_scvi found (can be used as integrated.scvi)")
    if 'X_umap' in adata.obsm.keys() or 'X_umap_scvi' in adata.obsm.keys():
        print("  ✓ UMAP found")
else:
    print("  None found")

