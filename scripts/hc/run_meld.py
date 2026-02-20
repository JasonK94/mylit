#!/usr/bin/env python3
"""
HC-only MELD: Per-cell condition likelihood for age and sex
Adapted from stroke_hc_v8_2/scripts/claude/15_meld.py (CSV-based approach)

Usage:
    conda activate scvi-tools
    Rscript scripts/hc/export_for_python.R  # export embeddings first
    python scripts/hc/run_meld.py

Output: /data/user3/sobj/hc_only_v1/meld/
"""

import os
import sys
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scanpy as sc
import meld
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from datetime import datetime

print(f"=== HC-only MELD Analysis ===")
print(f"Time: {datetime.now()}\n")

BASE = "/data/user3/sobj/hc_only_v1"
OUT = os.path.join(BASE, "meld")
os.makedirs(OUT, exist_ok=True)

# ---- Load from CSV exports (proven pattern from stroke MELD) ----
emb_path = os.path.join(OUT, "scvi_embeddings.csv")
umap_path = os.path.join(OUT, "umap_embeddings.csv")
meta_path = os.path.join(OUT, "metadata.csv")

if not os.path.exists(emb_path):
    print(f"ERROR: {emb_path} not found. Run export_for_python.R first:")
    print(f"  Rscript scripts/hc/export_for_python.R")
    sys.exit(1)

print("Loading from CSV exports...")
scvi_emb = pd.read_csv(emb_path, index_col=0)
meta = pd.read_csv(meta_path, index_col=0)
print(f"  Cells: {len(meta)}, scVI dims: {scvi_emb.shape[1]}")

# Create AnnData
adata = sc.AnnData(X=np.zeros((len(meta), 1)), obs=meta)
adata.obsm['X_scVI'] = scvi_emb.values

if os.path.exists(umap_path):
    umap_emb = pd.read_csv(umap_path, index_col=0)
    adata.obsm['X_umap'] = umap_emb.values

# ---- Build KNN graph ----
print("\nBuilding KNN graph from scVI embeddings...")
sc.pp.neighbors(adata, use_rep='X_scVI', n_neighbors=15, method='umap')

# ---- MELD: Sex (M vs F) ----
print("\n--- MELD: Sex (M vs F) ---")
sample_labels = adata.obs['name'].astype(str).values
sample_sex_map = dict(zip(adata.obs['name'].astype(str), adata.obs['sex'].astype(str)))

meld_op = meld.MELD(beta=67)
meld_op.fit(adata.obsp['connectivities'] + adata.obsp['connectivities'].T)

sample_densities = meld_op.transform(sample_labels)
print(f"  Sample densities shape: {sample_densities.shape}")

unique_samples = np.unique(sample_labels)
m_samples = [s for s in unique_samples if sample_sex_map.get(s) == 'M']
f_samples = [s for s in unique_samples if sample_sex_map.get(s) == 'F']
print(f"  M samples: {len(m_samples)}, F samples: {len(f_samples)}")

m_cols = [list(unique_samples).index(s) for s in m_samples]
f_cols = [list(unique_samples).index(s) for s in f_samples]

if isinstance(sample_densities, pd.DataFrame):
    m_density = sample_densities.iloc[:, m_cols].mean(axis=1).values
    f_density = sample_densities.iloc[:, f_cols].mean(axis=1).values
else:
    m_density = np.mean(sample_densities[:, m_cols], axis=1)
    f_density = np.mean(sample_densities[:, f_cols], axis=1)

total = m_density + f_density
total[total == 0] = 1
adata.obs['meld_male_likelihood'] = m_density / total

print(f"  Male likelihood: mean={adata.obs['meld_male_likelihood'].mean():.3f}, "
      f"std={adata.obs['meld_male_likelihood'].std():.3f}")

# ---- MELD: Age (Young vs Old) ----
print("\n--- MELD: Age (Young vs Old) ---")
sample_age_map = dict(zip(adata.obs['name'].astype(str), adata.obs['age_group'].astype(str)))

young_samples = [s for s in unique_samples if sample_age_map.get(s) == 'Young']
old_samples = [s for s in unique_samples if sample_age_map.get(s) == 'Old']
middle_samples = [s for s in unique_samples if sample_age_map.get(s) == 'Middle']
print(f"  Young: {len(young_samples)}, Middle: {len(middle_samples)}, Old: {len(old_samples)}")

young_cols = [list(unique_samples).index(s) for s in young_samples]
old_cols = [list(unique_samples).index(s) for s in old_samples]
middle_cols = [list(unique_samples).index(s) for s in middle_samples]

if isinstance(sample_densities, pd.DataFrame):
    young_density = sample_densities.iloc[:, young_cols].mean(axis=1).values
    old_density = sample_densities.iloc[:, old_cols].mean(axis=1).values
    middle_density = sample_densities.iloc[:, middle_cols].mean(axis=1).values
else:
    young_density = np.mean(sample_densities[:, young_cols], axis=1)
    old_density = np.mean(sample_densities[:, old_cols], axis=1)
    middle_density = np.mean(sample_densities[:, middle_cols], axis=1)

total_age = young_density + old_density + middle_density
total_age[total_age == 0] = 1
adata.obs['meld_old_likelihood'] = old_density / total_age
adata.obs['meld_young_likelihood'] = young_density / total_age
adata.obs['meld_age_gradient'] = (old_density - young_density) / total_age

print(f"  Old likelihood: mean={adata.obs['meld_old_likelihood'].mean():.3f}")
print(f"  Young likelihood: mean={adata.obs['meld_young_likelihood'].mean():.3f}")
print(f"  Age gradient: mean={adata.obs['meld_age_gradient'].mean():.3f}")

# ---- Per cell type summary ----
print("\n--- MELD Summary by Cell Type ---")
meld_by_ct = adata.obs.groupby('anno1').agg(
    mean_male=('meld_male_likelihood', 'mean'),
    mean_old=('meld_old_likelihood', 'mean'),
    mean_age_gradient=('meld_age_gradient', 'mean'),
    n_cells=('meld_male_likelihood', 'count')
).sort_values('mean_age_gradient')
meld_by_ct.to_csv(os.path.join(OUT, "meld_summary_by_celltype.csv"))
print("\n  Age gradient by cell type:")
print(meld_by_ct[['mean_age_gradient', 'mean_male', 'n_cells']].to_string())

# ---- Save scores ----
meld_cols = [c for c in adata.obs.columns if c.startswith('meld_')]
adata.obs[meld_cols + ['anno1', 'anno2', 'name']].to_csv(os.path.join(OUT, "meld_scores_all.csv"))

# ---- Plots ----
print("\nGenerating plots...")
if 'X_umap' in adata.obsm:
    # Age MELD
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    for ag, color in [('Young', 'blue'), ('Middle', 'grey'), ('Old', 'red')]:
        mask = adata.obs['age_group'] == ag
        axes[0].scatter(adata.obsm['X_umap'][mask, 0], adata.obsm['X_umap'][mask, 1],
                       s=0.1, alpha=0.2, label=ag, c=color)
    axes[0].set_title('Age Group'); axes[0].legend(markerscale=10)

    sc_plot = axes[1].scatter(adata.obsm['X_umap'][:, 0], adata.obsm['X_umap'][:, 1],
                              c=adata.obs['meld_age_gradient'], cmap='RdBu_r',
                              s=0.1, alpha=0.3, vmin=-0.5, vmax=0.5)
    plt.colorbar(sc_plot, ax=axes[1], label='Age gradient (Old-Young)')
    axes[1].set_title('MELD: Age Gradient')

    ct_meld = adata.obs.groupby('anno1')['meld_age_gradient'].mean().sort_values()
    colors = ['steelblue' if v < 0 else 'firebrick' for v in ct_meld.values]
    axes[2].barh(range(len(ct_meld)), ct_meld.values, color=colors)
    axes[2].set_yticks(range(len(ct_meld)))
    axes[2].set_yticklabels(ct_meld.index, fontsize=7)
    axes[2].set_xlabel('Mean Age Gradient'); axes[2].set_title('MELD by Cell Type (Age)')
    axes[2].axvline(0, color='grey', linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "meld_age_umap.png"), dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()
    print("  Saved: meld_age_umap.png")

    # Sex MELD
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    for sx, color in [('F', 'coral'), ('M', 'steelblue')]:
        mask = adata.obs['sex'] == sx
        axes[0].scatter(adata.obsm['X_umap'][mask, 0], adata.obsm['X_umap'][mask, 1],
                       s=0.1, alpha=0.2, label=sx, c=color)
    axes[0].set_title('Sex'); axes[0].legend(markerscale=10)

    sc_plot = axes[1].scatter(adata.obsm['X_umap'][:, 0], adata.obsm['X_umap'][:, 1],
                              c=adata.obs['meld_male_likelihood'], cmap='RdBu_r',
                              s=0.1, alpha=0.3, vmin=0.2, vmax=0.8)
    plt.colorbar(sc_plot, ax=axes[1], label='P(Male)')
    axes[1].set_title('MELD: Male Likelihood')

    ct_meld_sex = adata.obs.groupby('anno1')['meld_male_likelihood'].mean().sort_values()
    axes[2].barh(range(len(ct_meld_sex)), ct_meld_sex.values - 0.5)
    axes[2].set_yticks(range(len(ct_meld_sex)))
    axes[2].set_yticklabels(ct_meld_sex.index, fontsize=7)
    axes[2].set_xlabel('Mean P(Male) - 0.5'); axes[2].set_title('MELD by Cell Type (Sex)')
    axes[2].axvline(0, color='grey', linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "meld_sex_umap.png"), dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()
    print("  Saved: meld_sex_umap.png")

print(f"\n=== MELD Complete: {datetime.now()} ===")
print(f"Output: {OUT}")
