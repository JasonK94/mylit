#!/usr/bin/env python3
"""
HC-only scCODA: Compositional differential abundance analysis
Adapted from stroke_hc_v8_2/scripts/claude/13_sccoda_run.py

Usage:
    conda activate scvi-tools
    Rscript scripts/hc/export_for_python.R  # export composition tables first
    python scripts/hc/run_sccoda.py

Output: /data/user3/sobj/hc_only_v1/sccoda/
"""

import os
import sys
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from datetime import datetime

print(f"=== HC-only scCODA Analysis ===")
print(f"Time: {datetime.now()}\n")

BASE = "/data/user3/sobj/hc_only_v1"
OUT = os.path.join(BASE, "sccoda")
os.makedirs(OUT, exist_ok=True)

# Try pertpy
try:
    import pertpy as pt
    print("Using pertpy for scCODA")
except ImportError:
    print("ERROR: pertpy not installed. Run: pip install pertpy")
    sys.exit(1)


def run_sccoda(comp_df, condition_col, cell_type_cols, label, reference_ct):
    """Run scCODA via pertpy (same pattern as stroke analysis)"""
    print(f"\n--- {label} ---")
    print(f"  Condition: {condition_col}")
    print(f"  Patients: {len(comp_df)}, Cell types: {len(cell_type_cols)}")
    print(f"  Reference: {reference_ct}")

    X = comp_df[cell_type_cols].values.astype(np.float64)
    obs = comp_df.drop(columns=cell_type_cols).reset_index(drop=True)
    var = pd.DataFrame(index=cell_type_cols)

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obs[condition_col] = adata.obs[condition_col].astype(str)
    print(f"  Condition levels: {adata.obs[condition_col].value_counts().to_dict()}")

    try:
        sccoda = pt.tl.Sccoda()
        sccoda.prepare(adata, formula=f"C({condition_col})", reference_cell_type=reference_ct)
        sccoda.run_nuts(adata, num_warmup=500, num_samples=10000, rng_key=42)

        credible = sccoda.credible_effects(adata)
        print(f"\n  Credible effects:")
        print(credible)
        credible.to_csv(os.path.join(OUT, f"sccoda_{label}_credible.csv"))

        try:
            summary = sccoda.summary(adata)
            summary.to_csv(os.path.join(OUT, f"sccoda_{label}_summary.csv"))
        except Exception:
            pass

        return adata, sccoda

    except Exception as e:
        print(f"  scCODA failed: {e}")
        # Fallback: Mann-Whitney U
        from scipy import stats
        from statsmodels.stats.multitest import multipletests

        results = []
        ref_vals = adata.obs[condition_col].unique()
        for ct in cell_type_cols:
            idx = cell_type_cols.index(ct)
            groups = [X[adata.obs[condition_col] == cond, idx] for cond in ref_vals]
            if len(groups) == 2:
                stat, pval = stats.mannwhitneyu(groups[0], groups[1], alternative='two-sided')
                fc = np.mean(groups[1]) / max(np.mean(groups[0]), 0.01)
            else:
                stat, pval = stats.kruskal(*groups)
                fc = np.nan
            results.append({'cell_type': ct, 'stat': stat, 'pvalue': pval, 'fold_change': fc})

        result_df = pd.DataFrame(results)
        result_df['fdr'] = multipletests(result_df['pvalue'], method='fdr_bh')[1]
        result_df = result_df.sort_values('pvalue')
        print(f"\n  Fallback (Mann-Whitney):")
        print(result_df.to_string())
        result_df.to_csv(os.path.join(OUT, f"sccoda_{label}_fallback.csv"), index=False)
        return adata, None


# ---- Load data ----
anno1_file = os.path.join(OUT, "composition_anno1.csv")
if not os.path.exists(anno1_file):
    print(f"ERROR: {anno1_file} not found. Run export_for_python.R first.")
    sys.exit(1)

comp_df = pd.read_csv(anno1_file)
meta_cols = ['name', 'GEM', 'age', 'sex', 'age_group']
ct_cols = [c for c in comp_df.columns if c not in meta_cols]
print(f"Loaded: {len(comp_df)} patients, {len(ct_cols)} cell types")

# Reference: most abundant cell type
ref_ct = comp_df[ct_cols].sum().idxmax()
print(f"Reference cell type: {ref_ct}")

# ---- scCODA: Sex (M vs F) ----
try:
    run_sccoda(comp_df, "sex", ct_cols, "sex_M_vs_F", ref_ct)
except Exception as e:
    print(f"Sex scCODA failed: {e}")
    import traceback; traceback.print_exc()

# ---- scCODA: Age group (Young vs Old) ----
try:
    yo_df = comp_df[comp_df['age_group'].isin(['Young', 'Old'])].copy()
    print(f"\nYoung: {(yo_df['age_group']=='Young').sum()}, Old: {(yo_df['age_group']=='Old').sum()}")
    run_sccoda(yo_df, "age_group", ct_cols, "age_Young_vs_Old", ref_ct)
except Exception as e:
    print(f"Age scCODA failed: {e}")
    import traceback; traceback.print_exc()

# ---- scCODA: Age group (3 groups) ----
try:
    run_sccoda(comp_df, "age_group", ct_cols, "age_3groups", ref_ct)
except Exception as e:
    print(f"3-group age scCODA failed: {e}")

print(f"\n=== scCODA Complete: {datetime.now()} ===")
print(f"Output: {OUT}")
