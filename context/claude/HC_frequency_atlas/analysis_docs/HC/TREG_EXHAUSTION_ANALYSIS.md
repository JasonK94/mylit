# Treg Exhaustion Markers × Age — Analysis Report

> Created: 2026-02-18 | Dataset: hc_only_v1 | Focus: CD4 Treg (8,200 cells, 96 patients)

## 1. Motivation

User-specified 23 exhaustion-related markers tested for age association in healthy CD4 Treg cells. The goal was to determine whether Treg exhaustion features change with normal aging.

## 2. Marker List & Availability

### Requested (23)
```
PTGFRN, NR4A3, EGR2, EGR3, NR4A1, FCEB2, ZBED2, OSBPL6, IL1B2, XIBP1,
INS3, IL1RN, DHBS2, ITPKA, ITGB8, GLT1D1, UBE2QL1, PTGS1, GEM, PAPSS2,
CD276, PLEKHA7, SLC16A14
```

### Found in Data (15)
PTGFRN, NR4A3, EGR2, EGR3, NR4A1, ZBED2, IL1RN, ITPKA, ITGB8, GLT1D1, PTGS1, GEM, PAPSS2, PLEKHA7, SLC16A14

### Missing (8)
| Gene | Possible Match | Note |
|---|---|---|
| FCEB2 | — | No match found |
| OSBPL6 | OSBPL1A–11 exist | Family member typo? |
| IL1B2 | IL1B | Trailing digit typo? |
| XIBP1 | XBP1? | Typo? |
| INS3 | GINS3 | Partial match |
| DHBS2 | — | No match |
| UBE2QL1 | — | No match |
| CD276 | — | Not expressed in this dataset |

## 3. Sparsity Analysis (Critical Finding)

**10/15 markers have >70% zero patients at pseudobulk level.**

| Category | Genes | % Zero |
|---|---|---|
| **Well-expressed** | NR4A3, NR4A1 | 0% |
| **Analyzable** | GEM, EGR2, PLEKHA7, EGR3 | 5–57% |
| **Sparse (unanalyzable)** | PTGFRN, ZBED2, IL1RN, ITPKA, ITGB8, GLT1D1, PTGS1, PAPSS2, SLC16A14 | 74–83% |

Patient-level pseudobulk uses raw count sums across median 78 Treg cells per patient. Despite this, most markers remain at zero because they are **not typically expressed in Treg**.

## 4. Age Association Results

### 4.1 Spearman Correlation (pseudobulk logCPM vs age)

**No marker reached nominal significance (p<0.05).**

| Gene | rho | p-value | padj | Direction |
|---|---|---|---|---|
| IL1RN | −0.176 | 0.087 | 0.54 | Down |
| NR4A1 | −0.157 | 0.126 | 0.54 | Down |
| GLT1D1 | −0.153 | 0.136 | 0.54 | Down |
| SLC16A14 | −0.142 | 0.168 | 0.54 | Down |
| GEM | −0.126 | 0.221 | 0.54 | Down |
| ... | | | | |

### 4.2 limma-voom (covariate-controlled, `~ age + sex`)

**1 marker nominally significant (FDR not significant):**

| Gene | logFC/yr | p-value | adj.P | Direction |
|---|---|---|---|---|
| **ITGB8** | +0.020 | **0.006** | 0.55 | Up with age |
| PTGFRN | +0.013 | 0.088 | 0.74 | Up |
| NR4A1 | −0.011 | 0.175 | 0.80 | Down |

### 4.3 Age Group Stratification (Young/Middle/Old)

Kruskal-Wallis test across 3 age groups: **All 15 markers NS (p>0.05).**

### 4.4 Sex Stratification

| Gene | Female rho (n=57) | Female p | Male rho (n=39) | Male p |
|---|---|---|---|---|
| **ITGB8** | −0.072 | 0.60 | **+0.404** | **0.011** |
| All others | NS | | NS | |

ITGB8 shows age association **only in males** — but this gene is 77% zero at pseudobulk, so the signal comes from very few non-zero observations.

### 4.5 Direction Concordance

6/15 markers are **discordant** between Spearman and limma:
- ZBED2, ITPKA, GLT1D1, PTGS1, PAPSS2, SLC16A14

This high discordance rate confirms the signal is **indistinguishable from noise**.

## 5. Treg Subclustering

### 7 Subclusters (scVI, res=0.5)

| Cluster | Cells | Top Markers | Character |
|---|---|---|---|
| 0 | 2,856 | FOXP3, RTKN2, IL2RA, CCR6 | Classic Treg |
| 1 | 2,447 | ARMH1, NOSIP, KLRB1, TSHZ2 | Cytotoxic-like |
| 2 | 748 | PLK1, DLGAP5, CDC20 | Proliferating |
| 3 | 621 | IGFBP7, NCAM1, NCR1 | NK-like / innate |
| 4 | 605 | FXYD2, IL10, SESTD1 | IL10+ suppressive |
| 5 | 553 | CD8A, CD8B, CCL5 | CD8 contaminant? |
| 6 | 370 | LINC02694, FTX, IKZF2 | lncRNA-hi |

### Per-Subcluster Age Association

54 tests (gene × subcluster), **1 nominally significant**:
- Cluster 1 × GEM: rho=−0.212, p=0.039, padj=0.79

After FDR correction: **zero significant associations**.

## 6. FGS Genome-Wide Ranking

None of the 15 exhaustion markers appeared in the genome-wide FGS top-200 (out of 18,182 genes). Their mean rank positions ranged from 6,400 to 12,100.

## 7. Conclusions

1. **These exhaustion markers have no meaningful age association in healthy Treg**
2. Most markers are **too sparse** in Treg for pseudobulk analysis (>70% zero)
3. The marker set likely derives from **CD8 exhaustion or disease contexts** (chronic antigen, tumor), not healthy Treg aging
4. Subclustering does not rescue the signal
5. Sex stratification reveals ITGB8 in males as suggestive but unreliable (sparse data)
6. **In contrast**, the FGS pipeline identified 84 age-associated genes in Treg at 3+ method consensus — the real age-signal involves ISG genes (PLAC8, EPSTI1), JAML, MSC, and DGKA

## 8. Output Files

```
/data/user3/sobj/hc_only_v1/treg_exhaustion_age/
├── exhaustion_markers_age_stats.csv           # Combined Spearman + limma stats
├── treg_subclustered.qs                       # 7-cluster Treg object (107 MB)
├── treg_subcluster_markers.csv                # FindAllMarkers per subcluster
├── subcluster_age_associations.csv            # Per-subcluster age Spearman
├── scatter_all_markers_vs_age_white.png       # All 15 markers scatter (white BG)
├── scatter_markers_sex_stratified.png         # Sex-stratified scatter
├── boxplot_markers_by_age_group.png           # Young/Middle/Old boxplots
├── heatmap_exhaustion_markers_age_white.png   # Z-scored heatmap by age
├── dotplot_markers_rho_white.png              # Spearman rho dot plot
├── dotplot_markers_by_subcluster.png          # Expression across subclusters
├── barplot_subcluster_by_age_group.png        # Subcluster composition × age
└── barplot_subcluster_by_sex.png              # Subcluster composition × sex
```
