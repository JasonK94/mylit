# 🎯 Cursor Refactored Analysis Templates

## Overview

I've created **two comprehensive, well-structured RMarkdown templates** based on your refactored function library. These replace your messy working files with clean, reproducible analysis workflows.

---

## What's New?

### 📊 **stroke_cursor.Rmd** - scRNA-seq Analysis Template
**For:** Stroke project (or any scRNA-seq project)

**Complete workflow includes:**
1. ✅ Data loading from multiple GEMs
2. ✅ SNP/HTO-based demultiplexing (using your refactored `demulti_utils`)
3. ✅ Quality control and filtering
4. ✅ Ambient RNA removal (SoupX)
5. ✅ Integration (Harmony or CCA)
6. ✅ Clustering and cell type annotation
7. ✅ **Differential expression** (pseudobulk, linear regression, LMM)
8. ✅ **Pathway enrichment** (GO, KEGG, GSEA)
9. ✅ **Advanced analyses:**
   - Trajectory inference (Slingshot)
   - Cell-cell communication (NicheNet)
   - Gene signature discovery
10. ✅ Composition analysis and visualization
11. ✅ Automated results export

---

### 🔬 **mibd_cursor.Rmd** - GeoMx Spatial Analysis Template
**For:** mIBD project (or any GeoMx spatial profiling)

**Complete workflow includes:**
1. ✅ GeoMx data loading from Excel
2. ✅ Metadata integration (clinical + spatial)
3. ✅ Quality control (sample and gene level)
4. ✅ Q3 normalization (GeoMx standard)
5. ✅ Exploratory analysis (PCA, UMAP)
6. ✅ **Differential expression:**
   - Simple Wilcoxon tests
   - Limma with patient pairing
   - Segment-specific analysis
   - **Linear mixed models** for complex designs
7. ✅ **Pathway enrichment** (GO, KEGG, GSEA)
8. ✅ **Gene signature analysis:**
   - Module scoring
   - Biomarker discovery
9. ✅ Advanced visualizations
10. ✅ Automated Excel export with multiple sheets

---

## Key Features

### 🎨 **Clean Structure**
- Numbered sections with clear headers
- Logical flow from data loading → QC → analysis → results
- Collapsible code chunks
- Table of contents for easy navigation

### 📚 **Uses All Your Refactored Functions**

#### From `differential_expression.R`:
```r
# Pseudobulk DE
run_pseudobulk_deg(...)

# Linear regression  
linear_seurat(...)

# Linear mixed models
create_analysis_config(...)
run_lmm_multiple_genes(...)
find_response_differential_genes(...)
```

#### From `pathway_enrichment.R`:
```r
# Comprehensive pathway analysis
myGO(
  deg_results = ...,
  run_go = TRUE,
  run_kegg = TRUE,
  run_gsea = TRUE
)
```

#### From `signature_discovery.R`:
```r
# Find biomarkers
find_gene_signature(
  data = sobj,
  method = "lasso"  # or tree_based, limma, nmf, etc.
)
```

#### From `GeoMx.R`:
```r
# GeoMx-specific functions
q3_normalize(...)
find_deg_geomx(...)
plot_deg_volcano(...)
plot_deg_heatmap(...)
```

#### From `demulti_utils.R`:
```r
# Demultiplexing
get_barcode_mapping(...)
summarize_demulti_results(...)
```

#### And many more from your refactored library!

---

## How to Use

### For Stroke Analysis (scRNA-seq):

1. **Open the template:**
   ```r
   # In RStudio
   file.edit("/data/kjc1/mylit/rmd/stroke_cursor.Rmd")
   ```

2. **Customize paths:**
   - Update `project_dir` to your project location
   - Adjust GEM paths if needed
   - Set `output_dir` for results

3. **Run step by step:**
   - Start with Section 1 (Setup)
   - Run each section sequentially
   - Check outputs before proceeding

4. **Adjust parameters:**
   - QC thresholds (nFeature, percent.mt)
   - Integration method (Harmony vs CCA)
   - Clustering resolution
   - Cell type annotations (based on your markers)

5. **Choose DE method:**
   - **Pseudobulk**: For group comparisons with replicates
   - **Linear Regression**: For continuous/ordinal variables
   - **LMM**: For repeated measures or complex designs

6. **Generate report:**
   - Click "Knit to HTML" for a complete report
   - All figures and tables will be included

---

### For mIBD Analysis (GeoMx):

1. **Open the template:**
   ```r
   file.edit("/data/kjc1/mylit/rmd/mibd_cursor.Rmd")
   ```

2. **Customize file paths:**
   - Update `data_file` and `clinical_file`
   - Set `output_dir`

3. **Adjust metadata column names:**
   - The template uses placeholders like `timepoint`, `segment`, `patho_id`
   - Update these to match your actual column names
   - Update the `create_analysis_config()` call for LMM

4. **Run QC and check distributions:**
   - Review QC plots
   - Adjust filtering thresholds if needed

5. **Choose appropriate DE method:**
   - **Simple Wilcoxon**: Quick initial analysis
   - **Limma with pairing**: Accounts for patient matching
   - **Segment-specific**: Analyzes each tissue compartment separately
   - **LMM**: For complex designs (drug × timepoint × response)

6. **Customize gene sets:**
   - Update the `gene_sets` list for module scoring
   - Add your disease-relevant pathways

7. **Generate comprehensive report:**
   - Creates Excel file with all results
   - Multiple sheets for different comparisons

---

## Comparison: Old vs New

### Before (Your Working Files):

**mibd_251001.Rmd:**
```r
# Messy working file
# Korean notes mixed with code
# Functions defined inline
# Hard to follow the flow
# Redundant code chunks
# No clear structure
```

**stroke_total_250703.Rmd:**
```r
# Multiple versions (250517, 250605, 250623, 250703...)
# Experimental code mixed with production code
# Hard to find specific analyses
# Lots of commented-out sections
```

### After (Cursor Templates):

**stroke_cursor.Rmd / mibd_cursor.Rmd:**
```r
✅ Clear numbered sections
✅ English documentation
✅ Uses refactored functions from myR
✅ Logical workflow from start to finish
✅ Best practices throughout
✅ Easy to customize
✅ Publication-ready reports
✅ Comprehensive output (RDS + Excel)
```

---

## Template Sections Explained

### Common Structure (Both Templates):

#### 1. **Setup**
- Library loading
- Path definitions
- Global options

#### 2. **Data Loading**
- Read raw data
- Create Seurat objects
- Initial metadata

#### 3. **Quality Control**
- Calculate QC metrics
- Visualize distributions
- Apply filters

#### 4. **Normalization**
- Method-appropriate normalization
- Validation plots

#### 5. **Exploratory Analysis**
- PCA/UMAP
- Batch effect checking
- Sample clustering

#### 6. **Differential Expression**
- **Multiple methods available:**
  - Pseudobulk (edgeR)
  - Linear regression
  - Linear mixed models
- Appropriate for your experimental design

#### 7. **Pathway Enrichment**
- GO analysis
- KEGG pathways
- GSEA
- Visualization

#### 8. **Advanced Analyses**
- Trajectory (scRNA-seq)
- CCI (scRNA-seq)
- Biomarker discovery
- Gene signatures

#### 9. **Visualization**
- Publication-quality plots
- Heatmaps
- Custom visualizations

#### 10. **Save Results**
- RDS objects
- Excel exports
- Figure files

---

## Customization Tips

### Adjust for Your Data:

#### **Metadata Column Names:**
```r
# Find this in Section 2 and update throughout:
sample_col = "your_sample_column"
group_col = "your_grouping_column"
cluster_col = "your_cluster_column"
```

#### **QC Thresholds:**
```r
# Section 3, adjust based on your data quality:
nFeature_min <- 200      # Minimum genes per cell
nFeature_max <- 6000     # Maximum genes (for doublet filtering)
percent_mt_max <- 20     # Maximum mitochondrial percentage
```

#### **Cell Type Markers:**
```r
# Section 4, update with your markers:
markers_to_plot <- list(
  Your_CellType1 = c("MARKER1", "MARKER2", "MARKER3"),
  Your_CellType2 = c("MARKER4", "MARKER5", "MARKER6")
)
```

#### **Gene Sets for Scoring:**
```r
# Section 8, define your pathways:
gene_sets <- list(
  Your_Pathway1 = c("GENE1", "GENE2", "GENE3"),
  Your_Pathway2 = c("GENE4", "GENE5", "GENE6")
)
```

---

## Decision Guide: Which DE Method?

### Use **Pseudobulk** when:
- ✅ You have biological replicates (multiple samples per condition)
- ✅ Simple group comparisons (treated vs control)
- ✅ Want the gold standard for scRNA-seq DE

```r
results <- run_pseudobulk_deg(
  object = sobj,
  sample_col = "sample_id",
  group_col = "condition",
  comparison = c("Treated", "Control"),
  mode = "per_cluster"
)
```

---

### Use **Linear Regression** when:
- ✅ Continuous predictors (age, time, severity score)
- ✅ Ordinal variables (disease stage: mild/moderate/severe)
- ✅ Need to control for multiple covariates

```r
results <- linear_seurat(
  sobj = sobj,
  regressor = "age",
  regressor.type = "continuous",
  covariates = c("sex", "batch")
)
```

---

### Use **Linear Mixed Models** when:
- ✅ Repeated measures (same patient, pre/post treatment)
- ✅ Hierarchical data (cells within patients)
- ✅ Complex designs (drug × timepoint × response)
- ✅ Need to account for random effects

```r
config <- create_analysis_config(
  patient = "PatientID",
  drug = "Treatment",
  timepoint = "Time",
  response = "Responder"
)

results <- run_lmm_multiple_genes(
  seurat_obj = sobj,
  genes = candidate_genes,
  config = config
)
```

---

## Example Workflows

### Workflow 1: Simple scRNA-seq Analysis

```r
# 1. Load data
sobj <- Load10X(...)

# 2. QC
sobj <- subset(sobj, subset = nFeature_RNA > 200 & percent.mt < 20)

# 3. Normalize & cluster
sobj <- NormalizeData(sobj) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:30) %>%
  FindClusters()

# 4. DE analysis
results <- run_pseudobulk_deg(
  sobj,
  sample_col = "sample",
  group_col = "condition",
  mode = "overall"
)

# 5. Pathway analysis
pathways <- myGO(results, run_go = TRUE, run_gsea = TRUE)
```

---

### Workflow 2: GeoMx with Patient Pairing

```r
# 1. Load & prepare
q3_norm <- q3_normalize(count_matrix)

# 2. Limma with patient pairing
design <- model.matrix(~ 0 + timepoint + patient_id, data = metadata)
fit <- lmFit(log2(q3_norm + 1), design)
# ... (see template for full code)

# 3. LMM for complex design
config <- create_analysis_config(
  patient = "patient_id",
  timepoint = "timepoint",
  response = "response_status"
)

lmm_results <- run_lmm_multiple_genes(sobj, genes, config)

# 4. Find response-associated genes
response_genes <- find_response_differential_genes(
  lmm_results$summary,
  config
)
```

---

## Output Files

### Both templates generate:

1. **RDS files:**
   - Processed Seurat object
   - Normalized data
   - DE results
   - Pathway results
   - All intermediate objects

2. **Excel files:**
   - Multiple sheets for different comparisons
   - Well-formatted with gene names, stats, etc.
   - Easy to share with collaborators

3. **Figures:**
   - UMAP plots
   - Volcano plots
   - Heatmaps
   - Pathway plots
   - All as PDF (high resolution)

4. **HTML report:**
   - Click "Knit to HTML" in RStudio
   - Complete analysis with all code, figures, tables
   - Share with collaborators
   - Publication supplement

---

## Best Practices

### ✅ DO:
1. **Run sections sequentially** - Don't skip QC!
2. **Save intermediate objects** - Use `saveRDS()` after major steps
3. **Document your changes** - Add notes in the template
4. **Check distributions** - Always visualize before filtering
5. **Use appropriate DE methods** - See decision guide above
6. **Adjust p-value thresholds** - Based on your multiple testing burden
7. **Validate key findings** - Check top genes make biological sense

### ❌ DON'T:
1. **Don't run all chunks at once** - Review each step
2. **Don't skip QC plots** - They reveal data quality issues
3. **Don't ignore warnings** - They often indicate real problems
4. **Don't use default thresholds blindly** - Adjust for your data
5. **Don't over-interpret marginal results** - Focus on robust findings

---

## Troubleshooting

### Common Issues:

#### **"Cannot find function X"**
```r
# Make sure myR is loaded:
library(devtools)
setwd("/data/kjc1/mylit/myR")
load_all()
```

#### **"Metadata column not found"**
```r
# Check your column names:
colnames(sobj@meta.data)

# Update template to match your names
sample_col = "your_actual_column_name"
```

#### **"Object X not found"**
```r
# Run sections sequentially
# Make sure previous sections completed successfully
```

#### **LMM won't converge**
```r
# Try:
# 1. Reduce number of genes
# 2. Check for collinearity in design
# 3. Simplify the model formula
```

#### **Out of memory**
```r
# For large datasets:
# 1. Process in batches
# 2. Use fewer genes for LMM
# 3. Increase R memory: R --max-mem-size=64G
```

---

## Next Steps

1. **Try the templates on your current projects:**
   - Stroke: Use `stroke_cursor.Rmd`
   - mIBD: Use `mibd_cursor.Rmd`

2. **Customize for your specific needs:**
   - Update metadata column names
   - Adjust QC thresholds
   - Add project-specific gene sets

3. **Generate reports:**
   - Knit to HTML for quick review
   - Knit to PDF for publication supplements

4. **Share with collaborators:**
   - The HTML reports are self-contained
   - Excel files are ready for non-R users

5. **Iterate:**
   - Keep the messy working files for exploration
   - Copy final analyses into the clean templates
   - Maintain both versions if needed

---

## Template Maintenance

These templates use your **refactored function library**, so:

✅ **When you update functions in myR:**
- The templates automatically use the new versions
- Just reload with `load_all()`

✅ **When you add new functions:**
- Add examples to the appropriate section
- Update this guide

✅ **When you find bugs:**
- Fix in the function library (myR/R/)
- Templates benefit from the fix automatically

---

## Summary

### You now have:
1. ✅ **Two comprehensive analysis templates**
2. ✅ **Clean, reproducible workflows**
3. ✅ **Integration with your refactored functions**
4. ✅ **Publication-ready reports**
5. ✅ **Easy customization for new projects**

### Benefits:
- ⚡ **Faster analysis** - No more copy-pasting messy code
- 🎯 **Consistent** - Same structure for all projects
- 📚 **Documented** - Clear explanations throughout
- 🔬 **Comprehensive** - All analyses in one place
- 🤝 **Shareable** - Easy for collaborators to understand

---

**Happy analyzing! 🚀**

If you need help customizing the templates for specific use cases, just ask!

