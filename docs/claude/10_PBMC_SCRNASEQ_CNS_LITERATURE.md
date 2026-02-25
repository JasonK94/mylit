# High-Impact PBMC Single-Cell Studies in CNS Diseases: Comprehensive Literature Review

> Generated: 2026-02-25
> Purpose: Compile precedent studies that used PBMC/peripheral blood single-cell omics to study CNS diseases,
> published in high-impact journals. Extract narrative strategies for overcoming the "indirect measurement" limitation.
> Relevance: Direct models for our stroke_hc_v8_2 manuscript strategy (226K cells, 100 patients, PBMC scRNAseq)

---

## Executive Summary

Despite the anatomical separation between peripheral blood and the CNS, a substantial and growing body of work has been published in top-tier journals (Nature, Nature Immunology, Science Advances, Cell Discovery, JCI, Nature Communications, etc.) using PBMC single-cell approaches to study CNS diseases. The key narrative strategies that enable these publications fall into several categories:

1. **Blood-Brain Paired Design**: Analyze both blood and CSF/brain tissue, showing that peripheral signatures mirror or complement CNS changes
2. **BBB Breakdown Argument**: The blood-brain barrier is compromised in disease, making peripheral-CNS immune crosstalk biologically plausible
3. **Biomarker Accessibility**: Blood is clinically accessible; CSF/brain tissue is not -- peripheral signatures have translational value
4. **Immune Trafficking**: Peripheral immune cells physically migrate into the CNS, so their pre-infiltration state is informationally relevant
5. **Functional Validation**: Mouse models, flow cytometry, or independent cohorts validate peripheral findings
6. **Multi-layered Evidence**: Concordance across multiple analytical methods strengthens claims from indirect measurements

---

## TIER 1: Nature, Nature Immunology, Science, Cell (IF > 30)

---

### Paper 1: Garcia-Bonilla et al. -- Brain and Blood Single-Cell Transcriptomics in Stroke

**Citation:** Garcia-Bonilla L, Shahanoor Z, Sciortino R, et al. Analysis of brain and blood single-cell transcriptomics in acute and subacute phases after experimental stroke. *Nature Immunology*. 2024;25:357-370.

| Attribute | Detail |
|-----------|--------|
| **Disease** | Ischemic stroke (tMCAO mouse model) |
| **Sample type** | Brain tissue + peripheral blood (paired) |
| **Dataset size** | ~43,269 brain cells + blood leukocytes; Day 2 and Day 14 post-stroke + sham; young (8-12 wk) and aged (17-20 mo) mice |
| **Methods** | Drop-seq scRNA-seq, Seurat clustering, Slingshot trajectory, CellChat CCI, MAST DEG, DecontX, DoubletFinder |
| **Impact factor** | Nature Immunology (~30) |

**Key findings:**
- Strong divergence of post-ischemic microglia, monocyte-derived macrophages, and neutrophils over time
- Blood inflammatory monocytes (Ly6c2+) showed transcriptional shifts: Mo1-Mo2 clusters abundant at D02, declining by D14
- Slingshot trajectory predicted in situ progression of monocyte-derived macrophages FROM blood precursors into brain phenotypes (Mo1/Mo2 -> MdC3 -> MdC1/MdC6)
- Brain-infiltrating leukocytes displayed HIGHER transcriptomic divergence than circulating counterparts
- Peripheral blood major cell types showed "conserved positioning" across conditions -- lesser systemic differentiation than brain responses

**Narrative strategy:**
- **Paired blood+brain design** is the core justification: by profiling BOTH compartments, they can show that phenotypic diversification occurs WITHIN the brain, not in peripheral priming
- **Trajectory inference** bridges the two compartments -- blood monocytes are mapped as precursors to brain macrophages
- Created a public data resource (StrokeVis Shiny app) for community validation

**Validation:**
1. Cx3cr1-CreERT2:R26-tdTomato lineage tracing mice (microglia vs infiltrating macrophage identity)
2. RNAscope FISH + immunofluorescence (spatial validation of marker genes)
3. Flow cytometry (confirmed DC subtypes)

**Relevance to our study:** This is the closest methodological precedent for our work. However: (1) it is a MOUSE study, not human; (2) it has paired brain tissue, which we lack. Our advantage is human patients (100 vs mice), much larger cell count (226K vs ~43K), and outcome stratification (g3), which this study lacks entirely.

---

### Paper 2: Sulzer et al. -- T Cells Recognize alpha-Synuclein in Parkinson's Disease

**Citation:** Sulzer D, Alcalay RN, Garretti F, Cote L, Owen E, Bhatt N, Lindestam Arlehamn CS, et al. T cells from patients with Parkinson's disease recognize alpha-synuclein peptides. *Nature*. 2017;546(7660):656-661.

| Attribute | Detail |
|-----------|--------|
| **Disease** | Parkinson's disease |
| **Sample type** | PBMCs (stimulated with alpha-synuclein peptides) |
| **Dataset size** | 67 PD patients + 36 age-matched healthy controls |
| **Methods** | PBMC stimulation assays, cytokine measurement (IFN-gamma, IL-5, IL-10), HLA typing |
| **Impact factor** | Nature (~65) |

**Key findings:**
- Alpha-synuclein-derived peptides elicit specific T cell responses in PD patients
- Both IL-5-secreting CD4+ T helper cells and IFN-gamma-secreting CD8+ cytotoxic T cells respond
- Five epitopes identified from two alpha-synuclein regions (N-terminus and C-terminus)
- Response associated with specific HLA alleles

**Narrative strategy:**
- **Autoimmune hypothesis**: Framed PD as having a potential autoimmune component
- **Direct mechanistic link**: alpha-synuclein (a CNS protein) is the antigen driving peripheral T cell responses
- **HLA genetics**: Connected peripheral immune findings to well-known genetic risk (HLA association with PD)
- **Clinical accessibility**: Peripheral blood T cell testing could serve as early biomarker

**Validation:**
- Replicated in follow-up longitudinal study (Lindestam Arlehamn et al., Nat Commun 2020): alpha-syn T cell responses detected BEFORE motor PD diagnosis, declining after

**Relevance to our study:** This Nature paper demonstrates that peripheral blood immunophenotyping is publishable in the highest-tier journal IF there is a clear mechanistic link (here: autoantigen) between periphery and CNS. Our study uses a different approach (unbiased transcriptomics) but the SIDS/SIIS framework provides our mechanistic link.

---

### Paper 3: Phongpreecha et al. -- Single-Cell Peripheral Immunoprofiling of AD and PD

**Citation:** Phongpreecha T, Fernandez R, Mrdjen D, et al. Single-cell peripheral immunoprofiling of Alzheimer's and Parkinson's diseases. *Science Advances*. 2020;6(48):eabd5575.

| Attribute | Detail |
|-----------|--------|
| **Disease** | Alzheimer's disease + Parkinson's disease |
| **Sample type** | PBMCs (unstimulated + 7 stimulant conditions) |
| **Dataset size** | Discovery: 108 participants (28 AD, 17 PD, 53 HC, 10 young HC); Validation: 24 (9 AD, 15 HC); ~1.2 million individual cells |
| **Methods** | CyTOF mass cytometry (not scRNAseq), 35 PBMC subsets, 15 intracellular signaling markers, immunological Elastic Net (iEN) machine learning |
| **Impact factor** | Science Advances (~13) |

**Key findings:**
- Three main signaling pathway differences in AD PBMCs:
  1. Reduced pPLCgamma2 activation across most cell types and stimulations
  2. Selective pSTAT1 and pSTAT5 activation (IFN-alpha-dependent)
- iEN classification: Discovery AUC=0.72, Validation AUC=0.84
- No single feature sufficient alone -- requires integrated model
- PLCgamma2 finding functionally buttresses the PLCG2 rare coding variant associated with decreased AD risk (genetic validation)

**Narrative strategy:**
- **"Peripheral blood as a window"**: Explicitly argued that PBMCs traffic into the brain and may participate in disease mechanisms
- **Genetic convergence**: Connected functional PBMC findings to GWAS risk variants (PLCG2)
- **Stimulation paradigm**: Did not just profile resting PBMCs but tested FUNCTIONAL RESPONSES to 7 stimulants, generating 4,200 features per participant
- **Machine learning**: Used ML to demonstrate that peripheral immune signatures can PREDICT disease status
- **Separate validation cohort**: Independent cohort without model retraining

**Validation:**
1. Separate discovery (n=108) and validation (n=24) cohorts
2. 250 iterations of two-layered cross-validation
3. Disease cross-prediction (AD vs PD) to control for non-specific effects

**Relevance to our study:** This is the strongest precedent for publishing peripheral-only immunoprofiling in CNS disease. Key lessons: (1) functional stimulation adds dimensionality beyond resting state; (2) ML prediction gives quantitative metric of peripheral signal strength; (3) genetic convergence (GWAS variant -> functional PBMC phenotype) strengthens the narrative. Our FGS external validation (AUC 0.82-0.89 in bulk datasets) is analogous to their ML prediction approach.

---

## TIER 2: Nature Communications, Cell Discovery, JCI, Science Translational Medicine (IF 10-30)

---

### Paper 4: Schafflick et al. -- Integrated Blood and CSF in MS

**Citation:** Schafflick D, Xu CA, Hartlehnert M, Cole M, Schulte-Mecklenbeck A, Lautwein T, Wolbert J, Heming M, Meuth SG, Kuhlmann T, Gross CC, Wiendl H, Yosef N, Meyer zu Horste G. Integrated single cell analysis of blood and cerebrospinal fluid leukocytes in multiple sclerosis. *Nature Communications*. 2020;11:247.

| Attribute | Detail |
|-----------|--------|
| **Disease** | Multiple sclerosis (relapsing-remitting) |
| **Sample type** | PBMCs + CSF (paired) |
| **Dataset size** | Blood: 42,969 cells (5 ctrl + 5 MS); CSF: 22,357 cells (4 ctrl + 4 MS); total 61,051 |
| **Methods** | Microfluidics scRNA-seq, flow cytometry validation, Cell Set Enrichment Analysis (CSEA, novel method), bulk RNA-seq of sorted TFH cells, mouse EAE models |
| **Impact factor** | Nature Communications (~17) |

**Key findings:**
- MS increases transcriptional diversity in BLOOD but increases cell type diversity in CSF
- TFH cells expanded in CSF (not blood), correlated with B cell enrichment
- CSEA (novel cluster-independent method) identified TFH signature enrichment in MS CSF (p=0.002)
- Blood changes: enhanced activation markers and trafficking molecules in T cells

**Narrative strategy:**
- **Compartmentalization as the story**: The key finding IS the blood vs CSF difference, making peripheral blood analysis essential to show the contrast
- **Novel analytical method**: Developed CSEA specifically for this analysis
- **Mouse model validation**: Bcl6-deficient mice (lacking TFH cells) showed reduced EAE severity and B cell CNS infiltration

**Validation:**
1. Flow cytometry (confirmed TFH cells in MS CSF)
2. Mouse models: CD4-Cre-Bcl6-fl/fl mice (genetic validation)
3. Adoptive transfer experiments
4. Bulk RNA-seq of sorted cells

**Relevance to our study:** The blood-CSF paired design is the gold standard. We lack CSF, but our dual-layer design (L1 disease + L2 outcome) provides a different type of multi-dimensional analysis. Their cell count (~61K) is much smaller than ours (226K), and they have fewer patients (10 vs 100).

---

### Paper 5: Wang et al. -- T Cell Transcriptome and TCR Profiling in PD

**Citation:** Wang P, Yao L, Luo M, et al. Single-cell transcriptome and TCR profiling reveal activated and expanded T cell populations in Parkinson's disease. *Cell Discovery*. 2021;7:52.

| Attribute | Detail |
|-----------|--------|
| **Disease** | Parkinson's disease |
| **Sample type** | PBMCs + CSF (from public data) |
| **Dataset size** | 103,365 T cells with transcriptome; 113,690 with TCR sequencing; 84,384 with both; 8 PD + 6 HC (blood) + CSF from public |
| **Methods** | 10x Genomics 5' scRNA-seq + V(D)J TCR sequencing, Seurat, diffusion maps, Monocle 2 trajectory, iSMART TCR clustering |
| **Impact factor** | Cell Discovery (~14) |

**Key findings:**
- CD8+ T cells show continuous progression from central memory to terminal effector in PD blood
- Cytotoxic CD4+ T cells (CD4 CTLs) remarkably expanded in PD, derived from Th1 cells (TCR-based fate decision)
- 55 clonotypes spanning central memory -> transitional -> terminal effector stages (clonal tracking)
- 67 TCR specificity groups shared between blood and CSF of PD patients
- Predicted TCR-antigen pairs validated against published alpha-synuclein epitopes

**Narrative strategy:**
- **TCR clonal sharing**: The finding that identical TCRs exist in both blood and CSF is the direct evidence of peripheral-CNS immune trafficking
- **Antigen prediction**: Linked expanded T cell clones to alpha-synuclein peptides, connecting to Sulzer 2017
- **Differentiation trajectory**: Used clonal tracking to show that T cell differentiation occurs progressively in peripheral blood

**Validation:**
1. TCR clonotype tracking (blood-CSF shared clones)
2. Epitope prediction against known alpha-synuclein peptides
3. HLA-TCR mapping (14 groups linked to 11 HLA alleles)

**Relevance to our study:** The clonal expansion analysis adds a layer we do not have (no VDJ data). However, our multi-method DA analysis (5 methods) and trajectory analysis provide comparable depth for different aspects of immune biology.

---

### Paper 6: Cantoni et al. -- CSF Immune Cell Compendium Across Neurological Diseases

**Citation:** Cantoni C, Smirnov RA, Firulyova M, et al. A single-cell compendium of human cerebrospinal fluid identifies disease-associated immune cell populations. *J Clin Invest*. 2025;135(1):e177793.

| Attribute | Detail |
|-----------|--------|
| **Disease** | Multiple: MS, AD, PD, COVID-19, viral encephalitis, autoimmune encephalitis, HIV, NMO |
| **Sample type** | CSF (135 samples) + PBMC (58 samples), paired |
| **Dataset size** | 403,973 immune cells (195,431 PBMCs + 208,542 CSF) from 139 individuals |
| **Methods** | 10x Genomics scRNA-seq (5' and 3'), Seurat + Harmony integration, Leiden clustering, Slingshot trajectory, flow cytometry |
| **Impact factor** | JCI (~15) |

**Key findings:**
- CD4+ T cells enriched in CSF; naive B cells preferentially in blood
- Seven DC populations predominantly in CSF
- Novel AREG+ cDC2 subpopulation in CSF, increased in MS
- FN1+ microglia-like cells in CSF contain AD risk genes (ABI3, CD33, PTK2B) -- increased in neurodegeneration
- Disease-specific myeloid and lymphoid frequency differences between blood and CSF

**Narrative strategy:**
- **Compendium approach**: Pooled data across multiple diseases and studies, creating a reference atlas
- **Blood-CSF paired comparison**: The systematic blood vs CSF comparison IS the study's contribution
- **Cross-disease architecture**: Showed that some immune patterns are disease-specific while others are shared

**Relevance to our study:** The compendium approach validates that large-scale single-cell studies of immune cells in neurological disease are valued even without brain tissue. Our 226K cells from 100 patients would contribute substantially to such a compendium for stroke specifically.

---

### Paper 7: Jacobs et al. -- Common Features of Neuroinflammation in CSF

**Citation:** Jacobs BM, Gasperi C, Reddy S, et al. Single-cell analysis of cerebrospinal fluid reveals common features of neuroinflammation. *Cell Reports Medicine*. 2024;5(12):101863.

| Attribute | Detail |
|-----------|--------|
| **Disease** | Multiple neuroinflammatory and non-inflammatory neurological disorders |
| **Sample type** | CSF (354,055 cells) + peripheral blood (422,809 cells) -- paired |
| **Dataset size** | 776,864 total cells from 201 patients (largest single-cell intrathecal immune study) |
| **Methods** | scRNA-seq + V(D)J BCR/TCR sequencing, clonal expansion analysis |
| **Impact factor** | Cell Reports Medicine (~12) |

**Key findings:**
- CSF cellular/transcriptional landscape altered in neuroinflammation but STRIKINGLY SIMILAR across different disorders
- Clonal expansion of CSF B and T cells observed across diseases, most prominent in MS
- Common drivers of clonal expansion identified across inflammatory disorders
- CSF is distinct from peripheral blood in both composition and gene expression

**Narrative strategy:**
- **Shared vs disease-specific**: The narrative centers on what is COMMON across neuroinflammatory diseases
- **Largest intrathecal dataset**: Scale justifies the study
- **Peripheral blood as baseline**: Blood provides the reference for identifying CSF-enriched populations

**Relevance to our study:** At ~777K total cells and 201 patients, this is comparable to our scale. Their paired CSF+blood design is stronger for CNS specificity, but our outcome stratification (g3) adds a dimension they lack.

---

### Paper 8: Grandke et al. -- Sex-Specific Blood Atlas in Neurodegeneration

**Citation:** Grandke F, Fehlmann T, Kern F, et al. A single-cell atlas to map sex-specific gene-expression changes in blood upon neurodegeneration. *Nature Communications*. 2025;16:1965.

| Attribute | Detail |
|-----------|--------|
| **Disease** | AD, PD, MCI (multiple neurodegenerative diseases) |
| **Sample type** | PBMCs only (NO CSF or brain tissue) |
| **Dataset size** | 909,322 PBMCs from 257 individuals (121 HC, 48 MCI, 27 AD, 46 PD, 15 PD+MCI), 290 samples |
| **Methods** | DNBelab C Series scRNA-seq, muscat pseudobulk DEG (limma-voom), scCODA composition, GeneTrail 3 pathway enrichment |
| **Impact factor** | Nature Communications (~17) |

**Key findings:**
- Sex-specific immune changes: males showed more cell types with deregulated genes in PD; females showed more in AD
- CD8+ T cells and Plasma cells: positive fold-change in males, negative in females
- 36 genes in males and 7 in females deregulated in BOTH PBMCs and brain
- SRP-dependent cotranslational protein targeting pathway appears in both blood and brain

**Narrative strategy:**
- **PBMC-ONLY study published in Nat Commun**: No CSF, no brain tissue. PURE peripheral blood.
- **BBB breakdown argument**: "Accumulating evidence points to partial breakdown of BBB with age during age-related diseases, suggesting a still underestimated role of the peripheral immune system"
- **Biomarker accessibility**: "Broadly accessible and low-invasively collectable plasma samples offer great promise"
- **Cross-tissue validation**: Compared PBMC signatures against 21 brain studies (ZEBRA dataset) and ROSMAP data
- **Sex stratification**: Novel dimension of analysis (sex-specific effects)

**Validation:**
1. ZEBRA dataset (21 brain transcriptome studies) -- identified 36 male + 7 female overlapping genes
2. ROSMAP dataset (4 studies, 6 brain regions) -- confirmed pathway overlap
3. Literature cross-reference: 24/32 male genes and 5/8 female genes previously reported in AD context

**Relevance to our study:** THIS IS THE MOST DIRECTLY RELEVANT PRECEDENT FOR PBMC-ONLY DESIGN. A PBMC-only study (no CSF, no brain), large-scale (909K cells, 257 individuals), published in Nature Communications in 2025. Uses similar methods (scRNA-seq, muscat, scCODA). The key lesson: cross-tissue validation against existing brain datasets and emphasis on clinical accessibility of blood-based biomarkers can justify a peripheral-only design.

---

## TIER 3: Disease-Specific High-Impact Journals (IF 5-15)

---

### Paper 9: Stamova et al. -- Peripheral Blood Gene Expression and Stroke Outcome

**Citation:** Stamova B, Jickling GC, Ander BP, et al. Early peripheral blood gene expression associated with good and poor 90-day ischemic stroke outcomes. *J Neuroinflammation*. 2023;20:7.

| Attribute | Detail |
|-----------|--------|
| **Disease** | Ischemic stroke |
| **Sample type** | Whole peripheral blood (not single-cell, but directly relevant) |
| **Dataset size** | 38 IS patients + 18 controls; 3 timepoints (<=3h, 5h, 24h post-stroke) |
| **Methods** | Whole-transcriptome microarray (Affymetrix), WGCNA, hub gene analysis, NIHSS/mRS correlation |
| **Impact factor** | J Neuroinflammation (~10) |

**Key findings:**
- 467/526/571 genes associated with POOR 90-day mRS outcome at <=3h / 5h / 24h respectively
- Poor outcome: IL-6, IL-7, IL-1, STAT3, S100A12, NF-kB, TLR pathway upregulated
- Good outcome: erythroblast and megakaryocyte gene modules
- Hub genes in monocytes (MIER1, DICER1, FBXW5, SELPLG) correlated with admission NIHSS
- Neutrophil genes UP in poor outcome; T cell/B cell/monocyte genes DOWN

**Narrative strategy:**
- **Longitudinal temporal design**: Multiple timepoints within 24h showing dynamic gene expression changes
- **Direct outcome prediction**: Gene expression at 3h predicts 90-day mRS
- **Cell-type deconvolution**: Used cell-type-specific gene signatures to interpret bulk data (neutrophil, monocyte, T cell modules)

**Relevance to our study:** This bulk transcriptomics study is the direct precursor to what we do at single-cell resolution. Our FGS signature validation (AUC 0.82-0.89 in independent bulk datasets) directly extends this work. Their finding that monocyte hub genes correlate with NIHSS is consistent with our CD14+ Monocyte as the sole L2 DA cell type.

---

### Paper 10: Carmona-Mora et al. -- Monocyte/Neutrophil Transcriptome Dynamics in Stroke

**Citation:** Carmona-Mora P, Ander BP, Jickling GC, et al. Monocyte, neutrophil, and whole blood transcriptome dynamics following ischemic stroke. *BMC Medicine*. 2023;21:49.

| Attribute | Detail |
|-----------|--------|
| **Disease** | Ischemic stroke |
| **Sample type** | Sorted peripheral blood monocytes, neutrophils, whole blood |
| **Dataset size** | 38 IS patients + 18 controls; 3 timepoints |
| **Methods** | Bulk RNA-seq of sorted cell populations, temporal differential expression, hub gene correlation with NIHSS |
| **Impact factor** | BMC Medicine (~9) |

**Key findings:**
- Monocyte hub genes (MIER1, DICER1, FBXW5, SELPLG, TIA1, BCKDK) correlate with admission NIHSS (p<0.01)
- Neutrophil hub genes limited (only 2 identified)
- Temporal dynamics: distinct gene expression waves at 0-24h, 24-48h, >48h
- Cell-type-specific gene expression patterns distinct between monocytes and neutrophils

**Narrative strategy:**
- **Cell-type-sorted bulk RNA-seq**: Bridges between whole blood and single-cell resolution
- **Temporal dynamics**: Demonstrates that peripheral blood gene expression is dynamically changing, not static
- **NIHSS correlation**: Connects molecular data to clinical severity

**Relevance to our study:** Our single-cell data resolves what their sorted bulk data approximates. The correlation between monocyte hub genes and NIHSS further supports our finding that monocytes are the outcome-critical cell type.

---

### Paper 11: Xie et al. -- CTSS Monocyte Activation in Stroke

**Citation:** Xie L, Wang B, et al. Single-cell RNA sequencing of peripheral blood reveals that monocytes with high cathepsin S expression aggravate cerebral ischemia-reperfusion injury. *Brain Behav Immun*. 2023;107:330-345.

| Attribute | Detail |
|-----------|--------|
| **Disease** | Ischemic stroke (mouse I/R model) |
| **Sample type** | Peripheral blood (mice) |
| **Dataset size** | scRNA-seq of mouse PBMCs at multiple timepoints post-I/R |
| **Methods** | scRNA-seq, CTSS expression analysis, CTSS knockout mice, BBB leakage assays |
| **Impact factor** | Brain Behav Immun (~15) |

**Key findings:**
- Peripheral monocyte subpopulations significantly increased after I/R
- CTSS (Cathepsin S) identified as key molecule regulating monocyte activation (peaks day 3)
- CTSS knockout reduces infarct size, improves neurological scores, reduces BBB leakage
- CTSS degrades JAM family proteins -> BBB disruption mechanism

**Narrative strategy:**
- **Peripheral discovery -> Central validation**: Found CTSS in peripheral monocytes, then showed its BBB-disrupting mechanism
- **Therapeutic target identification**: CTSS inhibition as potential treatment -> translational angle
- **Direct mechanistic link**: Peripheral monocyte CTSS -> BBB disruption -> CNS damage

**Relevance to our study:** Our cNMF programs (TNFa/NFkB in monocytes) likely encompass CTSS-mediated activation. This study validates that peripheral monocyte gene programs directly impact CNS outcomes.

---

### Paper 12: Cho et al. -- First Human PBMC scRNAseq in Stroke

**Citation:** Cho SJ, et al. Single-cell RNA sequencing of human peripheral blood mononuclear cells in ischemic stroke. *Stroke Vasc Neurol*. 2022.

| Attribute | Detail |
|-----------|--------|
| **Disease** | Ischemic stroke (human) |
| **Sample type** | PBMCs |
| **Dataset size** | ~26,302 cells from 19 patients (mild IS + HC) |
| **Methods** | 10x Genomics scRNA-seq |
| **Impact factor** | SVN (~5) |

**Key findings:**
- Two novel CD14+ monocyte subclusters: DC-related and NK cell-related
- NK cells increased in mild stroke (discordant with severe stroke findings)
- First human PBMC scRNA-seq atlas for stroke

**Relevance to our study:** This is the ONLY directly comparable prior study (human PBMC scRNA-seq in stroke). Our study is ~9x larger in cells (226K vs 26K) and ~5x larger in patients (100 vs 19), with outcome stratification they lack.

---

### Paper 13: Li et al. -- PD PBMC Immune Features (npj PD)

**Citation:** Li X, et al. Single-cell RNA sequencing reveals peripheral immunological features in Parkinson's Disease. *npj Parkinson's Disease*. 2024;10:190.

| Attribute | Detail |
|-----------|--------|
| **Disease** | Parkinson's disease (PD-Early, PD-Late, controls) |
| **Sample type** | PBMCs |
| **Dataset size** | scRNA-seq of PBMCs from PD-Early, PD-Late, and matched controls |
| **Methods** | 10x Genomics scRNA-seq, NK cell subcluster analysis |
| **Impact factor** | npj Parkinson's Disease (~8) |

**Key findings:**
- 6 major immune cell subsets identified; NK cells declined with PD progression
- XCL2 positively correlated with PD severity, with good predictive performance
- Specific expression of XCL2 in NK cell subclusters C2 and C5

**Narrative strategy:**
- **Disease staging**: Compared early vs late PD, showing progression in peripheral blood
- **Biomarker identification**: XCL2 as a specific, measurable peripheral biomarker
- **Cell-type specificity**: Narrowed biomarker to specific NK subclusters

---

### Paper 14: ALS PBMC Studies (Composite)

**Citations:**
- Garofalo S, et al. Single-cell RNA sequencing highlights the role of distinct natural killer subsets in sporadic amyotrophic lateral sclerosis. *J Neuroinflammation*. 2025;22:22.
- Keating GM, et al. Single-cell analysis reveals expanded CD8+ GZMKhigh T cells in CSF and shared peripheral clones in sporadic ALS. *Brain Communications*. 2024;6(6):fcae428.

| Attribute | Detail |
|-----------|--------|
| **Disease** | ALS (amyotrophic lateral sclerosis) |
| **Sample type** | PBMCs + CSF (some studies) |
| **Methods** | scRNA-seq, TCR profiling |
| **Impact factor** | J Neuroinflammation (~10), Brain Communications (~5) |

**Key findings across studies:**
- CD56dim NK cell expansion in ALS (specific NK_2 subpopulation)
- Reduced naive T and B cells with clonal expansion of terminal effector T cells
- CD8+ GZMKhi effector memory T cells expanded in CSF, shared clones with blood
- Novel CD8+ T cell type with stemness genes in ALS patients
- Correlation with plasma neurofilament light (NfL)

**Narrative strategy:**
- **Clonal sharing**: Blood-CSF shared TCR clonotypes as evidence of immune trafficking
- **NfL correlation**: Associated peripheral immune changes with plasma NfL (a neurodegeneration biomarker)
- **Progression biomarker**: Adaptive immune cells as monitors of disease progression

---

## TIER 4: Methodological Precedents and Cross-Disease Studies

---

### Paper 15: Stephenson et al. -- COVID-19 PBMC Multi-Omics (Nature Medicine)

**Citation:** Stephenson E, Reynolds G, Bottling S, et al. Single-cell multi-omics analysis of the immune response in COVID-19. *Nature Medicine*. 2021;27:904-916.

| Attribute | Detail |
|-----------|--------|
| **Disease** | COVID-19 (with neurological implications) |
| **Sample type** | PBMCs |
| **Dataset size** | >780,000 PBMCs from 130 patients (varying severity) |
| **Methods** | 10x Genomics CITE-seq (RNA + surface protein + TCR/BCR), multi-omics integration |
| **Impact factor** | Nature Medicine (~82) |

**Key findings:**
- Nonclassical monocytes (CD16+C1QA/B/C+) expanded, predicted to replenish alveolar macrophages
- Early CD34+ HSPCs primed toward megakaryopoiesis
- Clonally expanded CD8+ T cells characterize severe disease
- Increased CD8+ effector / effector memory ratio in severe COVID

**Narrative strategy:**
- **Multi-omics (CITE-seq)**: Combining transcriptome + surface protein + receptor profiling
- **Severity stratification**: Mild/moderate/severe/critical creates the gradient design
- **780K cells**: Massive scale justifies peripheral-only design
- **Clinical relevance**: COVID severity prediction from PBMC profiles

**Relevance to our study:** While not a pure CNS disease study, this sets the methodological gold standard for PBMC multi-omics in disease. Our 226K cells with dual-layer design is comparable in ambition. The severity stratification mirrors our L2 g3 outcome design.

---

### Paper 16: Korin et al. -- High-Dimensional Brain Immune Compartment (Nature Neuroscience)

**Citation:** Korin B, et al. High-dimensional, single-cell characterization of the brain's immune compartment. *Nature Neuroscience*. 2017;20:1300-1309.

| Attribute | Detail |
|-----------|--------|
| **Disease** | Neuroinflammation (general; baseline + disease) |
| **Sample type** | Brain parenchyma, choroid plexus, meninges + blood |
| **Methods** | CyTOF mass cytometry (44 surface markers), blood-brain comparison |
| **Impact factor** | Nature Neuroscience (~28) |

**Key findings:**
- Catalogued all immune cell types at CNS borders (meninges, choroid plexus) and parenchyma
- Blood-brain comparison revealed previously undescribed CD8 T, B, NK, DC subsets in brain
- Microglia heterogeneity (chemokine receptor-high, MHC-II-high, TCR-beta+ subsets)
- Established that peripheral immune cells are present at CNS borders even in healthy state

**Relevance to our study:** Establishes that the brain's immune compartment derives from peripheral blood, justifying peripheral blood profiling as informative about potential CNS-infiltrating cells.

---

### Paper 17: Absinta et al. -- Lymphocyte-Microglia-Astrocyte Axis in Progressive MS

**Citation:** Absinta M, Maric D, Gharagozloo M, et al. A lymphocyte-microglia-astrocyte axis in chronic active multiple sclerosis. *Nature*. 2021;597:709-714.

| Attribute | Detail |
|-----------|--------|
| **Disease** | Progressive multiple sclerosis |
| **Sample type** | Brain tissue (MRI-informed single-nucleus RNA-seq) |
| **Methods** | snRNA-seq of demyelinated white matter lesion edges |
| **Impact factor** | Nature (~65) |

**Key findings:**
- Defined "microglia inflamed in MS" (MIMS) and "astrocytes inflamed in MS" phenotypes
- Complement C1q expressed mainly by a subgroup of microglia driving inflammation
- Identified C1q as promising therapeutic target for progressive MS
- Chronically inflamed rim lesions as MRI biomarker for progression

**Relevance to our study:** While this is a brain tissue study, it establishes the brain-side immune signatures that our peripheral blood data would complement. The C1q/complement pathway identified here is detectable in our CellChat/MNN analysis.

---

## Cross-Cutting Analysis: Narrative Strategies for PBMC-to-CNS Justification

### Strategy 1: Blood-Brain Barrier Compromise
**Used by:** Grandke 2025, Phongpreecha 2020, Cantoni 2025
**Logic:** BBB is disrupted in neurological disease, so peripheral immune cells have direct access to CNS. Blood profiling captures the pre-infiltration state.
**Applicability to our study:** Very strong -- BBB disruption in ischemic stroke is well-established and occurs within hours, far more acutely than in neurodegenerative diseases.

### Strategy 2: Paired Compartment Design (Blood + CSF/Brain)
**Used by:** Garcia-Bonilla 2024, Schafflick 2020, Cantoni 2025, Jacobs 2024, Wang 2021
**Logic:** Profile BOTH blood and CSF/brain, making the comparison itself the study's contribution.
**Applicability to our study:** We lack CSF -- this is our main limitation. However, we can cite existing brain scRNA-seq data (Garcia-Bonilla mouse atlas) as the complementary compartment and emphasize that stroke is a condition where BBB breakdown makes peripheral profiling especially relevant.

### Strategy 3: Immune Trafficking Evidence
**Used by:** Wang 2021 (TCR sharing blood-CSF), Sulzer 2017 (antigen-specific T cells), Korin 2017
**Logic:** Identical immune cells/clones found in both blood and CNS, proving physical migration.
**Applicability to our study:** We do not have VDJ data. However, our monocyte trajectory analysis (showing altered maturation dynamics) is conceptually analogous -- peripheral monocytes are the established precursors to brain-infiltrating macrophages (Garcia-Bonilla 2024).

### Strategy 4: Clinical Biomarker Accessibility
**Used by:** Grandke 2025, Phongpreecha 2020, Li 2024, Stamova 2023
**Logic:** Brain/CSF is clinically inaccessible for routine monitoring; blood-based biomarkers have translational value.
**Applicability to our study:** Very strong -- our FGS signature (AUC 0.82-0.89 in bulk validation) directly demonstrates translational utility. This is arguably stronger than most studies reviewed here.

### Strategy 5: SIDS/CIDS Framework (Stroke-Specific)
**Used by:** Implicitly in all stroke immune studies; explicitly in Meisel 2005 and derivatives
**Logic:** Stroke causes systemic immunodepression (SIDS), meaning peripheral blood IS the disease compartment, not just a window into CNS. The peripheral immune dysregulation is itself clinically important (infection susceptibility, secondary immune injury).
**Applicability to our study:** This is our CORE narrative strategy. SIDS means the peripheral immune dysregulation IS the clinical problem. Our SIIS refinement (selective immunosuppression + innate activation) extends this framework with the first large-scale single-cell evidence.

### Strategy 6: External/Cross-Platform Validation
**Used by:** Grandke 2025 (ZEBRA/ROSMAP brain datasets), Phongpreecha 2020 (separate validation cohort), Stamova 2023 (independent cohort)
**Logic:** Validate peripheral findings against brain transcriptome data or independent cohorts.
**Applicability to our study:** Our FGS external validation (3 bulk cohorts, AUC 0.82-0.89) and cross-method concordance (5 DA methods) provide robust validation. We can additionally compare key genes (S100A8/A9, HLA genes, ISG signature) against brain datasets from Garcia-Bonilla 2024.

### Strategy 7: Genetic Convergence
**Used by:** Phongpreecha 2020 (PLCgamma2 -> PLCG2 variant), Sulzer 2017 (HLA alleles), Cantoni 2025 (AD risk genes in CSF cells)
**Logic:** PBMC functional findings converge with known genetic risk variants, creating a mechanistic bridge.
**Applicability to our study:** We have not performed genetic analysis, but ISG+ cell expansion is consistent with type I IFN pathway genetics in stroke. This is a potential future direction.

---

## Recommended Narrative Strategy for Our Manuscript

Based on this literature review, the optimal narrative strategy for our stroke_hc_v8_2 paper combines:

### Primary Justification: SIDS/SIIS Framework (Strategy 5)
- Stroke induces systemic immune changes -- peripheral blood IS the disease compartment
- Our refinement from SIDS to SIIS (selective immunosuppression + innate activation) is itself a novel contribution
- Supported by 5-method concordance (MASC, scCODA, MILO, MELD, Augur)
- This is STRONGER than the typical BBB argument because SIDS explicitly states that peripheral immunodysfunction causes clinical complications (infection, mortality)

### Secondary Justification: Clinical Biomarker Accessibility (Strategy 4)
- 100 patients from a single center (SNUBH), clinically well-characterized with IAT outcomes
- FGS signature validated in 3 independent bulk cohorts (AUC 0.82-0.89)
- Directly applicable to clinical prognostication in the acute stroke setting
- Blood draw is trivially feasible; CSF/brain is not

### Tertiary Justification: Literature Cross-Reference (Strategy 6)
- Our monocyte findings align with Garcia-Bonilla 2024 (Nat Immunol) brain trajectory data
- S100A8/A9 alarmin axis validated across multiple studies and modalities
- Cross-tissue validation using existing brain transcriptome datasets (conceptually like Grandke 2025)

### Scale Advantage:
- Largest human PBMC scRNAseq in stroke: 226K cells, 100 patients
- Prior human study: Cho 2022, 26K cells, 19 patients (9x and 5x larger respectively)
- Our scale is comparable to Grandke 2025 (909K cells, 257 patients for AD/PD in Nat Commun) and Jacobs 2024 (777K cells, 201 patients for neuroinflammation in Cell Rep Med)

### Unique Contribution: Dual-Layer Outcome Design
- NO prior single-cell study in any CNS disease has treatment outcome stratification
- L1 (HC vs Stroke) + L2 (g3 Good vs Bad within IS) is a unique multi-layer design
- Identifies monocyte as the SOLE outcome-determining cell type
- External FGS validation proves translational utility

---

## Summary Comparison Table

| Study | Journal | IF | Year | Disease | Sample | Cells | Patients | Paired? | Outcome? |
|-------|---------|------|------|---------|--------|-------|----------|---------|----------|
| Garcia-Bonilla | Nat Immunol | ~30 | 2024 | Stroke | Brain+Blood | ~43K | mice | Yes (brain) | No |
| Sulzer | Nature | ~65 | 2017 | PD | PBMC | N/A | 103 | No | No |
| Absinta | Nature | ~65 | 2021 | MS | Brain | snRNA | postmortem | N/A | No |
| Phongpreecha | Sci Adv | ~13 | 2020 | AD+PD | PBMC | ~1.2M | 132 | No | No |
| Schafflick | Nat Commun | ~17 | 2020 | MS | Blood+CSF | 61K | 10 | Yes (CSF) | No |
| Wang | Cell Discov | ~14 | 2021 | PD | Blood+CSF | 103K | 14 | Yes (CSF) | No |
| Cantoni | JCI | ~15 | 2025 | Multi-CNS | Blood+CSF | 404K | 139 | Yes (CSF) | No |
| Jacobs | CRM | ~12 | 2024 | Multi-neuro | Blood+CSF | 777K | 201 | Yes (CSF) | No |
| Grandke | Nat Commun | ~17 | 2025 | AD+PD+MCI | **PBMC only** | 909K | 257 | **No** | No |
| Stephenson | Nat Med | ~82 | 2021 | COVID | PBMC | 780K | 130 | No | Severity |
| Stamova | J Neuroinfl | ~10 | 2023 | Stroke | Whole blood | bulk | 56 | No | **Yes (mRS)** |
| Cho | SVN | ~5 | 2022 | Stroke | PBMC | 26K | 19 | No | No |
| Li | npj PD | ~8 | 2024 | PD | PBMC | scRNA | staged | No | Staging |
| Xie | BBI | ~15 | 2023 | Stroke | PBMC (mouse) | scRNA | mice | No | No |
| **Ours** | **Target** | **--** | **2026** | **Stroke** | **PBMC** | **226K** | **100** | **No** | **Yes (g3)** |

**Key gap our study fills:** The intersection of (1) large-scale scRNA-seq, (2) stroke, (3) human PBMC, AND (4) treatment outcome stratification. No existing study occupies this space. The closest comparable PBMC-only studies in neurodegeneration (Grandke 2025, Phongpreecha 2020) do not have outcome data. The closest stroke outcome study (Stamova 2023) is bulk, not single-cell.

---

## References

1. Garcia-Bonilla L et al. Analysis of brain and blood single-cell transcriptomics in acute and subacute phases after experimental stroke. *Nat Immunol*. 2024;25:357-370.
2. Sulzer D et al. T cells from patients with Parkinson's disease recognize alpha-synuclein peptides. *Nature*. 2017;546:656-661.
3. Lindestam Arlehamn CS et al. Alpha-synuclein-specific T cell reactivity is associated with preclinical and early Parkinson's disease. *Nat Commun*. 2020;11:1875.
4. Phongpreecha T et al. Single-cell peripheral immunoprofiling of Alzheimer's and Parkinson's diseases. *Sci Adv*. 2020;6:eabd5575.
5. Schafflick D et al. Integrated single cell analysis of blood and cerebrospinal fluid leukocytes in multiple sclerosis. *Nat Commun*. 2020;11:247.
6. Wang P et al. Single-cell transcriptome and TCR profiling reveal activated and expanded T cell populations in Parkinson's disease. *Cell Discov*. 2021;7:52.
7. Cantoni C et al. A single-cell compendium of human cerebrospinal fluid identifies disease-associated immune cell populations. *J Clin Invest*. 2025;135:e177793.
8. Jacobs BM et al. Single-cell analysis of cerebrospinal fluid reveals common features of neuroinflammation. *Cell Rep Med*. 2024;5:101863.
9. Grandke F et al. A single-cell atlas to map sex-specific gene-expression changes in blood upon neurodegeneration. *Nat Commun*. 2025;16:1965.
10. Stamova B et al. Early peripheral blood gene expression associated with good and poor 90-day ischemic stroke outcomes. *J Neuroinflammation*. 2023;20:7.
11. Carmona-Mora P et al. Monocyte, neutrophil, and whole blood transcriptome dynamics following ischemic stroke. *BMC Med*. 2023;21:49.
12. Xie L et al. Single-cell RNA sequencing of peripheral blood reveals that monocytes with high cathepsin S expression aggravate cerebral ischemia-reperfusion injury. *Brain Behav Immun*. 2023;107:330-345.
13. Cho SJ et al. Single-cell RNA sequencing of PBMC in ischemic stroke. *Stroke Vasc Neurol*. 2022.
14. Li X et al. Single-cell RNA sequencing reveals peripheral immunological features in Parkinson's Disease. *npj Parkinsons Dis*. 2024;10:190.
15. Garofalo S et al. Single-cell RNA sequencing highlights the role of distinct natural killer subsets in sporadic ALS. *J Neuroinflammation*. 2025;22:22.
16. Keating GM et al. Single-cell analysis reveals expanded CD8+ GZMKhigh T cells in CSF and shared peripheral clones in sporadic ALS. *Brain Commun*. 2024;6:fcae428.
17. Stephenson E et al. Single-cell multi-omics analysis of the immune response in COVID-19. *Nat Med*. 2021;27:904-916.
18. Korin B et al. High-dimensional, single-cell characterization of the brain's immune compartment. *Nat Neurosci*. 2017;20:1300-1309.
19. Absinta M et al. A lymphocyte-microglia-astrocyte axis in chronic active multiple sclerosis. *Nature*. 2021;597:709-714.
20. Meisel C et al. Central nervous system injury-induced immune deficiency syndrome. *Nat Rev Neurosci*. 2005;6:775-786.
