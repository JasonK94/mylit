# Stroke scRNAseq & Immune Profiling Landmark Papers

## 1. Landmark scRNAseq Studies

### 1.1 Brain Tissue scRNAseq

#### Anrather Lab — Brain & Blood Longitudinal scRNAseq (Nature Immunology 2024)
- **Citation:** Garcia-Bonilla L, Shahanoor Z, Sciortino R, et al. "Brain and blood single-cell transcriptomics in acute and subacute phases after experimental stroke." *Nat Immunol* 25:357-370 (2024). DOI: 10.1038/s41590-023-01711-x
- **Design:** Longitudinal scRNAseq of brain + blood; young + aged mice; sham, Day 2, Day 14 post-MCAO; 43,269 brain cells + blood leukocytes
- **Key findings:**
  - Brain-infiltrating leukocytes show **higher transcriptomic divergence** than circulating counterparts
  - Local tissue microenvironment (not peripheral priming) determines phenotypes
  - MdC (monocyte-derived cells) increased 2% → 16% after stroke
  - 8 distinct microglia subclusters with DAM signature (Apoe, Cst7, Clec7a)
  - 6 neutrophil clusters: continuous de novo recruitment from blood
  - MdC trajectory: blood monocytes → brain macrophages with repair phenotype (Igf1) at D14
  - **ISG signature (Ifit3, Isg15, Rsad2)** across multiple cell types
  - Only 9% of pathways shared between human and mouse in acute phase; 47% in subacute
- **Biomarkers:** Spp1, Cst7, Apoe, Lgals3 (microglia); Lrg1 (endothelial); Igf1 (repair macrophages)
- **Database:** https://anratherlab.shinyapps.io/strokevis/
- **Relevance:** Most comprehensive brain+blood single-cell atlas in stroke

#### Zheng et al. — Ischemic Mouse Brain scRNAseq (J Cereb Blood Flow Metab 2022)
- **Citation:** Zheng K, et al. "Single-cell RNA-seq reveals the transcriptional landscape in ischemic stroke." *J Cereb Blood Flow Metab* 42(1):56-73 (2022). PMC8721774
- **Key findings:**
  - Key transcription factors: **STAT1, JUN, MYC**
  - Key chemokines: **CXCL9, CXCL13**
  - ISG (interferon-stimulated gene) clusters prominently upregulated
  - CellChat: total interactions +104%, strength +95% in stroke vs intact
  - **LGALS9-CD44** identified as critical signaling pathway

#### Integrating Spatial + scRNAseq (Science Translational Medicine 2023)
- **Citation:** [Multiple groups]. "Integrating spatial and single-cell transcriptomics to characterize the molecular and cellular architecture of the ischemic mouse brain." *Sci Transl Med* (2023).
- **Key findings:** Spatial context of immune cell infiltration, microglial activation zones

#### Yang et al. — Multi-compartment Myeloid Ontology (Advanced Science 2025)
- **Citation:** Yang S, et al. "Single-Cell Transcriptomes of Immune Cells from Multiple Compartments Redefine the Ontology of Myeloid Subtypes Post-Stroke." *Adv Sci* 12(13):2408722 (2025). DOI: 10.1002/advs.202408722
- **Design:** Mouse p-MCAO; scRNAseq of brain, calvarial bone marrow (CBM), femur bone marrow (FBM), PBMCs
- **Key findings:**
  - Neutrophils: 0.7% → 19.6% (dramatic increase)
  - Monocytes/macrophages: 8.1% → 19.2%
  - B cells: 43.2% → 10.9% (dramatic decrease)
  - Novel CD14+ neutrophil subpopulation resembling CBM neutrophils
  - Two distinct origins: disease inflammatory macrophages (from circulating monocytes, skull-derived) vs disease-associated microglia (from pre-existing homeostatic microglia)
- **Relevance:** First multi-compartment survey; skull bone marrow as source of infiltrating myeloid cells

### 1.2 Peripheral Blood (PBMC) scRNAseq — Relatively Rare

#### Cho et al. — Circulating Immune Cell Landscape in Mild IS (Stroke Vasc Neurol 2022)
- **Citation:** Cho YE, Lee H, Bae HR, et al. "Circulating immune cell landscape in patients who had mild ischaemic stroke." *Stroke Vasc Neurol* 7:319-327 (2022). DOI: 10.1136/svn-2021-001224
- **Design:** 10x Genomics scRNAseq of PBMCs; 10 acute IS (within 24h) + 9 matched HC; 12,156 stroke + 14,146 control cells
- **Key findings:**
  - CD14+ monocyte subtype cluster: **↓ in stroke** (3.8% vs 9.4% controls)
  - Two novel CD14+ monocyte subclusters: DC-related and NK cell-related
  - NK cell cluster: **↑ in stroke** (11.6% vs 6.2%, p=0.022)
  - NFkB signaling enhanced in NK cell population
  - CD4+ T cells: 14.3% stroke vs 11.8% control (trend ↑)
  - CD8+ T cells: 12.7% both groups (no change)
- **Biomarkers:** CX3CR1, GIMAP7 (NK cell activation); NFkB pathway
- **Significance:** First human PBMC scRNAseq in ischemic stroke

#### CyTOF Profiling — Location-Specific Immune Changes (Frontiers Immunology 2026)
- **Citation:** [Authors]. "CyTOF profiling identifies location-specific peripheral immune checkpoint and immune cell subset in mild ischemic stroke." *Front Immunol* (2026).
- **Key findings:**
  - **Subcortical infarction**: myeloid-dominant inflammatory program, persistently elevated classical monocytes, delayed checkpoint engagement
  - **Cortical infarction**: adaptive immune regulation, earlier modulation of CD4+ and CD8+ T-cell subsets

#### Decoding Immune Cell Dynamics — scRNAseq Analysis (PMC 2025)
- **Citation:** [Authors]. "Decoding immune cell dynamics in ischemic stroke: insights from single-cell RNA sequencing analysis." PMC12037566 (2025).
- **Key findings:** Comprehensive characterization of immune cell transcriptomic shifts

#### Xie et al. — CTSS-High Monocytes Aggravate I/R Injury (Brain Behav Immun 2023)
- **Citation:** Xie L, Zhang S, Huang L, et al. "Single-cell RNA sequencing of peripheral blood reveals that monocytes with high cathepsin S expression aggravate cerebral ischemia-reperfusion injury." *Brain Behav Immun* 107:330-344 (2023). DOI: 10.1016/j.bbi.2022.11.001
- **Design:** Mouse MCAO, scRNAseq of peripheral blood
- **Key findings:**
  - Peripheral monocyte subpopulations significantly ↑ after I/R
  - **CTSS (Cathepsin S)** identified as key molecule regulating monocyte activation (peaks day 3)
  - CTSS KO: reduced infarct, improved neurological scores, reduced BBB leakage
  - CTSS degrades JAM family proteins → BBB disruption
- **Biomarkers:** CTSS as therapeutic target; JAM proteins
- **Significance:** Specific monocyte activation pathway linking peripheral blood to BBB damage

#### ISC 2025 Abstract — PBMC scRNAseq in Early Neurological Deterioration
- **Citation:** "Circulating immune cell landscape in ischemic stroke patients with early neurological deterioration." *Stroke* 56(Suppl 1):TP307 (2025).
- **Design:** 5 END + 5 non-END patients; 20 PBMC samples (before/after END)
- **Key findings:**
  - Mononuclear macrophage proportion ↑ in END patients
  - DEG count significantly higher in END group
- **Significance:** First scRNAseq linking peripheral immune changes to early deterioration

### 1.3 Bulk Transcriptomics of Sorted Populations

#### Carmona-Mora et al. — Monocyte/Neutrophil Transcriptome Dynamics (BMC Med 2023)
- **Citation:** Carmona-Mora P, et al. "Monocyte, neutrophil, and whole blood transcriptome dynamics following ischemic stroke." *BMC Med* 21:65 (2023). DOI: 10.1186/s12916-023-02766-1
- **Design:** 38 IS patients + 18 controls; sorted cell populations; 3 timepoints (0-24h, 24-48h, >48h)
- **Key findings:**
  - Genes generally **upregulated in neutrophils** and **downregulated in monocytes**
  - Unique temporal patterns for each cell type and stroke etiology (CE, LVO, SAO)
  - Interleukin signaling enrichment varies by timepoint and etiology

#### Amini et al. — Early Gene Expression Predicting 90-Day Outcome (J Neuroinflammation 2023)
- **Citation:** Amini H, et al. "Early peripheral blood gene expression associated with good and poor 90-day ischemic stroke outcomes." *J Neuroinflammation* 20:7 (2023). DOI: 10.1186/s12974-022-02680-y
- **Design:** CLEAR trial patients; blood at ≤3h, 5h, 24h after stroke
- **Key findings:**
  - 467 genes at ≤3h, 526 at 5h, 571 at 24h → poor 90-day outcome
  - **Neutrophil-specific genes** enriched in poor outcome modules
  - **T cell-specific genes** enriched in good outcome modules
- **Significance:** Cell-type-specific blood gene expression within hours predicts long-term outcome

### 1.4 SIDS Mechanism Papers

#### Prass et al. — Sympathetic Mechanism (J Exp Med 2003)
- **Citation:** Prass K, et al. "Stroke-induced immunodeficiency promotes spontaneous bacterial infections and is mediated by sympathetic activation." *J Exp Med* 198:725-736 (2003).
- **Key findings:**
  - Catecholamine-mediated lymphocyte dysfunction is the key factor
  - SNS overactivation within hours drives immunodepression
  - Beta-receptor blockade prevents impaired T cell IFN-gamma → reduces lung bacteria
- **Significance:** Seminal paper establishing SNS mechanism of stroke immunosuppression

#### Meisel et al. — CIDS Framework (Nature Rev Neurosci 2005)
- **Citation:** Meisel C, Schwab JM, et al. "Central nervous system injury-induced immune deficiency syndrome." *Nat Rev Neurosci* 6:775-786 (2005). DOI: 10.1038/nrn1765
- **Key findings:**
  - Three axes: SNS, HPA, vagus nerve
  - Manifestations: lymphopenia, dysfunctional innate+adaptive immunity
  - Indicators: CD4+ T lymphopenia, ↑ NLR, ↓ monocyte HLA-DR
- **Significance:** Definitive CIDS framework; >1000 citations

#### Liesz et al. — Tregs as Cerebroprotective (Nature Medicine 2009)
- **Citation:** Liesz A, et al. "Regulatory T cells are key cerebroprotective immunomodulators in acute experimental stroke." *Nat Med* 15:192-199 (2009). DOI: 10.1038/nm.1927
- **Key findings:**
  - Treg depletion → enlarged infarct, worse neurological function
  - IL-10 signaling essential for Treg cerebroprotective effect
  - Treg-derived IL-10 is the main mediator
- **Significance:** Seminal Treg neuroprotection paper; reproduced >20 times

### 1.5 PBMC scRNAseq Gap in the Literature

> **Important:** PBMC scRNAseq in stroke patients remains **relatively rare**. Most scRNAseq studies focus on mouse brain tissue (MCAO model). Our study (100 patients, HC vs IS, PBMC scRNAseq with g3 stratification) fills a significant gap.

---

## 2. Flow Cytometry & Immune Cell Composition Studies

### 2.1 Vogelgesang et al. — Lymphocyte Subsets & Infection (Stroke 2008)
- **Citation:** Vogelgesang A, et al. "Analysis of lymphocyte subsets in patients with stroke and their influence on infection after stroke." *Stroke* 39(1):237-241 (2008).
- **Key findings:**
  - Dramatic loss of CD4+ and CD8+ T cells within 12 hours
  - Loss of CD4+ T-lymphocytes is **persistent** in patients developing infection
  - T helper cell functional defects contribute to immunosuppression

### 2.2 CD4+, CD8+, Double-Negative T Cells (Frontiers Cell Neurosci 2025)
- **Citation:** [Authors]. "Circulating CD4+, CD8+, and double-negative T cells in ischemic stroke and stroke-associated infection." *Front Cell Neurosci* (2025).
- **Design:** Prospective, 52 patients; flow cytometry at days 1, 3, 10, 90
- **Key findings:** Temporal profiling of T cell recovery trajectory

### 2.3 Monocyte Subtypes & Clinical Course (J Cereb Blood Flow Metab 2009)
- **Citation:** Urra X, et al. "Monocyte subtypes predict clinical course and prognosis in acute stroke." *J Cereb Blood Flow Metab* (2009).
- **Key findings:**
  | Subset | Markers | Change Post-Stroke | Clinical Significance |
  |--------|---------|-------------------|----------------------|
  | Classical | CD14++CD16- | Increased (acute/subacute) | 80-90% of circulating; infiltrate via CCL2-CCR2 |
  | Intermediate | CD14++CD16+ | Increased (acute) | Associated with CNS tissue damage |
  | Non-classical | CD14+CD16++ | Decreased | Associated with infection susceptibility |

### 2.4 CD163 as Monocyte Activation Marker (2021)
- **Citation:** [Authors]. PMC8268853 (2021).
- **Key findings:** CD163 proposed as biomarker of monocyte activation in IS

---

## 3. Multi-Omics Integration Studies

### 3.1 Machine Learning + scRNA-seq for IS Biomarkers (Scientific Reports 2024)
- **Citation:** [Authors]. "Integrating machine learning and single-cell transcriptomic analysis to identify potential biomarkers and analyze immune features of ischemic stroke." *Sci Rep* (2024).
- **Key findings:**
  - Five core mRNAs: **S100A6, ANXA3, NCF4, CAPG, ARPC1B**
  - Key microRNA: miR-298-5p
  - Combined ML + scRNAseq approach

### 3.2 Key Biomarkers via scRNAseq + WGCNA (Aging 2023)
- **Citation:** Tao J, Xie X, Luo M, Sun Q. "Identification of key biomarkers in ischemic stroke: single-cell sequencing and weighted co-expression network analysis." *Aging* (2023).
- **Key findings:**
  - **MRPS11, MRPS12** (mitochondrial ribosomal proteins) — both down-regulated in IS
  - MRPS12 expression decreased during pre-B cell CD34 differentiation

### 3.3 Multi-Omics Review (Ageing Research Reviews 2022)
- **Citation:** [Authors]. "Multi-omics research strategies in ischemic stroke: A multidimensional perspective." *Ageing Res Rev* (2022).
- **Key points:** Framework for integrating genomics, transcriptomics, proteomics, metabolomics

### 3.4 Single-Cell Omics Review (European J Pharmacology 2025)
- **Citation:** [Authors]. "Single-cell omics: moving towards a new era in ischemic stroke research." *Eur J Pharmacol* (2025).

### 3.5 Multilevel Omics for Stroke (Nature Reviews Neurology 2020)
- **Citation:** Montaner J, et al. "Multilevel omics for the discovery of biomarkers and therapeutic targets for stroke." *Nat Rev Neurol* 16:247-264 (2020).
- **Key points:** Authoritative framework — integromics needed to bridge discovery and clinical biomarkers

---

## 4. Key Biomarkers from scRNAseq Studies

### Genes/Proteins Consistently Identified

| Gene/Protein | Context | Direction in Stroke | Cell Type |
|-------------|---------|-------------------|-----------|
| STAT1 | TF hub | ↑ | Multiple |
| JUN | TF hub | ↑ | Multiple |
| MYC | TF hub | ↑ | Multiple |
| CXCL9 | Chemokine | ↑ | Myeloid |
| CXCL13 | Chemokine | ↑ | B cell zone |
| LGALS9 | Ligand (CellChat) | ↑ | Microglia/Macrophage |
| S100A6 | Calcium binding | ↑ | Myeloid |
| ANXA3 | Annexin | ↑ | Myeloid |
| NCF4 | NADPH oxidase | ↑ | Neutrophil/Monocyte |
| ISG15 | ISG marker | ↑ | ISG-hi cells |
| S100A8/A9 | Alarmin/DAMP | ↑ | Neutrophil/Monocyte |

### Cell Type Composition Changes (Consensus from Literature)

| Cell Type | Change in Stroke | Evidence Level | Mechanism |
|-----------|-----------------|----------------|-----------|
| **Classical Monocyte** | ↑↑ | Strong (flow + scRNAseq) | Bone marrow mobilization |
| **Intermediate Monocyte** | ↑ | Moderate (flow) | CNS damage marker |
| **Non-classical Monocyte** | ↓ | Moderate (flow) | Infection susceptibility |
| **NK Cells** | ↑ (blood), ↓ (spleen) | Moderate (scRNAseq) | Redistribution |
| **B Cells** | ↓↓ | Strong (flow + scRNAseq) | SIDS/apoptosis |
| **CD4+ T Cells** | ↓ | Strong (flow) | SIDS/lymphopenia |
| **CD8+ T Cells** | ↓ then recover | Moderate (flow) | Transient lymphopenia |
| **Tregs** | Variable | Conflicting | Neuroprotective role |
| **Neutrophils** | ↑↑ | Very strong (clinical) | First responders |
| **Gamma-delta T** | ↑ (brain) | Moderate (mouse) | IL-17 production |

---

## 5. Relevance to Our Study (stroke_hc_v8_2)

### Direct Validations Possible
1. **Monocyte expansion** — our MASC/MILO/scCODA all show CD14+ Mono DA → literature-consistent
2. **B cell decrease** — validate against SIDS literature (43.2% → 10.9%)
3. **NK cell changes** — check our data against the ↑ in mild stroke
4. **CD4/CD8 ratio** — compare with Vogelgesang temporal data
5. **ISG signature** — our FGS likely captures ISG-hi expansion

### Novel Contributions
1. **PBMC scRNAseq with 226K cells** — one of the largest human PBMC stroke datasets
2. **g3 stratification (IAT outcome)** — no prior scRNAseq study examines EVT response at single-cell level
3. **CellChat v2 with condition comparison** — can validate LGALS9-CD44, IL18 pathways
4. **Multi-method DA concordance** — MASC × MILO × scCODA 3-method approach is rigorous
5. **Cross-validation with UK Biobank proteomics** — bridge scRNAseq and circulating protein levels
