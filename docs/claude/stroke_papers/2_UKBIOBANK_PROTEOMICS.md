# UK Biobank Proteomics & Stroke Biomarkers

## 1. UKB-PPP Foundation

### Sun BB et al. — The UKB-PPP Atlas (Nature 2023)
- **Citation:** Sun BB, Chiou J, Traylor M, et al. "Plasma proteomic associations with genetics and health in the UK Biobank." *Nature* 622:329-338 (2023).
- **Design:** 54,219 UKB participants, Olink Explore platform, 2,923 proteins
- **Key findings:**
  - 14,287 primary pQTLs (81% previously undescribed)
  - ~83% of proteins have at least one pQTL; ~67% have cis-associations
  - >860 trans-pQTL loci with interacting partners
- **Relevance:** Master pQTL atlas for all subsequent UKB stroke proteomics studies

---

## 2. Proteomics-Based Stroke Risk Prediction

### 2.1 Zhong et al. — IS Protein Risk Score (Stroke 2025)
- **Citation:** Zhong X, et al. "Large-Scale Plasma Proteomics Profiles for Predicting Ischemic Stroke Risk in the General Population." *Stroke* 56(2) (2025).
- **Design:** 53,029 UKB-PPP participants, LASSO regression, 2,911 proteins
- **Key findings:**
  - **17 proteins** selected by LASSO (13 positive, 4 negative association)
  - **Top 5: GDF15, PLAUR (suPAR), NT-proBNP, IGFBP4, BCAN**
  - Protein risk score C-statistic = **0.765** vs age+sex (0.720), clinical RF (0.753), PRS (0.730)
  - PLAUR, MDK, GDF15 central in protein interaction network
- **Proteins:** GDF15, PLAUR, NT-proBNP, IGFBP4, BCAN, MDK + 11 others

### 2.2 Ho FK et al. — CVD & Dementia Prediction (Circulation 2025)
- **Citation:** Ho FK, Mark PB, et al. "A Proteomics-Based Approach for Prediction of Different Cardiovascular Diseases and Dementia." *Circulation* 151(6) (2025).
- **Design:** 51,859 UKB participants, median f/u 13.6 years, 4,857 MACE events
- **Key findings:**
  - 222 proteins predict across all outcomes; 86 for primary MACE
  - **Top associations:** NT-proBNP (HR 1.68/SD), proADM (HR 1.60), GDF-15 (HR 1.47), WFDC2 (HR 1.46), IGFBP4 (HR 1.41)
  - Protein model: NRI +0.09, C-stat +0.051 vs PREVENT score

### 2.3 Johansen et al. — Embolic vs. Thrombotic Stroke (Ann Neurol 2025)
- **Citation:** Johansen MC, et al. "Proteomics and the Risk of Incident Embolic and Thrombotic Stroke." *Annals of Neurology* 98(2) (2025).
- **Design:** 10,929 participants, ~20y f/u; 168 EIS, 459 TIS events
- **Key findings:**
  - 20 mid-life proteins associated with either EIS or TIS
  - **NPPB (NT-proBNP): EIS-specific, NOT associated with TIS** (p-diff = 9.14E-7)
  - EIS: cardiac dysfunction proteins (NPPB, SPINK4, OMG, NCAN)
  - TIS: inflammation proteins (interleukins)
- **Implication:** Stroke subtypes have fundamentally different proteomic signatures

### 2.4 Lind et al. — Plasma Protein Profiling (JAHA 2021)
- **Citation:** Lind L, et al. "Large-Scale Plasma Protein Profiling of Incident Myocardial Infarction, Ischemic Stroke, and Heart Failure." *JAHA* 10:e023330 (2021).
- **Design:** 742 proteins in 826 males (ULSAM cohort), age 70
- **Key findings:**
  - 13 proteins associated with incident IS
  - **GDF-15, WFDC2, KIM-1: associated with ALL THREE CVD outcomes** (MI, IS, HF)
  - LASSO panel of 11 proteins: +5.0% discrimination improvement

---

## 3. Mendelian Randomization (MR) Studies

### 3.1 Chong M et al. — Plasma Proteome for Stroke Targets (Nature Comms 2022)
- **Citation:** Chong M, Sjaarda J, et al. "Systematic Mendelian randomization using the human plasma proteome to discover potential therapeutic targets for stroke." *Nat Commun* 13:6143 (2022).
- **Design:** Two-sample MR, 308 plasma proteins, cis-pQTLs as instruments
- **Key findings:**
  - **6 causal proteins: TFPI, TMPRSS5, CD40, MMP12, IL6RA, CD6**
  - TFPI, CD40, IL6RA, MMP12 → lower stroke risk
  - TMPRSS5, CD6 → higher stroke risk
  - Colocalization confirmed for 4/6 proteins
  - IL6RA = tocilizumab target; CD40 = immunotherapy target

### 3.2 Chong M et al. — Novel Drug Targets by Subtype (Circulation 2019)
- **Citation:** Chong M, et al. "Novel Drug Targets for Ischemic Stroke Identified Through Mendelian Randomization Analysis of the Blood Proteome." *Circulation* 140(10):819-830 (2019).
- **Design:** 653 proteins, 3 IS subtypes (LAA, CE, SAO)
- **Key findings:**
  - **SCARA5**: novel target for cardioembolic stroke
  - Side-effect profiling for drug development prioritization

### 3.3 Cross-Ancestry MR — East Asian + European (Stroke 2025)
- **Citation:** [CKB/UKB consortium]. "Proteome-Wide Genetic Study in East Asians and Europeans Identified Multiple Therapeutic Targets for Ischemic Stroke." *Stroke* 56 (2025).
- **Design:** 2,923 Olink proteins; IS: 22,664 (East Asian) + 62,100 (European)
- **Key findings:**
  - **10 causal proteins: FGF5, TMPRSS5, FURIN, F11, ALDH2, ABO (positive); GRK5, KIAA0319, PROCR, MMP12 (inverse)**
  - Directionally consistent across ancestries
  - **F11 (Factor XI)**: active clinical trials for factor XI inhibitors (abelacimab, osocimab)
  - FURIN, FGF5: blood pressure pathway links

### 3.4 MVP/UKB — Stroke Progression Proteins (Stroke 2024)
- **Citation:** [MVP/UKB]. "Protein Identification for Stroke Progression via Mendelian Randomization in Million Veteran Program and UK Biobank." *Stroke* 55 (2024).
- **Design:** GWAS for subsequent MACE (n=51,929) after first stroke; MR of 1,463 proteins
- **Key findings:**
  - **CCL27, TNFRSF14: associated with SUBSEQUENT MACE but NOT initial stroke**
  - Different protein drivers for recurrence vs. first event
  - Inflammatory mechanisms for secondary prevention

### 3.5 Li et al. — CSF + Plasma MR (Int J Surgery 2024)
- **Citation:** Li C, et al. "Promising therapeutic targets for ischemic stroke identified from plasma and cerebrospinal fluid proteomes: a multicenter Mendelian randomization study." *Int J Surg* 110(2) (2024).
- **Key findings:**
  - Plasma: **MMP12, SWAP70** causally associated with IS
  - CSF: **prekallikrein** causally associated
  - Lanadelumab (anti-prekallikrein) as potential repurposed drug

### 3.6 Cross-Biobank Protein-CVD (Nature Cardiovasc Res 2024)
- **Citation:** [UKB/CKB consortium]. "Measured and genetically predicted protein levels and cardiovascular diseases in UK Biobank and China Kadoorie Biobank." *Nat Cardiovasc Res* 3 (2024).
- **Design:** 2,919 proteins in 52,164 UKB + validation in CKB
- **Key findings:**
  - 126 proteins associated with ALL THREE CVD outcomes (MI, IS, HF)
  - 118 replicated in CKB (cross-ancestry)
  - MR + coloc: 47 causal; 18 strong colocalization
  - **FGF5, PROCR, FURIN** = top drug targets

---

## 4. Olink Proteomics in Stroke

### 4.1 Li et al. — Immune Response Proteins (J Proteome Res 2024)
- **Citation:** Li Y, et al. "Serum Olink Proteomics-Based Identification of Protein Biomarkers Associated with the Immune Response in Ischemic Stroke." *J Proteome Res* 23(3) (2024).
- **Design:** 92 immune response proteins by Olink; 88 IS + 88 controls
- **Key findings:**
  - 59/92 differentially expressed
  - **6 biomarkers: MASP1, STC1, HCLS1, CLEC4D, PTH1R, PIK3AP1**
  - Diagnostic AUC: 0.962 (internal), 0.857 (external)
  - Pathways: C-type lectin receptor, PI3K-AKT, B-cell receptor signaling

### 4.2 Remon et al. — IS vs. ICH Discrimination (Biomarker Res 2025)
- **Citation:** Remon A, et al. "Plasma proteomics uncovers divergent molecular signatures in ischemic stroke and intracerebral hemorrhage." *Biomarker Res* 13 (2025).
- **Design:** 388 patients (344 IS, 44 ICH); 3,072 Olink proteins; within 6h of onset
- **Key findings:**
  - 67 significant after FDR
  - **GFAP: top discriminator (AUC 0.887)**
  - BCAN (0.820), SNAP25 (0.797), SPOCK1 (0.786) for ICH vs IS

---

## 5. Convergent Protein Biomarkers (Cross-Study Summary)

### Tier 1: Replicated in ≥3 Studies

| Protein | Full Name | Studies | Role |
|---------|-----------|--------|------|
| **GDF-15** | Growth Differentiation Factor 15 | Zhong, Ho, Lind, Gadd | Stress/inflammation; top predictor |
| **NT-proBNP** | N-terminal pro-BNP | Zhong, Ho, Johansen, Lind | Cardiac dysfunction; embolic stroke-specific |
| **MMP12** | Matrix Metalloproteinase 12 | Chong 2022, Li 2024, Cross-ancestry 2025 | MR-confirmed causal; inflammation |
| **IGFBP4** | IGF Binding Protein 4 | Zhong, Ho | Growth factor regulation |

### Tier 2: MR-Confirmed Causal

| Protein | Causal Direction | Therapeutic Potential |
|---------|-----------------|---------------------|
| **TFPI** | Protective | Coagulation pathway |
| **IL6RA** | Protective | Tocilizumab (existing drug) |
| **CD40** | Protective | Immunotherapy target |
| **F11** | Risk factor | Factor XI inhibitors in trials |
| **FURIN** | Risk factor | Blood pressure pathway |
| **FGF5** | Risk factor | Vascular biology |
| **PROCR** | Protective | Anticoagulant pathway |

### Tier 3: Novel/Emerging

| Protein | Context | Significance |
|---------|---------|-------------|
| **PLAUR/suPAR** | Risk score | Fibrinolysis/inflammation |
| **CLEC4D** | Diagnostic | Innate immunity; C-type lectin |
| **MASP1** | Diagnostic | Complement activation |
| **GFAP** | Differential Dx | IS vs ICH (AUC 0.887) |
| **CCL27/TNFRSF14** | Recurrence | Stroke progression-specific |

---

## 6. Relevance to Our Study

### scRNAseq ↔ Proteomics Cross-Validation Opportunities

| Proteomic Finding | Our scRNAseq Data Point | Validation Strategy |
|------------------|------------------------|-------------------|
| GDF15 ↑ in stroke risk | GDF15 expression in monocytes/macrophages | DEG L1 (HC vs IS), FGS |
| MMP12 causal | MMP12 in myeloid cells | DEG consensus, cell-type specific |
| IL6RA protective | IL6R expression across cell types | CellChat pathway analysis |
| CD40 protective | CD40 on B cells, DCs | CellChat, DEG |
| PLAUR/suPAR risk | PLAUR in monocytes | DEG, FGS gene overlap |
| CLEC4D diagnostic | CLEC4D in myeloid | FindAllMarkers check |
| NT-proBNP embolic | NPPB not in PBMCs (cardiac) | Cannot validate directly |
| F11 risk | F11 not typically in PBMCs | Cannot validate directly |

### Integration Strategy
1. **Gene-level overlap**: Compare DEG gene lists with MR-confirmed causal proteins
2. **Pathway-level**: Map CellChat pathways to protein interaction networks
3. **Cell-type resolution**: Which cells express the proteomic biomarkers?
4. **g3 stratification**: Do MR-confirmed proteins differ between IAT responders vs non?
