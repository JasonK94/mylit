# Key Biomarkers Summary & Connection to Our Data

## 1. Master Biomarker Table (Cross-Study Convergence)

### Tier 1: Most Replicated Protein Biomarkers (≥3 independent studies)

| Protein | Full Name | Evidence | Mechanism | Our Data Link |
|---------|-----------|----------|-----------|---------------|
| **GDF-15** | Growth Differentiation Factor 15 | 4+ UKB studies; top LASSO coefficient | Mitochondrial stress, inflammation | Check DEG in monocytes/macrophages |
| **NT-proBNP** | N-terminal pro-BNP | 4+ studies; embolic stroke-specific | Cardiac dysfunction | Cardiac origin — no direct PBMC validation |
| **MMP12** | Matrix Metalloproteinase 12 | 3 MR studies confirm causal | ECM remodeling, inflammation | Check DEG in myeloid cells |
| **IGFBP4** | IGF Binding Protein 4 | Zhong 2025, Ho 2025 | Growth factor regulation | Check expression across cell types |
| **IL6/IL6RA** | Interleukin-6 / Receptor | MR-confirmed causal; tocilizumab target | Central inflammatory mediator | IL6R expression in CellChat pathways |

### Tier 2: MR-Confirmed Causal Proteins (Drug Targets)

| Protein | Direction | Drug Potential | Status |
|---------|-----------|---------------|--------|
| **F11** (Factor XI) | Risk ↑ | Factor XI inhibitors | Phase 3 trials (abelacimab) |
| **IL6RA** | Protective | Tocilizumab | FDA-approved (other indications) |
| **CD40** | Protective | Anti-CD40 antibodies | Immunotherapy target |
| **TFPI** | Protective | Coagulation modulator | Preclinical |
| **FURIN** | Risk ↑ | Blood pressure pathway | Mechanistic insight |
| **FGF5** | Risk ↑ | Vascular biology | Mechanistic insight |
| **PROCR** | Protective | Anticoagulant | Activated protein C pathway |
| **SCARA5** | Novel target | Cardioembolic stroke-specific | Early discovery |

### Tier 3: scRNAseq-Identified Gene Biomarkers

| Gene | Cell Type | Direction in Stroke | Evidence |
|------|-----------|-------------------|----------|
| **S100A8/A9** | Neutrophil/Monocyte | ↑↑ | DAMP; multiple scRNAseq + proteomics |
| **S100A6** | Myeloid | ↑ | ML + scRNAseq (2024) |
| **ANXA3** | Myeloid | ↑ | ML + scRNAseq (2024) |
| **NCF4** | Neutrophil/Monocyte | ↑ | NADPH oxidase; ML biomarker |
| **STAT1** | Multiple | ↑ | TF hub in scRNAseq |
| **ISG15** | ISG-hi cells | ↑ | Interferon response signature |
| **CXCL9/CXCL13** | Myeloid/B zone | ↑ | Chemokine hubs |
| **LGALS9** | Macrophage | ↑ | CellChat: LGALS9-CD44 pathway |
| **CTSS** | CD14+ Monocytes | ↑ (peaks D3) | Cathepsin S; BBB disruption via JAM degradation |
| **SPP1** | Microglia/MdC | ↑ | Osteopontin; activation marker |
| **IGF1** | Repair macrophages | ↑ (D14) | MdC repair phase |
| **CX3CR1** | NK cells | ↑ | NK activation in mild stroke |
| **IFIT3/RSAD2** | Multiple | ↑ | ISG signature across cell types |

### Tier 4: Clinical Blood Biomarkers (Established)

| Biomarker | Formula/Measurement | Best Cutoff | AUC | Use Case |
|-----------|-------------------|-------------|-----|----------|
| **NLR** | Neutrophil/Lymphocyte | 4.30 (24h post-EVT) | ~0.70 | Prognosis, EVT outcome |
| **SII** | Platelet×Neutrophil/Lymphocyte | 652.73 | 0.698 | 90-day adverse outcome |
| **MHR** | Monocyte/HDL | 0.51 | ~0.65 | Death, disability |
| **PLR** | Platelet/Lymphocyte | 127.01 | 0.643 | 3-month outcome |
| **CRP** | C-reactive protein | 40.58 mg/L (24h) | 0.812 | Post-EVT prognosis |
| **IL-6** | Interleukin-6 | — | — | Infarct volume, NIHSS |

---

## 2. Cell Type Composition Changes (Literature Consensus)

### Stroke vs. Healthy Control (L1 comparison)

| Cell Type | Direction | Strength | Mechanism | Our MASC/MILO/scCODA |
|-----------|-----------|----------|-----------|---------------------|
| **Classical Monocytes (CD14+)** | ↑↑ | Very strong | BM mobilization, CCL2-CCR2 | Check L1 DA |
| **Intermediate Monocytes** | ↑ | Moderate | CNS damage marker | Check L1 DA |
| **Non-classical Monocytes** | ↓ | Moderate | SIDS, infection risk | Check L1 DA |
| **Neutrophils** | ↑↑ | Very strong | First responders | Check L1 DA |
| **B Cells** | ↓↓ | Strong | SIDS, apoptosis (43→11%) | Check L1 DA |
| **CD4+ T Cells** | ↓ | Strong | SIDS, lymphopenia | Check L1 DA |
| **CD8+ T Cells** | ↓→recover | Moderate | Transient lymphopenia | Check L1 DA |
| **NK Cells** | ↑ (mild stroke) | Moderate | Redistribution | Check L1 DA |
| **Tregs** | Variable | Conflicting | Neuroprotective | Check within-CD4 FGS |
| **Gamma-delta T** | ↑ (brain) | Mouse only | IL-17, not in PBMCs | N/A |

### IAT Responder vs Non-responder (L2 comparison)

| Finding | Literature | Our Data |
|---------|-----------|----------|
| CD14+ Monocyte is key | Strong support (monocyte subtypes predict outcome) | **Only DA cell type in L2** (MILO+scCODA) |
| NLR predicts EVT outcome | Meta-analysis confirmed | Can compute from our data |
| Post-procedure inflammation | CRP, IL-6 elevation | DEG L2 NEBULA results |
| Monocyte HLA-DR | ↓ HLA-DR = immunodepression | Check HLA-DR expression in L2 |

---

## 3. CCI Pathways to Validate

### From Literature → Check in Our CellChat

| Pathway | Source → Target | Expected in Stroke | Our Data |
|---------|----------------|-------------------|----------|
| LGALS9-CD44 | Macrophage → multiple | ↑↑ | Check L1 CellChat |
| IL18 signaling | NK/T → multiple | ↑ | Check L1 CellChat |
| XCL1 | NK/T → DC | ↑ | Check L1 CellChat |
| CCL2-CCR2 | Tissue → Monocyte | ↑ (recruitment) | Check L1 CellChat |
| MIF signaling | Multiple | Sex-dimorphic (HC data) | Check L1/L2 CellChat |
| Complement (C3/C5) | Multiple | ↑ | Check CellChat pathways |

### ML-Predicted Cytokine Predictors of Infarct Volume
CCR4, IFNA2, IL-9, CXCL3, IL-7, CCL4, IL-5, CCR3, TNF-alpha, IL-27
→ Check overlap with our CellChat identified pathways

---

## 4. FGS × Proteomics Gene Overlap

### Genes to Check in Our FGS Results

| Gene | Proteomic Context | Expected in FGS? |
|------|------------------|-------------------|
| GDF15 | Top UKB predictor | May appear in monocyte FGS |
| MMP12 | MR-confirmed causal | Check myeloid FGS |
| PLAUR | Top LASSO coefficient | Check monocyte FGS |
| CLEC4D | Diagnostic biomarker | Check myeloid FGS |
| MASP1 | Complement; diagnostic | Check B cell/myeloid FGS |
| IL6 | Central mediator | Check monocyte FGS |
| S100A8/A9 | DAMP/alarmin | Likely in inflammatory monocyte FGS |
| CD40 | MR-confirmed causal | Check B cell/DC FGS |

---

## 5. Novel Contribution Framing

### What Literature DOES NOT Have (= Our Novel Contributions)

1. **Large-scale PBMC scRNAseq in stroke** (226K cells, 100 patients)
   - Most studies: mouse brain tissue or small PBMC cohorts
   - Ours: largest human PBMC stroke scRNAseq to date

2. **IAT response at single-cell resolution**
   - NO prior scRNAseq study examines EVT/IAT responder vs non-responder
   - g3 stratification is entirely novel at single-cell level

3. **Multi-method DA concordance** (MASC × MILO × scCODA × MELD)
   - Literature typically uses 1 method; we use 3-4 methods

4. **scRNAseq ↔ UK Biobank proteomics cross-validation**
   - Bridge single-cell gene expression and population-level circulating proteins
   - Validate MR-confirmed proteins at cell-type resolution

5. **HC aging atlas (549K cells, 96 patients) as reference**
   - Independent healthy aging reference for age-confounding correction
   - Age-related immune changes (from hc_only_v1) distinguish aging from stroke effects

6. **Dual-layer design with HC vs IS + g3 stratification**
   - Separates disease effect (L1) from outcome effect (L2)
   - g3 analysis within IS only → removes disease confounding

---

## 6. Suggested Paper Framing

### Title Ideas
- "Single-cell transcriptomic atlas of peripheral immune landscape in ischemic stroke reveals monocyte-centric mechanisms of thrombectomy response"
- "Multi-method single-cell analysis of peripheral blood in ischemic stroke: cell-type composition, gene expression, and cell-cell communication predict thrombectomy outcome"

### Key Narrative Arc
1. **Introduction:** SIDS is well-known but cell-type resolution lacking in human PBMCs
2. **Methods:** 226K cells, HC vs IS, g3 stratification, multi-method framework
3. **Results L1:** Cell composition changes (validate against literature: monocyte↑, B cell↓, etc.)
4. **Results L2:** CD14+ Monocyte as ONLY outcome-associated cell type → novel finding
5. **CCI:** Condition-specific pathway rewiring (validate LGALS9-CD44, complement)
6. **Proteomics integration:** MR-confirmed proteins expressed in specific cell types
7. **Discussion:** Monocyte-centric model of EVT response; therapeutic implications (IL6RA, CD40)
