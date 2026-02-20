# Stroke Immunity, Cell-Cell Interaction, and Treatment Response

## 1. Innate Immunity in Stroke

### 1.1 Neutrophils

**Temporal dynamics:** First immune cells to arrive at ischemic brain; detectable in microvessels within 1 hour; peak at days 1-3; majority of immune cells in ischemic hemisphere by day 3.

**Pathogenic mechanisms:**
- BBB damage via ROS and MMPs
- Pro-inflammatory cytokines (IL-1beta, TNF-alpha)
- **NETs (Neutrophil Extracellular Traps)**: Present in virtually ALL ischemic stroke thrombi; impair tPA-induced thrombolysis; DNase 1 improves lysis of NET-rich thrombi
- Neutrophil-microglia crosstalk exacerbates ischemia-reperfusion injury (Jin et al., *Adv Sci* 2025)

**Key references:**
- [Role of Immune Cells Migrating to the Ischemic Brain](https://www.ahajournals.org/doi/10.1161/STROKEAHA.118.021474) — Stroke
- [Neutrophil Mobilization Triggers Microglial Change](https://advanced.onlinelibrary.wiley.com/doi/10.1002/advs.202503722) — Adv Sci 2025
- [NET Content Impairs tPA-Induced Thrombolysis](https://www.ahajournals.org/doi/10.1161/STROKEAHA.117.019896) — Stroke

### 1.2 Monocytes / Macrophages

**Monocyte subtypes:**

| Subset | Markers | Change | Role |
|--------|---------|--------|------|
| Classical | CD14++CD16- | ↑ acute/subacute | 80-90% circulating; CCL2-CCR2 axis infiltration |
| Intermediate | CD14++CD16+ | ↑ acute | CNS tissue damage marker |
| Non-classical | CD14+CD16++ | ↓ | Infection susceptibility |

**Polarization dynamics:**
- Ly6Chi/CCR2+ monocytes infiltrate CNS → M1-like macrophages/Tip-DCs
- Week 1: M1 → M2 transition begins (alternative polarization)
- Monocyte infiltration: detectable 24h, peaks day 4, persists weeks
- **MAFB** accelerates DAMP clearance through MSR1 (Shichita et al., *Nat Med*)

**Key references:**
- [Monocytes in Ischemic Stroke Pathobiology](https://pmc.ncbi.nlm.nih.gov/articles/PMC4761876/) — Frontiers
- [MAFB Prevents Excess Inflammation After Stroke](https://www.nature.com/articles/nm.4312) — Nature Medicine

### 1.3 DAMPs and Danger Signaling

**HMGB1:**
- Released from dying cells as key DAMP
- Recognized by TLR2, TLR4, RAGE
- Activates NF-kB → inflammatory cascade
- Local concentration correlates with infarct volume

**S100A8/A9 (Calprotectin):**
- Released by neutrophils/monocytes
- Associates with leukocyte infiltration via collateral pathways
- Clearance via MSR1/MARCO scavenger receptors

---

## 2. Adaptive Immunity in Stroke

### 2.1 T Cells

**Acute (days 1-3):** Antigen-independent reaction; related to infarct volume

**Subacute (days 3-7+):** Antigen-dependent; neuroantigen-specific T cell responses (e.g., MOG35-55) → Th1/Th17 exacerbation

**Th1/Th2 shift:** Stroke induces shift from Th1 (inflammatory) → Th2 (humoral anti-inflammatory) — key feature of SIDS

**Regulatory T cells (Tregs):**
- CD4+CD25+Foxp3+ → neuroprotective (acute) + brain repair (chronic)
- IL-10 is essential for Treg neuroprotective effects
- Adoptive transfer of 2M Tregs: acute protection + long-term recovery; therapeutic window up to 24h
- **Treg exhaustion markers** not significant by age (our HC data confirms this)

**Gamma-delta T cells:**
- Pathogenic via IL-17 secretion
- Activated by IL-1beta, IL-23, IL-18
- Promote M1 microglial polarization
- Recognize DAMPs via TLRs

**Key references:**
- [Immune Pathways in Ischemic Stroke](https://www.ahajournals.org/doi/10.1161/CIRCRESAHA.121.319994) — Circ Res 2021
- [Regulatory T Lymphocytes as Therapy](https://pmc.ncbi.nlm.nih.gov/articles/PMC10239790/)
- [Gamma-delta T Cell in Cerebral Ischemia](https://www.frontiersin.org/journals/neurology/articles/10.3389/fneur.2022.842212/full)

### 2.2 B Cells

**Dual role:**
- Acute: B-cell overshoot can be protective (via IL-10)
- Chronic: B cell clusters in stroke lesions (weeks-months) → **cognitive decline** and **post-stroke dementia**

**Autoimmune contribution:**
- Stroke triggers antibodies against CNS self-antigens (e.g., myelin basic protein)
- **B cell ablation prevents post-stroke cognitive impairment** (mice)
- Potentially treatable with FDA-approved B cell-targeting drugs

**Key references:**
- [B-Lymphocyte-Mediated Delayed Cognitive Impairment](https://www.jneurosci.org/content/35/5/2133) — J Neurosci 2015
- [Updates on B-cells in Ischemic Stroke](https://www.frontiersin.org/journals/cellular-neuroscience/articles/10.3389/fncel.2024.1340756/full)

### 2.3 NK Cells

**Mechanisms of brain injury:**
- Necrosis of neural cells via IFN-gamma
- BBB damage in response to IP-10
- Perforin-mediated neuronal killing
- NKG2D upregulated in ischemic hemisphere

**Peripheral changes:**
- Spleen atrophy → peripheral NK cell depletion
- Anti-NKG2D antibody alleviates stroke outcome

**Key references:**
- [Ischemic Neurons Recruit NK Cells](https://www.pnas.org/doi/10.1073/pnas.1315943111) — PNAS
- [NKG2D Signaling in Cerebral Ischemia](https://www.ahajournals.org/doi/10.1161/JAHA.122.029529) — JAHA

---

## 3. Cell-Cell Communication (CCI) in Stroke

### 3.1 Key Signaling Pathways from Literature

| Pathway | Source Cell | Target Cell | Role |
|---------|-----------|-------------|------|
| **LGALS9-CD44** | Microglia/Macrophage | Multiple | Critical post-ischemic signaling |
| **IL18** | NK/T cells | Multiple | Secretory axis |
| **XCL1** | NK/T cells | DCs | Secretory axis |
| **CCL2-CCR2** | Damaged tissue | Classical monocytes | Monocyte recruitment |
| **CX3CL1-CX3CR1** | Neurons | Microglia | Neuroprotection/damage |
| **HMGB1-TLR2/4/RAGE** | Dying cells | Immune cells | DAMP signaling |

### 3.2 ML-Predicted Cytokine Predictors of Infarct Volume
CCR4, IFNA2, IL-9, CXCL3, IL-7, CCL4, IL-5, CCR3, TNF-alpha, IL-27

### 3.3 Relevance to Our CellChat Results
- Our L1 CellChat identifies condition-specific pathway changes → validate against LGALS9-CD44, IL18
- Our rankNet shows pathway strength differences → compare with +104%/+95% from mouse data
- g3 stratification adds EVT-response dimension not in literature

---

## 4. Key Signaling Pathways

### 4.1 NLRP3 Inflammasome
- **Two-signal model:** (1) NF-kB/MAPK priming → NLRP3/pro-IL-1beta; (2) Assembly trigger → caspase-1 → mature IL-1beta/IL-18
- Triggers: K+ efflux, ROS, mitochondrial dysfunction, ATP, DAMPs
- NLRP3 deficiency attenuates ischemia-reperfusion injury

### 4.2 JAK-STAT
- STAT3 → NLRP3 promoter acetylation
- JAK2/STAT3 inhibition ameliorates neuroinflammation

### 4.3 NF-kB
- Central hub: HMGB1 → TLR2/4/RAGE → NF-kB → pro-IL-1beta, NLRP3
- Connects DAMPs, inflammasome, and cytokine production

### 4.4 Complement System

| Phase | Timing | C3aR Activity |
|-------|--------|---------------|
| Acute | Minutes-hours | Immune cell recruitment, BBB permeability, edema |
| Subacute | Hours-days | Modulate immune cells, secondary injury |
| Chronic | Days-weeks | Prolonged neuroinflammation, synaptic pruning |

- C3 → C3a (microglia/astrocyte activation) + C3b (opsonization)
- C5 → C5a (neuronal apoptosis) + C5b (MAC)
- IgM-dependent complement → NLRP3 in endothelial cells

---

## 5. Stroke-Induced Immunodepression Syndrome (SIDS)

### 5.1 Core Features

| Feature | Description |
|---------|-------------|
| Lymphopenia | Rapid, persistent decrease (blood, spleen, thymus) |
| Th1→Th2 shift | Loss of cell-mediated immunity |
| NK depletion | ↓ numbers and function |
| Monocyte deactivation | ↓ HLA-DR expression |
| Splenic atrophy | Apoptotic loss of lymphocytes |
| IgG deficiency | Acquired immunoglobulin G deficiency |

### 5.2 Mechanisms — Three Neuroimmune Axes
1. **SNS (Sympathetic):** Catecholamine-mediated immunosuppression
2. **HPA axis:** Cortisol-mediated
3. **Vagus nerve:** Cholinergic anti-inflammatory pathway

### 5.3 Clinical Consequences
- Pneumonia: 12%, UTI: 16% (most common post-stroke infections)
- Pneumonia has greatest effect on neurological recovery
- SIDS patients show less NIHSS improvement
- NLR within 24h predicts pneumonia (but NOT UTI)
- Stroke location matters: insular cortex + large infarcts → more immunodepression

**Key references:**
- [Stroke-Induced Immunosuppression: Prevention and Prediction](https://jneuroinflammation.biomedcentral.com/articles/10.1186/s12974-021-02177-0)
- [SIDS as Predictor of Functional Outcome](https://www.nature.com/articles/s41598-024-58562-1)

---

## 6. EVT/IAT Response Biomarkers

### 6.1 Laboratory Biomarkers

| Biomarker | Threshold/Finding | Predictive Value |
|-----------|-------------------|-----------------|
| **NLR (24h)** | >4.30 | Poor prognosis post-EVT |
| **NLR (3-7d f/u)** | Lower = better | Independent predictor of favorable outcome |
| **CRP (24h)** | >40.58 mg/L | AUC 0.812 for poor prognosis |
| **SII** | >652.73 | OR 2.915 for 90-day adverse outcome |
| **PLR** | >127.01 | AUC 0.643 (lower than NLR/SII) |

**Important:** Post-procedure NLR is MORE informative than admission NLR

### 6.2 Imaging Biomarkers
- mCCT (microvascular cerebral circulation time): predicts malignant edema
- Ischemic core volume on CTP
- CBV index: low ASPECTS + high CBV → may still benefit from EVT
- Pulse pressure >57.39 mmHg → predicts symptomatic ICH

### 6.3 Our g3 Variable Context
- g3 = IAT response outcome (responder=1 vs nonresponder=2)
- **L2 analysis** (IS only, g3=1 vs g3=2) directly maps to EVT response
- Our finding: **CD14+ Monocyte is the ONLY DA cell type in L2** (MILO + scCODA concordant)
  - This aligns with monocyte's critical role in post-EVT inflammation
  - Post-procedure monocyte activation may determine outcome

---

## 7. Well-Known Stroke Immune Biomarkers

### 7.1 NLR (Neutrophil-to-Lymphocyte Ratio)
- Meta-analysis (41 articles, 27,124 individuals): independently predicts prognosis
- EVT: 24h NLR >4.30 predicts poor outcome
- Follow-up NLR (3-7d) more informative than admission
- Predicts intracranial hemorrhage after EVT

### 7.2 MHR (Monocyte-to-HDL Ratio)
- Highest vs lowest quartile: OR 1.97 for death, 3.34 for disability
- Best cutoff: 0.51 (sensitivity 62.3%, specificity 66.6%)
- Predicts stroke-associated pneumonia

### 7.3 SII (Systemic Immune-Inflammation Index)
- Formula: Platelets × Neutrophils / Lymphocytes
- High SII: OR 2.915 for 90-day adverse outcomes
- Better reflects inflammation-immunosuppression interplay

### 7.4 Other Blood Biomarkers

| Biomarker | Association | Context |
|-----------|-------------|---------|
| IL-6 | ↑ with infarct volume & NIHSS | Acute phase |
| CRP | Poor prognosis at 1 & 3 months | Widely available |
| Fibrinogen | Independent predictor of poor prognosis | Coagulation + inflammation |
| Copeptin | Predicts infection within 5 days | Stress response |
| Procalcitonin | Predicts stroke-associated pneumonia | Infection-specific |

---

## 8. Temporal Dynamics Summary

| Time | Key Events |
|------|------------|
| **Min-Hours** | DAMPs (HMGB1, S100A8/A9); complement activation; microglial activation; neutrophil arrival |
| **Hours-Day 1** | Neutrophil infiltration begins; monocyte infiltration begins; BBB breakdown; NF-kB/NLRP3 |
| **Days 1-3** | Neutrophil peak; gamma-delta T cells (IL-17); NK-mediated neuronal killing; SIDS onset |
| **Days 3-7** | M1→M2 transition; monocyte peak (day 4); antigen-dependent T cell response; Treg expansion |
| **Week 1-2** | Alternative macrophage polarization; adaptive immunity established |
| **Weeks-Months** | B cell clusters → autoantibodies → cognitive decline; immune dysfunction persists ~6 weeks |

---

## 9. Key Review Papers

1. **Iadecola C, Buckwalter MS, Anrather J.** "Immune responses to stroke: mechanisms, modulation, and therapeutic potential." *J Clin Invest* 130(6):2777-2788 (2020). [PMC7260029](https://pmc.ncbi.nlm.nih.gov/articles/PMC7260029/)

2. **Iadecola C, Anrather J.** "The immunology of stroke: from mechanisms to translation." *Nature Medicine* 17:796-808 (2011). [PMC3137275](https://pmc.ncbi.nlm.nih.gov/articles/PMC3137275/)

3. **Shi K et al.** "Immune Pathways in Etiology, Acute Phase, and Chronic Sequelae of Ischemic Stroke." *Circulation Research* 130(4):527-543 (2022). [Link](https://www.ahajournals.org/doi/10.1161/CIRCRESAHA.121.319994)

4. **Meisel C et al.** "Stroke-induced immunodepression: experimental evidence and clinical relevance." (2005).

5. **Doyle KP et al.** "B-Lymphocyte-Mediated Delayed Cognitive Impairment following Stroke." *J Neurosci* 35(5):2133 (2015). [Link](https://www.jneurosci.org/content/35/5/2133)

6. **Montaner J et al.** "Multilevel omics for the discovery of biomarkers and therapeutic targets for stroke." *Nat Rev Neurol* 16:247-264 (2020). [Link](https://www.nature.com/articles/s41582-020-0350-6)

7. **Simats A et al.** "Stroke Proteomics: From Discovery to Diagnostic and Therapeutic Applications." *Circ Res* 130 (2022). [Link](https://www.ahajournals.org/doi/10.1161/CIRCRESAHA.122.320110)
