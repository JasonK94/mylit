# Prognosis Factor Discovery Strategy — stroke_hc_v8_2

> Created: 2026-02-25
> Purpose: PBMC scRNAseq에서 ischemic stroke prognosis factor를 어떻게 도출할 것인가

---

## 1. The Central Challenge: PBMC → CNS Gap

Ischemic stroke는 CNS event이지만, PBMC를 분석한다.
이 indirectness를 극복하는 것이 논문의 핵심 narrative.

### 1.1 Why PBMC Works for Stroke (Literature-Based Justification)

| Mechanism | Evidence | Our Data Connection |
|-----------|----------|-------------------|
| **BBB Breakdown** | Stroke 시 수시간 내 BBB 파괴, 말초 면역세포 CNS 직접 침투 | Monocyte = brain macrophage precursor (Garcia-Bonilla 2024) |
| **DAMP Signaling** | 죽은 뉴런 → HMGB1, S100A8/A9 → 말초 TLR 활성화 | S100A8/A9 = top MOFA+ weight, top MNN interaction, top gene dynamics |
| **SIDS/SIIS** | SNS/HPA axis → lymphopenia + monocyte deactivation | 5-method concordance: lymphopenia + selective innate activation |
| **Trained Immunity** | Stroke → BM HSPC 후성유전 재프로그래밍 → 장기 면역 변화 (Simats 2024, Cell) | Monocyte-centric outcome = trained immunity 반영 |
| **Clinical Consequence** | 말초 면역 상태 → post-stroke infection, cardiac complication 직접 결정 | FGS AUC 0.82-0.89 in bulk datasets |

### 1.2 Reframing: PBMC Is Not "Indirect" — It IS the Disease

SIDS framework에서 말초 면역 dysregulation은:
- CNS의 "window"가 아니라
- **독립적 clinical consequence의 원인** (infection, cardiac dysfunction, mortality)
- 따라서 PBMC 분석 = 직접적 outcome-relevant analysis

---

## 2. Competitive Landscape: Where We Stand

### 2.1 No Existing Study Fills Our Niche

```
Large-scale scRNAseq  ×  Human  ×  Stroke/PBMC  ×  Treatment Outcome
                          ↑                          ↑
                   (mouse만 있음)              (아무도 안 함)
```

| Study | Cells | Patients | scRNAseq? | Outcome? | Gap vs Us |
|-------|-------|----------|-----------|----------|-----------|
| Cho 2022 (SVN) | 26K | 19 | Yes | No | 9x cells, 5x patients, no outcome |
| Garcia-Bonilla 2024 (Nat Immunol) | 43K | mice | Yes | No | Mouse, brain tissue |
| Grandke 2025 (Nat Commun) | 909K | 257 | Yes | No | AD/PD (not stroke), no outcome |
| Stamova/Amini 2023 (J Neuroinfl) | bulk | 56 | No | Yes (mRS) | Bulk, not single-cell |
| **Ours** | **226K** | **100** | **Yes** | **Yes (g3 IAT)** | **All four criteria met** |

### 2.2 Comparable PBMC-Only Precedent

**Grandke et al. 2025 (Nat Commun)** — 가장 직접적 선례:
- PBMC only, no CSF, no brain tissue
- 909K cells, 257 individuals (AD, PD, MCI)
- muscat pseudobulk DEG + scCODA (우리와 동일 방법론)
- BBB breakdown + brain dataset cross-validation (ZEBRA) 논거

→ **우리가 이 precedent을 따라가되, outcome stratification + stroke specificity로 차별화**

---

## 3. Prognosis Discovery: 3-Stage Pipeline

### Stage 1: Compositional — WHO is different? (Cell type level)

**Question**: Good vs Bad outcome에서 어떤 세포가 비율적으로 다른가?

| Method | L2 Finding | Concordance |
|--------|-----------|-------------|
| scCODA (Bayesian) | CD14+ Mono **only** credible | ✅ |
| MILO (neighbourhood) | 19% DA nhoods, CD14+ Mono dominant | ✅ |
| MELD (continuous) | P(Bad) highest: CD14+ Mono 0.792 | ✅ |
| Augur (expression-based) | Near random (0.50-0.62) except Mono | ✅ |
| MASC (mixed-effects) | CD14+ Mono significant | ✅ |

**Conclusion**: **5-method concordance → CD14+ Monocyte is the sole compositional prognosis factor**

### Stage 2: Transcriptional — WHAT is different? (Gene program level)

**Question**: Monocyte 안에서 어떤 유전자/프로그램이 다른가?

| Analysis | L2 Key Finding |
|----------|---------------|
| DEG (NEBULA L2) | 89 DEGs, CD14+ Mono에 집중 |
| cNMF GEPs | TNFa/NFkB program (GEP4/5) → Bad enriched |
| FGS (n=200) | TOP25_DOWN AUC 0.82-0.89 (external) |
| Gene dynamics (GAMM) | 10 genes × 3 methods (Mono = only compartment with g3 effect) |
| CellChat (CCI) | ALL pathways ↑ in Bad; IL1 = Bad-only pathway |
| MNN v3 | S100A8-CD36 top interaction |
| LIANA | PF4-LDLR lost in Bad; SNCA-LAG3 gained |

**Conclusion**: **Maladaptive TNFa/NFkB + alarmin (S100A/A9) hyperactivation in Bad outcome monocytes**

### Stage 3: Subcellular — WHICH subpopulation? (Subcluster level)

**Question**: 어떤 monocyte subcluster가 Bad outcome을 drive하는가?

| Compartment | Bad-Enriched SC | OR | Biological Identity |
|------------|----------------|-----|-------------------|
| **Monocyte** | SC2 | 403 | Inflammatory classical |
| **Monocyte** | SC3 | 193 | TNFa/NFkB activated |
| CD4 T | SC3 | 18.6 | Th1 activated |
| CD8 T | SC8 | 20.6 | Terminal effector |

**Conclusion**: **Mono SC2/SC3 = pathogenic monocyte subpopulations (OR > 100)**

---

## 4. Paper Narrative Arc

### 4.1 Main Story (5 beats)

```
Beat 1: SIIS Establishment (L1)
  "Stroke induces SELECTIVE immunosuppression — not global shutdown
   but a rewiring from adaptive toward innate inflammatory pathways"
  → Figure 1 (overview) + Figure 2 (DA concordance)

Beat 2: Monocyte as Central Orchestrator (L1→L2 bridge)
  "Across all 19 analyses, monocyte emerges as THE cell type that
   (a) is most transcriptomically altered in stroke,
   (b) is the ONLY type with compositional difference by outcome,
   (c) is the ONLY compartment with trajectory-level g3 effect"
  → Figure 3 (DEG/FGS) + transition to L2

Beat 3: Outcome Determination is Monocyte-Intrinsic (L2)
  "Good vs Bad IAT outcome is not determined by global immune
   composition but by monocyte-specific gene programs:
   - Bad = maladaptive hyperactivation (TNFa/NFkB, IL1, all CCI ↑)
   - Good = controlled inflammatory response"
  → Figure 4 (CCI) + Figure 5 (trajectory/gene dynamics)

Beat 4: Pathogenic Subpopulations (L2 deep-dive)
  "Subclustering identifies specific monocyte substates (SC2/SC3,
   OR > 100) that are massively enriched in bad outcome —
   these express TNFa/NFkB + alarmin programs"
  → Figure 6 (subclustering + cNMF)

Beat 5: Clinical Translation (Validation)
  "A 25-gene signature derived from this single-cell outcome
   comparison achieves AUC 0.82-0.89 in 3 independent bulk
   cohorts — demonstrating that PBMC scRNAseq findings
   translate to clinically actionable biomarkers"
  → Figure 3 Panel D (validation ROC curves)
```

### 4.2 The "Indirect" Problem — How We Address It

**In Introduction:**
- Acknowledge upfront: "Ischemic stroke is a CNS event; PBMC profiling captures the systemic, not local, immune response"
- Immediately reframe: "However, peripheral immune dysregulation is not merely a mirror of CNS events — it is itself a major determinant of clinical outcome through..."
  - Post-stroke infection (Meisel 2005)
  - Cardiac complications (Simats 2024)
  - Immune cell trafficking (Garcia-Bonilla 2024)
  - DAMP-mediated peripheral activation (S100A8/A9, TLR pathway)

**In Discussion:**
- Cite Grandke 2025 (Nat Commun): PBMC-only → Nat Commun, no CSF/brain, with BBB argument
- Cite Garcia-Bonilla 2024: Blood monocyte = brain macrophage precursor
- Cite Simats 2024: Trained immunity in BM — stroke의 long-term effect는 말초에서 시작
- Our unique strength: outcome stratification (아무도 안 했음)

---

## 5. Additional Analyses to Strengthen

### Priority 1: Brain Cross-Reference (Medium effort, High impact)

**Goal**: Garcia-Bonilla 2024 mouse brain scRNAseq data → 우리 PBMC monocyte program genes overlap 확인

**Method**:
1. GEO에서 Garcia-Bonilla data 다운 (public)
2. Mouse brain macrophage DEG list 추출
3. Ortholog mapping → human gene name
4. 우리 L2 monocyte DEG / cNMF GEP genes와 overlap 분석
5. Fisher exact test for enrichment

**Expected output**: "X% of our outcome-associated monocyte genes are also differentially expressed in brain-infiltrating macrophages post-stroke"

### Priority 2: Clinical Variable Correlation (Low effort, High impact)

**Goal**: FGS score ↔ NIHSS, mRS, DWI volume, sampling timepoint 등 임상 변수 상관

**Requires**: 임상 데이터 접근 (현재 g3만 available)
- NIHSS admission score
- mRS 90-day score
- Onset-to-sampling time
- DWI lesion volume
- Recanalization status (TICI)

### Priority 3: Drug Target Mapping (Low effort, Medium impact)

**Goal**: TNFa/NFkB program genes → druggable targets

**Method**:
1. Mono GEP4/GEP5 top 50 genes
2. DrugBank / ChEMBL / TTD mapping
3. Existing clinical trials for identified targets in stroke

### Priority 4: Deconvolution-Based Validation (Medium effort, Medium impact)

**Goal**: CIBERSORTx로 우리 scRNAseq → signature matrix → bulk dataset에서 monocyte fraction 추정

**Method**:
1. Build signature matrix from our 21 cell types
2. Apply to GSE16561, GSE22255, GSE58294
3. Estimated monocyte fraction ↔ stroke status/outcome 검증

---

## 6. Figure–Narrative Mapping (Updated)

| Figure | Narrative Beat | Key Message |
|--------|---------------|-------------|
| Fig 1 | Overview | 226K cells, 100 pts, dual-layer design |
| Fig 2 | SIIS (L1) | 5-method DA concordance → selective immunosuppression |
| Fig 3 | DEG + FGS | Monocyte-dominant DEG, g3-derived signature, AUC 0.89 validation |
| Fig 4 | CCI | Global CCI collapse + selective RESISTIN/MIF activation; Bad = all ↑ |
| Fig 5 | Trajectory | Mono = only g3 trajectory effect; 10 genes × 3 methods |
| Fig 6 | Subcluster + Programs | Mono SC2/SC3 (OR>100 Bad); TNFa/NFkB cNMF program |

---

## 7. Key Sentences for Abstract

> "We present the largest human PBMC single-cell RNA-seq atlas of ischemic stroke
> (226,855 cells, 100 patients), uniquely incorporating thrombectomy treatment
> outcome stratification. Multi-method analysis reveals that stroke induces
> selective immunosuppression — global adaptive immune communication collapse
> coexists with paradoxical innate inflammatory pathway activation
> (RESISTIN, MIF, S100A8/A9). Among IAT-treated patients, prognosis is
> determined exclusively by monocyte biology: CD14+ monocytes are the sole
> compositionally altered cell type, with maladaptive TNFa/NFkB
> hyperactivation distinguishing poor from good outcomes. A 25-gene
> peripheral blood signature derived from these findings achieves
> AUC 0.82-0.89 in three independent cohorts."

---

## References (Key Precedents)

1. Garcia-Bonilla et al. Brain+blood scRNAseq in stroke. *Nat Immunol* 2024.
2. Grandke et al. PBMC-only atlas in neurodegeneration. *Nat Commun* 2025.
3. Phongpreecha et al. PBMC immunoprofiling in AD/PD. *Sci Adv* 2020.
4. Simats et al. Trained immunity after brain injury. *Cell* 2024.
5. Stamova/Amini et al. Blood gene expression & stroke outcome. *J Neuroinflammation* 2023.
6. Cho et al. First human PBMC scRNAseq in stroke. *SVN* 2022.
7. Meisel et al. CIDS framework. *Nat Rev Neurosci* 2005.
8. Stephenson et al. COVID PBMC 780K cells. *Nat Med* 2021.
