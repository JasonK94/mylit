# Paper Strategy — stroke_hc_v8_2

> Created: 2026-02-21 | 목적: 논문 작성을 위한 남은 작업 정의 및 에이전트 할당

---

## 현재 상태 요약 (Current State)

### 완료된 분석 (19 Phases)
모든 핵심 computational 분석이 완료됨:

| Category | Methods | Status |
|----------|---------|--------|
| **Integration** | scVI, scANVI | ✅ |
| **Frequency/DA** | MASC, scCODA, MILO, MELD, Augur, propeller | ✅ |
| **DEG** | muscat-edgeR, muscat-DESeq2, NEBULA, consensus | ✅ |
| **CCI** | CellChat v2, MNN v3, LIANA | ✅ |
| **Unbiased** | FGS (whole + 15 cell types), cNMF (16 types + downstream), MOFA+ | ✅ |
| **Trajectory** | Slingshot, Monocle3, GAMM, Lamian, pseudobulk GAMM | ✅ |
| **Validation** | External validation (3 bulk datasets), scANVI annotation QC | ✅ |
| **Subclustering** | Mono (12 SC), CD4_T (17 SC), CD8_T (16 SC) + frequency | ✅ |

### 핵심 발견 (Key Findings)

1. **SIIS (Stroke-Induced Immunosuppression)**: CCI -37%, 그러나 RESISTIN/MIF/CypA 선택적 활성화
2. **L2 (g3)**: CD14+ Mono가 유일한 compositionally DA cell type → g3 effect는 expression-level
3. **Trajectory**: Mono가 유일한 g3 effect compartment (10 genes by 3-method consensus)
4. **External validation**: FGS TOP25_DOWN AUC 0.82-0.89 (bulk 3 datasets)
5. **cNMF**: TNFa/NFkB program (Inflammatory Mono GEP5, CD14+ Mono GEP4)이 핵심

---

## 남은 작업 (Remaining Tasks)

### A. 생물학적 해석 (Biological Interpretation)

#### A1. Subcluster Biological Annotation
- **대상**: Mono 12 SC, CD4_T 17 SC, CD8_T 16 SC
- **작업**: markers_all.csv 기반으로 각 subcluster의 생물학적 identity 명명
- **출력**: annotation table (SC → biological name → key markers → literature reference)
- **참고**: stroke_papers/3_STROKE_IMMUNITY_CCI.md (monocyte subtypes: classical/intermediate/non-classical)

#### A2. cNMF GEP Biological Naming
- **대상**: 16 cell types × 8-11 GEPs = 164 programs
- **작업**: GSEA annotation + top genes 기반으로 각 GEP에 생물학적 이름 부여
- **출력**: GEP annotation master table (cell_type, gep, biological_name, top_pathway, condition_association)
- **핵심**: condition-associated GEPs 중 stroke immunology와 연결되는 것 식별

#### A3. Literature Cross-Reference
- **작업**: 우리 데이터 key findings vs 기존 문헌 체계적 비교
- **참고**: stroke_papers/ 5개 문서 (40+ papers, 83 references)
- **출력**: validation/novelty table
  - Validated: 문헌과 일치하는 발견 (e.g., Mono expansion, B cell decrease)
  - Novel: 우리만의 발견 (e.g., g3 stratification, 226K cells, cross-layer concordance)

### B. 시각화 (Figure Generation)

#### B1. Paper Figures v3
- **목표**: 실제 논문에 들어갈 최종 figure panel 생성
- **Figure 1**: Study overview + integration quality
  - UMAP (anno1, split by cohort)
  - DotPlot (canonical markers)
  - Stacked bar (anno2, sorted by FGS score)
  - Patient demographics table
- **Figure 2**: Differential Abundance
  - 3-method concordance heatmap (MASC × scCODA × MILO)
  - MELD UMAP (P(Stroke), P(Bad))
  - scCODA credible effects dotplot
- **Figure 3**: DEG & Gene Signature
  - FGS gene importance + AUC
  - Cross-layer overlap (UpSet plot)
  - Pathway enrichment (HALLMARK/KEGG)
  - External validation ROC curves
- **Figure 4**: CCI
  - CellChat rankNet (information flow)
  - CellChat bubble (top L-R pairs)
  - MNN circos (key interactions)
  - LIANA rank comparison
- **Figure 5**: Trajectory & Gene Dynamics
  - Pseudotime UMAP (3 compartments)
  - Pseudotime density (HC vs Stroke, Good vs Bad)
  - Top gene dynamics curves (10 consensus genes)
  - Effect size comparison (ABC vs Lamian vs pseudobulk)
- **Figure 6**: Subclustering & Gene Programs
  - Subcluster UMAP (Mono/CD4/CD8)
  - Subcluster frequency (HC vs Stroke)
  - cNMF top condition-associated GEPs
  - cNMF × FGS overlap

#### B2. CCI Consolidated Summary
- MNN v3 + CellChat v2 + LIANA 통합 요약 figure
- Pathway-level 통합: 같은 pathway가 여러 method에서 나오는지
- Sender/receiver analysis summary

### C. 문서화 (Documentation)

#### C1. Results Narrative
- 각 figure에 대한 results section 초안 (영어)
- Key finding → statistical evidence → biological interpretation

#### C2. Methods Draft
- 각 분석의 methods section (parameter, tool version, statistical test)

---

## 에이전트 할당 (Agent Task Assignments)

### Agent 1: Subcluster Biological Annotation
```
Input:  subclustering/{Mono,CD4_T,CD8_T}/markers_all.csv
        stroke_papers/3_STROKE_IMMUNITY_CCI.md (monocyte/T cell subtypes)
Task:   Read markers, cross-reference with known biology, assign names
Output: subclustering/subcluster_annotation.csv
        subclustering/annotation_summary.md
```

### Agent 2: cNMF GEP Naming + Summary
```
Input:  unbiased/cnmf/downstream/{all_gsea_results.csv, gep_annotations_hallmark.csv, gep_condition_association.csv}
        unbiased/cnmf/downstream/gsea/{celltype}_top_genes.csv
Task:   Assign biological names to GEPs, identify condition-relevant programs
Output: unbiased/cnmf/downstream/gep_master_annotation.csv
        unbiased/cnmf/downstream/gep_summary.md
```

### Agent 3: Paper Figures v3
```
Input:  All results across reorganized directories
Task:   Generate publication-ready multi-panel figures
Output: figures/v3_paper/{Figure1-6}.png
Script: scripts/descriptive/generate_paper_figures_v3.R
```

### Agent 4: CCI Consolidated Summary
```
Input:  cci/cellchat/, cci/mnn/, cci/liana/, cci/cci_deg_crossref/
Task:   Cross-method integration, pathway concordance, final summary
Output: cci/consolidated_summary.csv
        cci/consolidated_summary.md
```

### Agent 5: Literature Cross-Reference
```
Input:  docs/claude/stroke_papers/ (5 docs, 83 refs)
        docs/claude/7_PROGRESS_STATUS.md (19 findings)
Task:   Systematic comparison of our findings vs literature
Output: docs/claude/9_LITERATURE_CROSSREF.md
```

---

## 우선순위

| Priority | Task | Agent | Reason |
|----------|------|-------|--------|
| 1 | Subcluster annotation | Agent 1 | Subclusters currently unnamed (SC0-16) → need biology |
| 1 | cNMF GEP naming | Agent 2 | GEPs currently unnamed (GEP1-11) → need biology |
| 2 | Literature cross-ref | Agent 5 | Needed for Discussion section |
| 3 | Paper figures v3 | Agent 3 | After annotations complete |
| 3 | CCI summary | Agent 4 | After annotations complete |

---

## 파일 위치 (Updated after reorganization)

```
stroke/ (= /data/user3/sobj/stroke_hc_v8_2/)
├── frequency/       # MASC, milo, meld, augur, sccoda, propeller, scanvi
├── cci/             # cellchat, mnn, liana, cci_deg_crossref, plots
├── unbiased/        # cnmf, fgs, mofa, summary/
├── deg/             # consensus, pathway, FindMarkers, external_validation, external_data
├── trajectory/      # v1/, v2/, v3/
├── subclustering/   # Monocyte, CD4_T, CD8_T, frequency/
├── pipeline/        # P0_PIPE, P1_pt_umap, P2_MASC, solo, integration, annotation_qc
├── figures/         # v1_initial/, v2_paper/
├── scripts/         # cci/, deg/, descriptive/, frequency/, pipeline/, trajectory/, unbiased/
└── logs/
```
