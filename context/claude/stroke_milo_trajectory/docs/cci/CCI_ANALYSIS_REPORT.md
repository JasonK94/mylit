# CCI Analysis Report — stroke_hc_v8_2

> Updated: 2026-02-18 | Cell-Cell Interaction 분석 결과, 해석, 문제점, 향후 계획

---

## 1. 분석 개요 (Overview)

### Dual-Layer Design
- **Layer 1 (L1)**: HC vs Stroke (cohort), anno2 (8 compartments) — 131K cells, 56 samples
- **Layer 2 (L2)**: g3 Good vs Bad (IS only), anno2 — 54K cells, 32 samples

### Methods
| Method | Level | L1 | L2 | Status |
|--------|-------|----|----|--------|
| CellChat | anno2 | ✅ | ✅ | Complete |
| MNN (MultiNicheNet) | anno2 | ✅ | ✅ | Complete (group_prioritization 0 rows) |
| MNN (MultiNicheNet) | anno1 (19 types) | ✅ | ✅ | Complete (group_prioritization 0 rows) |

---

## 2. CellChat 결과 (Results)

### Method 2A: Sample-wise Merged
개별 환자(sample) 단위로 CellChat 실행 후 condition별 merge → 비교.

### L1: HC vs Stroke (SIIS Discovery)

**Summary Statistics**:
| Metric | HC (n=20) | Stroke (n=36) | Difference |
|--------|-----------|---------------|------------|
| Total Interaction Count (mean) | 114.75 | 72.39 | **-42.36 (-37%)** |
| Total Interaction Strength (mean) | 7.35 | 5.00 | **-2.34 (-32%)** |
| Cell Types | 7 | 7 | — |

**핵심 발견: SIIS (Stroke-Induced Immunosuppression)**

1. **Global CCI Collapse**: Stroke에서 총 상호작용 37% 감소
2. **DC Hub Disruption**: DC가 HC에서 주요 CCI hub → Stroke에서 incoming/outgoing 모두 급감
3. **Innate→Adaptive Failure**: Mono/DC → Tc/Bc 경로 단절이 가장 두드러짐
4. **Pathway-level Changes**:
   - MIF: HC에서 strong, Stroke에서 유지 (pro-inflammatory sustain)
   - GALECTIN: HC에서 active → Stroke에서 감소
   - ANNEXIN: 양 조건 모두 존재하나 구조 변화
   - IL16: HC에서만 활성 (CD4 Tc migration signal)

### L2: g3 Good vs Bad (Outcome Prediction)

**Summary Statistics**:
| Metric | g3 Good (n=11) | g3 Bad (n=21) | Difference |
|--------|----------------|---------------|------------|
| Total Interaction Count (mean) | 69.18 | 76.24 | **+7.06 (+10%)** |
| Total Interaction Strength (mean) | 4.71 | 5.19 | **+0.48 (+10%)** |
| Cell Types | 7 | 7 | — |

**핵심 발견: Bad Outcome에서의 면역 과활성화**

1. **Paradoxical CCI Increase**: SIIS (면역억제) 위에 Bad outcome에서 +10% CCI 증가
2. **Mono→DC/Bc Over-activation**: Monocyte에서 DC, B cell로의 signal 과잉 (불균형적 항원제시 활성화)
3. **NKc Functional Suppression**: NK cell의 outgoing CCI 감소 → 세포독성 기능 저하
4. **Pathway-level**:
   - RESISTIN: Bad에서만 활성 (metabolic inflammation marker)
   - MIF: Bad > Good (pro-inflammatory amplification)
   - GALECTIN, ANNEXIN: 구조적 차이보다 양적 차이

### CellChat 생성 플롯 목록

#### Comparison Plots (Script: `run_cellchat_comparison.R`)
```
L1 (cellchat_L1_cohort/):
  01-08: circle_diff, heatmap_diff, heatmap_sidebyside, boxplot (count/weight) — 16 PNG + PDF

L2 (cellchat_L2_g3/):
  01-08: same structure — 16 PNG + PDF
```

#### Native CellChat Plots (Script: `plot_cci_comprehensive.R`)
```
L1 pathways/:
  pathway_{PATHWAY}_{condition}.png — circle plots (14 PNG)
  contrib_{PATHWAY}_{condition}.png — L-R contribution (14 PNG)
  {PATHWAY}_{condition}.pdf — native circle (5 PDF)
  Total: 33 files

L2 pathways/:
  pathway_{PATHWAY}_{condition}.png — circle plots (16 PNG)
  contrib_{PATHWAY}_{condition}.png — L-R contribution (16 PNG)
  {PATHWAY}_{condition}.pdf — native circle (5 PDF)
  Total: 37 files
```

#### Bubble Plots (from merged objects)
```
11_bubble_{condition}.png, 12_bubble_comparison.png
```

---

## 3. MNN (MultiNicheNet) 결과

### anno2 Level
- L1 + L2 실행 완료
- **문제**: `group_prioritization_tbl` 0 rows — anno2가 너무 coarse하여 DE→LR prioritization filter 통과 불가
- 기본적인 DE volcano, ligand activity plots 생성됨

### anno1 Level (19 Cell Types)

#### L1: HC vs Stroke

**DE Genes Summary (p_adj < 0.05)**:
| Cell Type | Total Genes | Sig | Up | Down |
|-----------|------------|-----|-----|------|
| **cDC1** | 1,751 | **594** | 216 | **378** |
| CD4_S100A8_CD14 | 1,004 | 142 | 10 | 132 |
| CD4_S100A8_CSF3R | 983 | 84 | 17 | 67 |
| MAIT | 915 | 61 | 6 | 55 |
| cDC2 | 1,293 | 55 | 2 | 53 |
| ISG+ Myeloid | 957 | 30 | 1 | 29 |
| CD8+ Trm | 818 | 29 | 1 | 28 |
| Treg | 863 | 24 | 1 | 23 |
| B_cell | 711 | 15 | 0 | 15 |
| Plasma_cell | 814 | 3 | 0 | 3 |

**Key Findings**:
1. **cDC1 dominates** with 594 significant genes — predominantly **downregulated** (378↓ vs 216↑) in Stroke
2. **Global downregulation pattern**: All cell types show more down than up → Stroke = transcriptional suppression
3. **CD4_S100A8 subtypes**: Strongly affected (142 + 84 sig genes) — 이 세포들은 myeloid-like CD4 T cells로 SIIS의 주요 marker
4. **Major cell types (Mono, NK, CD4 Naive, CD8 Cytotoxic) show 0 sig genes**: MNN의 DE filter가 이들에게 너무 strict

**Ligand Activities**: 1.2M rows in `ligand_activities_target_de_tbl`
- Top ligands: cDC1, cDC2로 수렴하는 activity가 가장 높음
- SIIS의 molecular mechanism: DC-targeting ligands의 activity가 지배적

#### L2: g3 Good vs Bad (IS Only)

**DE Genes Summary (p_adj < 0.05)**:
| Cell Type | Total Genes | Sig | Up | Down |
|-----------|------------|-----|-----|------|
| **ISG+ T_cell** | 520 | **20** | 8 | 12 |
| ISG+ Myeloid | 763 | 13 | 4 | 9 |
| B_cell | 491 | 6 | 0 | 6 |
| CD4_S100A8_CSF3R | 528 | 4 | 0 | 4 |
| MAIT | 516 | 3 | 2 | 1 |
| Plasma_cell | 530 | 3 | 0 | 3 |
| CD8+ Trm | 419 | 2 | 0 | 2 |

**Key Findings**:
1. **ISG+ cells 최다 DE**: ISG+ T_cell (20) + ISG+ Myeloid (13) → interferon response가 g3의 핵심 차이
2. **전반적으로 sparse**: L2 sample 수 적고 (32 patients) g3 변수 noise → 약한 signal
3. **Ligand activities**: 525K rows, ISG-related ligands가 상위권

### MNN 생성 플롯 목록

```
mnn_L1_cohort_anno1/:
  01_de_summary_bar.png          — Cell type별 sig DE gene 수
  02_volcano_grid.png            — Top 9 cell type volcanos
  03_ligand_activity_grid.png    — Top 6 receiver의 ligand activity bars
  04_top_ligand_receiver.png     — Top 40 ligand→receiver pairs
  05_circos_ligand_receiver.{pdf,png} — Chord diagram (ligand→receiver)
  06_ligand_receiver_heatmap.png — Heatmap (receiver × ligand activity)
  + de_summary_by_celltype.csv, la_de_top500.csv

mnn_L2_g3_anno1/:
  (same structure)
```

---

## 4. 기술적 문제점 (Technical Issues)

### 4.1 MNN group_prioritization_tbl Empty (0 Rows)

**원인**: MultiNicheNet의 `generate_prioritization_tables()`가 DE genes → LR network → sender-receiver pair로 이어지는 cascading filter를 적용. anno1에서 19개 cell type에 대해 sender×receiver 조합이 폭발적으로 많아지면서, 각 조합 내 DE→LR overlap이 임계값 미달.

**영향**: 표준 MNN 시각화 (prioritized LR pairs, top interactions) 생성 불가. `plot_multinichenet_circos()` 함수도 작동하지 않음.

**해결 방법 (적용됨)**: `ligand_activities_target_de_tbl` (1.2M rows)를 직접 활용하여 커스텀 시각화 생성.

**향후 계획**: 더 관대한 threshold로 MNN 재실행 고려:
- `fraction_cutoff`: 0.05 → 0.01
- `logFC_threshold`: 0.10 → 0.05
- `p_val_threshold`: 0.05 → 0.10

### 4.2 CellChat mergeCellChat() Failure

**문제**: Method 2A에서 condition-level merged 객체의 `@data` slot이 빈 상태 (ncol=0). `mergeCellChat()`이 이를 처리하지 못함.

**원인**: Sample-wise CellChat 실행 후 condition별 aggregate는 `@data.signaling`과 `@net`만 보존하고 원본 expression matrix (`@data`)는 미보존.

**해결**: 개별 sample CellChat 객체를 직접 로드하여 condition별 대표 sample (cell 수 최대)을 선택, 네이티브 함수에 전달.

### 4.3 rankNet / diff_scatter Failure

**문제**: "number of items to replace is not a multiple of replacement length"

**원인**: 대표 sample 간 cell type 구성이 다름 (HC sample에 특정 cell type 존재하지만 Stroke sample에 없거나 그 반대). CellChat의 `rankNet()`, `netAnalysis_signalingChanges_scatter()`는 두 객체의 cell type이 동일해야 함.

**상태**: 미해결. 해결하려면 공통 cell type만 가진 sample pair를 선별하거나, missing cell type에 대해 0-padding 필요.

### 4.4 MNN Wrapper Duplicate Parameter

**파일**: `myR/R/cci_multinichenet_wrapper.R`
**문제**: `p_val_thresh = 0.05` 파라미터가 line 52-53에서 중복 선언.
**영향**: R에서 duplicate formal argument error. 직접 source 불가.
**상태**: 미수정 (wrapper 대신 직접 코드 작성으로 우회).

### 4.5 Patient Name Mapping

**문제**: CellChat sample 디렉토리명이 R의 `make.names()` 처리된 환자명 사용.
- HC: "h.2460083" 등 → 그대로
- Stroke: "10001480" 등 → "X10001480" (숫자 시작이므로 X prefix)

**해결**: `make.names(patient_name)`으로 mapping 생성하여 56/56 (L1), 32/32 (L2) 완전 매칭.

---

## 5. 통합 해석 (Integrated Interpretation)

### CellChat + MNN Cross-validation

| Finding | CellChat Evidence | MNN Evidence |
|---------|-------------------|--------------|
| **SIIS (L1)** | -37% total CCI, DC hub collapse | cDC1: 594 DE genes (378↓), global downregulation |
| **DC centrality loss** | DC incoming/outgoing 급감 | cDC1 + cDC2 = 649 sig DE genes |
| **Myeloid-like CD4** | — (anno2 level) | CD4_S100A8: 226 sig genes (mostly down) |
| **g3 Bad = CCI increase** | +10% total CCI | ISG+ cells: 33 sig DE genes (top) |
| **ISG role in outcome** | — | ISG+ T_cell (20), ISG+ Myeloid (13) = g3의 핵심 |
| **NKc suppression** | NK outgoing 감소 | NK_cell: 0 sig MNN genes (too filtered) |

### Biological Model

```
HC (Baseline)
    ↓ Stroke (SIIS)
    └── Total CCI -37%
    └── DC hub collapse (594 genes down)
    └── Innate→Adaptive communication failure
    └── MIF, GALECTIN pathways disrupted
        ↓ Bad Outcome (g3=2)
        └── Paradoxical +10% CCI (maladaptive over-activation)
        └── Mono→DC/Bc over-signaling
        └── ISG+ cells: interferon response aberration
        └── RESISTIN pathway activation (metabolic inflammation)
        └── NKc functional suppression
```

---

## 6. 향후 계획 (Future Plans)

### 6.1 MNN Permissive Re-run (Priority: Medium)
- 더 관대한 DE filter로 MNN 재실행
- 목표: `group_prioritization_tbl`에 실제 prioritized pairs 확보
- anno1 level 유지, threshold 완화

### 6.2 CellChat rankNet Fix (Priority: Low)
- 공통 cell type subset으로 대표 sample pair 재선택
- 또는 cell type padding 구현

### 6.3 Sender-Receiver Network Visualization (Priority: Medium)
- MNN `ligand_activities_target_de_tbl`에서 sender 정보 추출
- Sender→Ligand→Receiver 3-layer network 시각화
- Condition-specific network comparison

### 6.4 CCI × DEG Cross-reference (Priority: High)
- MNN DE genes와 DEG Consensus (L1) / NEBULA (L2) 결과 교차 비교
- 공통 DE genes → CCI에서 ligand/receptor로 작용하는 유전자 식별
- FGS CMGI genes 중 CCI 관련 유전자 (MT-CO1, HLA-DQA2 등) 심층 분석

---

## 7. 결과물 위치 (Output Locations)

```
/data/user3/sobj/stroke_hc_v8_2/cci/
├── cellchat/
│   ├── L1_cohort_anno2/          # Per-sample CellChat objects
│   └── L2_g3_anno2/              # Per-sample CellChat objects
├── mnn/
│   ├── L1_cohort_anno1/          # multinichenet_results.qs
│   ├── L1_cohort_anno2/
│   ├── L2_g3_anno1/
│   └── L2_g3_anno2/
└── plots/
    ├── cellchat_L1_cohort/       # 16 comparison PNG + 33 pathway files
    │   └── pathways/
    ├── cellchat_L2_g3/           # 16 comparison PNG + 37 pathway files
    │   └── pathways/
    ├── mnn_L1_cohort_anno1/      # 6 PNG + 1 PDF + CSVs
    ├── mnn_L1_cohort_anno2/      # 11 files
    ├── mnn_L2_g3_anno1/          # 6 PNG + 1 PDF + CSVs
    └── mnn_L2_g3_anno2/          # (minimal)
```

### Scripts
```
/data/user3/sobj/stroke_hc_v8_2/scripts/
├── plot_cci_comprehensive.R      # CellChat native plots (bubble, river, pathway circles)
└── plot_mnn_anno1.R              # MNN anno1 DE volcanos, circos, heatmaps
```
