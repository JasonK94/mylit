# Agent: CCI Re-implementation & Multi-analysis

> **Agent session ID**: `b27e086f-9707-4dc4-8d70-b3ed35221e70`
> **Last synced**: 2026-02-18 14:18 KST
> **Role**: CCI (CellChat/MNN) 감사 + 재구현, trajectory effect sizes, MILO, 추가분석 기획

---

## What This Agent Did

### 1. CCI CellChat v2 재구현 (핵심)
- **문제**: CellChat v2 `mergeCellChat()`이 `@data` slot을 drop → rankNet, bubble 등 비교 함수 불가
- **해결**: Seurat에서 condition별 직접 `createCellChat()` → @data 보존
- **결과**: L1 59 PNG + L2 59 PNG (rankNet, bubble, heatmap, scatter, pathway chord 모두 정상)
- **핵심 스크립트**: `scripts/stroke/reimplement_cellchat_v2.R`
- **Output**: `/data/user3/sobj/stroke_hc_v8_2/cci/plots/cellchat_{L1_cohort,L2_g3}_v2/`

### 2. MNN (MultiNicheNet) 감사
- **결론**: 근본적으로 비호환 — `lr_target_prior_cor`에 celltype당 sample ≥5 필요, 불충족
- **Permissive thresholds 재실행**: 여전히 0 rows → 확정 abandoned
- **대안**: LIANA (multi-method consensus CCI) 검토 필요

### 3. Trajectory Effect Sizes (4 phases)
- Phase 1-2: ABC/rABC/RMISE + comparison figures
- Phase 3: Lamian validation (chisq method)
- Phase 4: Pseudobulk GAMM sensitivity
- **핵심 발견**: 10 genes confirmed by all 3 methods (mono/cohort); mono is ONLY compartment with g3 effect

### 4. MILO Neighbourhood DA
- L1: 11,739 nhoods (63% sig), L2: 4,664 nhoods (19% sig)
- 3-method concordance (MILO x MASC x scCODA): 14/19 DA by ≥2 methods

### 5. Documentation & Cross-Agent Context
- CLAUDE.md에 cross-agent hints 추가 (cNMF, MELD, Augur, MOFA, propeller 구현 가이드)
- CCI 전용 문서 정리: `docs/claude/cci/`
- Progress tracker 지속 업데이트

---

## Directory Structure

```
agent_cci_reimplement/
├── README.md                   ← 이 파일
├── CLAUDE.md                   ← 프로젝트 루트 CLAUDE.md 스냅샷
├── memory/
│   └── MEMORY.md               ← Claude auto-memory
├── conversations/
│   └── b27e086f-*.jsonl        ← 이 agent의 대화 로그 (26MB)
├── docs/
│   ├── cci/                    ← CCI 분석 전용
│   │   ├── CCI_ANALYSIS_REPORT.md
│   │   └── CCI_AUDIT_AND_REIMPLEMENTATION.md
│   ├── stroke/                 ← stroke_hc_v8_2 핵심 문서
│   │   ├── 1_STROKE_HC_V8_2_CONTEXT.md
│   │   ├── 4_PROJECT_GOALS.md
│   │   └── 7_PROGRESS_STATUS.md
│   ├── trajectory/
│   │   ├── 1_TRAJECTORY_ANALYSIS_STRATEGY.md
│   │   ├── 2_DIFFERENTIAL_DYNAMICS_METHODS.md
│   │   ├── 3_CURRENT_STATE_AND_PLAN.md
│   │   └── 4_TRAJECTORY_RESULTS_INTEGRATED.md
│   └── HC/
│       ├── HC_ANALYSIS_CONTEXT.md
│       ├── HC_NOVEL_ANALYSIS_PLAN.md
│       ├── HC_RESULTS_INTERIM.md
│       ├── FGS_AGE_RESULTS.md
│       └── TREG_EXHAUSTION_ANALYSIS.md
└── scripts/
    ├── stroke/                 ← stroke analysis scripts (16 files)
    │   ├── reimplement_cellchat_v2.R      ← ★ 핵심
    │   ├── fix_cellchat_heatmaps.R
    │   ├── fix_cellchat_v2_final.R
    │   ├── reimplement_cellchat_comparison.R
    │   ├── reimplement_mnn_permissive.sh
    │   ├── compute_trajectory_effect_sizes.R
    │   ├── run_lamian_validation.R
    │   ├── run_pseudobulk_gamm.R
    │   ├── run_milo.R
    │   ├── run_gene_dynamics_v3.R
    │   ├── run_descriptive_figures.R
    │   ├── run_descriptive_part2.R
    │   ├── run_pathway_crosslayer.R
    │   ├── run_crosslayer_concordance.R
    │   ├── plot_cci_comprehensive.R
    │   ├── plot_mnn_anno1.R
    │   └── plot_fgs_results.R
    └── hc/                     ← HC analysis scripts (12 files)
        ├── run_fgs_continuous_v2.R
        ├── run_fgs_by_compartment.R
        ├── run_masc_hc_downsampled.R
        ├── run_masc_hc.R
        ├── run_milo_age.R
        ├── run_meld.py
        ├── run_sccoda.py
        ├── run_treg_exhaustion_deep.R
        ├── validate_aging_signatures.R
        ├── export_for_python.R
        ├── interpret_cellchat_age.R
        └── interpret_mnn.R
```

## Known Issues Discovered

```
CellChat v2 mergeCellChat → @data = 0 cells (by design)
CellChat centrality lost during merge → compute on individual, then lift→merge
ComplexHeatmap + png() → must explicitly draw() HeatmapList
MNN lr_target_prior_cor needs ≥5 samples/celltype/condition → always NULL
Lamian sample = patient_name (NOT GEM)
MASC 549K → all NA (must downsample to ~120K)
```

## Recovery Instructions

새 세션에서 이 agent의 context 복원:
1. `CLAUDE.md` 읽기 → 전체 프로젝트 구조
2. `memory/MEMORY.md` 읽기 → 누적 지식
3. `docs/stroke/7_PROGRESS_STATUS.md` 읽기 → 현재 진행 상황
4. `docs/cci/CCI_AUDIT_AND_REIMPLEMENTATION.md` 읽기 → CCI 상세 상태

## Sync Command

```bash
# 이 agent의 context 재동기화 (로컬 → 원격)
BASE="context/claude/agent_cci_reimplement"
cp ~/.claude/projects/-data-user3-git-repo-mylit/memory/MEMORY.md $BASE/memory/
cp CLAUDE.md $BASE/
cp ~/.claude/projects/-data-user3-git-repo-mylit/b27e086f-*.jsonl $BASE/conversations/
cp docs/claude/cci/*.md $BASE/docs/cci/
cp docs/claude/trajectory/*.md $BASE/docs/trajectory/
cp docs/claude/HC/*.md $BASE/docs/HC/
cp docs/claude/{1,4,7}_*.md $BASE/docs/stroke/
```
