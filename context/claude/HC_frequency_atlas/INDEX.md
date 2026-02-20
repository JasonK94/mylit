# Claude Context Backup — HC Frequency Atlas Agent

> Agent ID: `claude_HC_frequency_atlas`
> Created: 2026-02-18 | Purpose: Remote backup of all Claude agent context
> Focus: HC-only PBMC atlas multi-method frequency analysis (MELD, MILO, MASC, scCODA, CellChat, MNN)
> Restore: Copy `memory/MEMORY.md` → `~/.claude/projects/-data-user3-git-repo-mylit/memory/MEMORY.md`
> Restore conversations: Copy `conversations/*.jsonl` → `~/.claude/projects/-data-user3-git-repo-mylit/`

---

## Directory Structure

```
context/claude_HC_frequency_atlas/
├── INDEX.md                          # This file — navigation guide
├── CLAUDE.md                         # Project-level instructions (copy of repo root)
├── memory/
│   └── MEMORY.md                     # Persistent agent memory (key patterns, fixes, state)
├── conversations/                    # Full conversation transcripts (JSONL)
│   ├── c210ff2f-...jsonl             # 35MB — HC multi-method frequency analysis (current session)
│   ├── b27e086f-...jsonl             # 26MB — Stroke analysis (trajectory, gene dynamics, MILO)
│   ├── 83aeb54f-...jsonl             # 7.4MB — CCI CellChat reimplementation + FGS
│   └── (12 smaller sessions)
├── analysis_docs/
│   ├── stroke/                       # Stroke-HC dataset context & results
│   │   ├── 1_STROKE_HC_V8_2_CONTEXT.md
│   │   ├── 2_gene_dynamics_approaches.md
│   │   ├── 3_REFACTORING_PLAN.md
│   │   ├── 4_PROJECT_GOALS.md        # 11-item analysis plan
│   │   ├── 5_ANALYSIS_REPORT.md
│   │   ├── 6_ADDITIONAL_ANALYSES_PLAN.md
│   │   └── 7_PROGRESS_STATUS.md      # Master progress tracker
│   ├── HC/                           # HC-only atlas context & results
│   │   ├── HC_ANALYSIS_CONTEXT.md    # Master context (549K cells, 96 patients)
│   │   ├── HC_NOVEL_ANALYSIS_PLAN.md # Multi-method frequency plan
│   │   ├── HC_RESULTS_INTERIM.md     # *** Comprehensive results (7 methods, ~430 lines) ***
│   │   ├── FGS_AGE_RESULTS.md        # FGS continuous results (6 methods)
│   │   └── TREG_EXHAUSTION_ANALYSIS.md
│   ├── cci/                          # Cell-cell interaction docs
│   │   ├── CCI_ANALYSIS_REPORT.md
│   │   └── CCI_AUDIT_AND_REIMPLEMENTATION.md
│   └── trajectory/                   # Trajectory analysis docs
│       ├── 1_TRAJECTORY_ANALYSIS_STRATEGY.md
│       ├── 2_DIFFERENTIAL_DYNAMICS_METHODS.md
│       ├── 3_CURRENT_STATE_AND_PLAN.md
│       └── 4_TRAJECTORY_RESULTS_INTEGRATED.md
└── scripts/
    ├── hc/                           # HC analysis scripts (12 files)
    │   ├── export_for_python.R, run_meld.py, run_sccoda.py
    │   ├── run_masc_hc.R, run_masc_hc_downsampled.R, run_milo_age.R
    │   ├── interpret_cellchat_age.R, interpret_mnn.R
    │   ├── validate_aging_signatures.R
    │   └── run_fgs_continuous_v2.R, run_fgs_by_compartment.R, run_treg_exhaustion_deep.R
    └── stroke/                       # Stroke analysis scripts (44 files)
        ├── run_cci_all.sh, reimplement_cellchat_v2.R
        ├── run_trajectory*.R, run_gene_dynamics_v3.R
        ├── run_deg_and_milo.sh, run_milo.R
        └── (see full listing in directory)
```

---

## Key Conversation Sessions

### 1. HC Multi-Method Frequency Analysis (c210ff2f, 35MB) — CURRENT
**Date**: 2026-02-18 | **Status**: All complete

Ran 7 independent frequency analysis methods on 549K HC PBMC cells (96 Korean donors, age 19-66):
- **MELD**: Per-cell age/sex density estimation
- **MILO**: Neighbourhood-level DA (continuous age)
- **MASC**: Mixed-effects logistic regression (downsampled to 115K)
- **scCODA**: Bayesian compositional DA (sex + age binary + 3-group)
- **CellChat**: Cell-cell interaction comparison (age + sex)
- **MNN**: MultiNicheNet CCI-aware DE interpretation
- **Literature validation**: 12 aging signatures cross-referenced

**Core finding**: 5-method concordance on 5 cell types (CD8 Naive↓, CD8 Effector↑, NK CD16hi↑, MAIT↓, CD8 TEM↑). Novel: inflammatory mono decline, MAIT as #2 declining, male-specific naive decline, MIF sex dimorphism.

### 2. Stroke Analysis (b27e086f, 26MB)
**Date**: 2026-02-18 | **Scope**: Trajectory, gene dynamics, MILO, MASC, DEG

Comprehensive stroke-HC analysis including:
- Trajectory v2-v3 (Slingshot + Monocle3, 3 compartments)
- Gene dynamics v3 (NB-GAMM, all compartments × conditions)
- Trajectory effect sizes (ABC/rABC/RMISE, Lamian validation, pseudobulk GAMM)
- MILO v2 (L1 + L2 DA testing)
- scCODA (L1 14/21 credible, L2 only CD14+ Mono)

### 3. CCI + FGS (83aeb54f, 7.4MB)
**Date**: 2026-02-18 | **Scope**: CellChat reimplementation, FGS pipeline

- CellChat v2 re-implementation from Seurat (fixed @data empty issue)
- FGS pipeline runs (n=50, n=100)
- Descriptive figure generation (63 files)

---

## Critical Files for Resuming Work

If starting a new session, the agent should read these files in order:
1. `CLAUDE.md` — Project instructions, conventions, environment setup
2. `memory/MEMORY.md` — Persistent patterns, known issues, completed work
3. `analysis_docs/stroke/7_PROGRESS_STATUS.md` — Master progress tracker
4. `analysis_docs/HC/HC_RESULTS_INTERIM.md` — HC comprehensive results
5. `analysis_docs/stroke/4_PROJECT_GOALS.md` — 11-item analysis plan

---

## Data Locations (not backed up — too large)

```
/data/user3/sobj/
├── stroke_hc_v8_2/         # Primary stroke-HC dataset
│   ├── 5_strokev8_clean.qs # 205K cells
│   ├── 5_1_hc_is.qs        # L1: 131K (HC+IS)
│   ├── 5_2_is_g3.qs        # L2: 54K (IS g3=1 vs g3=2)
│   └── (results in subdirs)
├── hc_only_v1/              # HC-only atlas
│   ├── 2_hc_annotated.qs   # 549K cells (was 5_hc_clean.qs)
│   ├── meld/               # MELD output
│   ├── milo_age/           # MILO age DA
│   ├── masc/               # MASC (downsampled)
│   ├── sccoda/             # scCODA
│   ├── cellchat_interpretation/  # CellChat age/sex
│   ├── aging_signatures/   # Literature validation
│   └── mnn_interpretation/ # MNN DE/ligand
```

---

## Pending Work (as of 2026-02-18)

### stroke_hc_v8_2
- FGS n=200 (running), within-cell-type FGS (14 pending)
- cNMF (K selection done for 2, NK factorizing)
- LIANA (CCI alternative), Subclustering, MELD, Augur, MOFA+, scANVI, External validation

### hc_only_v1
- DEG consensus (age/sex) — per-cell-type pseudobulk
- Trajectory analysis — monocyte → DC differentiation
- scANVI, MrVI interpretation, Subclustering
- External validation (Oelen et al. 2022)

---

## Environment Notes

- R 4.3.2 at `/opt/R/4.3.2/`
- Always set `.libPaths()` before `library()` calls (renv + user lib)
- Python: `conda activate scvi-tools` for MELD/scCODA/cNMF
- GPU: RTX 4090 (scVI)
- CPU affinity: `taskset -c 0-15` or `16-31` for long jobs
