# Claude Session Index — stroke_hc_v8_2 + hc_only_v1

> Backup date: 2026-02-18
> Agent folder: `stroke_milo_trajectory`
> Purpose: Disaster recovery — if local `.claude/` is lost, this backup preserves
>          full conversation history, memory, and project context documents.

## Active Sessions (by size, most recent first)

### 1. `c210ff2f` (35M, 25,277 lines) — **Main analysis session**
- **Period**: Feb 13–18, 2026 (multi-continuation)
- **Scope**: Full project from initial codebase survey to MILO completion
- **Key work**:
  - Codebase survey (myR package, docs, scripts)
  - HC-only dataset: MrVI, MASC, MILO, DEG Consensus, FGS (6 methods), CCI (CellChat+MNN), Treg exhaustion, MELD, scCODA, Literature validation
  - HC multi-method frequency: 5-method concordance (CD8 Naive, CD8 Eff, NK CD16hi, MAIT, CD8 TEM)
  - stroke_hc_v8_2: Descriptive figures (63), FGS (n=50/100), MILO (L1+L2)
  - Context backup (this file)

### 2. `b27e086f` (26M, 8,488 lines) — **Stroke analysis core**
- **Period**: Feb 12–18, 2026 (multi-continuation)
- **Scope**: Core stroke_hc_v8_2 analysis pipeline
- **Key work**:
  - Annotation QC, MASC (6 runs), DEG Consensus (L1: 3 methods, L2: NEBULA)
  - CellChat v2 re-implementation (rankNet, bubble, scatter)
  - MNN investigation + abandonment (0 prioritized interactions)
  - Trajectory v2 (Slingshot + Monocle3, scVI UMAP)
  - Gene dynamics v3 (batch-corrected GAMM)
  - scCODA (L1: 14/21 credible, L2: CD14+ Mono only)
  - FGS pipeline, pathway enrichment (HALLMARK/KEGG/GOBP)
  - Cross-layer concordance (134 genes)

### 3. `83aeb54f` (7.4M, 2,695 lines) — **Trajectory effect sizes**
- **Period**: Feb 17–18, 2026 (multi-continuation)
- **Scope**: Systematic trajectory analysis with validation
- **Key work**:
  - Phase 1-2: ABC/rABC/RMISE effect sizes
  - Phase 3: Lamian validation (chisq method)
  - Phase 4: Pseudobulk GAMM sensitivity analysis
  - Phase 5: Biological interpretation document
  - Key finding: 10 genes confirmed by all 3 methods (mono/cohort)
  - Mono = ONLY compartment with g3 effect

### 4. Small sessions (< 16K each)
| Session ID | Size | Date | Content |
|------------|------|------|---------|
| `98b05255` | 13K | Feb 18 | Short follow-up |
| `0f83b884` | 8.8K | Feb 18 | Short follow-up |
| `c7a61481` | 8.2K | Feb 18 | Short follow-up |
| `4df200be` | 5.1K | Feb 14 | HC dataset work |
| `8e60e572` | 5.1K | Feb 13 | HC dataset work |
| `c01ffd58` | 4.9K | Feb 14 | HC dataset work |
| `c205f327` | 3.7K | Feb 18 | Brief session |
| `5e191901` | 2.0K | Feb 13 | Brief session |
| `e4572bb6` | 1.7K | Feb 18 | Brief session |
| `da356c66` | 1.5K | Feb 13 | Brief session |
| `dafa6922` | 1.1K | Feb 13 | Brief session |
| `e65e0106` | 957B | Feb 13 | Brief session |

## Session Continuation Chain

```
c210ff2f (initial) → continuation → continuation → ... → current
83aeb54f (trajectory focus) → continuation → continuation → current
b27e086f (stroke core) → continuation → ... → current
```

## Recovery Instructions

If local `.claude/` is lost:
1. Copy `conversations/*.jsonl` back to `~/.claude/projects/-data-user3-git-repo-mylit/`
2. Copy `memory/MEMORY.md` to `~/.claude/projects/-data-user3-git-repo-mylit/memory/`
3. The `docs/` folder is already in the git repo at `docs/claude/`
4. `CLAUDE.md` is already at project root
5. Start new session — MEMORY.md will auto-load and provide full context
