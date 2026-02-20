# Claude Context Backup — stroke_milo_trajectory

> Created: 2026-02-18 14:17 KST
> Agent: Claude Opus 4.6 (claude-opus-4-6)
> Project: mylit (scRNAseq analysis)

## Purpose

Disaster recovery backup of all Claude agent context.
If `~/.claude/projects/-data-user3-git-repo-mylit/` is lost or corrupted,
this backup contains everything needed to resume work.

## Directory Structure

```
stroke_milo_trajectory/
├── README.md                  ← This file
├── CLAUDE.md                  ← Project instructions (319 lines)
├── conversations/             ← All JSONL conversation transcripts (68M total)
│   ├── c210ff2f-*.jsonl       ← Main session (35M, HC + stroke + MILO)
│   ├── b27e086f-*.jsonl       ← Stroke core session (26M)
│   ├── 83aeb54f-*.jsonl       ← Trajectory session (7.4M)
│   └── (12 smaller sessions)
├── memory/
│   └── MEMORY.md              ← Persistent agent memory (auto-loaded)
├── docs/                      ← All docs/claude/ context documents
│   ├── 1_STROKE_HC_V8_2_CONTEXT.md
│   ├── 4_PROJECT_GOALS.md     ← Master analysis plan (11 items)
│   ├── 7_PROGRESS_STATUS.md   ← Progress tracker
│   ├── HC/                    ← HC-only analysis context
│   ├── cci/                   ← CCI analysis reports
│   └── trajectory/            ← Trajectory analysis docs
└── session_meta/
    └── SESSION_INDEX.md       ← Session ID → purpose mapping
```

## Current Analysis State (as of backup)

### stroke_hc_v8_2 (primary dataset)
- 205,277 cells / 100 patients (HC + IS)

| Analysis | Status |
|----------|--------|
| Pre-processing + Annotation | DONE |
| Descriptive Figures (63 files) | DONE |
| MASC (6 comparisons) | DONE |
| scCODA (L1: 14/21, L2: CD14+ Mono) | DONE |
| MILO (L1: 11K nhoods, L2: 4.6K) | DONE |
| DEG Consensus (L1: 3 methods, L2: NEBULA) | DONE |
| CellChat v2 (L1+L2) | DONE |
| Trajectory v2 (3 compartments) | DONE |
| Gene Dynamics v3 + effect sizes + Lamian | DONE |
| FGS n=50/100 | DONE |
| FGS n=200 | RUNNING |
| FGS within-celltype (CD14 Mono) | RUNNING |
| cNMF v2 | RUNNING |
| MELD / Augur / External validation | PENDING |

### hc_only_v1 (HC-only atlas)
- 549,340 cells / 96 patients
- Multi-method frequency analysis COMPLETE (5-method concordance)
- FGS, CCI, DEG Consensus all DONE

## Restoring from Backup

```bash
# 1. Restore conversation history
cp context/claude/stroke_milo_trajectory/conversations/*.jsonl \
   ~/.claude/projects/-data-user3-git-repo-mylit/

# 2. Restore memory
mkdir -p ~/.claude/projects/-data-user3-git-repo-mylit/memory/
cp context/claude/stroke_milo_trajectory/memory/MEMORY.md \
   ~/.claude/projects/-data-user3-git-repo-mylit/memory/

# 3. docs/claude/ and CLAUDE.md are already in git repo
```
