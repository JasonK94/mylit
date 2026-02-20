# Codebase Refactoring Plan — mylit Repository

## 1. Current State Assessment

### Architecture
The repo uses **git worktrees** — 13 branches, each checked out in `/data/user3/git_repo/_wt/{module}/`. Each worktree has its own **complete copy of myR/** plus module-specific scripts.

```
main (mylit/)           ← myR + docs + config
├── _wt/analysis/       ← branch: analysis
├── _wt/cci/            ← branch: mnn
├── _wt/cellchat/       ← branch: cc
├── _wt/deg-consensus/  ← branch: deg-consensus-dev
├── _wt/fgs/            ← branch: fgs-refactoring
├── _wt/lds/            ← branch: lds
├── _wt/masc/           ← branch: masc
├── _wt/milo/           ← branch: milo
├── _wt/pipe/           ← branch: pipe
├── _wt/plots/          ← branch: plots-dev
├── _wt/pseudotime/     ← branch: pseudotime-dev
├── _wt/pt.umap/        ← branch: pt.umap
└── _wt/ds/             ← branch: ds
```

### Problem Summary

| Issue | Scale | Impact |
|-------|-------|--------|
| **myR duplicated across 13 worktrees** | ~60 R files × 13 copies | Fixes in main don't propagate; bugs recur |
| **All copies have diverged** | Every file's md5 differs from main | No single source of truth |
| **Module-specific code mixed into myR** | e.g., CCI functions in cellchat/myR/, FGS in fgs/myR/ | Hard to find authoritative version |
| **Test/dev cruft** | test_to_delete.R (6.5K LOC), test_claude.R (2.1K) | Clutters package |
| **~900 test/run scripts** across worktrees | Many one-off experiments | Hard to find actual entry points |
| **88+ documentation files** | Some duplicated between main and worktrees | Inconsistent guidance |

### Root Cause
The worktree pattern was originally useful for isolated module development. But since each worktree carries a full myR copy, and development happened independently on each branch, the codebase fractured into 13 semi-independent copies.

---

## 2. Proposed Architecture

### 2.1 Target State

```
mylit/
├── myR/                          # SINGLE source of truth for R package
│   ├── R/
│   │   ├── core/                 # Core utilities (seurat, plots, validation)
│   │   ├── analysis/             # Analysis functions (DEG, trajectory, etc.)
│   │   ├── cci/                  # Cell-cell interaction
│   │   ├── fgs/                  # Feature gene signature
│   │   ├── deg_consensus/        # DEG consensus pipeline
│   │   └── legacy/               # Deprecated (tagged for removal)
│   ├── DESCRIPTION
│   ├── NAMESPACE
│   └── man/
├── scripts/                      # ALL runnable scripts, organized by module
│   ├── cellchat/
│   │   └── run_cellchat_cli.R
│   ├── mnn/
│   │   └── run_multinichenet.R
│   ├── fgs/
│   │   ├── run_fgs_pipeline.R
│   │   └── run_fgs_continuous.R  # NEW: continuous variable fork
│   ├── masc/
│   │   └── run_masc.R
│   ├── deg/
│   │   └── run_deg_consensus_cli.R
│   ├── milo/
│   │   └── run_milo.R
│   ├── pipe/
│   │   └── run_pipeline.sh
│   └── pseudotime/
│       └── run_trajectory.R
├── mypy_re/                      # Python pipeline (scVI, Solo, etc.)
├── config/                       # Manifests + configs
├── docs/                         # Documentation (single source)
└── tests/                        # Formal test suite
```

### 2.2 Key Principles

1. **One myR**: Eliminate all worktree copies. Scripts reference a single installed or `source()`-d myR.
2. **Scripts as thin CLIs**: Each `scripts/{module}/run_*.R` is a thin CLI wrapper that sources myR functions and parses arguments.
3. **No worktrees for code isolation**: Use feature branches with regular merges instead.
4. **Module-specific code stays in myR**: CCI code → `myR/R/cci/`, FGS → `myR/R/fgs/`, etc.
5. **Git Repo symlink unnecessary**: Once scripts are in main, worktrees can be retired.

---

## 3. Migration Strategy

### Phase 1: Consolidate myR (HIGH priority)
**Goal**: Merge all unique functions from worktree myR copies into main myR.

#### 3.1.1 Audit each worktree's myR
For each worktree, diff against main:
```bash
diff -rq /data/user3/git_repo/mylit/myR/R/ /data/user3/git_repo/_wt/{module}/myR/R/
```
- Identify files that exist only in worktrees (new functions)
- Identify files where worktree version is newer/better
- Keep the most recent, most complete version of each function

#### 3.1.2 Merge priority (by module activity)
1. **milo** — recently fixed (`calcNhoodDistance` bugfix in main)
2. **deg-consensus** — actively used, critical analysis
3. **fgs** — active development (continuous fork coming)
4. **cellchat** — active CCI analysis
5. **cci/mnn** — active CCI analysis
6. **pseudotime** — active trajectory analysis
7. **masc** — completed, stable
8. Others: lds, pipe, plots, pt.umap, analysis, ds — lower priority

#### 3.1.3 Strategy per module
```
For each module M:
  1. git diff main..{branch} -- myR/R/ > /tmp/diff_M.patch
  2. Review changes: genuine improvements vs divergence artifacts
  3. Cherry-pick genuine improvements into main's myR
  4. Test: library(myR) on main with the merged code
  5. Update scripts/{M}/ to use main myR
```

### Phase 2: Consolidate Scripts (MEDIUM priority)
**Goal**: Move all CLI scripts from worktrees to `scripts/` in main.

#### Current vs Target
| Current | Target |
|---------|--------|
| `_wt/cellchat/scripts/cellchat/run_cellchat_cli.R` | `scripts/cellchat/run_cellchat_cli.R` |
| `_wt/cci/scripts/cci/mnn/run_multinichenet.R` | `scripts/mnn/run_multinichenet.R` |
| `_wt/fgs/scripts/fgs/run_fgs_pipeline.R` | `scripts/fgs/run_fgs_pipeline.R` |
| `_wt/masc/scripts/masc/run_masc.R` | `scripts/masc/run_masc.R` |
| `_wt/deg-consensus/scripts/consensus/run_deg_consensus_cli.R` | `scripts/deg/run_deg_consensus_cli.R` |

#### Steps
1. Copy authoritative CLI scripts to `scripts/{module}/`
2. Update `source()` paths to reference main myR
3. Test each CLI with the same arguments used in production
4. Update CLAUDE.md and docs with new script locations

### Phase 3: Clean Up (LOW priority)
1. **Remove test cruft**: Delete `test_to_delete.R`, `test_claude.R`, `test_cursor.R`
2. **Archive worktrees**: Once scripts work from main, mark worktrees as archived
3. **Consolidate docs**: Merge worktree docs into main `docs/`
4. **NAMESPACE cleanup**: Export pseudotime functions properly
5. **Remove project-specific code**: `IBD_delta_reduction.R`, `IBD_PCLMM_LDA.R` → separate project or `legacy/`

---

## 4. Risk Assessment

### Risks
| Risk | Mitigation |
|------|-----------|
| **Breaking active analysis** | Keep worktrees functional during migration; don't delete until verified |
| **Merge conflicts in myR** | Do module-by-module, test each merge |
| **CCI scripts depend on worktree myR** | Update `source()` paths; test with actual data |
| **FGS/TML code is most divergent** | Special attention; may need manual reconciliation |

### Non-Risks
- Config/manifests: Already centralized in main
- Python code (mypy_re): Separate directory, no worktree issues
- Data files: External to repo, not affected

---

## 5. Effort Estimate

| Phase | Scope | Effort |
|-------|-------|--------|
| Phase 1.1: Audit all worktree diffs | 13 worktrees × 60 files | 2-3 hours |
| Phase 1.2: Merge top-6 modules into main myR | 6 modules, ~50 files with conflicts | 4-6 hours |
| Phase 1.3: Merge remaining 7 modules | Lower priority, less divergence | 2-3 hours |
| Phase 2: Move CLI scripts | 8-10 scripts to relocate & test | 2-3 hours |
| Phase 3: Cleanup | Test files, docs, NAMESPACE | 1-2 hours |
| **Total** | | **~15-20 hours** |

This is a substantial but not overwhelming effort. Can be done incrementally over several sessions, starting with the most actively used modules.

---

## 6. Immediate Quick Wins (can do now)

1. **Delete `myR/R/test_to_delete.R`** (6.5K LOC of cruft)
2. **Delete `myR/R/test_claude.R`** and other test_*.R files
3. **Copy FGS continuous script** to `scripts/fgs/run_fgs_continuous.R`
4. **Update CLAUDE.md** with current script locations from worktrees
5. **Export pseudotime functions** in NAMESPACE (fixes the constant `source()` requirement)

---

## 7. Decision Points for User

1. **Do you want to fully retire worktrees?** Or keep them as read-only archives?
2. **Priority**: Which modules are most critical to consolidate first?
3. **Testing**: Do you have any formal test suite, or should we create one during refactoring?
4. **Git history**: Preserve branch history (merge) or start fresh on main?
5. **Package installation**: Should myR be installable via `devtools::install()` and used via `library(myR)`, or continue with `source()`?
