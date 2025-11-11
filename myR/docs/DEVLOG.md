# Developer Log (DEVLOG)

This log captures the narrative context, decisions, and reasoning behind major development milestones of `myR`.

## 2025-05-09 — Project Kickoff (`1d9a350`)
- **Author**: kjc17  
- **Summary**: Seeded the repository with demultiplexing, integration, and Seurat workflows.
- **Details**:
  - Added `demuxalot`, `soupx`, `scVI` scripts plus seven-step preprocessing R pipelines.
  - Created initial `myR/` layout with DESCRIPTION/NAMESPACE scaffolding and large method exports.
- **Next Steps**:
  - Validate load/install flow for the newly imported modules.
  - Curate large artifacts to keep git history manageable.

## 2025-05-16 — Bulk Feature Ingestion (`338a7f1`)
- **Author**: kjc17  
- **Summary**: Imported core analysis modules and supporting notebooks.
- **Details**:
  - Checked in `CCI.R`, `pseudobulk_deg.R`, `pseudotime.R`, `legacy.R`, and supporting tests.
  - Versioned NicheNet outputs, plots, and Rmd reports used to verify early pipelines.
- **Next Steps**:
  - Organize modules into thematic subdirectories.
  - Remove oversized binary outputs from version control.

## 2025-05-22 — Recovery & Merge (`71c090f`, `b70456c`)
- **Author**: kjc17  
- **Summary**: Recovered overwritten work and reconciled divergent branches.
- **Details**:
  - Restored removed files, reintroduced `demulti_utils`, and merged working copies.
  - Cleaned git history by clearing cached large files and fixing `.gitignore`.
- **Next Steps**:
  - Stabilize utility scripts and start documenting public interfaces.

## 2025-10-10 — Formula Utility Improvements (`9f04cc9`, `c3afe5a`, `3f44856`)
- **Author**: kjc_server1  
- **Summary**: Harmonized formula-based helpers across analysis utilities.
- **Details**:
  - Updated FLSG/RLMG helpers to accept `formula_str` consistently.
  - Patched `fmsg`, refreshed `test.Rmd`, and laid groundwork for revised testing strategy.
- **Next Steps**:
  - Propagate the formula pattern into downstream test harnesses.
  - Document expectation changes for collaborators.

## 2025-10-31 — Pseudobulk & Testing Overhaul (`881589c`, `25b205b`, `e785399`, `39a6159`)
- **Author**: kjc_server1  
- **Summary**: Promoted stable routines and expanded the analytics test suite.
- **Details**:
  - Moved curated functions into `test_working.R`, refactored pseudobulk logic (`pb_deg_cursor`), and migrated Seurat usage from `slot` to `layer`.
  - Added TestLISI, PERMANOVA, PCP, PTMFL utilities; declared `NMF` dependency explicitly.
- **Next Steps**:
  - Wrap pseudobulk refactor into a dedicated module for easier maintenance.
  - Capture examples for the new QA utilities.

## 2025-11-05 — Infrastructure Tidying (`497eab8`, `80260af`, `80f7941`)
- **Author**: kjc_server1  
- **Summary**: Reorganized utility scripts and pruned unstable experiments.
- **Details**:
  - Relocated base R helpers, added MUSCAT/NEBULA runners, and trimmed MAST code.
  - Updated `.gitignore`; stopped tracking `projects/` and `trash/` outputs to streamline the repo.
- **Next Steps**:
  - Keep experimental notebooks out of git to avoid future churn.
  - Align documentation with the new directory structure.

## 2025-11-09 — Test Harness Consolidation (`8670f8c`)
- **Author**: kjc_server1  
- **Summary**: Folded day-to-day experimentation into the curated test suite.
- **Details**:
  - Synced deletions with ignore rules and ensured runnable entry points in `test_working.R`.
- **Next Steps**:
  - Begin drafting formal documentation for test workflows.

## 2025-11-10 — Signature v5.2 & Milo Opus6 (`2d9c50e`, `18f1ad1`, `b260bb1`, `a360837`, `7bf2904`)
- **Author**: kjc_server1  
- **Summary**: Rebuilt signature scoring and introduced Milo differential abundance analysis.
- **Details**:
  - Delivered signature workflow v5.2 with validation utilities and refreshed docs.
  - Added `milo_opus6`, improved meta-learner parallel defaults, captured the `opus6` integration snapshot, and expanded test coverage (`test_gpt.R`).
- **Next Steps**:
  - Publish usage examples for Milo and signature pipelines.
  - Formalize versioning to capture these feature milestones.

## 2025-11-10 — Documentation Formalization (`cf6d13e`, `ac3cb7b`, `18c2810`)
- **Author**: kjc_server1  
- **Summary**: Ran `cinit` and relocated documentation to repository root.
- **Details**:
  - Generated context, bootstrap instructions, and dual-language docs.
  - Moved `DEVLOG.md` and `CHANGELOG.md` into place to anchor ongoing documentation.
- **Next Steps**:
  - Keep context files synchronized with evolving package functionality.
  - Schedule regular updates to DEVLOG/CHANGELOG alongside major commits.

---

### Upcoming Focus
- Backfill post–signature v5.2 updates into `function_analysis*.md`.
- Decide on semantic version tags capturing Milo/signature milestones before the next release.
