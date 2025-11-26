# Documentation Management Context

This file describes the strategy for managing core documentation in the `myR` project.
**Agents should read this to understand how to handle `guide.md`, `README.md`, and other core docs.**

## 1. Core Philosophy: Git Subtree Management

*   **Independent Repository**: `docs-main/` is managed as a Git Subtree pointing to a separate repository (`mylit-docs`)
*   **Multi-Worktree Editing**: Documentation can be edited from any worktree/branch
*   **Conflict Minimization**: Independent versioning reduces merge conflicts when multiple worktrees edit simultaneously
*   **Synchronization**: Changes are synced using `git subtree push/pull` commands

## 2. Directory Structure (`docs-main/`)

The `docs-main/` folder is a Git Subtree managed independently:

```
docs-main/
├── guide.md                    # Main project guide (English)
├── guide_KR.md                # Main project guide (Korean)
├── README.md                   # Root README
├── integrated/                 # Integrated Guides
├── rules/                      # Documentation Rules
├── dev/                        # DEVLOG, CHANGELOG, context
└── config/                     # Shared scripts like vars_config.R
```

## 3. Handling Documentation Edits

### Scenario A: Editing in `main` branch (Recommended)
1.  Edit files directly in `docs-main/`.
2.  Commit changes: `git add docs-main/ && git commit -m "..."`.
3.  Push to main repo: `git push origin main`.
4.  Push to docs repo: `./scripts/docs-sync.sh push main`.
5.  **Result**: Official documentation is updated and synced.

### Scenario B: Editing in Feature Branch/Worktree
1.  Edit files directly in `docs-main/`.
2.  Commit changes to current branch.
3.  Push to docs repo: `./scripts/docs-sync.sh push main`.
4.  **Action**: Other worktrees can pull latest: `./scripts/docs-sync.sh pull main`.

### Scenario C: Syncing Feature Branch with Latest Docs
1.  Run `./scripts/docs-sync.sh pull main` in the feature branch/worktree.
2.  **Result**: Latest docs from remote are merged into current branch.

## 4. Technical Mechanisms

*   **Git Subtree**: `docs-main/` is a subtree pointing to `docs-main-origin` remote
*   **Sync Scripts**: `scripts/docs-sync.sh` provides push/pull commands
*   **Remote Repository**: Separate GitHub repository for independent versioning
*   **Worktree Support**: Each worktree can have its own `docs-main-origin` remote

## 5. Agent Instructions

*   **When reading**: Access files directly in `docs-main/` directory.
*   **When writing**:
    *   Edit files in `docs-main/` directory (works from any branch/worktree).
    *   After committing, use `./scripts/docs-sync.sh push main` to sync to docs repository.
    *   Acknowledge that this change will be available to other worktrees after they pull.
    *   If editing from a feature branch, consider also updating main branch for consistency.

## 6. Quick Reference

```bash
# Push local docs changes to remote
./scripts/docs-sync.sh push main

# Pull latest docs from remote
./scripts/docs-sync.sh pull main

# Check docs remote
git remote -v | grep docs-main
```

## 7. Related Documentation

- `README_SUBTREE.md` - Detailed subtree documentation (this directory)
- `SETUP_INSTRUCTIONS.md` - Initial setup guide
- `WORKTREE_SETUP.md` - Worktree configuration
- `SUBTREE_STRATEGY_PLAN.md` - Strategy planning document

---
*Last Updated: 2025-11-25*
*Strategy: Git Subtree with independent repository*

