# Documentation Subtree Management

This directory (`docs-main/`) is managed as a **Git Subtree**, allowing it to be edited from any worktree/branch while maintaining independent version control and minimizing merge conflicts.

## Structure

```
docs-main/
├── guide.md                    # Main project guide (English)
├── guide_KR.md                # Main project guide (Korean)
├── README.md                   # Root README
├── README_SUBTREE.md          # This file
├── CONTEXT_DOCUMENTATION.md   # Agent context documentation
├── integrated/
│   ├── INTEGRATED_GUIDE.md     # Comprehensive module guide (English)
│   └── INTEGRATED_GUIDE_KR.md  # Comprehensive module guide (Korean)
├── rules/
│   └── DOCS_ORGANIZATION_RULE.md  # Documentation organization rules
├── dev/
│   ├── DEVLOG.md              # Development log (English)
│   ├── DEVLOG_Korean.md       # Development log (Korean)
│   ├── CHANGELOG.md           # Changelog (English)
│   ├── CHANGELOG_Korean.md    # Changelog (Korean)
│   ├── context.md             # Project context (English)
│   ├── context_Korean.md      # Project context (Korean)
│   └── README_Korean.md       # Package README (Korean)
└── config/
    └── vars_config.R          # Variable configuration script
```

## Purpose

This directory is managed as a **Git Subtree** to:
- **Prevent merge conflicts**: Independent repository reduces conflicts when multiple worktrees edit simultaneously
- **Enable worktree editing**: Can be edited from any worktree/branch
- **Independent versioning**: Has its own git history in a separate repository
- **Easy synchronization**: Push/pull changes independently using `git subtree` commands

## How It Works

`docs-main/` is a Git Subtree pointing to a separate repository (`mylit-docs`):
- **Remote repository**: `git@github.com:JasonK94/mylit-docs.git` (or your configured remote)
- **Remote name**: `docs-main-origin`
- **Integration**: Embedded in main repository but can be synced independently

## Usage

### Editing from Any Worktree

```bash
# From any worktree
cd /path/to/worktree

# Edit files
vim docs-main/guide.md

# Commit changes (to current branch)
git add docs-main/
git commit -m "Update guide.md"

# Push to docs repository
./scripts/docs-sync.sh push main
```

### Getting Latest Changes

```bash
# Pull latest docs from remote
./scripts/docs-sync.sh pull main
```

### Main Branch Only (Safer Approach)

```bash
# Switch to main branch
cd /home/user3/data_user3/git_repo/mylit
git checkout main

# Edit and commit
vim docs-main/guide.md
git add docs-main/
git commit -m "Update guide.md"
git push origin main

# Push to docs repository
./scripts/docs-sync.sh push main
```

## Setup

If you haven't set up the subtree yet, see:
- `SETUP_INSTRUCTIONS.md` - Initial setup guide
- `WORKTREE_SETUP.md` - Worktree configuration guide
- `SUBTREE_STRATEGY_PLAN.md` - Detailed strategy documentation

## Conflict Resolution

If multiple worktrees edit simultaneously:

1. Pull latest version first:
   ```bash
   ./scripts/docs-sync.sh pull main
   ```

2. Resolve conflicts manually if needed

3. Push resolved changes:
   ```bash
   git add docs-main/
   git commit -m "Resolve conflicts"
   ./scripts/docs-sync.sh push main
   ```

## Notes

- Always pull before pushing to avoid conflicts
- Small, frequent commits are better than large changes
- The subtree preserves full git history
- Each worktree needs the `docs-main-origin` remote configured (see `WORKTREE_SETUP.md`)

