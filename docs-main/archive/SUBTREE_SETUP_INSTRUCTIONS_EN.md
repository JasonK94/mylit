# Git Subtree Setup Instructions

## Step 1: Create GitHub Repository

1. Go to https://github.com/new
2. Repository name: `mylit-docs` (or your preferred name)
3. Description: "Documentation repository for mylit project"
4. Set to **Private** (recommended) or Public
5. **DO NOT** initialize with README, .gitignore, or license
6. Click "Create repository"

## Step 2: Note the Repository URL

After creating, you'll see a URL like:
- SSH: `git@github.com:JasonK94/mylit-docs.git`
- HTTPS: `https://github.com/JasonK94/mylit-docs.git`

**Use SSH format** (same as your main repo): `git@github.com:JasonK94/mylit-docs.git`

## Step 3: Run Migration Script

After creating the repository, run:
```bash
cd /home/user3/data_user3/git_repo/mylit
./scripts/setup-docs-subtree.sh <repository-url>
```

Example:
```bash
./scripts/setup-docs-subtree.sh git@github.com:JasonK94/mylit-docs.git
```

## What the Script Does

1. Creates temporary directory with docs-main/ contents
2. Initializes new git repository
3. Adds all files and commits with preserved history
4. Pushes to GitHub
5. Adds subtree to main repository
6. Configures remotes
7. Cleans up temporary files

## Verification

After setup, verify:
```bash
# Check subtree remote
git remote -v | grep docs-main

# Check docs-main/ is tracked
git ls-files docs-main/

# Test pull
./scripts/docs-sync.sh pull main
```

