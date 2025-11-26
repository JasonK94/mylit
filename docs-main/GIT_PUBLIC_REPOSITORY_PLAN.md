# Git Public Repository Migration Plan

## Overview

This document outlines the strategy for creating a public version of the repository while maintaining the private development repository.

---

## Current State

- **Private Repository**: Contains full development history, worktrees, and potentially sensitive information
- **Worktrees**: Multiple worktrees (`_wt/`) for parallel development
- **Branches**: Multiple feature branches (`deg-consensus-dev`, `pseudotime-dev`, etc.)

---

## Goal

1. Create a `main-public` branch that is clean and ready for public release
2. Convert current public remote repository to private
3. Create a new public remote repository
4. Maintain synchronization between private and public repositories

---

## Step-by-Step Plan

### Phase 1: Preparation

#### 1.1 Audit Current Repository

**Actions:**
- [x] Review codebase for sensitive information (see `SECURITY_REVIEW_REPORT.md`)
- [ ] Create checklist of files to clean
- [ ] Document current remote repository URLs

**Commands:**
```bash
# Check current remotes
git remote -v

# Check current branch
git branch --show-current

# List all branches
git branch -a
```

#### 1.2 Backup Current State

**Actions:**
- [ ] Create backup tag before major changes
- [ ] Document current remote repository settings

**Commands:**
```bash
# Create backup tag
git tag backup-before-public-prep-$(date +%Y%m%d)

# Push tag to remote (optional, for safety)
git push origin backup-before-public-prep-$(date +%Y%m%d)
```

---

### Phase 2: Create main-public Branch

#### 2.1 Create and Switch to main-public Branch

**Actions:**
- [ ] Create new branch from current main/master
- [ ] Switch to the new branch

**Commands:**
```bash
# Ensure you're on main/master branch
git checkout main  # or master, depending on your default branch

# Create and switch to main-public branch
git checkout -b main-public

# Verify you're on the new branch
git branch --show-current
```

#### 2.2 Clean Sensitive Information

**Actions:**
- [ ] Remove/replace absolute paths in documentation
- [ ] Update configuration files to use environment variables
- [ ] Replace user-specific paths with generic placeholders
- [ ] Review and clean all files listed in `SECURITY_REVIEW_REPORT.md`

**Process:**
1. Use search and replace for common patterns:
   ```bash
   # Find all occurrences (review before replacing)
   grep -r "/home/user3" docs-main/ scripts/ myR/
   grep -r "/data/user3" docs-main/ scripts/ myR/
   ```

2. Replace with generic paths or environment variables
3. Test that code still works with new paths
4. Commit changes incrementally

**Example Commits:**
```bash
# Clean README files
git add docs-main/README.md docs-main/README_KR.md
git commit -m "docs: Remove absolute paths from README files"

# Clean guide files
git add docs-main/guide.md docs-main/guide_KR.md
git commit -m "docs: Replace absolute paths with generic examples in guides"

# Update configuration
git add scripts/vars_config.R docs-main/config/vars_config.R
git commit -m "config: Make data paths configurable via environment variables"

# Continue for other files...
```

#### 2.3 Create Configuration System

**Actions:**
- [ ] Create `config/config.example.R` template
- [ ] Update code to use configuration system
- [ ] Document configuration in README

**Files to Create:**
- `config/config.example.R` - Example configuration
- `config/README.md` - Configuration documentation

#### 2.4 Final Review

**Actions:**
- [ ] Review all changes
- [ ] Test that examples work with generic paths
- [ ] Ensure no sensitive information remains
- [ ] Update documentation to reflect changes

**Commands:**
```bash
# Review changes
git log --oneline main..main-public

# Check for any remaining sensitive paths
grep -r "user3\|/home/user3\|/data/user3" --exclude-dir=.git

# Review diff
git diff main..main-public --stat
```

---

### Phase 3: Repository Management

#### 3.1 Make Current Remote Private

**Actions:**
- [ ] Access repository settings (GitHub/GitLab/etc.)
- [ ] Change repository visibility to Private
- [ ] Verify access restrictions

**Note:** This should be done on the hosting platform (GitHub, GitLab, etc.)

#### 3.2 Create New Public Repository

**Options:**

**Option A: Fork on Same Platform**
- Create a new public repository
- Add as new remote
- Push `main-public` branch

**Option B: Create Separate Repository**
- Create new repository on hosting platform
- Initialize as public
- Add as remote

**Commands:**
```bash
# Add new public remote (replace URL with actual public repo URL)
git remote add public <new-public-repo-url>

# Verify remotes
git remote -v

# Push main-public branch to public repository
git push public main-public:main

# Set upstream for future pushes
git push -u public main-public:main
```

#### 3.3 Update Remote Names (Optional but Recommended)

**Actions:**
- [ ] Rename current remote to `private` or `origin-private`
- [ ] Keep `public` remote for public repository

**Commands:**
```bash
# Rename origin to origin-private (if desired)
git remote rename origin origin-private

# Or keep origin as private and add public
git remote add public <public-repo-url>
```

---

### Phase 4: Maintenance Strategy

#### 4.1 Regular Sync Process

**Workflow for Syncing Private → Public:**

1. **Work on private repository** (main branch)
2. **Periodically sync to public:**
   ```bash
   # Switch to main-public branch
   git checkout main-public
   
   # Merge from main (this will bring in new changes)
   git merge main
   
   # Clean up any sensitive information that came from main
   # (Review and remove absolute paths, etc.)
   
   # Commit cleanup
   git add .
   git commit -m "chore: Clean sensitive info from main-public"
   
   # Push to public
   git push public main-public:main
   ```

#### 4.2 Automated Cleanup Script (Optional)

Create a script to automate path sanitization:

```bash
#!/bin/bash
# scripts/cleanup_for_public.sh

# Replace common patterns
find . -type f -name "*.md" -o -name "*.R" | while read file; do
  sed -i 's|/home/user3/data_user3/git_repo|path/to/repo|g' "$file"
  sed -i 's|/data/user3/sobj|path/to/data|g' "$file"
  sed -i 's|/home/user3/GJC_KDW_250721|path/to/workdir|g' "$file"
done
```

**Usage:**
```bash
# Before merging main → main-public
./scripts/cleanup_for_public.sh
git add .
git commit -m "chore: Sanitize paths for public release"
```

#### 4.3 Branch Strategy

**Private Repository:**
- `main` - Main development branch (private)
- `deg-consensus-dev` - Feature branch
- `pseudotime-dev` - Feature branch
- Other feature branches...

**Public Repository:**
- `main` - Public release branch (synced from `main-public`)
- Optionally: `dev` - Development preview (if desired)

**Workflow:**
1. Develop on private `main` or feature branches
2. Merge feature branches to private `main`
3. Periodically sync to `main-public` with cleanup
4. Push `main-public` to public `main`

---

## Checklist

### Pre-Migration
- [ ] Review `SECURITY_REVIEW_REPORT.md`
- [ ] Backup current repository state
- [ ] Document current remote URLs
- [ ] List all sensitive files to clean

### Branch Creation
- [ ] Create `main-public` branch
- [ ] Clean all sensitive information
- [ ] Update configuration system
- [ ] Test examples with generic paths
- [ ] Review all changes

### Repository Setup
- [ ] Make current remote private
- [ ] Create new public repository
- [ ] Add public remote
- [ ] Push `main-public` to public repository
- [ ] Verify public repository is accessible

### Documentation
- [ ] Update README with setup instructions
- [ ] Document configuration system
- [ ] Add contribution guidelines (if applicable)
- [ ] Create setup guide for new users

### Maintenance
- [ ] Set up sync workflow
- [ ] Create cleanup script (optional)
- [ ] Document maintenance process
- [ ] Schedule regular sync reviews

---

## Risk Mitigation

### Potential Issues

1. **Accidentally pushing sensitive info to public**
   - **Mitigation**: Always review changes before pushing
   - **Mitigation**: Use automated cleanup script
   - **Mitigation**: Keep private and public remotes clearly separated

2. **Breaking changes in public version**
   - **Mitigation**: Test examples before pushing
   - **Mitigation**: Use tags for releases
   - **Mitigation**: Maintain changelog

3. **Sync conflicts**
   - **Mitigation**: Regular, small syncs rather than large batches
   - **Mitigation**: Document merge process
   - **Mitigation**: Keep `main-public` as a dedicated branch

### Rollback Plan

If issues occur:

```bash
# Revert to previous state
git checkout main-public
git reset --hard <previous-commit-hash>

# Force push (use with caution)
git push public main-public:main --force
```

---

## Timeline Estimate

- **Phase 1 (Preparation)**: 1-2 hours
- **Phase 2 (Branch Creation & Cleanup)**: 4-8 hours (depending on codebase size)
- **Phase 3 (Repository Setup)**: 30 minutes - 1 hour
- **Phase 4 (Maintenance Setup)**: 1-2 hours

**Total**: Approximately 1-2 days of focused work

---

## Post-Migration Tasks

1. **Monitor public repository**
   - Check for issues reported by users
   - Review pull requests (if applicable)
   - Update documentation based on feedback

2. **Regular maintenance**
   - Weekly or bi-weekly sync from private to public
   - Review and clean sensitive information
   - Update public documentation

3. **Version management**
   - Tag releases in public repository
   - Maintain changelog
   - Document version differences if needed

---

## Notes

- Keep `main-public` branch local and sync to public repository
- Never merge from public to private (to avoid bringing in external changes)
- Consider using GitHub Actions or similar for automated cleanup (future enhancement)
- Document any differences between private and public versions

---

## Questions to Resolve

Before proceeding, clarify:

1. **Repository Hosting**: Which platform? (GitHub, GitLab, etc.)
2. **Access Control**: Who should have access to private repository?
3. **Public License**: What license for public repository?
4. **Contributions**: Will you accept contributions to public repository?
5. **Sync Frequency**: How often to sync private → public?
6. **Worktree Strategy**: How to handle worktrees in public repository? (May not be relevant for public users)

---

## References

- `SECURITY_REVIEW_REPORT.md` - Detailed security review
- `docs-main/README.md` - Public-facing README
- `docs-main/README_KR.md` - Korean README

