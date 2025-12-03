# Security Review & Codebase Cleanup Report

## Executive Summary

This report identifies sensitive information, hardcoded paths, and user-specific data that should be removed or generalized before making the repository public. It also provides recommendations for improving file naming conventions and project structure for better user experience.

---

## 1. Sensitive Information Found

### 1.1 Absolute Paths with User Information

**Critical Issues:**
- `/home/user3/` - User-specific home directory paths
- `/data/user3/` - User-specific data directory paths
- `/home/user3/data_user3/git_repo/` - Full repository paths with user info
- `/home/user3/GJC_KDW_250721` - Project-specific working directory

**Files Affected:**
- `docs-main/README_KR.md` - ✅ Fixed (replaced with generic paths)
- `docs-main/README.md` - ✅ Fixed (replaced with generic paths)
- `docs-main/guide.md` - Contains multiple absolute paths
- `docs-main/guide_KR.md` - Contains multiple absolute paths
- `scripts/vars_config.R` - Contains hardcoded data paths
- `docs-main/config/vars_config.R` - Contains hardcoded data paths
- `myR/R/CCI.R` - Contains default paths like `/data/user3/git_repo/human`
- Multiple documentation files in `docs/` directories
- Multiple test scripts in `scripts/` directories

**Recommendation:**
- Replace all absolute paths with relative paths or environment variables
- Use `here::here()` or `getwd()` for relative path resolution
- Create configuration files that users can customize
- Use `file.path()` with user-configurable base directories

### 1.2 Dataset-Specific File Names

**Found References:**
- `IS6_sex_added_251110.qs` - Specific dataset identifier
- `IS6_sex_added_251110_ds2500.qs` - Downsampled version
- `IS_scvi_251107_ds2500.qs` - Another dataset variant
- `data_seurat_251104.qs` - GeoMx dataset
- `IBD_biologics.rds` - IBD dataset

**Files Affected:**
- `scripts/vars_config.R` - Hardcoded dataset paths
- `docs-main/guide.md` - Example dataset references
- Multiple documentation files

**Recommendation:**
- Keep dataset names in example code but make paths configurable
- Add comments explaining these are example datasets
- Consider using placeholder names like `example_dataset.qs`

### 1.3 Project-Specific Directory Names

**Found:**
- `GJC_KDW_250721` - Appears to be a project identifier with initials and date
- `_wt/` - Worktree directory (may be acceptable if documented)

**Recommendation:**
- `GJC_KDW_250721` should be replaced with generic examples or environment variables
- `_wt/` is acceptable if properly documented as a worktree convention

### 1.4 Internal Server Paths

**Found:**
- `/data/user3/git_repo/human` - NicheNet data directory
- `/data/user3/git_repo/mouse` - NicheNet data directory

**Files Affected:**
- `myR/R/CCI.R` - Default paths for NicheNet data
- `docs/cci/CCI_INTEGRATED_GUIDE.md` - Documentation references

**Recommendation:**
- Make these configurable with fallback to current working directory
- Document that users should set their own data directories
- Use environment variables or configuration files

---

## 2. File Naming Improvements

### 2.1 Inconsistent Naming Conventions

**Current Issues:**

1. **Mixed Case in File Names:**
   - `CCI.R` vs `cci/` directory (should be consistent)
   - `run_cci_analysis.R` (good) vs `CCI.R` (inconsistent)

2. **Abbreviations:**
   - `pt.umap` - Uses dot notation (uncommon in file systems)
   - `deg-consensus-dev` - Uses hyphen (good)
   - `pseudotime-dev` - Uses hyphen (good)

3. **Test Files:**
   - `test_cci.R` (good)
   - `test_cci_interactive.R` (good)
   - `test_cci_minimal.R` (good)
   - `test_cci_simple.R` (good)
   - `test_cci_user_code.R` (good)
   - `test_cci_actual.R` (unclear naming)

**Recommendations:**

1. **Standardize Module Names:**
   - Consider renaming `pt.umap` to `pt_umap` or `patient_umap` for better filesystem compatibility
   - Keep hyphenated names for branches (`deg-consensus-dev`)

2. **Function File Naming:**
   - Use lowercase with underscores: `cci_analysis.R` instead of `CCI.R`
   - Or use PascalCase consistently: `CCI.R`, `DEGConsensus.R`, etc.
   - Current mix of `CCI.R` and `run_cci_analysis.R` is confusing

3. **Test File Organization:**
   - Group related tests: `test_cci_*.R` is good
   - Consider subdirectories: `tests/unit/`, `tests/integration/`
   - Archive old test files: `tests/archive/`

### 2.2 Documentation File Naming

**Current Structure (Good):**
- `*_INTEGRATED_GUIDE.md` - Clear naming
- `*_INTEGRATED_GUIDE_KR.md` - Clear localization
- `README.md` - Standard
- `TEST_INSTRUCTIONS.md` - Clear purpose

**Recommendations:**
- Current documentation naming is generally good
- Consider adding version numbers for major guides: `CCI_INTEGRATED_GUIDE_v1.md`
- Archive old versions: `archive/CCI_INTEGRATED_GUIDE_v0.md`

### 2.3 Directory Structure Improvements

**Current Issues:**

1. **Mixed Directory Patterns:**
   ```
   myR/R/cci/          # lowercase
   myR/R/analysis/     # lowercase
   myR/R/deg_consensus/ # lowercase with underscore
   ```

2. **Worktree References:**
   - `_wt/` directory name is not self-explanatory
   - Consider `worktrees/` or document it clearly

**Recommendations:**

1. **Standardize Directory Names:**
   - Use lowercase with underscores: `deg_consensus` instead of `deg-consensus` (for directories)
   - Keep hyphens for branch names only

2. **Document Worktree Structure:**
   - Add `WORKTREE_README.md` explaining the `_wt/` convention
   - Or rename to `worktrees/` for clarity

---

## 3. Code Quality Issues

### 3.1 Hardcoded Configuration

**Files Requiring Updates:**

1. **`scripts/vars_config.R`** and **`docs-main/config/vars_config.R`:**
   ```r
   base_path <- "/data/user3/sobj"  # Should be configurable
   ```

2. **`myR/R/CCI.R`:**
   ```r
   default_paths <- c(
     file.path("/data/user3/git_repo", species),
     file.path(getwd(), paste0("nichenet_data_", species))
   )
   ```

**Recommendation:**
- Create `config/config.R` or use environment variables
- Provide example configuration files: `config/config.example.R`
- Document required environment variables

### 3.2 Example Code in Documentation

**Current State:**
- Many documentation files contain hardcoded paths in examples
- Examples are not runnable by external users

**Recommendation:**
- Replace with generic placeholders: `path/to/your/data.qs`
- Add setup instructions for users to configure their paths
- Provide working examples that use relative paths or `here::here()`

---

## 4. Git Branch Strategy for Public Repository

### 4.1 Proposed Branch Structure

**Current State:**
- Main repository (private)
- Multiple worktrees for development

**Proposed Structure:**

1. **`main-public`** branch:
   - Clean, public-ready version
   - No sensitive information
   - Generic examples only
   - Well-documented

2. **`main`** branch (private):
   - Original private repository
   - Contains all development history
   - May contain sensitive paths (acceptable for private repo)

3. **Development branches:**
   - Keep existing branch structure
   - `deg-consensus-dev`, `pseudotime-dev`, etc.

### 4.2 Migration Plan

**Step 1: Create `main-public` branch**
```bash
git checkout -b main-public
# Clean up sensitive information
# Remove hardcoded paths
# Update documentation
git commit -m "Prepare public release"
```

**Step 2: Update Remote Repository**
```bash
# Make current remote private
# (Change repository settings on GitHub/GitLab)

# Create new public repository
# Fork or create new repo for public version
git remote add public <new-public-repo-url>
git push public main-public:main
```

**Step 3: Maintain Sync**
- Regularly merge `main` → `main-public` (with cleanup)
- Use automated scripts to sanitize paths before merging
- Keep sensitive information only in private `main` branch

### 4.3 Files to Exclude from Public Repository

**Consider `.gitignore` additions:**
- User-specific configuration files
- Local test data paths
- Internal documentation with sensitive information

**Consider separate branches for:**
- Internal documentation (`docs-internal/`)
- Development logs with sensitive paths (`devlogs-internal/`)

---

## 5. Priority Actions

### High Priority (Before Public Release)

1. ✅ **Replace absolute paths in README files** - DONE
2. **Replace paths in `guide.md` and `guide_KR.md`**
3. **Update `scripts/vars_config.R` to use configurable paths**
4. **Update `myR/R/CCI.R` default paths**
5. **Review and clean all documentation files in `docs/`**
6. **Review and clean all test scripts in `scripts/`**

### Medium Priority (Improvements)

1. **Standardize file naming conventions**
2. **Create configuration system for paths**
3. **Add example configuration files**
4. **Document worktree structure clearly**
5. **Archive old test files**

### Low Priority (Nice to Have)

1. **Rename `pt.umap` to `pt_umap` or `patient_umap`**
2. **Rename `_wt/` to `worktrees/`**
3. **Add version numbers to major guides**
4. **Reorganize test file structure**

---

## 6. Recommended File Changes Summary

### Files Requiring Immediate Attention

| File | Issue | Priority |
|-----|-------|----------|
| `docs-main/guide.md` | Multiple absolute paths | High |
| `docs-main/guide_KR.md` | Multiple absolute paths | High |
| `scripts/vars_config.R` | Hardcoded data paths | High |
| `docs-main/config/vars_config.R` | Hardcoded data paths | High |
| `myR/R/CCI.R` | Default server paths | High |
| `docs/cci/CCI_INTEGRATED_GUIDE.md` | Example paths | Medium |
| `docs/cci/CCI_INTEGRATED_GUIDE_KR.md` | Example paths | Medium |
| All `scripts/*/test_*.R` files | Hardcoded paths | Medium |
| All `docs/*/archive/*.md` files | Historical paths | Low |

### Files with Good Practices (No Changes Needed)

- `docs-main/README.md` - ✅ Already cleaned
- `docs-main/README_KR.md` - ✅ Already cleaned
- Most documentation structure is well-organized

---

## 7. Configuration System Proposal

### 7.1 Recommended Approach

Create a configuration system:

```r
# config/config.R
get_config <- function() {
  list(
    data_dir = Sys.getenv("MYR_DATA_DIR", unset = file.path(getwd(), "data")),
    worktree_dir = Sys.getenv("MYR_WORKTREE_DIR", unset = file.path(getwd(), "worktrees")),
    nichenet_data_dir = Sys.getenv("MYR_NICHENET_DIR", unset = file.path(getwd(), "nichenet_data"))
  )
}
```

### 7.2 Example Configuration File

```r
# config/config.example.R
# Copy this to config/config.R and customize

MYR_CONFIG <- list(
  data_dir = "/path/to/your/data",
  worktree_dir = "/path/to/worktrees",
  nichenet_data_dir = "/path/to/nichenet/data"
)
```

---

## Conclusion

The codebase contains significant amounts of user-specific and environment-specific information that should be removed or made configurable before public release. The main concerns are:

1. **Absolute paths** throughout documentation and code
2. **Hardcoded configuration** in multiple files
3. **Dataset-specific references** that may not be relevant to external users

The file naming conventions are generally good but could benefit from standardization. The proposed `main-public` branch strategy will allow maintaining a clean public version while preserving the full development history in the private repository.

**Next Steps:**
1. Review and approve this report
2. Create `main-public` branch
3. Systematically clean files according to priority
4. Set up configuration system
5. Test public repository setup
6. Create public repository and sync

