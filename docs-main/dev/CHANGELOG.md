# Changelog

All notable changes to `myR` are documented here. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/) and targets [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Signature workflow v5.2 with validation utilities and refreshed documentation (`2d9c50e`, `39a6159`).
- Milo differential abundance pipeline (`milo_opus6`) and meta-level gene importance diagnostics (`18f1ad1`, `b260bb1`, `a360837`, `7bf2904`).
- Expanded analysis test suite: TestLISI, PERMANOVA, PCP, PTMFL, plus safer parallel defaults for meta learners (`25b205b`, `5863176`).
- Comprehensive pseudobulk DEG module `pb_deg.R` and migration guide (`881589c` lineage, incorporating Claude branch refinements).
- Project bootstrap documentation (`context.md`, `NEXT_STEPS.md`) and bilingual records added during initialization (`cf6d13e`, `18c2810`).

### Changed
- Migrated Seurat usage from `slot` to `layer` conventions to align with v5 (`881589c`, `39a6159`).
- Reorganized R scripts, quarantined deprecated modules, and standardized MUSCAT/NEBULA utilities (`497eab8`, `bc19660`).
- Hardened `.gitignore` and repository hygiene to exclude large project outputs (`80f7941`, `9774761`).

### Removed
- Cleansed legacy `projects/` notebooks and `trash/` artifacts from version control (`80f7941`).
- Pruned unstable MUSCAT/MAST experiments, leaving only verified routines (`80260af`).

### Fixed
- Declared `NMF` dependency in DESCRIPTION to resolve `nmf()` availability issues (`e785399`).
- Improved PERMANOVA metadata robustness and rowname handling in testing utilities (`9cb9d3a`).

