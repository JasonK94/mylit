# SNUBH Server Environment Reproduction Guide

> Exported: 2026-02-18 | Source: SNUBH server (user3)
> Purpose: 다른 서버에서 동일한 R + Python 분석 환경을 재현하기 위한 파일들

---

## Files

| File | Description | Size |
|------|-------------|------|
| `renv.lock` | R renv lockfile (1112 packages, R 4.3.2) | 1.1MB |
| `r_installed_packages.csv` | All R packages with versions and library paths | — |
| `r_session_info.txt` | R version, key package versions, library paths | — |
| `conda_scvi-tools.yml` | **Primary** conda env: scVI, MELD, scCODA, cNMF, scanpy | 387 deps |
| `pip_scvi-tools.txt` | pip freeze of scvi-tools env | 336 pkgs |
| `conda_scvi-env.yml` | Legacy scvi env (older versions) | — |
| `conda_cellbender-env.yml` | CellBender (ambient RNA removal) | — |
| `conda_demux_env.yml` | Souporcell + Solo (demultiplexing) | — |
| `system_info.txt` | OS, GPU, paths, notes | — |

---

## Reproduction Steps

### 1. R Environment (renv)

```bash
# Install R 4.3.2
# Then:
mkdir -p ~/project && cd ~/project
cp renv.lock .

# In R:
install.packages("renv")
renv::restore()  # Installs all 1112 packages from lockfile
```

**Important**: renv.lock은 Bioconductor 3.18 + CRAN snapshot 기반. R 4.3.x 필수.

### 2. Python Environment (conda)

```bash
# Primary analysis env (scVI, MELD, scCODA, cNMF)
conda env create -f conda_scvi-tools.yml

# Or from pip (if conda fails due to platform differences):
conda create -n scvi-tools python=3.10
conda activate scvi-tools
pip install -r pip_scvi-tools.txt
```

**Key Python packages**:
- `scvi-tools` — scVI, scANVI integration
- `pertpy` — scCODA compositional analysis
- `meld` — condition density estimation
- `cnmf` — consensus NMF gene programs
- `scanpy`, `anndata` — single-cell data handling
- `cellrank` — trajectory inference (optional)

### 3. Optional Envs

```bash
# CellBender (GPU required)
conda env create -f conda_cellbender-env.yml

# Demultiplexing (Souporcell + Solo)
conda env create -f conda_demux_env.yml
```

---

## Key R Packages (versions)

| Package | Version | Purpose |
|---------|---------|---------|
| Seurat | 5.3.1 | Core scRNAseq framework |
| SeuratObject | 5.2.0 | Data structures |
| qs | 0.27.3 | Fast serialization (.qs files) |
| CellChat | 2.2.0 | Cell-cell interaction |
| multinichenetr | 2.1.0 | MultiNicheNet CCI |
| nichenetr | 2.2.1.1 | NicheNet ligand-receptor |
| muscat | 1.16.0 | Pseudobulk DEG |
| nebula | 1.5.6 | NBMM DEG (large datasets) |
| MAST | 1.28.0 | Single-cell DEG |
| DESeq2 | 1.42.1 | Differential expression |
| edgeR | 4.0.16 | Differential expression |
| miloR | 1.10.0 | Neighbourhood DA |
| slingshot | 2.10.0 | Trajectory inference |
| monocle3 | 1.4.26 | Trajectory + pseudotime |
| lme4 | 1.1.37 | Mixed-effects models (MASC) |
| speckle | 1.2.0 | propeller compositional DA |
| ComplexHeatmap | 2.18.0 | Heatmap visualization |
| MOFA2 | 1.12.1 | Multi-omics factor analysis |
| NMF | 0.27 | Non-negative matrix factorization |
| glmGamPoi | 1.14.3 | Fast GLM for scRNAseq |
| Augur | 1.0.3 | Cell type prioritization |

---

## Notes

1. **monocle3**는 renv 경로에만 설치됨. `.libPaths()` 설정 필수:
   ```r
   .libPaths(c("path/to/renv/library/R-4.3/x86_64-pc-linux-gnu", .libPaths()))
   ```

2. **S4Vectors** 충돌: `as.data.frame` override 방지를 위해 Suggests에 배치

3. **dplyr::first()** masking: Bioconductor 패키지 로드 후 dplyr를 마지막에 로드

4. **GPU**: scVI/scANVI/CellBender는 CUDA GPU 필요 (RTX 4090 사용)

5. **Platform**: `conda_*.yml`은 linux-64 기준. macOS/ARM에서는 `pip_scvi-tools.txt`로 재현 권장
