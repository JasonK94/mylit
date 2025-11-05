# Install required R packages for debugging
# Run this in the myr_cursor conda environment

# Install BiocManager if not available
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install CRAN packages
cran_packages <- c(
  "Seurat", "dplyr", "qs", "Matrix", "limma", "edgeR", "data.table"
)

# Install Bioconductor packages
bioc_packages <- c(
  "SingleCellExperiment", "SummarizedExperiment", "muscat", 
  "MAST", "nebula", "S4Vectors"
)

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install CRAN packages
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing CRAN package: ", pkg)
    install.packages(pkg)
  } else {
    message("Package already installed: ", pkg)
  }
}

# Install Bioconductor packages
BiocManager::install(bioc_packages, ask = FALSE, update = FALSE)

message("\nPackage installation complete!")

