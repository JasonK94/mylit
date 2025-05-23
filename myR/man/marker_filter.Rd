% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/markers.R
\name{marker_filter}
\alias{marker_filter}
\title{Filter out unwanted genes from marker results}
\usage{
marker_filter(markers, filter = c("rb", "mt", "hb", "AC", "ENSG", "LINC"))
}
\arguments{
\item{markers}{A data frame from Seurat's FindMarkers or FindAllMarkers results}

\item{filter}{Character vector specifying which gene types to filter out.
Options: "rb" (ribosomal), "mt" (mitochondrial), "hb" (hemoglobin),
"AC" (AC/AL genes), "ENSG" (ENSG genes), "LINC" (LINC genes)}
}
\value{
A filtered data frame with unwanted genes removed
}
\description{
This function removes genes matching specific patterns (ribosomal, mitochondrial,
hemoglobin, etc.) from marker gene results.
}
\examples{
\dontrun{
# Remove ribosomal and mitochondrial genes
markers_filtered <- marker_filter(markers, filter = c("rb", "mt"))

# Remove all unwanted gene types
markers_clean <- marker_filter(markers, filter = c("rb", "mt", "hb", "AC", "ENSG", "LINC"))
}
}
