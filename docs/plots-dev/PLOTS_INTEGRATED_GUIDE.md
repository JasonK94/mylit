# Plotting Module Integrated Guide

This document is the integrated guide for the Plotting module. (Work in progress)

## 1. Introduction

### Purpose
Provides various visualization functions based on Seurat objects and data frames to effectively represent analysis results.

## 2. Functions

The main functions currently under development/refactoring are as follows (see `archive/function_issues.md` for detailed analysis):

*   **upset_gene_lists**: Visualization of intersections between gene lists (Upset plot).
*   **vln_p**: Violin plot with statistical significance (p-value) display.
*   **cmb**: Proportional Bar Graph.
*   **acmb**: Absolute Count Bar Graph.
*   **cml**: Cumulative Line Graph.
*   **cdf**: Cumulative Distribution Function plot.
*   **cdf_multi**: Multi-variable CDF plot.

## 3. Development Status

*   Refactoring is in progress to ensure input flexibility (supporting both Seurat/data.frame) and resolve package dependency issues.

