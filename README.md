# GSE65682 Galectin Analysis (MSCA-PF Preliminary Data)

**Project**: Glyco-Sep – Decoding and Targeting the Glyco-Immune Landscape in Sepsis  
**Author**: Luis Castillo Montanez  
**Date**: 19 May 2026  
**Repository**: https://github.com/luiscastillo718-cell/GSE65682-Galectin-Analysis-MSCA

## Overview
Reproducible re-analysis of GSE65682 (MARS cohort) focusing on LGALS1, LGALS3, and LGALS9 transcript expression across SRS/Mars endotypes and 28-day survival.

## Key Files
- `GSE65682_Galectin_Analysis_V13.R` – complete reproducible script  
- `MSCA_Individual_Galectin_Analysis.pdf` – multi-panel figure (KM curves, endotype boxplots, Cox HRs, correlation matrix)  
- `GSE65682_Robust_Galectin_Immune_Analysis.R` – latest robust ANCOVA + FDR analysis of galectin–immune cell associations  
- `Galectin_Immune_Adjusted_Robust.csv` – statistical results table  
- `Robust_Galectin_Immune_Heatmap.pdf` – publication-ready heatmap  
- `session_info.txt` – exact R/Bioconductor versions used

## How to Reproduce
1. Clone the repository  
2. Place `GSE65682_series_matrix.txt.gz` in the working directory  
3. Run the `.R` scripts in R 4.5.1 / Bioconductor 3.21

All code and outputs are provided under the MIT license for full transparency and compliance with Horizon Europe FAIR data principles.

---

**Latest Addition (19 May 2026)**: Robust galectin–immune cell association analysis using ANCOVA models adjusted for age and SRS endotype, with Benjamini–Hochberg FDR correction. This exploratory analysis provides preliminary evidence of differential galectin–immune subset relationships and supports the need for cell-type-resolved multi-omics in Glyco-Sep.
How to update it: