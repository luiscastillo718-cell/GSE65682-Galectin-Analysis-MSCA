# GSE65682 Galectin Analysis (MSCA-PF Preliminary Data)

**Project**: Glyco-Sep – Decoding and Targeting the Glyco-Immune Landscape in Sepsis  
**Author**: Luis Castillo Montanez  
**Date**: 16 August 2026  
**Repository**: https://github.com/luiscastillo718-cell/GSE65682-Galectin-Analysis-MSCA

## Overview
Reproducible re-analysis of the MARS consortium dataset GSE65682 focusing on LGALS1, LGALS3 and LGALS9 transcript expression in human sepsis.  
The repository documents a progressive series of analyses that strengthen the rationale for the Glyco-Sep project.

---

## Analysis Workflow (Logical Order)

### 1. Original Endotype Analysis
Basic characterisation of LGALS1, LGALS3 and LGALS9 expression across the four Mars endotypes, including survival associations.

**Main files**
- `GSE65682_Galectin_Analysis_V13.R`
- `MSCA_Individual_Galectin_Analysis.pdf`

### 2. Robust Association with Inferred Immune Cell Fractions (ANCOVA)
Exploratory analysis testing whether galectin expression is associated with specific immune compartments after adjustment for age and endotype (Benjamini–Hochberg FDR correction).

**Main files**
- `GSE65682_Robust_Galectin_Immune_Analysis.R`
- `Robust_Galectin_Immune_Heatmap.pdf`
- `Galectin_Immune_Adjusted_Robust.csv`

### 3. Composition-Adjusted Endotype Effects (Monocyte + Age + Sex)
Addresses the critical question: *Do endotype differences in galectin expression persist after adjusting for monocyte abundance, age and sex?*

- Linear models: `LGALS ~ Mars endotype + Monocyte_score + Age + Sex`
- Estimated marginal means for all four endotypes
- Pairwise contrasts (Tukey-adjusted)

**Monocyte score genes**  
`CD14`, `CD68`, `CSF1R`, `LYZ`, `MAFB`, `FCGR1A`, `ITGAM`, `S100A8`, `S100A9`

**Key finding**: LGALS3 remains significantly higher in Mars1 after full adjustment.

**Main files**
- `GSE65682_Adjusted_Means_Significance.R`
- `Galectin_Adjusted_Means_All_Endotypes.pdf`
- `Galectin_Adjusted_Means_with_Significance.pdf`
- `Forest_Endotype_Effect_Full_Adjustment.pdf`
- `Galectin_Pairwise_Contrasts_Adjusted.csv`
- `GSE65682_galectin_monocyte_sex_master.csv`

### 4. Glycosylation-Related Transcriptional Modules
First “glyco” signal at the transcriptomic level. Curated modules of glycosylation machinery genes were evaluated for differences across Mars endotypes and correlation with LGALS1/3/9.

**Module gene composition**

| Module | Genes included |
|--------|----------------|
| **N-glycan branching** | MGAT1, MGAT2, MGAT3, MGAT4A, MGAT4B |
| **LacNAc / poly-LacNAc extension** | B3GNT2, B3GNT3, B3GNT5, B3GNT8, GCNT1, GCNT2 |
| **Core-1 O-glycosylation** | C1GALT1, C1GALT1C1 |
| **Sialylation** | ST3GAL1, ST3GAL4, ST3GAL6, ST6GAL1, ST6GALNAC1, ST6GALNAC2, ST6GALNAC4 |
| **Fucosylation** | FUT4, FUT7, FUT8 |

**Main files**
- `GSE65682_Glyco_Modules_Analysis.R`
- `Glyco_Modules_by_Endotype_Heatmap.pdf`
- `Glyco_Modules_vs_Galectins_Correlation.pdf`
- `Glyco_Modules_Endotype_Stats.csv`
- `Glyco_Modules_vs_Galectins_Correlation.csv`

---

## How to Reproduce
1. Clone the repository
2. Place `GSE65682_series_matrix.txt.gz` in the working directory
3. Run the R scripts in R 4.5.1 / Bioconductor 3.21

All code and outputs are provided under the MIT license for full transparency and compliance with Horizon Europe FAIR data principles.
