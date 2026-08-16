# GSE65682 Galectin Analysis (MSCA-PF Preliminary Data)

**Project**: Glyco-Sep – Decoding and Targeting the Glyco-Immune Landscape in Sepsis  
**Author**: Luis Castillo Montanez  
**Date**: 16 August 2026  
**Repository**: https://github.com/luiscastillo718-cell/GSE65682-Galectin-Analysis-MSCA

## Overview
Reproducible re-analysis of the MARS consortium dataset GSE65682 focusing on LGALS1, LGALS3 and LGALS9.  
The analyses address two key scientific questions raised during proposal development:

1. Do endotype differences in galectin expression persist after adjustment for monocyte fraction, age and sex?
2. Are there coordinated glycosylation-related transcriptional programmes that co-vary with the galectins?

## Analyses included

### 1. Endotype-selective galectin expression (composition-adjusted)
- Linear models: `LGALS ~ Mars endotype + Monocyte_score + Age + Sex`
- Estimated marginal means for all four Mars endotypes
- Pairwise contrasts (Tukey-adjusted)
- Key finding: LGALS3 remains significantly higher in Mars1 after full adjustment

### 2. Glycosylation-related transcriptional modules
Curated panel of glycosylation machinery genes grouped into five modules:
- N-glycan branching (MGATs)
- LacNAc / poly-LacNAc extension (B3GNTs, GCNTs)
- Core-1 O-glycosylation (C1GALT1, C1GALT1C1)
- Sialylation (ST3GALs, ST6GALs, ST6GALNACs)
- Fucosylation (FUTs)

**Results**:
- All five modules differ highly significantly across Mars endotypes
- Multiple modules show significant Spearman correlations with LGALS1, LGALS3 and LGALS9 (BH-adjusted)

### 3. Earlier analyses (still present)
- Original endotype boxplots and survival analyses
- Robust ANCOVA of galectin–immune cell associations

## Key Files
| File | Description |
|------|-------------|
| `GSE65682_Galectin_Analysis_V13.R` | Original multi-panel analysis |
| `GSE65682_Robust_Galectin_Immune_Analysis.R` | ANCOVA + FDR immune association analysis |
| `GSE65682_Glyco_Modules_Analysis.R` | Glycosylation modules + correlations with galectins |
| `Galectin_Adjusted_Means_*.pdf / .csv` | Fully adjusted expression across endotypes |
| `Glyco_Modules_by_Endotype_Heatmap.pdf` | Module scores across Mars endotypes |
| `Glyco_Modules_vs_Galectins_Correlation.pdf` | Correlation heatmap (modules × galectins) |
| `session_info.txt` | Exact R / Bioconductor versions |

## How to Reproduce
1. Clone the repository
2. Place `GSE65682_series_matrix.txt.gz` in the working directory
3. Run the R scripts in R 4.5.1 / Bioconductor 3.21

All code and outputs are provided under the MIT license for full transparency and compliance with Horizon Europe FAIR data principles.
