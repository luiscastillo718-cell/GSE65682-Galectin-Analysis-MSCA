# ==============================================================================
# Script: GSE65682_Robust_Galectin_Immune_Analysis.R
# Purpose: Robust analysis of Galectin-Immune cell associations
# ==============================================================================

# 1. Load Libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4) # Kept if needed, though we use lm() for robustness

# 2. Load Data (Assumes master_data is already in environment)
# master_data <- readRDS("Processed_Caches/GSE65682_survival_deconvoluted_master.rds")

# 3. Define Variables
cell_cols <- c("B_Cells", "Monocytes_Myel", "Neutrophils", "NK_Cells", 
               "T_Cells_CD4", "T_Cells_CD8")
gene_cols <- c("LGALS1", "LGALS3", "LGALS9")

# 4. Define Robust Statistical Function (ANCOVA-based)
get_robust_stats_fixed <- function(df) {
  # Skip if insufficient data
  if (sd(df$Proportion, na.rm = TRUE) == 0) return(data.frame(Correlation = NA, P.Value = NA))
  
  # Model as fixed effects (ANCOVA approach)
  model <- lm(Proportion ~ Expression + Age + SRS_Endotype, data = df)
  
  s <- summary(model)
  coefs <- s$coefficients
  
  # Extract results for 'Expression' (the gene of interest)
  if ("Expression" %in% rownames(coefs)) {
    est <- coefs["Expression", "Estimate"]
    p_val <- coefs["Expression", "Pr(>|t|)"]
    return(data.frame(Correlation = est, P.Value = p_val))
  } else {
    return(data.frame(Correlation = NA, P.Value = NA))
  }
}

# 5. Execute Analysis Pipeline
robust_stats <- master_data %>%
  # Reshape data
  pivot_longer(cols = all_of(cell_cols), names_to = "cell_type", values_to = "Proportion") %>%
  pivot_longer(cols = all_of(gene_cols), names_to = "Gene", values_to = "Expression") %>%
  # Group and Calculate
  group_by(cell_type, Gene) %>%
  do(get_robust_stats_fixed(.)) %>%
  ungroup() %>%
  # Apply FDR correction
  mutate(FDR = p.adjust(P.Value, method = "BH"),
         sig = case_when(
           FDR < 0.05 ~ "*",
           TRUE       ~ ""
         ))

# 6. Export Results
write.csv(robust_stats, "Galectin_Immune_Adjusted_Robust.csv", row.names = FALSE)

# 7. Generate Visualization
ggplot(robust_stats, aes(x = Gene, y = cell_type, fill = Correlation)) +
  geom_tile() +
  geom_text(aes(label = sig), size = 6) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  theme_minimal() +
  labs(title = "Robust Galectin-Immune Association",
       subtitle = "Coefficient Estimates (Adjusted for Age and Endotype)",
       x = "Galectin Gene", y = "Cell Type")

# Save the plot
ggsave("Robust_Galectin_Immune_Heatmap.pdf", width = 8, height = 6)