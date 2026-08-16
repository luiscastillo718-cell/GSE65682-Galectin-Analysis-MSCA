# ==============================================================================
# GSE65682 – Adjusted Galectin Expression with Key Significance Annotations
# Fully adjusted for Monocyte fraction + Age + Sex
# Standalone script (new R session)
# ==============================================================================

# 1. Setup
setwd("/Users/lcastillo6/Desktop/Personal/Caixa Research Institute/Human sepsis re-analysis/GSE65682_Galectin_Analysis_MSCA")
dir.create("GSE65682_Analysis_Outputs", showWarnings = FALSE)
set.seed(12345)

if (!require("pacman")) install.packages("pacman")
pacman::p_load(GEOquery, tidyverse, ggplot2, emmeans, ggsignif)

# 2. Load data
cat("Loading GSE65682 series matrix...\n")
gse <- getGEO(filename = "GSE65682_series_matrix.txt.gz", GSEMatrix = TRUE)

if (is.list(gse) && length(gse) >= 1) {
  eset <- gse[[1]]
} else {
  eset <- gse
}

expr_matrix <- exprs(eset)
p_data <- pData(eset)
cat("Data loaded. Samples:", ncol(expr_matrix), "\n")

# 3. Build final dataset
exp_clean <- data.frame(
  Sample = colnames(expr_matrix),
  LGALS1 = expr_matrix["11743721_at", ],
  LGALS3 = expr_matrix["11725937_a_at", ],
  LGALS9 = colMeans(expr_matrix[c("11716406_x_at", "11746558_x_at"), ])
)

final_df <- p_data %>%
  transmute(
    Sample       = geo_accession,
    SRS_Endotype = `endotype_class:ch1`,
    Age          = as.numeric(as.character(`age:ch1`)),
    Sex          = `gender:ch1`
  ) %>%
  mutate(Sex = factor(Sex, levels = c("female", "male"))) %>%
  inner_join(exp_clean, by = "Sample") %>%
  filter(!is.na(SRS_Endotype), SRS_Endotype != "Healthy", !is.na(Sex), !is.na(Age)) %>%
  mutate(SRS_Endotype = factor(SRS_Endotype, levels = c("Mars1", "Mars2", "Mars3", "Mars4")))

cat("Final cohort:", nrow(final_df), "samples\n")

# 4. Monocyte score
monocyte_markers <- c("CD14", "CD68", "CSF1R", "LYZ", "MAFB", "FCGR1A", "ITGAM", "S100A8", "S100A9")
fdata <- fData(eset)
symbol_col <- grep("Symbol|symbol|Gene.Symbol|GENE_SYMBOL", colnames(fdata), value = TRUE, ignore.case = TRUE)[1]

probe_map <- fdata %>%
  rownames_to_column("PROBEID") %>%
  dplyr::select(PROBEID, Symbol = all_of(symbol_col)) %>%
  filter(Symbol %in% monocyte_markers)

monocyte_probes <- intersect(probe_map$PROBEID, rownames(expr_matrix))
final_df$Monocyte_score <- colMeans(expr_matrix[monocyte_probes, final_df$Sample, drop = FALSE])
cat("Monocyte probes used:", length(monocyte_probes), "\n")

# 5. Fit fully adjusted models
model_LGALS1 <- lm(LGALS1 ~ SRS_Endotype + Monocyte_score + Age + Sex, data = final_df)
model_LGALS3 <- lm(LGALS3 ~ SRS_Endotype + Monocyte_score + Age + Sex, data = final_df)
model_LGALS9 <- lm(LGALS9 ~ SRS_Endotype + Monocyte_score + Age + Sex, data = final_df)

# 6. Estimated Marginal Means
emm_df <- bind_rows(
  as.data.frame(emmeans(model_LGALS1, ~ SRS_Endotype)) %>% mutate(Gene = "LGALS1"),
  as.data.frame(emmeans(model_LGALS3, ~ SRS_Endotype)) %>% mutate(Gene = "LGALS3"),
  as.data.frame(emmeans(model_LGALS9, ~ SRS_Endotype)) %>% mutate(Gene = "LGALS9")
)

print(emm_df)

# 7. Key pairwise contrasts (the ones you requested)
get_contrast_p <- function(model, g1, g2) {
  emm <- emmeans(model, ~ SRS_Endotype)
  cont <- as.data.frame(contrast(emm, method = "pairwise", adjust = "tukey"))
  idx <- which(cont$contrast == paste(g1, "-", g2) | cont$contrast == paste(g2, "-", g1))
  if (length(idx) == 0) return(NA_real_)
  cont$p.value[idx[1]]
}

key_contrasts <- tribble(
  ~Gene,    ~group1,  ~group2,  ~model,
  "LGALS1", "Mars2",  "Mars1",  model_LGALS1,
  "LGALS1", "Mars2",  "Mars3",  model_LGALS1,
  "LGALS1", "Mars2",  "Mars4",  model_LGALS1,
  "LGALS3", "Mars1",  "Mars2",  model_LGALS3,
  "LGALS3", "Mars1",  "Mars3",  model_LGALS3,
  "LGALS3", "Mars1",  "Mars4",  model_LGALS3,
  "LGALS9", "Mars4",  "Mars1",  model_LGALS9,
  "LGALS9", "Mars4",  "Mars2",  model_LGALS9,
  "LGALS9", "Mars4",  "Mars3",  model_LGALS9
)

key_contrasts <- key_contrasts %>%
  mutate(
    p.adj = pmap_dbl(list(model, group1, group2), get_contrast_p),
    label = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    )
  )

print(key_contrasts)

# 8. Set y-positions for the brackets (tuned for free scales)
y_positions <- emm_df %>%
  group_by(Gene) %>%
  summarise(max_y = max(upper.CL, na.rm = TRUE)) %>%
  mutate(
    y1 = max_y * 1.08,
    y2 = max_y * 1.16,
    y3 = max_y * 1.24
  )

key_contrasts <- key_contrasts %>%
  left_join(y_positions, by = "Gene") %>%
  group_by(Gene) %>%
  mutate(y.position = c(y1[1], y2[1], y3[1])) %>%
  ungroup()

# 9. Final plot with significance
p_final <- ggplot(emm_df, aes(x = SRS_Endotype, y = emmean, colour = SRS_Endotype)) +
  geom_point(size = 3.8) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.18, linewidth = 0.9) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 3) +
  scale_colour_manual(values = c(
    "Mars1" = "#E64B35",
    "Mars2" = "#4DBBD5",
    "Mars3" = "#00A087",
    "Mars4" = "#3C5488"
  )) +
  geom_signif(
    data = key_contrasts,
    aes(xmin = group1, xmax = group2, annotations = label, y_position = y.position),
    manual = TRUE,
    tip_length = 0.01,
    textsize = 4.2,
    vjust = -0.15,
    inherit.aes = FALSE
  ) +
  labs(
    title = "Adjusted Galectin Expression Across Mars Endotypes",
    subtitle = "Estimated marginal means ± 95% CI (adjusted for monocyte fraction, age and sex)\nKey pairwise contrasts indicated",
    x = "Mars Endotype",
    y = "Adjusted log₂ Expression"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 13),
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 10.5, hjust = 0.5, colour = "grey30"),
    axis.text.x = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold")
  )

print(p_final)

ggsave("GSE65682_Analysis_Outputs/Galectin_Adjusted_Means_with_Significance.pdf",
       p_final, width = 12, height = 5.5, dpi = 300)

# 10. Save supporting tables
write_csv(emm_df, "GSE65682_Analysis_Outputs/Galectin_Adjusted_Means_All_Endotypes.csv")
write_csv(key_contrasts, "GSE65682_Analysis_Outputs/Key_Pairwise_Contrasts.csv")

cat("\nSUCCESS\n")
cat("Plot saved: GSE65682_Analysis_Outputs/Galectin_Adjusted_Means_with_Significance.pdf\n")