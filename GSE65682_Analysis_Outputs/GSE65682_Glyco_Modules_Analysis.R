# ==============================================================================
# GSE65682 – Glycosylation-related transcriptional modules
# 1. Differences across Mars endotypes
# 2. Correlation with LGALS1, LGALS3, LGALS9
# Full standalone script (new R session)
# ==============================================================================

# 1. Setup
setwd("/Users/lcastillo6/Desktop/Personal/Caixa Research Institute/Human sepsis re-analysis/GSE65682_Galectin_Analysis_MSCA")
dir.create("GSE65682_Analysis_Outputs", showWarnings = FALSE)
set.seed(12345)

if (!require("pacman")) install.packages("pacman")
pacman::p_load(GEOquery, tidyverse, pheatmap, RColorBrewer)

# 2. Load series matrix (robust)
cat("Loading GSE65682 series matrix...\n")
gse <- getGEO(filename = "GSE65682_series_matrix.txt.gz", GSEMatrix = TRUE)

if (is.list(gse) && length(gse) >= 1) {
  eset <- gse[[1]]
} else if (inherits(gse, "ExpressionSet")) {
  eset <- gse
} else {
  stop("Unexpected object returned by getGEO")
}

expr_matrix <- exprs(eset)
p_data <- pData(eset)
cat("Data loaded. Samples:", ncol(expr_matrix), "\n")

# 3. Build final_df with galectins + endotype
exp_clean <- data.frame(
  Sample = colnames(expr_matrix),
  LGALS1 = expr_matrix["11743721_at", ],
  LGALS3 = expr_matrix["11725937_a_at", ],
  LGALS9 = colMeans(expr_matrix[c("11716406_x_at", "11746558_x_at"), ])
)

final_df <- p_data %>%
  transmute(
    Sample = geo_accession,
    SRS_Endotype = `endotype_class:ch1`
  ) %>%
  inner_join(exp_clean, by = "Sample") %>%
  filter(!is.na(SRS_Endotype), SRS_Endotype != "Healthy") %>%
  mutate(SRS_Endotype = factor(SRS_Endotype, levels = c("Mars1", "Mars2", "Mars3", "Mars4")))

cat("Final cohort:", nrow(final_df), "samples\n")
print(table(final_df$SRS_Endotype))

# 4. Curated glycosylation gene panel
glyco_genes <- list(
  N_branching  = c("MGAT1", "MGAT2", "MGAT3", "MGAT4A", "MGAT4B", "MGAT5"),
  LacNAc       = c("B3GNT2", "B3GNT3", "B3GNT5", "B3GNT8", "GCNT1", "GCNT2"),
  Core1_O      = c("C1GALT1", "C1GALT1C1"),
  Sialylation  = c("ST3GAL1", "ST3GAL4", "ST3GAL6", "ST6GAL1", "ST6GALNAC1", "ST6GALNAC2", "ST6GALNAC4"),
  Fucosylation = c("FUT4", "FUT7", "FUT8", "FUT9", "FUT11")
)

all_glyco <- unique(unlist(glyco_genes))

# 5. Map genes to probes
fdata <- fData(eset)
symbol_col <- grep("Symbol|symbol|Gene.Symbol|GENE_SYMBOL", colnames(fdata), value = TRUE, ignore.case = TRUE)[1]

probe_map <- fdata %>%
  rownames_to_column("PROBEID") %>%
  dplyr::select(PROBEID, Symbol = all_of(symbol_col)) %>%
  filter(Symbol %in% all_glyco)

available_genes <- unique(probe_map$Symbol)
cat("Glycosylation genes found on array:", length(available_genes), "\n")
print(sort(available_genes))

# 6. Average multiple probes per gene
glyco_expr <- sapply(available_genes, function(g) {
  probes <- probe_map$PROBEID[probe_map$Symbol == g]
  if (length(probes) == 1) {
    expr_matrix[probes, final_df$Sample]
  } else {
    colMeans(expr_matrix[probes, final_df$Sample, drop = FALSE])
  }
}) %>% as.data.frame() %>% rownames_to_column("Sample")

final_df <- final_df %>% inner_join(glyco_expr, by = "Sample")

# 7. Create module scores
for (mod in names(glyco_genes)) {
  genes_in_mod <- intersect(glyco_genes[[mod]], available_genes)
  if (length(genes_in_mod) >= 1) {
    final_df[[paste0(mod, "_score")]] <- rowMeans(final_df[, genes_in_mod, drop = FALSE])
    cat(mod, "score created with", length(genes_in_mod), "genes:", paste(genes_in_mod, collapse = ", "), "\n")
  }
}

module_cols <- grep("_score$", colnames(final_df), value = TRUE)

# 8. Statistical test: modules across endotypes
module_stats <- map_dfr(module_cols, function(m) {
  kw <- kruskal.test(as.formula(paste(m, "~ SRS_Endotype")), data = final_df)
  data.frame(Module = m, Kruskal_p = kw$p.value)
}) %>%
  mutate(Kruskal_p_adj = p.adjust(Kruskal_p, method = "BH"))

cat("\nModule differences across Mars endotypes:\n")
print(module_stats)

# 9. Spearman correlation: modules vs galectins
cor_results <- expand_grid(
  Module = module_cols,
  Galectin = c("LGALS1", "LGALS3", "LGALS9")
) %>%
  mutate(
    test = map2(Module, Galectin, ~cor.test(final_df[[.x]], final_df[[.y]], method = "spearman")),
    rho  = map_dbl(test, "estimate"),
    p    = map_dbl(test, "p.value")
  ) %>%
  select(-test) %>%
  mutate(
    p_adj = p.adjust(p, method = "BH"),
    Significance = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ ""
    )
  )

cat("\nCorrelations with galectins:\n")
print(cor_results)

# 10. Heatmap 1 – Modules across endotypes (z-scored) – FIXED
module_means <- final_df %>%
  filter(!is.na(SRS_Endotype)) %>%
  group_by(SRS_Endotype) %>%
  summarise(across(all_of(module_cols), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
  column_to_rownames("SRS_Endotype") %>%
  as.matrix() %>%
  t()

module_means <- module_means[complete.cases(module_means), ]
module_means_z <- t(scale(t(module_means)))

pdf("GSE65682_Analysis_Outputs/Glyco_Modules_by_Endotype_Heatmap.pdf", width = 6.5, height = 5)
pheatmap(module_means_z,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         main = "Glycosylation-related transcriptional modules\nacross Mars endotypes (z-scored)",
         fontsize = 12,
         border_color = NA)
dev.off()

# 11. Heatmap 2 – Modules vs Galectins
cor_mat <- cor(final_df[, c(module_cols, "LGALS1", "LGALS3", "LGALS9")],
               method = "spearman", use = "pairwise.complete.obs")
cor_plot_mat <- cor_mat[module_cols, c("LGALS1", "LGALS3", "LGALS9")]

sig_mat <- matrix(cor_results$Significance,
                  nrow = length(module_cols), byrow = FALSE,
                  dimnames = list(module_cols, c("LGALS1", "LGALS3", "LGALS9")))

pdf("GSE65682_Analysis_Outputs/Glyco_Modules_vs_Galectins_Correlation.pdf", width = 6, height = 5)
pheatmap(cor_plot_mat,
         display_numbers = sig_mat,
         number_color = "black",
         fontsize_number = 13,
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         breaks = seq(-0.5, 0.5, length.out = 101),
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         main = "Glycosylation-related modules vs Galectins\n(Spearman ρ + BH-adjusted significance)",
         fontsize = 12,
         border_color = NA)
dev.off()

# 12. Save all tables
write_csv(module_stats, "GSE65682_Analysis_Outputs/Glyco_Modules_Endotype_Stats.csv")
write_csv(cor_results, "GSE65682_Analysis_Outputs/Glyco_Modules_vs_Galectins_Correlation.csv")
write_csv(final_df %>% select(Sample, SRS_Endotype, LGALS1, LGALS3, LGALS9, all_of(module_cols)),
          "GSE65682_Analysis_Outputs/GSE65682_glyco_modules_master.csv")

cat("\n✅ Analysis complete. All files saved in GSE65682_Analysis_Outputs/\n")