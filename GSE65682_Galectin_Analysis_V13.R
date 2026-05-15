# ==============================================================================
# MSCA-GRADE ANALYSIS V13.0: Galectin_Analysis_MSCA
# ==============================================================================

# 1. ENVIRONMENT & REPRODUCIBILITY
set.seed(12345)
if (!require("pacman")) install.packages("pacman")
pacman::p_load(GEOquery, limma, tidyverse, ggpubr, survival, 
               survminer, broom, RColorBrewer, rstatix, patchwork, ggcorrplot)

# 2. ROBUST DATA LOADING (Local Matrix Access)
# Using local file to ensure reproducibility independent of GEO servers
cat("Step 1: Loading local data...\n")
gse <- getGEO(filename = "GSE65682_series_matrix.txt.gz", GSEMatrix = TRUE)
p_data <- pData(gse)
expr_matrix <- exprs(gse)

# 3. PROBE SELECTION & DOCUMENTATION
# LGALS1: 11743721_at | LGALS3: 11725937_a_at 
# LGALS9: Averaging 11716406_x_at & 11746558_x_at (High-coverage verification)
cat("Step 2: Processing verified probes...\n")
exp_clean <- data.frame(
  Sample = colnames(expr_matrix),
  LGALS1 = expr_matrix["11743721_at", ],
  LGALS3 = expr_matrix["11725937_a_at", ],
  LGALS9 = colMeans(expr_matrix[c("11716406_x_at", "11746558_x_at"), ])
)

# 4. PHENOTYPE CLEANING & ATTRITION REPORTING
cat("Step 3: Defining Sepsis cohort and reporting attrition...\n")
n_initial <- nrow(p_data)

final_df <- p_data %>%
  select(Sample = geo_accession, SRS_Endotype = `endotype_class:ch1`, 
         Age = `age:ch1`, Status = `mortality_event_28days:ch1`, Time = `time_to_event_28days:ch1`) %>%
  mutate(across(c(Age, Status, Time), ~as.numeric(as.character(.)))) %>%
  inner_join(exp_clean, by = "Sample") %>%
  # Restriction: Sepsis patients only (exclude Healthy and missing outcomes)
  filter(!is.na(Time), !is.na(SRS_Endotype), SRS_Endotype != "Healthy") %>%
  mutate(SRS_Endotype = factor(SRS_Endotype, levels = c("Mars1", "Mars2", "Mars3", "Mars4")))

n_final <- nrow(final_df)
cat(paste0("CONSORT: Initial N=", n_initial, " | Final Sepsis N=", n_final, "\n"))

# 5. GENERATING PANELS

# A. KM SURVIVAL SUITE (Clinical Impact)
create_km <- function(marker, title) {
  df_temp <- final_df
  df_temp$Grp <- ifelse(df_temp[[marker]] > median(df_temp[[marker]]), "High", "Low")
  fit <- survfit(Surv(Time, Status) ~ Grp, data = df_temp)
  ggsurvplot(fit, data = df_temp, pval = TRUE, palette = "npg", 
             title = title, ggtheme = theme_pubr())$plot
}
p_km1 <- create_km("LGALS1", "LGALS1 Survival")
p_km2 <- create_km("LGALS3", "LGALS3 Survival")
p_km3 <- create_km("LGALS9", "LGALS9 Survival")

# B. BOXPLOTS WITH PAIRWISE WILCOXON (Biological Discovery)
plot_df <- final_df %>%
  pivot_longer(cols = c(LGALS1, LGALS3, LGALS9), names_to = "Marker", values_to = "Expression")

stat.test <- plot_df %>%
  group_by(Marker) %>%
  wilcox_test(Expression ~ SRS_Endotype) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance() %>%
  add_y_position()

dummy_limits <- data.frame(
  Marker = c("LGALS1", "LGALS1", "LGALS3", "LGALS3", "LGALS9", "LGALS9"),
  Expression = c(0, 14, 0, 14, 0, 7), # Custom Y-scales as requested
  SRS_Endotype = factor("Mars1", levels = levels(final_df$SRS_Endotype))
)

p1 <- ggboxplot(plot_df, x = "SRS_Endotype", y = "Expression", fill = "SRS_Endotype", 
                palette = "npg", facet.by = "Marker", scales = "free_y", 
                add = "jitter", add.params = list(alpha = 0.2)) +
  geom_blank(data = dummy_limits) +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", hide.ns = TRUE) +
  stat_compare_means(method = "kruskal.test", label.y.npc = "top") +
  labs(title = "Galectin Expression by MARS Endotype", y = "Log2 Expression") +
  theme_pubr() + theme(legend.position = "none")

# C. PAIRED HAZARD RATIOS (Statistical Rigor)
cox_models <- list(
  "LGALS1 (Standard)"    = coxph(Surv(Time, Status) ~ LGALS1 + Age, data = final_df),
  "LGALS1 (Sensitivity)" = coxph(Surv(Time, Status) ~ LGALS1 + Age + SRS_Endotype, data = final_df),
  "LGALS3 (Standard)"    = coxph(Surv(Time, Status) ~ LGALS3 + Age, data = final_df),
  "LGALS3 (Sensitivity)" = coxph(Surv(Time, Status) ~ LGALS3 + Age + SRS_Endotype, data = final_df),
  "LGALS9 (Standard)"    = coxph(Surv(Time, Status) ~ LGALS9 + Age, data = final_df),
  "LGALS9 (Sensitivity)" = coxph(Surv(Time, Status) ~ LGALS9 + Age + SRS_Endotype, data = final_df)
)

cox_res <- map_dfr(cox_models, tidy, conf.int = TRUE, .id = "Model") %>%
  filter(term %in% c("LGALS1", "LGALS3", "LGALS9")) %>%
  mutate(MarkerGroup = str_extract(Model, "LGALS\\d+"),
         Model = factor(Model, levels = rev(names(cox_models))))

p2 <- ggplot(cox_res, aes(x = estimate, y = Model, color = MarkerGroup)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  geom_point(size = 4) +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), width = 0.3) +
  scale_color_manual(values = c("LGALS1" = "#E74C3C", "LGALS3" = "#2E86C1", "LGALS9" = "#27AE60")) +
  labs(title = "Hazard Ratios (Paired Sensitivity)", x = "Log Hazard Ratio (95% CI)") +
  theme_pubr() + theme(legend.position = "none")

# D. CORRELATION WITH P-VALUES
p3 <- ggcorrplot(cor(final_df[, c("LGALS1", "LGALS3", "LGALS9")]), type = "lower", lab = TRUE, 
                 p.mat = cor_pmat(final_df[, c("LGALS1", "LGALS3", "LGALS9")]), 
                 insig = "blank", title = "Correlation Matrix")

# 6. ASSEMBLY & EXPORT
final_layout <- (p_km1 | p_km2 | p_km3) / 
  (p1) / 
  (p2 | p3) +
  plot_annotation(title = "MSCA ANALYSIS: GALECTIN PROGNOSTIC SIGNATURE",
                  subtitle = "Clinical Outcomes -> Endotype Profile -> Multivariable Robustness")

ggsave("MSCA_Final_Clinical_First.pdf", final_layout, width = 16, height = 20)
writeLines(capture.output(sessionInfo()), "session_info.txt")
cat("\nSUCCESS: 'MSCA_Final_Clinical_First.pdf' and session info generated.\n")