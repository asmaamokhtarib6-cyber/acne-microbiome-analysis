## =============================================================================
## Sun et al. 2025: Skin Microbiome Analysis (Short vs. Long Acne Duration)
## Mirrors the methodology used for the Fadilah 26-sample age-comparison project
## Run this in Google Colab (R runtime). Upload level-6.csv and metadata
## exported from QIIME2's taxa-barplot.qzv before running.
## =============================================================================

## ---------------------------------------------------------------------------
## 1. SETUP
## ---------------------------------------------------------------------------
required_pkgs <- c("dplyr", "tidyr", "stringr", "ggplot2", "forcats",
                    "effectsize", "vegan", "pheatmap")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[,"Package"])]
if (length(new_pkgs) > 0) install.packages(new_pkgs, quiet = TRUE)

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(forcats)
library(effectsize)
library(vegan)      # adonis2 (PERMANOVA), diversity indices
library(pheatmap)

cat("R version:", R.version.string, "\n")
cat("All packages loaded successfully.\n")


## ---------------------------------------------------------------------------
## 2. DATA IMPORT
##
## Upload "level-6.csv" (genus-level relative/absolute counts, exported
## directly from QIIME2's taxa-barplot.qzv at Level 6) to the Colab session
## before running this cell.
## ---------------------------------------------------------------------------
raw <- read.csv("level-6.csv", check.names = FALSE)

# Identify metadata columns vs. taxon count columns
meta_cols <- c("index", "SampleNumber", "SubjectID", "Gender",
                "DurationGroup", "Location", "Group")
taxon_cols <- setdiff(colnames(raw), meta_cols)

cat("Samples:", nrow(raw), "\n")
cat("Taxon columns (genus-level, full lineage):", length(taxon_cols), "\n")

metadata <- raw %>% select(all_of(meta_cols)) %>% rename(sample_id = index)
counts   <- raw %>% select(all_of(c("index", taxon_cols))) %>% rename(sample_id = index)


## ---------------------------------------------------------------------------
## 3. CLEAN TAXON LABELS: extract a short genus name from the full lineage
## string (e.g. "d__Bacteria;p__...;g__Cutibacterium" -> "Cutibacterium").
## Where genus is unresolved ("g__" empty), fall back to the family name.
## ---------------------------------------------------------------------------
extract_genus <- function(lineage) {
  parts <- str_split(lineage, ";")[[1]]
  last  <- str_remove(parts[length(parts)], "^g__")
  if (last == "" || last == "__") {
    fam <- str_remove(parts[length(parts) - 1], "^f__")
    return(paste0(fam, " (unresolved genus)"))
  }
  return(last)
}
genus_labels <- sapply(taxon_cols, extract_genus)
names(genus_labels) <- taxon_cols


## ---------------------------------------------------------------------------
## 4. CONTAMINANT CHECK
## Same logic as the Fadilah project: flag any genus present at high
## abundance in effectively ALL samples that is implausible as a skin
## commensal (soil/water/reagent-associated genera are the classic
## "kitome" contamination signature).
## ---------------------------------------------------------------------------
rel_check <- counts %>%
  mutate(across(all_of(taxon_cols), ~ . / rowSums(across(all_of(taxon_cols))) * 100))

contam_candidates <- data.frame(
  taxon = taxon_cols,
  genus = genus_labels[taxon_cols],
  pct_samples_present = colSums(rel_check[taxon_cols] > 0) / nrow(rel_check) * 100,
  mean_rel_abund = colMeans(rel_check[taxon_cols])
) %>%
  filter(pct_samples_present > 90, mean_rel_abund > 1) %>%
  arrange(desc(mean_rel_abund))

cat("\n=== Taxa present in >90% of samples with >1% mean abundance ===\n")
print(contam_candidates)
# REVIEW THIS TABLE: Cutibacterium and Staphylococcus are expected skin
# commensals and should NOT be removed. Only remove a genus here if it is
# implausible as a skin organism (e.g. soil/water/reagent-associated, as
# Lysobacter was in the Fadilah dataset). No such genus was identified in
# this dataset at the time this script was written, so no removal was applied.


## ---------------------------------------------------------------------------
## 5. RELATIVE ABUNDANCE TABLE
## ---------------------------------------------------------------------------
rel_abund <- counts %>%
  mutate(across(all_of(taxon_cols), ~ . / rowSums(across(all_of(taxon_cols))) * 100)) %>%
  left_join(metadata, by = "sample_id")

# Long format for plotting / summaries, with clean genus names
rel_long <- rel_abund %>%
  pivot_longer(cols = all_of(taxon_cols), names_to = "taxon", values_to = "rel_pct") %>%
  mutate(genus = genus_labels[taxon])

top15_genera <- rel_long %>%
  group_by(genus) %>%
  summarise(mean_pct = mean(rel_pct)) %>%
  arrange(desc(mean_pct)) %>%
  slice_head(n = 15)

cat("\n=== Top 15 genera by mean relative abundance (all 70 samples) ===\n")
print(top15_genera)


## ---------------------------------------------------------------------------
## 6. PRIMARY COMPARISON: SHORT vs. LONG ACNE DURATION
## Healthy controls are excluded from this comparison. See Section 9 for a
## separate acne-vs-healthy exploratory check.
## ---------------------------------------------------------------------------
sl_data <- rel_abund %>% filter(DurationGroup %in% c("short", "long"))
cat("\nGroup sizes: short =", sum(sl_data$DurationGroup == "short"),
    " long =", sum(sl_data$DurationGroup == "long"), "\n")

top15_cols <- names(genus_labels)[genus_labels %in% top15_genera$genus]

diff_abund_results <- data.frame()
for (col in top15_cols) {
  short_vals <- sl_data[[col]][sl_data$DurationGroup == "short"]
  long_vals  <- sl_data[[col]][sl_data$DurationGroup == "long"]

  wt <- wilcox.test(short_vals, long_vals)
  eff <- rank_biserial(short_vals, long_vals)

  diff_abund_results <- rbind(diff_abund_results, data.frame(
    genus       = genus_labels[[col]],
    mean_short  = mean(short_vals),
    mean_long   = mean(long_vals),
    raw_p       = wt$p.value,
    effect_r    = eff$r_rank_biserial
  ))
}

diff_abund_results$padj_BH <- p.adjust(diff_abund_results$raw_p, method = "BH")
diff_abund_results <- diff_abund_results %>% arrange(raw_p)

cat("\n=== Genus-level differential abundance: short vs. long duration ===\n")
print(diff_abund_results, digits = 4)
write.csv(diff_abund_results, "sun2025_genus_diff_abundance.csv", row.names = FALSE)


## ---------------------------------------------------------------------------
## 7. ALPHA DIVERSITY (Shannon, Simpson, Observed richness)
## ---------------------------------------------------------------------------
count_matrix <- counts %>% select(all_of(taxon_cols)) %>% as.matrix()
rownames(count_matrix) <- counts$sample_id

alpha_div <- data.frame(
  sample_id = counts$sample_id,
  Shannon   = diversity(count_matrix, index = "shannon"),
  Simpson   = diversity(count_matrix, index = "simpson"),
  Observed  = specnumber(count_matrix)
) %>% left_join(metadata, by = "sample_id")

alpha_sl <- alpha_div %>% filter(DurationGroup %in% c("short", "long"))

cat("\n=== Alpha diversity: short vs. long duration ===\n")
for (idx in c("Shannon", "Simpson", "Observed")) {
  s <- alpha_sl[[idx]][alpha_sl$DurationGroup == "short"]
  l <- alpha_sl[[idx]][alpha_sl$DurationGroup == "long"]
  wt <- wilcox.test(s, l)
  cat(sprintf("%-10s mean_short=%.3f  mean_long=%.3f  p=%.4f\n",
              idx, mean(s), mean(l), wt$p.value))
}

# Boxplot
ggplot(alpha_sl, aes(x = DurationGroup, y = Shannon, fill = DurationGroup)) +
  geom_boxplot() +
  labs(title = "Shannon Diversity: Short vs. Long Acne Duration",
       x = "Duration Group", y = "Shannon Index") +
  theme_minimal()
ggsave("sun2025_alpha_diversity_boxplot.png", width = 6, height = 5)


## ---------------------------------------------------------------------------
## 8. BETA DIVERSITY: PERMANOVA (Bray-Curtis) + PCoA
## ---------------------------------------------------------------------------
sl_ids <- sl_data$sample_id
sl_matrix <- count_matrix[sl_ids, ]
bc_dist <- vegdist(sl_matrix, method = "bray")

permanova_result <- adonis2(bc_dist ~ DurationGroup, data = sl_data, permutations = 999)
cat("\n=== PERMANOVA (Bray-Curtis, short vs. long, 999 permutations) ===\n")
print(permanova_result)

# PCoA ordination
pcoa <- cmdscale(bc_dist, k = 2, eig = TRUE)
pcoa_df <- data.frame(PC1 = pcoa$points[,1], PC2 = pcoa$points[,2],
                        sample_id = sl_ids) %>%
  left_join(metadata, by = "sample_id")

var_explained <- round(100 * pcoa$eig / sum(pcoa$eig[pcoa$eig > 0]), 1)

ggplot(pcoa_df, aes(x = PC1, y = PC2, color = DurationGroup)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.68) +
  labs(title = "PCoA (Bray-Curtis): Short vs. Long Acne Duration",
       x = paste0("Axis 1 (", var_explained[1], "%)"),
       y = paste0("Axis 2 (", var_explained[2], "%)")) +
  theme_minimal()
ggsave("sun2025_pcoa_plot.png", width = 6, height = 5)


## ---------------------------------------------------------------------------
## 9. EXPLORATORY CHECK 1: LOCATION (skin pore "in" vs. skin surface "out")
## Not the primary hypothesis, but it's the paper's own primary axis, so
## worth checking for completeness with the same BH-correction discipline.
## ---------------------------------------------------------------------------
loc_data <- rel_abund  # all 70 samples (includes health_out as "out")
loc_results <- data.frame()
for (col in top15_cols) {
  in_vals  <- loc_data[[col]][loc_data$Location == "in"]
  out_vals <- loc_data[[col]][loc_data$Location == "out"]
  wt <- wilcox.test(in_vals, out_vals)
  loc_results <- rbind(loc_results, data.frame(
    genus = genus_labels[[col]],
    mean_in = mean(in_vals), mean_out = mean(out_vals),
    raw_p = wt$p.value
  ))
}
loc_results$padj_BH <- p.adjust(loc_results$raw_p, method = "BH")
loc_results <- loc_results %>% arrange(raw_p)
cat("\n=== Exploratory: genus-level differences, skin pore (in) vs. surface (out) ===\n")
print(loc_results, digits = 4)
write.csv(loc_results, "sun2025_genus_diff_location.csv", row.names = FALSE)


## ---------------------------------------------------------------------------
## 10. EXPLORATORY CHECK 2: ACNE (short_out + long_out) vs. HEALTHY (health_out)
## Compares surface samples only, to isolate disease status from location.
## ---------------------------------------------------------------------------
acne_health <- rel_abund %>%
  filter(Location == "out") %>%
  mutate(status = ifelse(DurationGroup == "health", "healthy", "acne"))

health_results <- data.frame()
for (col in top15_cols) {
  acne_vals   <- acne_health[[col]][acne_health$status == "acne"]
  health_vals <- acne_health[[col]][acne_health$status == "healthy"]
  wt <- wilcox.test(acne_vals, health_vals)
  health_results <- rbind(health_results, data.frame(
    genus = genus_labels[[col]],
    mean_acne = mean(acne_vals), mean_healthy = mean(health_vals),
    raw_p = wt$p.value
  ))
}
health_results$padj_BH <- p.adjust(health_results$raw_p, method = "BH")
health_results <- health_results %>% arrange(raw_p)
cat("\n=== Exploratory: genus-level differences, acne (surface) vs. healthy (surface) ===\n")
print(health_results, digits = 4)
write.csv(health_results, "sun2025_genus_diff_acne_vs_healthy.csv", row.names = FALSE)


## ---------------------------------------------------------------------------
## 11. HIERARCHICAL CLUSTERING HEATMAP (top 12 genera, all 70 samples)
## ---------------------------------------------------------------------------
top12_cols <- top15_cols[1:12]
heatmap_matrix <- rel_abund %>%
  select(sample_id, all_of(top12_cols)) %>%
  column_to_rownames("sample_id") %>%
  as.matrix()
colnames(heatmap_matrix) <- genus_labels[top12_cols]

annotation_row <- metadata %>%
  select(sample_id, Group) %>%
  column_to_rownames("sample_id")

pheatmap(
  log10(heatmap_matrix + 0.01),
  annotation_row = annotation_row,
  clustering_distance_rows = "correlation",
  clustering_method = "average",
  filename = "sun2025_heatmap.png",
  width = 8, height = 12
)


## ---------------------------------------------------------------------------
## 12. MASTER SUMMARY TABLE
## ---------------------------------------------------------------------------
summary_table <- data.frame(
  Analysis = c("Genus-level (Cutibacterium)", "Alpha diversity (Shannon)",
               "Alpha diversity (Simpson)", "Beta diversity (PERMANOVA)"),
  Comparison = rep("Short vs. Long duration", 4),
  Key_Statistic = c(
    paste0("r = ", round(diff_abund_results$effect_r[diff_abund_results$genus == "Cutibacterium"], 3)),
    "Wilcoxon", "Wilcoxon",
    paste0("R2 = ", round(permanova_result$R2[1], 3))
  ),
  p_value = c(
    round(diff_abund_results$padj_BH[diff_abund_results$genus == "Cutibacterium"], 4),
    round(wilcox.test(alpha_sl$Shannon[alpha_sl$DurationGroup=="short"],
                        alpha_sl$Shannon[alpha_sl$DurationGroup=="long"])$p.value, 4),
    round(wilcox.test(alpha_sl$Simpson[alpha_sl$DurationGroup=="short"],
                        alpha_sl$Simpson[alpha_sl$DurationGroup=="long"])$p.value, 4),
    round(permanova_result$`Pr(>F)`[1], 4)
  )
)
cat("\n=== MASTER SUMMARY TABLE ===\n")
print(summary_table)
write.csv(summary_table, "sun2025_master_summary.csv", row.names = FALSE)

cat("\nAll outputs saved: sun2025_genus_diff_abundance.csv, ",
    "sun2025_alpha_diversity_boxplot.png, sun2025_pcoa_plot.png, ",
    "sun2025_genus_diff_location.csv, sun2025_genus_diff_acne_vs_healthy.csv, ",
    "sun2025_heatmap.png, sun2025_master_summary.csv\n", sep="")
