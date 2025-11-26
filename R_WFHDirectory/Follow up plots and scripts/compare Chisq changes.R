#install.packages("dplyr")
#install.packages("tidyr")

library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(viridisLite)

#split tag
df_all <- df_all %>%
  separate(tag, into = c("source","type","metric"), sep = "_", remove = FALSE)

# Compute percentage difference per chromosome per group
df_all <- df_all %>%
  mutate(pct_diff = 100 * (obs - expected) / expected)

# Summarize duplicates by chromosome & source
df_summarized <- df_all %>%
  group_by(chromosome, source) %>%
  summarise(pct_diff = mean(pct_diff), .groups = "drop")  # mean handles multiple entries

# Pivot wider so each group/source is a numeric column
df_wide <- df_summarized %>%
  pivot_wider(names_from = source, values_from = pct_diff)

# Compute per-chromosome differences between groups
df_wide <- df_wide %>%
  mutate(
    diff_High_GEO     = Highstringency - GEO,
    diff_High_Encode  = Highstringency - Encode,
    diff_Low_GEO      = Lowstringency - GEO,
    diff_Low_Encode   = Lowstringency - Encode,
    diff_High_Low     = Highstringency - Lowstringency
  )

# Compute genome-wide correlations across all chromosomes
group_correlations <- df_wide %>%
  select(-chromosome) %>%
  cor(use = "pairwise.complete.obs")

# Optionally, add per-chromosome differences back to original df_all
df_all_final <- df_all %>%
  left_join(df_wide %>%
              select(chromosome, diff_High_GEO, diff_High_Encode, diff_Low_GEO, diff_Low_Encode, diff_High_Low),
            by = "chromosome")

#Save data
write.csv(df_all_final, "percentage_differences_and_deltas.csv", row.names = FALSE)

# ---- Plot genome-wide correlations ----
# ---- Genome-wide correlations (only original groups) ----
original_sources <- c("GEO", "Encode", "Highstringency", "Lowstringency")

group_correlations <- df_wide %>%
  select(all_of(original_sources)) %>%
  cor(use = "pairwise.complete.obs", method = "pearson")

# Convert correlation matrix to long format for plotting
cor_df <- as.data.frame(group_correlations) %>%
  tibble::rownames_to_column("Group1") %>%
  pivot_longer(-Group1, names_to = "Group2", values_to = "Correlation") %>%
  filter(Group1 < Group2)  # keep only unique pairs


cols24 <- viridis(6, option = "plasma")
cor_df$Comparison <- paste(cor_df$Group1, "vs", cor_df$Group2)

# Plot genome-wide correlations
tiff("PctChange_Pergroup.tiff", width = 2000, height = 2000, res = 300)

ggplot(cor_df, aes(x = interaction(Group1, Group2, sep = " vs "), 
                   y = Correlation, 
                   fill = Comparison)) +
  geom_col() +
  scale_fill_manual(values = cols24) +
  labs(title = "Genome-wide Correlations (Original Groups Only)",
       x = "Group Comparison",
       y = "Correlation (Pearson)") +
  theme_minimal() +
  ylim(0, 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

#Plot
tiff("PctChange_PerChromosome.tiff", width = 3508, height = 2480, res = 300) 
ggplot(df_all, aes(x = chromosome, y = pct_diff, fill = source)) +
  geom_col(position = "dodge") + 
  ylim(-100, 400) +
  labs(title = "Percentage Difference per Chromosome by Group",
       x = "Chromosome",
       y = "Percentage Difference (%)") +
  theme_minimal(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


# ---- Compute correlation significance ----
library(Hmisc)  # for rcorr
library(openxlsx)  # for writing Excel files

# Select original groups
data_for_corr <- df_wide %>% select(all_of(original_sources))

# Compute correlation matrix and p-values
corr_res <- rcorr(as.matrix(data_for_corr), type = "pearson")

# Extract correlation coefficients and p-values
cor_matrix <- corr_res$r
p_matrix <- corr_res$P

# Convert matrices to long format
cor_long <- as.data.frame(cor_matrix) %>%
  tibble::rownames_to_column("Group1") %>%
  pivot_longer(-Group1, names_to = "Group2", values_to = "Correlation") %>%
  filter(Group1 < Group2)

p_long <- as.data.frame(p_matrix) %>%
  tibble::rownames_to_column("Group1") %>%
  pivot_longer(-Group1, names_to = "Group2", values_to = "P_value") %>%
  filter(Group1 < Group2)

# Combine correlation and p-value
cor_results <- cor_long %>%
  left_join(p_long, by = c("Group1", "Group2"))

# Export to Excel
write.xlsx(cor_results, "GenomeWide_Correlations_with_Pvalues.xlsx", rowNames = FALSE)
