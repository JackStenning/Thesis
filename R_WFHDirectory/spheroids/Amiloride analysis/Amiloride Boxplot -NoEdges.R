# Load necessary libraries
library(ggplot2)
library(car)
library(multcompView)
library(dplyr)
library(readxl)
library(tidyverse)
library(ggpubr)

### Early Treatment
#Set working directory
setwd("C:/Users/jps558/OneDrive - University of York/Desktop/R_WorkingDirectory/R_WFHDirectory") 

#Path data
Path_E = "EarlyAmilorideNoedge.xlsx"

#Get Data
raw_data_E <- read_excel(Path_E, sheet = "Sheet1")

#Create order
desired_order <- c("Vehicle control", "No Amiloride", "1 nM Amiloride", "10 nM Amiloride", "100 nM Amiloride", "1 µM Amiloride", "10 µM Amiloride", "100 µM Amiloride")

# Reshape the data
long_data <- raw_data_E %>%
  pivot_longer(
    cols = starts_with("Day"),
    names_to = "Day",
    values_to = "Size"
  ) %>%
  mutate(
    Day = str_remove(Day, "Day"),                 # Remove "Day" prefix
    Day = factor(Day, levels = unique(Day)),      # Ensure day order is preserved
    Condition = as.factor(Condition),     # Factor for grouping
    Size = as.numeric(Size)                       # Ensure Size is numeric
  )

# Quick check: print structure
# str(long_data)

long_data <- long_data %>%
  mutate(
    Condition = factor(Condition, levels = desired_order)
  )

# Plot with cleaner layout
ggplot(long_data, aes(x = Day, y = Size, fill = Condition)) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    outlier.shape = NA, width = 0.6
  ) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
    alpha = 0.4, size = 1, shape = 21
  ) +
  labs(
    title = "(no edge) Spheroid Size Over Time by Amiloride Concentration after treatment on day 2",
    x = "Day",
    y = "Spheroid Size (µm)",
    fill = "Concentration"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 11),
    legend.position = "right",
    plot.title = element_text(size = 14, face = "bold")
  )

#########Plot Day 7-13

# Convert Day to numeric for filtering
long_data_filtered <- long_data %>%
  mutate(Day_num = as.numeric(as.character(Day))) %>%
  filter(Day_num >= 7 & Day_num <= 13)

ggplot(long_data_filtered, aes(x = Day, y = Size, fill = Condition)) +
  geom_boxplot(position = position_dodge(0.75), outlier.shape = NA, width = 0.6) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
              alpha = 0.4, size = 1, , shape = 21) +
  labs(
    title = "(no edge) Spheroid Size from Day 7–13 by Amiloride Concentration after treatment on day 2",
    x = "Day",
    y = "Spheroid Size (µm)",
    fill = "Concentration"
  ) +
  theme_minimal()

#ANOVA test
# Prepare empty list for sheets
tukey_sheets <- list()
all_significant <- data.frame()

# Loop through days 7 to 13
for (d in 7:13) {
  df_day <- long_data_filtered %>% filter(Day_num == d)
  if (nrow(df_day) < 2) next
  
  # Run ANOVA
  aov_result <- aov(Size ~ Condition, data = df_day)
  aov_summary <- tidy(aov_result)
  f_value <- aov_summary$statistic[1]
  aov_p <- aov_summary$p.value[1]
  
  # Run Tukey HSD
  tukey <- TukeyHSD(aov_result)
  tukey_df <- as.data.frame(tukey$Condition)
  tukey_df$comparison <- rownames(tukey_df)
  
  # Process into clean format
  tukey_df <- tukey_df %>%
    separate(comparison, into = c("group1", "group2"), sep = "-") %>%
    mutate(
      Day = d,
      F_value = round(f_value, 3),
      ANOVA_p = signif(aov_p, 3),
      p_label = case_when(
        `p adj` < 0.001 ~ "***",
        `p adj` < 0.01  ~ "**",
        `p adj` < 0.05  ~ "*",
        TRUE            ~ "ns"
      )
    ) %>%
    select(Day, group1, group2, diff, lwr, upr, `p adj`, p_label, F_value, ANOVA_p)
  
  # Add to list of sheets
  tukey_sheets[[paste0("Day_", d)]] <- tukey_df
  
  # Add significant only
  sig_df <- tukey_df %>% filter(`p adj` < 0.05)
  all_significant <- bind_rows(all_significant, sig_df)
}

# Add a sheet for only significant comparisons
tukey_sheets[["Significant_Only"]] <- all_significant

# Export to Excel file
write_xlsx(tukey_sheets, "early_noedges_tukey_anova_results_by_day.xlsx")

### Late Treatment

#Path Data 
Path_L = "LateAmilorideNoedge.xlsx"

#Get data
raw_data_L <- read_excel(Path_L, sheet = "Sheet1")

# Reshape the data
long_data_L <- raw_data_L %>%
  pivot_longer(
    cols = starts_with("Day"),
    names_to = "Day",
    values_to = "Size"
  ) %>%
  mutate(
    Day = str_remove(Day, "Day"),                 # Remove "Day" prefix
    Day = factor(Day, levels = unique(Day)),      # Ensure day order is preserved
    Condition = as.factor(Condition),     # Factor for grouping
    Size = as.numeric(Size)                       # Ensure Size is numeric
  )

# Quick check: print structure
# str(long_data)

long_data_L <- long_data_L %>%
  mutate(
    Condition = factor(Condition, levels = desired_order)
  )

# Plot with cleaner layout
ggplot(long_data_L, aes(x = Day, y = Size, fill = Condition)) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    outlier.shape = NA, width = 0.6
  ) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
    alpha = 0.4, size = 1, shape = 21
  ) +
  labs(
    title = "(no edge) Spheroid Size Over Time by Amiloride Concentration after treatment on day 9",
    x = "Day",
    y = "Spheroid Size (µm)",
    fill = "Concentration"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 11),
    legend.position = "right",
    plot.title = element_text(size = 14, face = "bold")
  )

#########Plot Day 7-13

# Convert Day to numeric for filtering
long_data_L_filtered <- long_data_L %>%
  mutate(Day_num = as.numeric(as.character(Day))) %>%
  filter(Day_num >= 7 & Day_num <= 13)

ggplot(long_data_L_filtered, aes(x = Day, y = Size, fill = Condition)) +
  geom_boxplot(position = position_dodge(0.75), outlier.shape = NA, width = 0.6) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
              alpha = 0.4, size = 1, , shape = 21) +
  labs(
    title = "(no edge) Spheroid Size from Day 7–13 by Amiloride Concentration after treatment on day 9",
    x = "Day",
    y = "Spheroid Size (µm)",
    fill = "Concentration"
  ) +
  theme_minimal()

#ANOVA test
# Prepare empty list for sheets
Ltukey_sheets <- list()
Lall_significant <- data.frame()

# Loop through days 7 to 13
for (d in 7:13) {
  Ldf_day <- long_data_L_filtered %>% filter(Day_num == d)
  if (nrow(Ldf_day) < 2) next
  
  # Run ANOVA
  Laov_result <- aov(Size ~ Condition, data = Ldf_day)
  Laov_summary <- tidy(Laov_result)
  Lf_value <- Laov_summary$statistic[1]
  Laov_p <- Laov_summary$p.value[1]
  
  # Run Tukey HSD
  Ltukey <- TukeyHSD(Laov_result)
  Ltukey_df <- as.data.frame(Ltukey$Condition)
  Ltukey_df$comparison <- rownames(Ltukey_df)
  
  # Process into clean format
  Ltukey_df <- Ltukey_df %>%
    separate(comparison, into = c("group1", "group2"), sep = "-") %>%
    mutate(
      Day = d,
      F_value = round(Lf_value, 3),
      ANOVA_p = signif(Laov_p, 3),
      p_label = case_when(
        `p adj` < 0.001 ~ "***",
        `p adj` < 0.01  ~ "**",
        `p adj` < 0.05  ~ "*",
        TRUE            ~ "ns"
      )
    ) %>%
    select(Day, group1, group2, diff, lwr, upr, `p adj`, p_label, F_value, ANOVA_p)
  
  # Add to list of sheets
  Ltukey_sheets[[paste0("Day_", d)]] <- Ltukey_df
  
  # Add significant only
  Lsig_df <- Ltukey_df %>% filter(`p adj` < 0.05)
  Lall_significant <- bind_rows(Lall_significant, Lsig_df)
}

# Add a sheet for only significant comparisons
Ltukey_sheets[["Significant_Only"]] <- Lall_significant

# Export to Excel file
write_xlsx(Ltukey_sheets, "late_noedge_tukey_anova_results_by_day.xlsx")

