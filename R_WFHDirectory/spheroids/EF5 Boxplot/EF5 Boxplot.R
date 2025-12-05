# Load necessary libraries
library(ggplot2)
library(car)
library(multcompView)
library(dplyr)
library(readxl)
library(tidyverse)
library(ggpubr)
library(tidyr)
library(readr)
library(writexl)
library(broom)


### Early Treatment
#Set working directory
setwd("C:/Users/jps558/OneDrive - University of York/Desktop/R_WorkingDirectory/R_WFHDirectory") 

#Path data
Path = "EF5 Sizes.xlsx"

#Get Data
raw_data <- read_excel(Path, sheet = "Sheet1")

# Clean and prepare
long_data <- raw_data %>%
  mutate(
    Day = factor(Day, levels = sort(unique(Day))),  # Ensure Day is a factor for ordering
    Condition = as.factor(Condition),
    Size = as.numeric(Size)                         # Ensure Size is numeric
  )

# Plot with cleaner layout
tiff("Average Spheroid Size Over Time boxplot (EF5 Conditions).tiff", width = 2481, height = 1749, res = 300)
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
    title = "Spheroid Size Over Time before treatment with EF5",
    x = "Day",
    y = "Spheroid Size (µm)",
    fill = "Condition"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 11),
    legend.position = "right",
    plot.title = element_text(size = 14, face = "bold")
  )
dev.off()

#Line plot
# Summary stats by Day and Condition
summary_data <- long_data %>%
  group_by(Day, Condition) %>%
  summarise(
    mean_size = mean(Size, na.rm = TRUE),
    sem = sd(Size, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Line plot with error bars
tiff("Average Spheroid Size Over Time (EF5 Conditions).tiff", width = 2481, height = 1749, res = 300)
ggplot(summary_data, aes(x = Day, y = mean_size, color = Condition, group = Condition)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_size - sem, ymax = mean_size + sem), width = 0.3) +
  labs(
    title = "Average Spheroid Size Over Time (EF5 Conditions)",
    x = "Day",
    y = "Mean Spheroid Size (µm)",
    color = "Condition"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  )
dev.off()

# Prepare output containers
tukey_results <- list()
all_significant <- data.frame()

# Loop over unique days in the dataset
for (d in sort(unique(long_data$Day))) {
  df_day <- long_data %>% filter(Day == d)
  if (n_distinct(df_day$Condition) < 2) next  # Skip if only one condition
  
  # ANOVA
  aov_result <- aov(Size ~ Condition, data = df_day)
  tidy_aov <- tidy(aov_result)
  
  f_value <- tidy_aov$statistic[1]
  aov_p <- tidy_aov$p.value[1]
  
  # Tukey HSD
  tukey <- TukeyHSD(aov_result)
  tukey_df <- as.data.frame(tukey$Condition)
  tukey_df$comparison <- rownames(tukey_df)
  
  # Clean and annotate
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
  
  # Save each day's full results
  tukey_results[[paste0("Day_", d)]] <- tukey_df
  
  # Save only significant comparisons
  all_significant <- bind_rows(all_significant, tukey_df %>% filter(`p adj` < 0.05))
}

# Add sheet with all significant results
tukey_results[["Significant_Only"]] <- all_significant

# Export to Excel
write_xlsx(tukey_results, "ef5_tukey_anova_results_by_day.xlsx")
