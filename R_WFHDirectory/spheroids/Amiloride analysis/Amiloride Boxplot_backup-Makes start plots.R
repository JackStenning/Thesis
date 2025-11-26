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

### Early Treatment
#Set working directory
setwd("//userfs/jps558/w2k/Desktop/R_WorkingDirectory") 

#Path data
Path_E = "EarlyAmiloride.xlsx"

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
    title = "Spheroid Size Over Time by Amiloride Concentration after treatment on day 2",
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
    title = "Spheroid Size from Day 7–13 by Amiloride Concentration after treatment on day 2",
    x = "Day",
    y = "Spheroid Size (µm)",
    fill = "Concentration"
  ) +
  theme_minimal()

#ANOVA test
# Create an empty list to collect Tukey results
tukey_list <- list()

# Initialize cumulative results dataframe
all_stats_df <- data.frame()

# Loop through each day and run ANOVA + Tukey
for (d in c(" 7", " 8", " 9", " 10", " 11", " 12", " 13")) {
  
  df_day <- long_data %>% filter(Day == d)
  
  if (nrow(df_day) < 2) next  # skip empty or invalid
  
  aov_result <- aov(Size ~ Condition, data = df_day)
  tukey <- TukeyHSD(aov_result)
  
  # Convert to data frame
  tukey_df <- as.data.frame(tukey$Condition)
  tukey_df$comparison <- rownames(tukey_df)
  
  # Split into groups
  tukey_df <- tukey_df %>%
    filter(`p adj` < 0.05) %>%
    separate(comparison, into = c("group1", "group2"), sep = "-") %>%
    mutate(
      Day = d,
      y.position = max(df_day$Size, na.rm = TRUE) + seq(5, 5 * n(), 5),
      p.label = case_when(
        `p adj` < 0.001 ~ "***",
        `p adj` < 0.01  ~ "**",
        `p adj` < 0.05  ~ "*"
      )
    )
  
  tukey_list[[d]] <- tukey_df
}

# Combine all Tukey results
stat_annotations <- bind_rows(tukey_list)

# Prepare annotations for faceted plot
stat_annotations_facet <- stat_annotations %>%
  mutate(
    xmin = factor(group1, levels = condition_levels),
    xmax = factor(group2, levels = condition_levels),
    x = group1,  # Any valid placeholder, will be overridden
    Day = as.factor(Day)  # for facetting
  )

# Ensure Day in filtered data is a factor too
long_data_filtered <- long_data_filtered %>%
  mutate(Day = as.factor(Day),
         Condition = factor(Condition, levels = condition_levels))

# Create the faceted boxplot
ggplot(long_data_filtered, aes(x = Condition, y = Size, fill = Condition)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(position = position_jitter(width = 0.1), alpha = 0.4, size = 1, shape = 21) +
  facet_wrap(~ Day, nrow = 2) +
  stat_pvalue_manual(
    stat_annotations_facet,
    label = "p.label",
    y.position = "y.position",
    xmin = "xmin",
    xmax = "xmax",
    inherit.aes = FALSE  # 🔑 disables ggplot's inherited aesthetics
  ) +
  labs(
    title = "Spheroid Size from Day 7–13 by Amiloride Concentration (Faceted by Day)",
    x = "Amiloride Concentration",
    y = "Spheroid Size (µm)",
    fill = "Concentration"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )


### Late Treatment

#Path Data 
Path_L = "LateAmiloride.xlsx"

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
    title = "Spheroid Size Over Time by Amiloride Concentration after treatment on day 9",
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
    title = "Spheroid Size from Day 7–13 by Amiloride Concentration after treatment on day 9",
    x = "Day",
    y = "Spheroid Size (µm)",
    fill = "Concentration"
  ) +
  theme_minimal()

#ANOVA
#ANOVA test
# Create an empty list to collect Tukey results
Ltukey_list <- list()

# Loop through each day and run ANOVA + Tukey
for (d in c(" 7", " 8", " 9", " 10", " 11", " 12", " 13")) {
  
  Ldf_day <- long_data_L %>% filter(Day == d)
  
  if (nrow(Ldf_day) < 2) next  # skip empty or invalid
  
  Laov_result <- aov(Size ~ Condition, data = Ldf_day)
  Ltukey <- TukeyHSD(Laov_result)
  
  # Convert to data frame
  Ltukey_df <- as.data.frame(Ltukey$Condition)
  Ltukey_df$comparison <- rownames(Ltukey_df)
  
  # Split into groups
  Ltukey_df <- Ltukey_df %>%
    filter(`p adj` < 0.05) %>%
    separate(comparison, into = c("group1", "group2"), sep = "-") %>%
    mutate(
      Day = d,
      y.position = max(Ldf_day$Size, na.rm = TRUE) + seq(5, 5 * n(), 5),
      p.label = case_when(
        `p adj` < 0.001 ~ "***",
        `p adj` < 0.01  ~ "**",
        `p adj` < 0.05  ~ "*"
      )
    )
  
  Ltukey_list[[d]] <- Ltukey_df
}

# Combine all Tukey results
Lstat_annotations <- bind_rows(Ltukey_list)
