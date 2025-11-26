# Load necessary libraries
library(ggplot2)
library(car)
library(multcompView)
library(dplyr)
library(readxl)
library(tidyverse)


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

# Filter for a specific day
e_day7_data <- long_data %>% filter(Day == " 7")
e_day8_data <- long_data %>% filter(Day == " 8")
e_day9_data <- long_data %>% filter(Day == " 9")
e_day10_data <- long_data %>% filter(Day == " 10")
e_day11_data <- long_data %>% filter(Day == " 11")
e_day12_data <- long_data %>% filter(Day == " 12")
e_day13_data <- long_data %>% filter(Day == " 13")

# Run one-way ANOVA
e_anova_day7 <- aov(Size ~ Condition, data = e_day7_data)
e_anova_day8 <- aov(Size ~ Condition, data = e_day8_data)
e_anova_day9 <- aov(Size ~ Condition, data = e_day9_data)
e_anova_day10 <- aov(Size ~ Condition, data = e_day10_data)
e_anova_day11 <- aov(Size ~ Condition, data = e_day11_data)
e_anova_day12 <- aov(Size ~ Condition, data = e_day12_data)
e_anova_day13 <- aov(Size ~ Condition, data = e_day13_data)

# View summary
summary(e_anova_day7)
summary(e_anova_day8)
summary(e_anova_day9)
summary(e_anova_day10)
summary(e_anova_day11)
summary(e_anova_day12)
summary(e_anova_day13)

#Levene's test - make sure these values are more than 0.05 to use anova
leveneTest(Size ~ Condition, data = e_day7_data)
leveneTest(Size ~ Condition, data = e_day8_data)
leveneTest(Size ~ Condition, data = e_day9_data)
leveneTest(Size ~ Condition, data = e_day10_data)
leveneTest(Size ~ Condition, data = e_day11_data)
leveneTest(Size ~ Condition, data = e_day12_data)
leveneTest(Size ~ Condition, data = e_day13_data)

#Tukey test

ed7_tukey_result <- TukeyHSD(e_anova_day7)
ed8_tukey_result <- TukeyHSD(e_anova_day8)
ed9_tukey_result <- TukeyHSD(e_anova_day9)
ed10_tukey_result <- TukeyHSD(e_anova_day10)
ed11_tukey_result <- TukeyHSD(e_anova_day11)
ed12_tukey_result <- TukeyHSD(e_anova_day12)
ed13_tukey_result <- TukeyHSD(e_anova_day13)

#Compact letter display
ed7_tukey_letters <- multcompLetters4(e_anova_day7, ed7_tukey_result)
ed8_tukey_letters <- multcompLetters4(e_anova_day8, ed8_tukey_result)
ed9_tukey_letters <- multcompLetters4(e_anova_day9, ed9_tukey_result)
ed10_tukey_letters <- multcompLetters4(e_anova_day10, ed10_tukey_result)
ed11_tukey_letters <- multcompLetters4(e_anova_day11, ed11_tukey_result)
ed12_tukey_letters <- multcompLetters4(e_anova_day12, ed12_tukey_result)
ed13_tukey_letters <- multcompLetters4(e_anova_day13, ed13_tukey_result)

ed7_letter_df <- as.data.frame.list(ed7_tukey_letters$Condition)
ed7_letter_df$Condition <- rownames(ed7_letter_df)
colnames(ed7_letter_df)[1] <- "letters"

ed8_letter_df <- as.data.frame.list(ed8_tukey_letters$Condition)
ed8_letter_df$Condition <- rownames(ed8_letter_df)
colnames(ed8_letter_df)[1] <- "letters"

ed9_letter_df <- as.data.frame.list(ed9_tukey_letters$Condition)
ed9_letter_df$Condition <- rownames(ed9_letter_df)
colnames(ed9_letter_df)[1] <- "letters"

ed10_letter_df <- as.data.frame.list(ed10_tukey_letters$Condition)
ed10_letter_df$Condition <- rownames(ed10_letter_df)
colnames(ed10_letter_df)[1] <- "letters"

ed11_letter_df <- as.data.frame.list(ed11_tukey_letters$Condition)
ed11_letter_df$Condition <- rownames(ed11_letter_df)
colnames(ed11_letter_df)[1] <- "letters"

ed12_letter_df <- as.data.frame.list(ed12_tukey_letters$Condition)
ed12_letter_df$Condition <- rownames(ed12_letter_df)
colnames(ed12_letter_df)[1] <- "letters"

ed13_letter_df <- as.data.frame.list(ed13_tukey_letters$Condition)
ed13_letter_df$Condition <- rownames(ed13_letter_df)
colnames(ed13_letter_df)[1] <- "letters"

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
