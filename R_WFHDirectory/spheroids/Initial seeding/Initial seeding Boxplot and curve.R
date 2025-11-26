# Load Libraries
library(ggplot2)
library(dplyr)
library(readxl)
library(tidyr)
library(ggpubr)
library(stringr)

# Set working directory
setwd("C:/Users/jps558/OneDrive - University of York/Desktop/R_WorkingDirectory/R_WFHDirectory/spheroids/Initial seeding") 

# Path to Excel file
Path <- "Initial Spheroid sizes full data.xlsx"

# Sheets of interest (Excel names)
sheets <- c("500 cellswell", "1000 cellswell", "5000 cellswell", "10000 cellswell")

# Create pretty labels (for plots only)
labels <- str_replace(sheets, "cellswell", "cells/well")

# Read and combine data safely
all_data <- lapply(sheets, function(s) {
  df <- read_excel(Path, sheet = s)
  
  # Convert all Day columns to numeric
  day_cols <- grep("^Day", names(df), value = TRUE)
  df[day_cols] <- lapply(df[day_cols], function(x) as.numeric(as.character(x)))
  
  df$Condition <- s   # store sheet name as condition
  return(df)
}) %>% bind_rows()

# Reshape to long format
long_data <- all_data %>%
  pivot_longer(
    cols = starts_with("Day"),
    names_to = "Day",
    values_to = "Size"
  ) %>%
  mutate(
    Day = str_remove(Day, "Day "),     
    Day = str_remove(Day, " size"),    
    Day = as.numeric(Day),
    Size = as.numeric(Size),
    Condition = factor(Condition, levels = sheets, labels = labels)  # pretty labels
  ) %>%
  filter(Day != 1)   # remove Day 1 (empty)

# -------- BOX PLOT --------
tiff("Spheroid Size Over Time by Seeding Density Boxplot.tiff", width = 2480, height = 1748, res = 300) 
ggplot(long_data, aes(x = factor(Day), y = Size, fill = Condition)) +
  geom_boxplot(position = position_dodge(width = 0.75), outlier.shape = NA, width = 0.6) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
              alpha = 0.4, size = 1, shape = 21) +
  labs(
    title = "Spheroid Diameter Over Time by Seeding Density",
    x = "Day",
    y = "Spheroid Diameter (µm)",
    fill = "Condition"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    plot.title = element_text(size = 15, face = "bold"),
    legend.position = "right"
  )
dev.off()

# -------- LINE PLOT --------
# Calculate summary stats
summary_data <- long_data %>%
  group_by(Day, Condition) %>%
  summarise(
    mean_size = mean(Size, na.rm = TRUE),
    sem = sd(Size, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

tiff("Mean Spheroid Growth Over Time by Seeding Density Curve.tiff", width = 2480, height = 1748, res = 300) 
ggplot(summary_data, aes(x = Day, y = mean_size, color = Condition, group = Condition)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_size - sem, ymax = mean_size + sem), width = 0.3) +
  labs(
    title = "Mean Spheroid Growth Over Time by Seeding Density",
    x = "Day",
    y = "Mean Spheroid Diameter (µm)",
    color = "Condition"
  ) +
  theme_minimal() +
  scale_y_continuous(limits = c(200, 1000)) +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    plot.title = element_text(size = 15, face = "bold")
  )
dev.off()
