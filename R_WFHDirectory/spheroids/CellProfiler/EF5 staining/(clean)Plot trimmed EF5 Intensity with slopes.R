# Load necessary packages
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(broom)
library(readxl)

# Set working directory
setwd("C:/Users/jps558/OneDrive - University of York/Desktop/R_WorkingDirectory/R_WFHDirectory/spheroids/CellProfiler/EF5 staining") 

# Read the Excel file (specify the sheet if needed)
df <- read_excel("EF5 Intensity from core_NoAmilorideSpheroidMaskAmilorideremoved.xlsx", 
                 sheet = 1)  # or use the sheet name instead of 1

# Select relevant columns
df_selected <- df %>%
  select(
    FileName_Brightfield,
    matches("^RadialDistribution_(FracAtD|MeanFrac)_Fluorescence_[1-6]of6$")
  )

# Pivot to long format
df_long <- df_selected %>%
  pivot_longer(
    cols = -FileName_Brightfield,
    names_to = c("Metric", "Bin"),
    names_pattern = "RadialDistribution_(.*)_Fluorescence_([1-6])of6",
    values_to = "Value"
  ) %>%
  mutate(Bin = as.integer(Bin))

# Extract sphere & treatment info
df_long <- df_long %>%
  mutate(
    SphereRaw = str_extract(FileName_Brightfield, "Sphere\\d+(_\\d+)?"),
    SphereNum = ifelse(str_detect(SphereRaw, "_"),
                       str_extract(SphereRaw, "(?<=_)\\d+"), "1"),
    ImageExportNum = str_extract(FileName_Brightfield, "Image Export-\\d+"),
    SphereID = paste0(SphereRaw, "_", ImageExportNum),
    
    IsDMSO = str_starts(FileName_Brightfield, "DMSO_"),
    PrimaryCompound = ifelse(IsDMSO, "DMSO", "EF5"),
    SecondaryCompound = ifelse(IsDMSO, "EF5", NA_character_),
    TreatmentGroup = ifelse(is.na(SecondaryCompound), PrimaryCompound,
                            paste(PrimaryCompound, SecondaryCompound, sep = "+")),
    CleanTreatmentGroup = TreatmentGroup,
    
    EF5DayMatch = str_extract(FileName_Brightfield, "EF5Day\\d+(-\\d+)?"),
    DayNum = as.integer(str_extract(EF5DayMatch, "(?<=EF5Day)\\d+")),
    Day = factor(paste0("Day ", DayNum), levels = c("Day 7", "Day 9", "Day 11", "Day 13"))
  )

# Assign a simple numeric ID to each sphere per metric/day
df_long <- df_long %>%
  group_by(Metric, Day) %>%
  mutate(SphereSimpleID = factor(match(SphereID, unique(SphereID)))) %>%
  ungroup()

# -----------------------------
# Parameters for saving plots
# -----------------------------
plot_width <- 6      # inches
plot_height <- 5     # inches
plot_dpi <- 300      # resolution

metrics <- c("FracAtD", "MeanFrac")
days <- c(7, 9, 11, 13)

# -----------------------------
# Individual-day plots
# -----------------------------
for (metric in metrics) {
  for (day in days) {
    p <- ggplot(df_long %>% filter(Metric == metric, DayNum == day),
                aes(x = Bin, y = Value, group = SphereSimpleID, color = SphereSimpleID)) +
      geom_line() +
      labs(title = paste("EF5 -", metric, "- Day", day),
           x = "Bin (Distance from Center)",
           y = metric) +
      theme_minimal() +
      theme()
    
    ggsave(filename = paste0("EF5_", metric, "_Day", day, ".tiff"),
           plot = p,
           width = plot_width,
           height = plot_height,
           dpi = plot_dpi)
  }
}

# -----------------------------
# Faceted plots across days
# -----------------------------
for (metric in metrics) {
  p <- ggplot(df_long %>% filter(Metric == metric),
              aes(x = Bin, y = Value, group = SphereSimpleID, color = as.factor(DayNum))) +
    geom_line(alpha = 0.7) +
    facet_wrap(~Day) +
    labs(title = paste("EF5 -", metric, "across Days 7, 9, 11, 13"),
         x = "Bin (Distance from Center)",
         y = metric,
         color = "Day") +
    theme_minimal() +
    theme(legend.position = "none")
  
  ggsave(filename = paste0("EF5_", metric, "_Faceted.tiff"),
         plot = p,
         width = plot_width * 1.5,
         height = plot_height,
         dpi = plot_dpi)
}
# -----------------------------
# Faceted plots across EF5 days + Vehicle with trend lines
# -----------------------------
for (metric in metrics) {
  p <- ggplot(df_long %>% filter(Metric == metric),
              aes(x = Bin, y = Value, group = SphereSimpleID, color = Day)) +
    geom_line(alpha = 0.6) +
    # add trend line per facet (linear regression fit)
    geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black", linetype = "dashed") +
    facet_wrap(~Day) +
    labs(title = paste("EF5 + Vehicle -", metric, "across Days (with trend lines)"),
         x = "Bin (Distance from Center)",
         y = metric,
         color = "Group") +
    theme_minimal()
  
  ggsave(filename = paste0("EF5_Vehicle_", metric, "_Faceted_withTrend.tiff"),
         plot = p,
         width = plot_width * 1.5,
         height = plot_height,
         dpi = plot_dpi)
}

# -----------------------------
# Slope calculation per sphere
# -----------------------------
slopes_df <- df_long %>%
  group_by(Metric, Day, SphereSimpleID) %>%
  do(tidy(lm(Value ~ Bin, data = .))) %>%
  filter(term == "Bin") %>%
  rename(Slope = estimate,
         Slope_SE = std.error,
         Slope_t = statistic,
         Slope_p = p.value) %>%
  ungroup()

# Summarize slopes across spheres
slopes_summary <- slopes_df %>%
  group_by(Metric, Day) %>%
  summarise(
    Mean_Slope = mean(Slope),
    SD_Slope = sd(Slope),
    N = n(),
    t_test_p = t.test(Slope)$p.value
  )

# Add significance labels
slopes_summary <- slopes_summary %>%
  mutate(SignifLabel = ifelse(t_test_p < 0.001, "***",
                              ifelse(t_test_p < 0.01, "**",
                                     ifelse(t_test_p < 0.05, "*", "ns"))))

# -----------------------------
# Slope boxplot with significance
# -----------------------------
p_slopes <- ggplot(slopes_df, aes(x = Day, y = Slope, fill = Metric)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, position = position_dodge(width = 0.8)) +
  geom_jitter(aes(color = Metric),
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(data = slopes_summary,
            aes(x = Day, y = max(slopes_df$Slope) + 0.05, label = SignifLabel, group = Metric),
            position = position_dodge(width = 0.8),
            vjust = 0) +
  labs(title = "Slope of EF5 Metric across Bins per Sphere",
       y = "Slope (Value ~ Bin)",
       x = "Day") +
  theme_minimal() +
  theme(legend.position = "top")

# Save as TIFF
ggsave(filename = "EF5_Slopes_Boxplot_withPoints.tiff",
       plot = p_slopes,
       width = 9,
       height = 5,
       dpi = 300)

