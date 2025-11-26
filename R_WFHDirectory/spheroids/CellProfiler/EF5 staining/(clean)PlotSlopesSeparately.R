# ================================================================
# Spheroid Fluorescence Intensity Analysis
# ================================================================

# -----------------------------
# Load packages
# -----------------------------
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(broom)
library(readxl)
library(purrr)

# -----------------------------
# Set working directory
# -----------------------------
setwd("C:/Users/jps558/OneDrive - University of York/Desktop/R_WorkingDirectory/R_WFHDirectory/spheroids/CellProfiler/EF5 staining") 

# -----------------------------
# Read data
# -----------------------------
df <- read_excel("EF5 Intensity from core_NoAmilorideSpheroidMaskAmilorideremoved.xlsx", sheet = 1)

# -----------------------------
# Select relevant columns
# -----------------------------
df_selected <- df %>%
  select(
    FileName_Brightfield,
    matches("^RadialDistribution_(FracAtD|MeanFrac)_Fluorescence_[1-6]of6$")
  )

# -----------------------------
# Pivot to long format
# -----------------------------
df_long <- df_selected %>%
  pivot_longer(
    cols = -FileName_Brightfield,
    names_to = c("Metric", "Bin"),
    names_pattern = "RadialDistribution_(.*)_Fluorescence_([1-6])of6",
    values_to = "Value"
  ) %>%
  mutate(Bin = as.integer(Bin))

# -----------------------------
# Extract sphere & treatment info
# -----------------------------
df_long <- df_long %>%
  mutate(
    SphereRaw = str_extract(FileName_Brightfield, "Sphere\\d+(_\\d+)?"),
    SphereNum = ifelse(str_detect(SphereRaw, "_"),
                       str_extract(SphereRaw, "(?<=_)\\d+"), "1"),
    ImageExportNum = str_extract(FileName_Brightfield, "Image Export-\\d+"),
    SphereID = paste0(SphereRaw, "_", ImageExportNum),
    
    # Identify Vehicle, DMSO, EF5
    IsVehicle = str_detect(FileName_Brightfield, "Vehicle"),
    IsDMSO = str_detect(FileName_Brightfield, "DMSO"),
    TreatmentGroup = case_when(
      IsVehicle ~ "Vehicle",
      IsDMSO ~ "DMSO",
      TRUE ~ "EF5"
    ),
    
    EF5DayMatch = str_extract(FileName_Brightfield, "EF5Day\\d+(-\\d+)?"),
    DayNum = as.integer(str_extract(EF5DayMatch, "(?<=EF5Day)\\d+")),
    Day = factor(paste0("Day ", DayNum), levels = c("Day 7", "Day 9", "Day 11", "Day 13"))
  )

# -----------------------------
# Assign numeric ID per sphere
# -----------------------------
df_long <- df_long %>%
  group_by(Metric, Day) %>%
  mutate(SphereSimpleID = factor(match(SphereID, unique(SphereID)))) %>%
  ungroup()

# -----------------------------
# Slope calculation per sphere
# -----------------------------
slopes_df <- df_long %>%
  group_by(Metric, Day, SphereSimpleID, TreatmentGroup) %>%
  do(tidy(lm(Value ~ Bin, data = .))) %>%
  filter(term == "Bin") %>%
  rename(Slope = estimate,
         Slope_SE = std.error,
         Slope_t = statistic,
         Slope_p = p.value) %>%
  ungroup()

# -----------------------------
# T-tests: EF5 & DMSO vs Vehicle
# -----------------------------
vehicle_name <- "Vehicle"
metrics <- unique(slopes_df$Metric)
days <- unique(slopes_df$Day)
treatments <- setdiff(unique(slopes_df$TreatmentGroup), vehicle_name)

comparison_results <- list()

for (metric in metrics) {
  for (day in days) {
    for (treat in treatments) {
      df_sub <- slopes_df %>%
        filter(Metric == metric,
               Day == day,
               TreatmentGroup %in% c(vehicle_name, treat))
      
      if(length(unique(df_sub$TreatmentGroup)) == 2) {
        ttest <- t.test(Slope ~ TreatmentGroup, data = df_sub)
        comparison_results <- append(comparison_results,
                                     list(data.frame(
                                       Metric = metric,
                                       Day = day,
                                       Treatment = treat,
                                       estimate_vehicle = ttest$estimate[vehicle_name],
                                       estimate_treatment = ttest$estimate[treat],
                                       estimate_diff = unname(diff(rev(ttest$estimate))),
                                       conf.low = ttest$conf.int[1],
                                       conf.high = ttest$conf.int[2],
                                       statistic = ttest$statistic,
                                       p.value = ttest$p.value
                                     )))
      }
    }
  }
}

comparison_results <- bind_rows(comparison_results)

# Save results
write.csv(comparison_results, "Slope_Comparison_vs_Vehicle.csv", row.names = FALSE)

# -----------------------------
# Plot slopes (separate plots per metric)
# -----------------------------
plot_width <- 6
plot_height <- 5
plot_dpi <- 300

for (metric in unique(slopes_df$Metric)) {
  p <- ggplot(slopes_df %>% filter(Metric == metric),
              aes(x = Day, y = Slope, fill = TreatmentGroup)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA,
                 position = position_dodge(width = 0.7)) +
    geom_jitter(aes(color = TreatmentGroup),
                position = position_jitterdodge(jitter.width = 0.15,
                                                dodge.width = 0.7),
                size = 2, alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = paste("Slope of", metric, "across Days"),
         y = "Slope (Value ~ Bin)",
         x = "Day") +
    theme_minimal() +
    theme(legend.position = "top")
  
  ggsave(filename = paste0("EF5_Slopes_", metric, ".tiff"),
         plot = p,
         width = plot_width,
         height = plot_height,
         dpi = plot_dpi)
}

