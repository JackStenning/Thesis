# ================================================================
# Spheroid Fluorescence Intensity Analysis: Slopes vs Pooled Vehicle
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
  select(FileName_Brightfield, matches("^RadialDistribution_(FracAtD|MeanFrac)_Fluorescence_[1-6]of6$"))

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
    SphereNum = ifelse(str_detect(SphereRaw, "_"), str_extract(SphereRaw, "(?<=_)\\d+"), "1"),
    ImageExportNum = str_extract(FileName_Brightfield, "Image Export-\\d+"),
    SphereID = paste0(SphereRaw, "_", ImageExportNum),
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
  rename(
    Slope = estimate,
    Slope_SE = std.error,
    Slope_t = statistic,
    Slope_p = p.value
  ) %>%
  ungroup()

# -----------------------------
# Pooled Vehicle t-tests
# -----------------------------
vehicle_name <- "Vehicle"
treatments <- setdiff(unique(slopes_df$TreatmentGroup), vehicle_name)

# Pooled Vehicle slopes
vehicle_df <- slopes_df %>% filter(TreatmentGroup == vehicle_name)

# EF5/DMSO slopes
ef5_df <- slopes_df %>% filter(TreatmentGroup %in% treatments)

# Unique combinations of Metric, Day, Treatment for EF5/DMSO
combos <- ef5_df %>%
  distinct(Metric, Day, TreatmentGroup) %>%
  rename(Treatment = TreatmentGroup)

# Run rowwise t-tests safely
comparison_results <- combos %>%
  rowwise() %>%
  mutate(ttest = list({
    cur_metric <- Metric
    cur_day <- Day
    cur_treat <- Treatment
    
    # EF5/DMSO subset for this Metric/Day/Treatment
    ef5_sub <- ef5_df %>% filter(Metric == cur_metric, Day == cur_day, TreatmentGroup == cur_treat)
    
    # Vehicle subset (pooled across all days)
    vehicle_sub <- vehicle_df %>% filter(Metric == cur_metric)
    
    # Combine for t-test
    df_sub <- bind_rows(
      ef5_sub %>% select(Slope, TreatmentGroup),
      vehicle_sub %>% select(Slope, TreatmentGroup)
    )
    
    # Factor order: Vehicle first
    df_sub$TreatmentGroup <- factor(df_sub$TreatmentGroup, levels = c(vehicle_name, cur_treat))
    
    # Run t-test
    t.test(Slope ~ TreatmentGroup, data = df_sub) %>% tidy()
  })) %>%
  unnest(ttest) %>%
  rename(
    estimate_vehicle = estimate1,
    estimate_treatment = estimate2
  ) %>%
  mutate(
    estimate_diff = estimate_treatment - estimate_vehicle
  ) %>%
  select(Metric, Day, Treatment, estimate_vehicle, estimate_treatment, estimate_diff, conf.low, conf.high, statistic, p.value)

# Save results
write.csv(comparison_results, "Slope_Comparison_vs_PooledVehicle.csv", row.names = FALSE)

# -----------------------------
# Plot slopes (optional)
# -----------------------------
plot_width <- 6
plot_height <- 5
plot_dpi <- 300

for (metric in unique(slopes_df$Metric)) {
  p <- ggplot(slopes_df %>% filter(Metric == metric), aes(x = Day, y = Slope, fill = TreatmentGroup)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA, position = position_dodge(width = 0.7)) +
    geom_jitter(aes(color = TreatmentGroup), position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7), size = 2, alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = paste("Slope of", metric, "across Days"), y = "Slope (Value ~ Bin)", x = "Day") +
    theme_minimal() +
    theme(legend.position = "top")
  
  ggsave(filename = paste0("EF5_Slopes_", metric, ".tiff"), plot = p, width = plot_width, height = plot_height, dpi = plot_dpi)
}
