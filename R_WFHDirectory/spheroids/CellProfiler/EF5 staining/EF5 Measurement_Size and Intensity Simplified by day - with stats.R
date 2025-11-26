# -------------------------
# Load packages
# -------------------------
library(tidyverse)
library(ggplot2)
library(stringr)
library(rstatix)
library(ggpubr)
library(writexl)

# -------------------------
# Set working directory
# -------------------------
setwd("C:/Users/jps558/OneDrive - University of York/Desktop/R_WorkingDirectory/R_WFHDirectory/spheroids/CellProfiler/EF5 staining") 

# -------------------------
# Load data
# -------------------------
high_signal <- read_csv("High Signal RegionsFull_HighSignal.csv")
whole_sphere <- read_csv("High Signal RegionsWholeSphereMerged.csv")
minus_high_signal <- read_csv("High Signal RegionsSpheroidMinusHigh.csv")

# -------------------------
# Select relevant columns and tag Region
# -------------------------
high_signal_selected <- high_signal %>%
  select(FileName_Brightfield, AreaShape_Area, AreaShape_MaxFeretDiameter,
         Intensity_MeanIntensity_Fluorescence, Intensity_MaxIntensity_Fluorescence) %>%
  mutate(Region = "HighSignal")

whole_sphere_selected <- whole_sphere %>%
  select(FileName_Brightfield, AreaShape_Area, AreaShape_MaxFeretDiameter,
         Intensity_MeanIntensity_Fluorescence, Intensity_MaxIntensity_Fluorescence) %>%
  mutate(Region = "WholeSphere")

minus_high_selected <- minus_high_signal %>%
  select(FileName_Brightfield, AreaShape_Area, AreaShape_MaxFeretDiameter,
         Intensity_MeanIntensity_Fluorescence, Intensity_MaxIntensity_Fluorescence) %>%
  mutate(Region = "SpheroidMinusHigh")

# -------------------------
# Combine data
# -------------------------
combined_data <- bind_rows(high_signal_selected, whole_sphere_selected, minus_high_selected)

# Convert pixels to µm
px_to_um <- 702 / 3120  # 0.225 µm/pixel
area_conversion <- px_to_um^2  # 0.050625

# -------------------------
# Parse identifiers
# -------------------------
combined_data <- combined_data %>%
  mutate(
    Treatment1 = case_when(
      str_detect(FileName_Brightfield, regex("Amiloride", ignore_case = TRUE)) ~ "Amiloride",
      str_detect(FileName_Brightfield, regex("DMSO", ignore_case = TRUE)) ~ "DMSO",
      TRUE ~ "None"
    ),
    Treatment2 = case_when(
      str_detect(FileName_Brightfield, regex("(?i)^vehicle")) ~ "Vehicle",
      str_detect(FileName_Brightfield, regex("(?i)EF5")) ~ "EF5",
      TRUE ~ "Other"
    ),
    Day = str_extract(FileName_Brightfield, regex("(?i)day[_-]?\\d+")) %>%
      str_remove_all("\\D+") %>%
      as.integer(),
    SphereID = str_extract(FileName_Brightfield, "Sphere\\d+_?\\d*[-_][^_]+_c\\d+"),
    Area_um2 = AreaShape_Area * area_conversion,
    MaxFeretDiameter_um = AreaShape_MaxFeretDiameter * px_to_um
  )

# -------------------------
# EF5-only dataset
# -------------------------
ef5_only <- combined_data %>%
  filter(Treatment1 == "None", Treatment2 == "EF5")

# -------------------------
# Vehicle dataset (prefer WholeSphere)
# -------------------------
vehicle_all <- combined_data %>%
  filter(Treatment1 == "None", Treatment2 == "Vehicle")
vehicle_whole <- vehicle_all %>% filter(Region == "WholeSphere")
if (nrow(vehicle_whole) == 0) vehicle_whole <- vehicle_all

# -------------------------
# Size over time with Vehicle
# -------------------------
size_timecourse_ef5 <- ef5_only %>%
  filter(Region %in% c("HighSignal", "WholeSphere")) %>%
  pivot_longer(
    cols = c(Area_um2, MaxFeretDiameter_um),
    names_to = "Measurement",
    values_to = "Value"
  ) %>%
  mutate(
    Measurement = recode(Measurement,
                         Area_um2 = "Area (µm²)",
                         MaxFeretDiameter_um = "Max Feret Diameter (µm)"),
    Region = recode(Region,
                    HighSignal = "High Signal Region",
                    WholeSphere = "Whole Sphere")
  )

vehicle_size_pooled <- vehicle_whole %>%
  pivot_longer(
    cols = c(Area_um2, MaxFeretDiameter_um),
    names_to = "Measurement",
    values_to = "Value"
  ) %>%
  mutate(
    Measurement = recode(Measurement,
                         Area_um2 = "Area (µm²)",
                         MaxFeretDiameter_um = "Max Feret Diameter (µm)"),
    Region = "Vehicle"
  )

size_timecourse_plot <- bind_rows(size_timecourse_ef5, vehicle_size_pooled) %>%
  mutate(Region = factor(Region, levels = c("High Signal Region", "Whole Sphere", "Vehicle")))

tiff("Spheroid Size Over Time (EF5 with Vehicle control added).tiff",
     width = 2780, height = 2048, res = 300)
ggplot(size_timecourse_plot, aes(x = factor(Day), y = Value, fill = Region)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  facet_wrap(~Measurement, scales = "free_y") +
  theme_minimal(base_size = 15) +   # increased font size
  labs(
    title = "Spheroid Size Over Time (EF5 with Vehicle control added)",
    x = "Day",
    y = NULL,
    fill = "Sample"
  ) +
  theme(
    axis.text.x = element_text(size = 14, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 14),
    plot.title = element_text(size = 17, face = "bold"),
    strip.text = element_text(size = 15),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 15)
  )
dev.off()

# -------------------------
# Intensity over time with Vehicle (Mean Intensity)
# -------------------------
intensity_timecourse_ef5 <- ef5_only %>%
  filter(Region %in% c("HighSignal", "WholeSphere")) %>%
  select(FileName_Brightfield, Region, Day, Intensity_MeanIntensity_Fluorescence) %>%
  pivot_longer(
    cols = c(Intensity_MeanIntensity_Fluorescence),
    names_to = "Measurement",
    values_to = "Value"
  ) %>%
  mutate(
    Measurement = "Mean Intensity",
    Region = recode(Region,
                    HighSignal = "High Signal Region",
                    WholeSphere = "Whole Sphere")
  )

vehicle_intensity_pooled <- vehicle_whole %>%
  select(FileName_Brightfield, Region, Day, Intensity_MeanIntensity_Fluorescence) %>%
  pivot_longer(
    cols = c(Intensity_MeanIntensity_Fluorescence),
    names_to = "Measurement",
    values_to = "Value"
  ) %>%
  mutate(
    Measurement = "Mean Intensity",
    Region = "Vehicle"
  )

intensity_timecourse_plot <- bind_rows(intensity_timecourse_ef5, vehicle_intensity_pooled) %>%
  mutate(Region = factor(Region, levels = c("High Signal Region", "Whole Sphere", "Vehicle")))

tiff("Fluorescence Intensity Over Time (High Signal vs Whole Sphere vs Vehicle).tiff",
     width = 2780, height = 2048, res = 300)
ggplot(intensity_timecourse_plot, aes(x = factor(Day), y = Value, fill = Region)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  facet_wrap(~Measurement, scales = "free_y") +
  theme_minimal(base_size = 15) +   # increased font size
  labs(
    title = "Fluorescence Intensity Over Time (High Signal vs Whole Sphere vs Vehicle)",
    x = "Day",
    y = "Mean Intensity",
    fill = "Sample"
  ) +
  theme(
    axis.text.x = element_text(size = 14, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 14),
    plot.title = element_text(size = 17, face = "bold"),
    strip.text = element_text(size = 15),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 15)
  )
dev.off()
# -------------------------
# Pairwise significance: Size
# -------------------------
pairwise_size_stats <- size_timecourse_plot %>%
  group_by(Measurement, Day) %>%
  group_modify(~{
    compare_means(
      Value ~ Region,
      data = .x,
      p.adjust.method = "BH"
    ) %>%
      filter(group1 == "High Signal Region" | group2 == "High Signal Region")
  }) %>%
  ungroup() %>%
  add_significance()

# -------------------------
# Pairwise significance: Intensity
# -------------------------
pairwise_intensity_stats <- intensity_timecourse_plot %>%
  group_by(Measurement, Day) %>%
  group_modify(~{
    compare_means(
      Value ~ Region,
      data = .x,
      p.adjust.method = "BH"
    ) %>%
      filter(group1 == "High Signal Region" | group2 == "High Signal Region")
  }) %>%
  ungroup() %>%
  add_significance()

# -------------------------
# Export to Excel
# -------------------------
write_xlsx(
  list(
    Size = pairwise_size_stats,
    Intensity = pairwise_intensity_stats
  ),
  path = "pairwise_stats.xlsx"
)

