# Load necessary packages
library(tidyverse)
library(ggplot2)
library(stringr)

#Set working directory
setwd("//userfs/jps558/w2k/Desktop/R_WorkingDirectory/CellProfiler/EF5 staining") 

#Load data
high_signal <- read_csv("High Signal RegionsFull_HighSignal.csv")
whole_sphere <- read_csv("High Signal RegionsWholeSphereMerged.csv")
minus_high_signal <- read_csv("High Signal RegionsSpheroidMinusHigh.csv")


#Select data
high_signal_selected <- high_signal %>%
  select(FileName_Brightfield, AreaShape_Area, AreaShape_MaxFeretDiameter, Intensity_IntegratedIntensity_Fluorescence, Intensity_MaxIntensity_Fluorescence) %>%
  mutate(Region = "HighSignal")

whole_sphere_selected <- whole_sphere %>%
  select(FileName_Brightfield, AreaShape_Area, AreaShape_MaxFeretDiameter, Intensity_IntegratedIntensity_Fluorescence, Intensity_MaxIntensity_Fluorescence) %>%
  mutate(Region = "WholeSphere")

minus_high_selected <- minus_high_signal %>%
  select(FileName_Brightfield, AreaShape_Area, AreaShape_MaxFeretDiameter, Intensity_IntegratedIntensity_Fluorescence, Intensity_MaxIntensity_Fluorescence) %>%
  mutate(Region = "SpheroidMinusHigh")

#Combine Data
combined_data <- bind_rows(high_signal_selected, whole_sphere_selected, minus_high_selected)

#Convert Pixels to µm
px_to_um <- 702 / 3120  # = 0.225 µm/pixel
area_conversion <- px_to_um^2  # = 0.050625

#Split out identifiers
combined_data <- combined_data %>%
  mutate(
    Treatment1 = case_when(
      str_detect(FileName_Brightfield, "Amiloride") ~ "Amiloride",
      str_detect(FileName_Brightfield, "DMSO") ~ "DMSO",
      TRUE ~ "None"
    ),
    Treatment2 = "EF5",
    Day = str_extract(FileName_Brightfield, "Day\\d+") %>% str_remove("Day") %>% as.integer(),
    SphereID = str_extract(FileName_Brightfield, "Sphere\\d+_?\\d*[-_][^_]+_c\\d+"),
    Area_um2 = AreaShape_Area * area_conversion,
    MaxFeretDiameter_um = AreaShape_MaxFeretDiameter * px_to_um
  )

#Filter only EF5 single treatment
ef5_only <- combined_data %>%
  filter(Treatment1 == "None")

#Make data Long
#Area and Diameter
size_long <- ef5_only %>%
  pivot_longer(
    cols = c(Area_um2, MaxFeretDiameter_um),
    names_to = "Measurement",
    values_to = "Value"
  ) %>%
  mutate(Measurement = recode(Measurement,
                              Area_um2 = "Area (µm²)",
                              MaxFeretDiameter_um = "Max Feret Diameter (µm)"))

#Intensity
intensity_long <- ef5_only %>%
  pivot_longer(
    cols = c(Intensity_IntegratedIntensity_Fluorescence, Intensity_MaxIntensity_Fluorescence),
    names_to = "Measurement",
    values_to = "Value"
  ) %>%
  mutate(Measurement = recode(Measurement,
                              IntegratedIntensity = "Integrated Intensity",
                              MaxIntensity = "Max Intensity"))

#Plot
ggplot(size_long, aes(x = Region, y = Value, fill = Region)) +
  geom_boxplot(outlier.shape = 21, outlier.color = "black", alpha = 0.8) +
  facet_wrap(~Measurement, scales = "free_y") +
  theme_minimal(base_size = 13) +
  labs(
    title = "Spheroid Size Comparison (EF5-only)",
    x = "Region",
    y = NULL
  )

############### Intensity plot ###############

#plot
ggplot(intensity_long, aes(x = Region, y = Value, fill = Region)) +
  geom_boxplot(outlier.shape = 21, outlier.color = "black", alpha = 0.8) +
  facet_wrap(~Measurement, scales = "free_y") +
  theme_minimal(base_size = 13) +
  labs(
    title = "Fluorescence Intensity Comparison (EF5-only)",
    x = "Region",
    y = "Intensity"
  )

#Plot by day

# Filter to regions of interest (excluding WholeSphere if you only want High vs MinusHigh)
timecourse_data <- ef5_only %>%
  filter(Region %in% c("HighSignal", "SpheroidMinusHigh")) %>%
  pivot_longer(
    cols = c(Intensity_IntegratedIntensity_Fluorescence, Intensity_MaxIntensity_Fluorescence),
    names_to = "Measurement",
    values_to = "Value"
  ) %>%
  mutate(
    Measurement = recode(Measurement,
                         IntegratedIntensity = "Integrated Intensity",
                         MaxIntensity = "Max Intensity"),
    Region = recode(Region,
                    HighSignal = "High Signal Region",
                    SpheroidMinusHigh = "Spheroid Minus High Signal")
  )

#Plot
ggplot(timecourse_data, aes(x = factor(Day), y = Value, fill = Region)) +
  geom_boxplot(outlier.shape = 21, outlier.color = "black", alpha = 0.8) +
  facet_wrap(~Measurement, scales = "free_y") +
  theme_minimal(base_size = 13) +
  labs(
    title = "Fluorescence Intensity Over Time (EF5-only Spheroids)",
    x = "Day",
    y = "Intensity",
    fill = "Region"
  )
