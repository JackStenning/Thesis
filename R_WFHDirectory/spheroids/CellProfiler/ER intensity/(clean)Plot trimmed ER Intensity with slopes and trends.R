
# -----------------------------
# Load necessary packages
# -----------------------------
library(readr)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(broom)
library(writexl)  # for Excel output

# -----------------------------
# Set working directory
# -----------------------------
setwd("C:/Users/jps558/OneDrive - University of York/Desktop/R_WorkingDirectory/R_WFHDirectory/spheroids/CellProfiler/ER intensity") 

# -----------------------------
# Read ER data
# -----------------------------
df <- read_excel("trimmedER Intensity from core_NoAmilorideCompleteSphere.xlsx", sheet = 1)

# -----------------------------
# Select relevant columns
# -----------------------------
df_selected <- df %>%
  select(
    FileName_Brightfield,
    matches("^RadialDistribution_(FracAtD|MeanFrac)_Fluorescence_[1-6]of6$")
  ) %>%
  mutate(across(matches("^RadialDistribution_.*_Fluorescence_[1-6]of6$"), as.numeric))

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
# Extract sphere & treatment info (with ordered Day factor)
# -----------------------------
df_long <- df_long %>%
  mutate(
    SphereRaw = str_extract(FileName_Brightfield, "Sphere\\d+(_\\d+)?"),
    ImageExportNum = str_extract(FileName_Brightfield, "Image Export-\\d+"),
    SphereID = paste0(SphereRaw, "_", ImageExportNum),
    
    # Identify treatment
    PrimaryCompound = ifelse(str_starts(FileName_Brightfield, "DMSO_"), "DMSO", "EF5"),
    SecondaryCompound = ifelse(str_starts(FileName_Brightfield, "Secondary"), "Secondary", NA_character_),
    TreatmentGroup = ifelse(is.na(SecondaryCompound), PrimaryCompound,
                            paste(PrimaryCompound, SecondaryCompound, sep = "+")),
    
    # Extract day only for EF5; secondary will be grouped as one "Secondary" day
    EF5DayMatch = str_extract(FileName_Brightfield, "EF5Day\\d+(-\\d+)?"),
    DayNum = as.integer(str_extract(EF5DayMatch, "(?<=EF5Day)\\d+")),
    Day = case_when(
      !is.na(DayNum) ~ paste0("Day ", DayNum),
      !is.na(SecondaryCompound) ~ "Secondary",
      TRUE ~ NA_character_
    ),
    # Ensure factor has the desired order
    Day = factor(Day, levels = c("Day 7", "Day 9", "Day 11", "Day 13", "Secondary"))
  )

# Assign numeric IDs per metric/day
df_long <- df_long %>%
  group_by(Metric, Day) %>%
  mutate(SphereSimpleID = factor(match(SphereID, unique(SphereID)))) %>%
  ungroup()

# ================================================================
# Slope analysis
# ================================================================

# -----------------------------
# Slope calculation per sphere (include Secondary as control)
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
# Pooled Secondary t-tests
# -----------------------------
control_name <- "Secondary"
treatments <- setdiff(unique(slopes_df$TreatmentGroup), control_name)

# Secondary slopes (pooled across all days)
secondary_df <- slopes_df %>% filter(TreatmentGroup == control_name)

# EF5/DMSO slopes
treat_df <- slopes_df %>% filter(TreatmentGroup %in% treatments)

# Unique combinations of Metric, Day, Treatment for EF5/DMSO
combos <- treat_df %>%
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
    treat_sub <- treat_df %>% 
      filter(Metric == cur_metric, Day == cur_day, TreatmentGroup == cur_treat)
    
    # Secondary subset (pooled across all days, same metric)
    secondary_sub <- secondary_df %>% 
      filter(Metric == cur_metric)
    
    # Combine for t-test
    df_sub <- bind_rows(
      treat_sub %>% select(Slope, TreatmentGroup),
      secondary_sub %>% select(Slope, TreatmentGroup)
    )
    
    # Factor order: Secondary first
    df_sub$TreatmentGroup <- factor(df_sub$TreatmentGroup, 
                                    levels = c(control_name, cur_treat))
    
    # Run t-test
    t.test(Slope ~ TreatmentGroup, data = df_sub) %>% tidy()
  })) %>%
  unnest(ttest) %>%
  rename(
    estimate_secondary = estimate1,
    estimate_treatment = estimate2
  ) %>%
  mutate(
    estimate_diff = estimate_treatment - estimate_secondary
  ) %>%
  select(Metric, Day, Treatment,
         estimate_secondary, estimate_treatment, estimate_diff,
         conf.low, conf.high, statistic, p.value)

# -----------------------------
# Save results to Excel
# -----------------------------
write_xlsx(comparison_results, "Slope_Comparison_vs_Secondary.xlsx")

# ================================================================
# Slope plots
# ================================================================
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
                size = 2, alpha = 0.7)
  