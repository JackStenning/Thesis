# Load necessary libraries
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)

# Read the Excel file (modify the file path accordingly)
setwd("//userfs/jps558/w2k/Desktop/R_WorkingDirectory") 

# Read in the data
file_path <- "Follow up plots and scripts/Puromycin relative survival.xlsx"
df <- read_excel(file_path)

# Reshape the data: Gather all relative survival columns into long format
df_long <- df %>%
  pivot_longer(
    cols = -1, # Keep the first column (Concentration) intact
    names_to = c("Day", "Replicate"),
    names_pattern = "Day (\\d+) RS rep (\\d+)",
    values_to = "Relative_Survival"
  ) %>%
  mutate(Day = as.numeric(Day))

# Calculate mean and standard deviation for each day and concentration
df_summary <- df_long %>%
  group_by(`Puromycin concentration`, Day) %>%
  summarise(
    Mean_RS = mean(Relative_Survival, na.rm = TRUE),
    SD_RS = sd(Relative_Survival, na.rm = TRUE),
    .groups = 'drop'
  )

# Plot the kill curve
ggplot(df_summary, aes(x = Day, y = Mean_RS, color = as.factor(`Puromycin concentration`), group = `Puromycin concentration`)) +
  geom_line(size = 1) +       # Draw trend lines
  geom_point(size = 2) +      # Show individual points
  geom_errorbar(aes(ymin = Mean_RS - SD_RS, ymax = Mean_RS + SD_RS), width = 0.2) +
  scale_x_continuous(breaks = 0:4) +  # Ensure Days are shown as whole numbers
  labs(title = "Puromycin Kill Curve Over 4 Days",
       x = "Days",
       y = "Relative Survival (%)",
       color = "Puromycin \nconcentration \n(µg/mL)") +
  theme_minimal()

# Box plot
ggplot(df_long, aes(x = as.factor(Day), y = Relative_Survival, fill = as.factor(`Puromycin concentration`))) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  # Box plot with transparency and no outliers
  geom_jitter(aes(color = as.factor(`Puromycin concentration`)), width = 0.2, size = 1.5, alpha = 0.6) + # Add jittered points
  labs(title = "Puromycin Kill Curve (Box Plot)",
       x = "Days",
       y = "Relative Survival (%)",
       fill = "Concentration (µg/mL)",
       color = "Concentration (µg/mL)") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2") +  # Nice color scheme
  scale_color_brewer(palette = "Dark2")  # Different colors for points
