##############################################################
# Load Libraries
##############################################################

library(ggplot2)
library(dplyr)
library(readxl)
library(stringr)
library(car)

# Install agricolae if missing
if (!require(agricolae)) install.packages("agricolae")
library(agricolae)

##############################################################
# Data Import
##############################################################

setwd("C:/Users/jps558/OneDrive - University of York/Desktop/R_WorkingDirectory/R_WFHDirectory/spheroids/Supplement test")

Path <- "B27+Glutamax Spheroid sizes.xlsx"

sheets <- c("500 cellswell + sup", "500 cellswell Red", "500 cellswell Clear",
            "1000 cellswell + sup", "1000 cellswell Red", "1000 cellswell Clear")

desired_order <- c(
  "500 cellswell Red",
  "1000 cellswell Red",
  "500 cellswell Clear",
  "1000 cellswell Clear",
  "500 cellswell + sup",
  "1000 cellswell + sup"
)

all_data <- lapply(sheets, function(s) {
  df <- read_excel(Path, sheet = s)
  
  df <- df[, colSums(!is.na(df)) > 0]
  
  size_col <- grep("Day\\s*7.*size", names(df), ignore.case = TRUE, value = TRUE)
  if(length(size_col) != 1){
    stop("Could not identify Day 7 size column in sheet: ", s)
  }
  
  names(df)[names(df) == size_col] <- "Day_7_size_um"
  df$Treatment <- s
  return(df)
}) %>% bind_rows()

all_data <- all_data %>%
  mutate(
    Day_7_size_um = as.numeric(Day_7_size_um),
    Treatment = factor(Treatment, levels = desired_order)
  )

##############################################################
# ANOVA + Tukey HSD + Significance Letters
##############################################################

anova_model <- aov(Day_7_size_um ~ Treatment, data = all_data)
print(summary(anova_model))

# Assumptions
print(shapiro.test(residuals(anova_model)))
print(leveneTest(Day_7_size_um ~ Treatment, data = all_data))

# Tukey HSD with agricolae
tukey_res <- HSD.test(anova_model, "Treatment", group = TRUE)
cld_table <- tukey_res$groups
cld_table$Treatment <- rownames(cld_table)
rownames(cld_table) <- NULL

cld_table <- cld_table %>%
  select(Treatment, groups) %>%
  rename(Letters = groups)

##############################################################
# Summary stats + merge letters
##############################################################

summary_data <- all_data %>%
  group_by(Treatment) %>%
  summarise(
    mean_size = mean(Day_7_size_um, na.rm = TRUE),
    sem = sd(Day_7_size_um, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

summary_data$Treatment <- factor(summary_data$Treatment, levels = desired_order)

summary_data <- merge(summary_data, cld_table, by = "Treatment", all.x = TRUE)

##############################################################
# Determine letter heights for plots
##############################################################

letter_height_bar <- max(summary_data$mean_size + summary_data$sem) + 20
letter_height_box <- max(all_data$Day_7_size_um, na.rm = TRUE) + 20

##############################################################
# BOX PLOT with significance letters
##############################################################

tiff("Spheroid_Size_Day7_by_Treatment_Boxplot_with_Significance.tiff",
     width = 2480, height = 1748, res = 300)

ggplot(all_data, aes(x = Treatment, y = Day_7_size_um, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(position = position_jitter(width = 0.15),
              alpha = 0.4, size = 1.5, shape = 21) +
  geom_text(data = summary_data,
            aes(x = Treatment, y = letter_height_box, label = Letters),
            size = 6) +
  labs(
    title = "Spheroid Diameter at Day 7 by Treatment",
    x = "Treatment",
    y = "Spheroid Diameter (µm)"
  ) +
  theme_minimal() +
  coord_cartesian(ylim = c(300, letter_height_box + 30)) +
  theme(
    axis.text.x = element_text(size = 12, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 13),
    plot.title = element_text(size = 15, face = "bold"),
    legend.position = "none"
  )

dev.off()

##############################################################
# BAR PLOT with significance letters
##############################################################

tiff("Mean_Spheroid_Size_Day7_with_Significance.tiff",
     width = 2480, height = 1748, res = 300)

ggplot(summary_data, aes(x = Treatment, y = mean_size, fill = Treatment)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = mean_size - sem, ymax = mean_size + sem),
                width = 0.3) +
  geom_text(aes(label = Letters, y = letter_height_bar), size = 6) +
  labs(
    title = "Mean Spheroid Diameter at Day 7 by Treatment",
    x = "Treatment",
    y = "Mean Spheroid Diameter (µm)"
  ) +
  theme_minimal() +
  coord_cartesian(ylim = c(300, letter_height_bar + 30)) +
  theme(
    axis.text.x = element_text(size = 12, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 13),
    plot.title = element_text(size = 15, face = "bold"),
    legend.position = "none"
  )

dev.off()

##############################################################
# END OF SCRIPT
##############################################################
