# Install packages if needed
# install.packages("readxl")
# install.packages("dplyr")
# install.packages("ggplot2")

install.packages(c("readxl","dplyr","ggplot2","stringr","janitor"))
library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)
library(janitor)


#Set working directory
setwd("C:/Users/jps558/OneDrive - University of York/Desktop/R_WorkingDirectory/R_WFHDirectory/Chapter 3") 


# ---- 1) Read the Excel file ----
df <- read_excel("Cell counting across conditions.xlsx", sheet = 1) %>% 
  clean_names()   # makes names like top_percent_red_1, condition_name, etc.

# ---- 2) find which column holds the condition text ----
cond_col <- if ("condition_name" %in% names(df)) {
  "condition_name"
} else if ("condition" %in% names(df)) {
  "condition"
} else {
  stop("Couldn't find a condition column (looked for 'condition_name' or 'condition').")
}

# ---- 3) identify the four '% Red' columns (based on 'percent' in cleaned names) ----
percent_cols <- names(df)[str_detect(names(df), "percent")]
if (length(percent_cols) < 1) {
  stop("No percent columns found. Check column names -- expect names containing 'percent' after clean_names().")
}

# (optional) print the percent columns found:
message("Using these percent columns for per-row averages: ", paste(percent_cols, collapse = ", "))

# ---- 4) compute per-row (replicate) average across the found % columns ----
df2 <- df %>%
  # ensure percent columns are numeric
  mutate(across(all_of(percent_cols), ~ as.numeric(.))) %>%
  rowwise() %>%
  mutate(rep_mean = mean(c_across(all_of(percent_cols)), na.rm = TRUE)) %>%
  ungroup()

# ---- 5) clean the condition label and get numeric dose for ordering ----
df2 <- df2 %>%
  mutate(
    condition_clean = str_remove(.data[[cond_col]], "\\s*-\\s*Rep.*"),           # remove "- Rep 1/2"
    condition_numeric = as.numeric(str_extract(condition_clean, "^[0-9.]+"))    # extract leading numeric (e.g. 0.35)
  )

# ---- 6) summarize across replicates: mean and SD of the replicate means ----
df_summary <- df2 %>%
  group_by(condition_clean, condition_numeric) %>%
  summarise(
    Mean = mean(rep_mean, na.rm = TRUE),
    SD   = sd(rep_mean, na.rm = TRUE),
    n_rep = sum(!is.na(rep_mean)),
    .groups = "drop"
  ) %>%
  mutate(
    ymin = pmax(Mean - SD, 0),    # avoid negative lower error
    ymax = Mean + SD
  )

# show the computed table
print(df_summary)

# ---- 7) Plot: bars ordered by numeric dose, with SD error bars ----
p <- ggplot(df_summary, aes(x = reorder(condition_clean, condition_numeric),
                            y = Mean,
                            fill = condition_clean)) +
  geom_col() +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2, color = "black") +
  theme_minimal(base_size = 14) +
  labs(title = "Average Co-transfection Efficiency (mean ± SD)",
       x = "Condition",
       y = "Mean percentage of co-transfected cells (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none")

# If you want to zoom the y-axis (without dropping data), uncomment and set your limits:
p <- p + coord_cartesian(ylim = c(0, 12))

print(p)

# ---- 8) (optional) Save plot as A5 @ 300 dpi ----
 ggsave("Average Co-transfection Efficiency (mean+- SD).png", plot = p,
        width = 210/25.4, height = 148/25.4)
