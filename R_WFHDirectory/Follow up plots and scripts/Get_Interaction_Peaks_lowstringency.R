#Load modules
library(dplyr)
library(tidyr)
library(ggplot2)

#Set working directory
setwd("C:/Users/jps558/OneDrive - University of York/Desktop/R_WorkingDirectory/R_WFHDirectory") 

#Load data
input_dir <- "overlap/Chia_cc"

files <- list.files(
  input_dir,
  pattern = "^basic_Chia_ChiapetRep[123]_peak_data_ER_WT4",
  full.names = TRUE
)

# Read and combine with replicate labels
all_data <- lapply(files, function(f) {
  df <- read.table(f, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  df$Replicate <- basename(f)  # add file name as replicate label
  return(df)
}) %>% bind_rows()

#Colnames
colnames(all_data)[1:27] <- c("chr_chia","start_chia","end_chia","name","score","strand",
                              "thickStart","thickEnd","itemRgb",
                              "blockCount","blockSizes","blockStarts", "Chr_CC",	"Start_CC",	"End_CC",	"Center",	"Experiment Insertions",	"Background insertions",	"Reference Insertions",	"pvalue Reference",	"pvalue Background",	"Fraction Experiment",	"TPH Experiment",	"Fraction background",	"TPH background",	"TPH background subtracted",	"pvalue_adj Reference")

# Subset relevant columns
subset_data <- all_data %>%
  select(chr_chia, start_chia, end_chia, name,
         Chr_CC, Start_CC, End_CC, Replicate)

# Split into anchor1 and anchor2 datasets
anchor1 <- subset_data %>%
  group_by(name, Replicate) %>%
  slice(1) %>%  # first row = anchor1
  ungroup() %>%
  rename(chr_chia_anchor1 = chr_chia,
         start_chia_anchor1 = start_chia,
         end_chia_anchor1 = end_chia,
         Chr_CC_anchor1 = Chr_CC,
         Start_CC_anchor1 = Start_CC,
         End_CC_anchor1 = End_CC)

anchor2 <- subset_data %>%
  group_by(name, Replicate) %>%
  slice(2:n()) %>%  # remaining rows = anchor2
  ungroup() %>%
  rename(chr_chia_anchor2 = chr_chia,
         start_chia_anchor2 = start_chia,
         end_chia_anchor2 = end_chia,
         Chr_CC_anchor2 = Chr_CC,
         Start_CC_anchor2 = Start_CC,
         End_CC_anchor2 = End_CC)

# 1. Keep all duplicates (all combinations)
final_df_all <- anchor1 %>%
  inner_join(anchor2, by = c("name", "Replicate"))

# 2. Keep only one row per ChIA-PET interaction
final_df_unique <- final_df_all %>%
  group_by(name, Replicate) %>%
  slice(1) %>%   # keep only the first combination
  ungroup()

# Output all combinations with duplicates
write.csv(final_df_all, "Chiapet_CC_highstringency_all_duplicates.csv", row.names = FALSE)

# Output only one row per ChIA-PET interaction
write.csv(final_df_unique, "Chiapet_CC_highstringency_unique.csv", row.names = FALSE)

## Total number of ER interactions from Fullwood 2007 = 6114
# Total number of ChIA-PET interactions
total_interactions <- 6114

# Count unique ChIA-PET interactions with at least one overlapping calling card
overlapping_interactions <- length(final_df_unique$name)

# Calculate percentage
percent_overlap <- (overlapping_interactions / total_interactions) * 100
percent_no_overlap <- 100 - percent_overlap

# Create data frame for plotting
plot_df <- data.frame(
  Status = c("Overlapping", "Not overlapping"),
  Percentage = c(percent_overlap, percent_no_overlap)
)

# Plot
tiff("Percentage of ChIA-PET Interactions Overlapping with low stringency Calling Card Peaks.tiff", width = 3510, height = 2481, res = 300)
ggplot(plot_df, aes(x = "", y = Percentage, fill = Status)) +
  geom_bar(stat = "identity", width = 0.5) +
  coord_flip() +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")), 
            position = position_stack(vjust = 0.5), size = 5) +
  labs(title = "Percentage of ChIA-PET Interactions Overlapping with low stringency aCalling Card Peaks",
       x = NULL, y = "Percentage") +
  theme_minimal() +
  scale_fill_manual(values = c("steelblue", "lightgray"))
dev.off()