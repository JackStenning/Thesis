# Load necessary library
library(ggplot2)
library(dplyr)
library(IRanges)

#Set working directory
setwd("//userfs/jps558/w2k/Desktop/R_WorkingDirectory") 


#Set column names
qbedcol_names = c("Chr",	"Start",	"End",	"Number of Insertions",	"Strand",	"Barcode")

pyCCcol_names = c("Chr",	"Start",	"End",	"Center",	"Experiment Insertions",	"Background insertions",	"Reference Insertions",	"pvalue Reference",	"pvalue Background",	"Fraction Experiment",	"TPH Experiment",	"Fraction background",	"TPH background",	"TPH background subtracted",	"pvalue_adj Reference")

# Load data
ER_peak2 = read.table("Normalising replicates/Peaks/peak_data_ER_WT2.bed", header = TRUE, sep = '\t', col.names = pyCCcol_names)
ER_peak4 = read.table("Normalising replicates/Peaks/peak_data_ER_WT4.bed", header = TRUE, sep = '\t', col.names = pyCCcol_names)

ER_Rep1_4 = read.table("Normalising replicates/Peaks/a_filtered_sortedER1_EKDL230008858.qbed_b_peakset4_.bed", header = TRUE, sep = '\t', col.names = qbedcol_names)
ER_Rep2_4 = read.table("Normalising replicates/Peaks/a_filtered_sortedER2_EKDL230008858.qbed_b_peakset4_.bed", header = TRUE, sep = '\t', col.names = qbedcol_names)
ER_Rep3_4 = read.table("Normalising replicates/Peaks/a_filtered_sortedER3_EKDL230008858.qbed_b_peakset4_.bed", header = TRUE, sep = '\t', col.names = qbedcol_names)
ER_Rep5_4 = read.table("Normalising replicates/Peaks/a_filtered_sortedER5_EKDL230008858.qbed_b_peakset4_.bed", header = TRUE, sep = '\t', col.names = qbedcol_names)
ER_Rep6_4 = read.table("Normalising replicates/Peaks/a_filtered_sortedER6_EKDL230008858.qbed_b_peakset4_.bed", header = TRUE, sep = '\t', col.names = qbedcol_names)

#Subset data

ER_peak2_subset <- ER_peak2[, c(1, 2, 3, 5)]

ER_peak4_subset <- ER_peak4[, c(1, 2, 3, 5)]

ER_Rep1_4_subset = ER_Rep1_4[, 1:4]
ER_Rep2_4_subset = ER_Rep2_4[, 1:4]
ER_Rep3_4_subset = ER_Rep3_4[, 1:4]
ER_Rep5_4_subset = ER_Rep5_4[, 1:4]
ER_Rep6_4_subset = ER_Rep6_4[, 1:4]

#Combine Replicate rows of rep 1 and 2
#rep1
ER_Rep1_4_subset_combined <- ER_Rep1_4_subset %>%
  group_by(Chr, Start, End) %>%
  summarize(
    Number.of.Insertions = sum(Number.of.Insertions, na.rm = TRUE)  # Sum insertions
  ) %>%
  ungroup()  # Remove grouping structure

# View the result
print(ER_Rep1_4_subset_combined)

#Rep2
ER_Rep2_4_subset_combined <- ER_Rep2_4_subset %>%
  group_by(Chr, Start, End) %>%
  summarize(
    Number.of.Insertions = sum(Number.of.Insertions, na.rm = TRUE)  # Sum insertions
  ) %>%
  ungroup()  # Remove grouping structure

# View the result
print(ER_Rep2_4_subset_combined)

# Create IRanges objects for peaks and replicates
ir_peak4 <- IRanges(start = ER_peak4_subset$Start, end = ER_peak4_subset$End)

ir_rep1 <- IRanges(start = ER_Rep1_4_subset_combined$Start, end = ER_Rep1_4_subset_combined$End)
ir_rep2 <- IRanges(start = ER_Rep2_4_subset_combined$Start, end = ER_Rep2_4_subset_combined$End)

# Find overlaps for replicate 1
overlap_peak4_rep1 <- findOverlaps(ir_peak4, ir_rep1, maxgap = 1000)

# Initialize a vector with 0 for Number.of.Insertions
rep1_insertions <- numeric(length(ir_peak4))
rep1_insertions[unique(queryHits(overlap_peak4_rep1))] <- tapply(
  ER_Rep1_4_subset_combined$Number.of.Insertions[subjectHits(overlap_peak4_rep1)],
  queryHits(overlap_peak4_rep1),
  sum
)

# Create summary for replicate 1
rep1_summary <- data.frame(
  Chr = ER_peak4_subset$Chr,
  Start = ER_peak4_subset$Start,
  End = ER_peak4_subset$End,
  Number.of.Insertions_Rep1 = rep1_insertions
)

# Find overlaps for replicate 2
overlap_peak4_rep2 <- findOverlaps(ir_peak4, ir_rep2, maxgap = 1000)

# Initialize a vector with 0 for Number.of.Insertions
rep2_insertions <- numeric(length(ir_peak4))
rep2_insertions[unique(queryHits(overlap_peak4_rep2))] <- tapply(
  ER_Rep2_4_subset_combined$Number.of.Insertions[subjectHits(overlap_peak4_rep2)],
  queryHits(overlap_peak4_rep2),
  sum
)

# Create summary for replicate 2
rep2_summary <- data.frame(
  Chr = ER_peak4_subset$Chr,
  Start = ER_peak4_subset$Start,
  End = ER_peak4_subset$End,
  Number.of.Insertions_Rep2 = rep2_insertions
)

# Merge the results from replicate 1 and replicate 2 with the peak data
final_data <- merge(rep1_summary, rep2_summary, by = c("Chr", "Start", "End"), all = TRUE)

# View the final data
print(final_data)

#Find Total insertions

peak4_Rep1_total = sum(final_data$Number.of.Insertions_Rep1)
peak4_Rep2_total = sum(final_data$Number.of.Insertions_Rep2)

# Filter out rows where any Number.of.Insertions columns have a 0
final_data_filtered <- final_data %>%
  filter(Number.of.Insertions_Rep1 != 0 & Number.of.Insertions_Rep2 != 0)

#calculate IPM
final_data_filtered$IPM_Rep1 = (final_data_filtered$Number.of.Insertions_Rep1/peak4_Rep1_total)*1e6
final_data_filtered$IPM_Rep2 = (final_data_filtered$Number.of.Insertions_Rep2/peak4_Rep2_total)*1e6

#Convert to log 10

final_data_filtered$log10_IPM_Rep1 = log10(final_data_filtered$IPM_Rep1)
final_data_filtered$log10_IPM_Rep2 = log10(final_data_filtered$IPM_Rep2)


#Find correlation
correlation1_4 <- cor(final_data_filtered$IPM_Rep1, final_data_filtered$IPM_Rep2, method = "pearson")
correlation1_4log <- cor(final_data_filtered$log10_IPM_Rep1, final_data_filtered$log10_IPM_Rep2, method = "pearson")


#Plot IPM
ggplot(final_data_filtered, aes(x = IPM_Rep1, y = IPM_Rep2)) +
  geom_point(color = "blue", size = 3) +  # Scatter plot points
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +  # Line of unity
  annotate("text", x = max(final_data_filtered$IPM_Rep1) * 0.5, 
           y = max(final_data_filtered$IPM_Rep2) * 0.9, 
           label = paste("Pearson r =", round(correlation1_4, 2)), 
           color = "black", size = 5) +  # Add correlation text
  labs(
    title = "Reproducibility of Insertions at MACCs peaks from Replicates 1 and 2 (IPM)",
    x = "Replicate 1 (IPM)",
    y = "Replicate 2 (IPM)"
  ) +
  theme_minimal()  # Minimal theme

#Plot log10(IPM)
ggplot(final_data_filtered, aes(x = log10_IPM_Rep1, y = log10_IPM_Rep2)) +
  geom_point(color = "blue", size = 3) +  # Scatter plot points
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +  # Line of unity
  annotate("text", x = max(final_data_filtered$log10_IPM_Rep1) * 0.5, 
           y = max(final_data_filtered$log10_IPM_Rep2) * 0.9, 
           label = paste("Pearson r =", round(correlation1_4log, 2)), 
           color = "black", size = 5) +  # Add correlation text
  labs(
    title = "Reproducibility of Insertions at MACCs peaks from Replicates 1 and 2 (log10 IPM)",
    x = "Replicate 1 (log10 IPM)",
    y = "Replicate 2 (log10 IPM)"
  ) +
  theme_light() +
  coord_cartesian(xlim = c(0, 5), ylim = c(0, 5), expand = FALSE) # Minimal theme

