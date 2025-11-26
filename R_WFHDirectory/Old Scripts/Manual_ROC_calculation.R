#Load packages
library(dplyr)

#Set working directory
setwd("//userfs/jps558/w2k/Desktop/R_WorkingDirectory") 


#Set column names
pyCCcol_names = c("Chr",	"Start",	"End",	"Center",	"Experiment Insertions",	"Background insertions",	"Reference Insertions",	"pvalue Reference",	"pvalue Background",	"Fraction Experiment",	"TPH Experiment",	"Fraction background",	"TPH background",	"TPH background subtracted",	"pvalue_adj Reference", "Number of Overlaps with ChIP-Seq Peaks")

##Read data
#Variables named <Window '-a'><Window'-b'>
ESR1_pyCC_2_Andy_q = read.table("results/1kbp_A_peak_data_ER_WT2.bed_B_Andy_ChIP_q0.05_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)

ESR1_pyCC_2_Encode2_q = read.table("results/1kbp_A_peak_data_ER_WT2.bed_B_Encode_ER_rep2_ChIP_q0.05_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)

ESR1_pyCC_4_Andy_q = read.table("results/1kbp_A_peak_data_ER_WT4.bed_B_Andy_ChIP_q0.05_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)

#Change overlap values >0 to 1
ESR1_pyCC_2_Andy_q <- ESR1_pyCC_2_Andy_q %>%
  mutate(Number.of.Overlaps.with.ChIP.Seq.Peaks = ifelse(Number.of.Overlaps.with.ChIP.Seq.Peaks > 0, 1, Number.of.Overlaps.with.ChIP.Seq.Peaks))

ESR1_pyCC_2_Encode2_q <- ESR1_pyCC_2_Encode2_q %>%
  mutate(Number.of.Overlaps.with.ChIP.Seq.Peaks = ifelse(Number.of.Overlaps.with.ChIP.Seq.Peaks > 0, 1, Number.of.Overlaps.with.ChIP.Seq.Peaks))

ESR1_pyCC_4_Andy_q <- ESR1_pyCC_4_Andy_q %>%
  mutate(Number.of.Overlaps.with.ChIP.Seq.Peaks = ifelse(Number.of.Overlaps.with.ChIP.Seq.Peaks > 0, 1, Number.of.Overlaps.with.ChIP.Seq.Peaks))

# Example data PyCC_2_Andy
true_labels <- ESR1_pyCC_2_Andy_q$Number.of.Overlaps.with.ChIP.Seq.Peaks  # True labels (1 for overlap, 0 for no overlap)
adjusted_p_values <- ESR1_pyCC_2_Andy_q$pvalue_adj.Reference  # Adjusted p-values

# Function to calculate TP, FP, TN, FN for a given threshold
calculate_metrics <- function(true_labels, p_values, threshold) {
  predicted_labels <- ifelse(p_values < threshold, 1, 0)  # Classify based on threshold
  
  TP <- sum(predicted_labels == 1 & true_labels == 1)  # True Positives
  FP <- sum(predicted_labels == 1 & true_labels == 0)  # False Positives
  TN <- sum(predicted_labels == 0 & true_labels == 0)  # True Negatives
  FN <- sum(predicted_labels == 0 & true_labels == 1)  # False Negatives
  
  return(c(TP = TP, FP = FP, TN = TN, FN = FN))
}

# Define thresholds to evaluate
thresholds <- seq(0, 1, by = 1e-3)  # Thresholds from 0 to 1 in steps of 0.1

# Initialize an empty list to store results
results <- data.frame()

# Iterate through thresholds and calculate metrics for each
for (threshold in thresholds) {
  metrics <- calculate_metrics(true_labels, adjusted_p_values, threshold)
  results <- rbind(results, c(threshold, metrics))  # Store the result for each threshold
}

# Name the columns for clarity
colnames(results) <- c("Threshold", "TP", "FP", "TN", "FN")

# Print the results
print(results)
