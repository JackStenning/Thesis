#Load packages

install.packages("pracma")
library(pracma)
library(dplyr)

#Set working directory
setwd("//userfs/jps558/w2k/Desktop/R_WorkingDirectory") 

#Set column names
pyCCcol_names = c("Chr",	"Start",	"End",	"Center",	"Experiment Insertions",	"Background insertions",	"Reference Insertions",	"pvalue Reference",	"pvalue Background",	"Fraction Experiment",	"TPH Experiment",	"Fraction background",	"TPH background",	"TPH background subtracted",	"pvalue_adj Reference", "Number of Overlaps with ChIP-Seq Peaks")

##Read data
#Variables named <Window '-a'><Window'-b'>
ESR1_pyCC_2_Andy_q = read.table("results/1kbp_A_peak_data_ER_WT2.bed_B_Andy_ChIP_q0.05_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)

#THE TWO BELOW HAVE NOT CURRENTLY BEEN PLOTTED
ESR1_pyCC_2_Encode2_q = read.table("results/1kbp_A_peak_data_ER_WT2.bed_B_Encode_ER_rep2_ChIP_q0.05_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)

ESR1_pyCC_4_Andy_q = read.table("results/1kbp_A_peak_data_ER_WT4.bed_B_Andy_ChIP_q0.05_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)

#Change overlap values >0 to 1
ESR1_pyCC_2_Andy_q <- ESR1_pyCC_2_Andy_q %>%
  mutate(Number.of.Overlaps.with.ChIP.Seq.Peaks = ifelse(Number.of.Overlaps.with.ChIP.Seq.Peaks > 0, 1, Number.of.Overlaps.with.ChIP.Seq.Peaks))

ESR1_pyCC_2_Encode2_q <- ESR1_pyCC_2_Encode2_q %>%
  mutate(Number.of.Overlaps.with.ChIP.Seq.Peaks = ifelse(Number.of.Overlaps.with.ChIP.Seq.Peaks > 0, 1, Number.of.Overlaps.with.ChIP.Seq.Peaks))

ESR1_pyCC_4_Andy_q <- ESR1_pyCC_4_Andy_q %>%
  mutate(Number.of.Overlaps.with.ChIP.Seq.Peaks = ifelse(Number.of.Overlaps.with.ChIP.Seq.Peaks > 0, 1, Number.of.Overlaps.with.ChIP.Seq.Peaks))

# Example data (you can replace this with your actual data)
true_labels <- ESR1_pyCC_2_Andy_q$Number.of.Overlaps.with.ChIP.Seq.Peaks # True labels (1 for overlap, 0 for no overlap)
adjusted_p_values <- ESR1_pyCC_2_Andy_q$pvalue_adj.Reference  # Adjusted p-values

# Function to calculate TP, FP, TN, FN for a given threshold
calculate_metrics <- function(true_labels, p_values, threshold) {
  # Set p-values that are essentially zero to a small value (e.g., 1e-100)
  p_values[p_values == 0] <- 1e-30  # Replace zero p-values with a very small value
  
  # Classify based on threshold
  predicted_labels <- ifelse(p_values < threshold, 1, 0)  # 1 if below threshold, else 0
  
  TP <- sum(predicted_labels == 1 & true_labels == 1)  # True Positives
  FP <- sum(predicted_labels == 1 & true_labels == 0)  # False Positives
  FN <- sum(predicted_labels == 0 & true_labels == 1)  # False Negatives
  TN <- sum(predicted_labels == 0 & true_labels == 0)  # True Negatives
  
  return(c(TP = TP, FP = FP, TN = TN, FN = FN))
}

# Logarithmic thresholds from very small to 1
log_thresholds <- 10^seq(-30, 0, by = 0.02)  # Logarithmic scale with steps of 0.2

# Initialize an empty data frame to store results
results <- data.frame()

# Iterate through thresholds and calculate metrics
for (threshold in log_thresholds) {
  metrics <- calculate_metrics(true_labels, adjusted_p_values, threshold)
  results <- rbind(results, c(threshold, metrics))  # Store results for each threshold
}

# Set column names
colnames(results) <- c("Threshold", "TP", "FP", "TN", "FN")

# Calculate True Positive Rate (TPR) and False Positive Rate (FPR)
results$TPR <- results$TP / (results$TP + results$FN)  # Sensitivity
results$FPR <- results$FP / (results$FP + results$TN)  # 1 - Specificity

# Plot ROC curve
plot(results$FPR, results$TPR, type = "l", col = "blue", xlab = "False Positive Rate", ylab = "True Positive Rate", main = "ROC Curve",  xaxs="i", yaxs="i")

# Optionally, add the random classifier line (AUC = 0.5)
abline(a = 0, b = 1, col = "red", lty = 2)

# Calculate AUC using the trapezoidal rule
auc_value <- trapz(results$FPR, results$TPR)

# Add AUC value as text to the plot
text(0.7, 0.2, paste("AUC =", round(auc_value, 3)), col = "black", cex = 1.2)
