#Install and load packages
#install.packages("ROCR")

library(ROCR)
library(dplyr)
library(ggplot2)
library(ggthemes)

#Set working directory
setwd("C:/Users/jps558/OneDrive - University of York/Desktop/R_WorkingDirectory/R_WFHDirectory") 


#Set column names
pyCCcol_names = c("Chr",	"Start",	"End",	"Center",	"Experiment Insertions",	"Background insertions",	"Reference Insertions",	"pvalue Reference",	"pvalue Background",	"Fraction Experiment",	"TPH Experiment",	"Fraction background",	"TPH background",	"TPH background subtracted",	"pvalue_adj Reference", "Number of Overlaps with ChIA-PET Peaks")
ChIP_names = c("Chr",	"Start",	"End",	"Name",	"Score",	"Strand",	"Signal Value",	"pvalue",	"qvalue",	"Peak", "Number of Overlaps with ChIA-PET Peaks")

##Read data
#Variables named <Window '-a'><Window'-b'>
#Encode Rep 2 was taken as it has more peaks than replicate 1 - even though there is less overlap with my data
ESR1_pyCC_2_anchor1 = read.table("overlap/cc_Chia/count_Chia_peak_data_ER_WT2_Chiapet_hg38_combined.bed_anchor1.txt", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_2_anchor2 = read.table("overlap/cc_Chia/count_Chia_peak_data_ER_WT2_Chiapet_hg38_combined.bed_anchor2.txt", header = TRUE, sep = "\t", col.names = pyCCcol_names)


ESR1_pyCC_4_anchor1 = read.table("overlap/cc_Chia/count_Chia_peak_data_ER_WT4_Chiapet_hg38_combined.bed_anchor1.txt", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_4_anchor2 = read.table("overlap/cc_Chia/count_Chia_peak_data_ER_WT4_Chiapet_hg38_combined.bed_anchor2.txt", header = TRUE, sep = "\t", col.names = pyCCcol_names)

#Counting the total reads

ESR1_pyCC_2_total=nrow(ESR1_pyCC_2_anchor1)

ESR1_pyCC_4_total=nrow(ESR1_pyCC_4_anchor1)

#Count the number of reads with no overlap
ESR1_pyCC_2_anchor1_overlap = as.numeric(length(which(ESR1_pyCC_2_anchor1$Number.of.Overlaps.with.ChIA.PET.Peaks==0)))
ESR1_pyCC_2_anchor2_overlap = as.numeric(length(which(ESR1_pyCC_2_anchor2$Number.of.Overlaps.with.ChIA.PET.Peaks==0)))



ESR1_pyCC_4_anchor1_overlap = as.numeric(length(which(ESR1_pyCC_4_anchor1$Number.of.Overlaps.with.ChIA.PET.Peaks==0)))
ESR1_pyCC_4_anchor2_overlap = as.numeric(length(which(ESR1_pyCC_4_anchor2$Number.of.Overlaps.with.ChIA.PET.Peaks==0)))

#Determine overlap percentage
Percentage_ESR1_pyCC_2_anchor1 =(1-(ESR1_pyCC_2_anchor1_overlap/ESR1_pyCC_2_total))*100
Percentage_ESR1_pyCC_2_anchor2 =(1-(ESR1_pyCC_2_anchor2_overlap/ESR1_pyCC_2_total))*100



Percentage_ESR1_pyCC_4_anchor1 =(1-(ESR1_pyCC_4_anchor1_overlap/ESR1_pyCC_4_total))*100
Percentage_ESR1_pyCC_4_anchor2 = (1-(ESR1_pyCC_4_anchor2_overlap/ESR1_pyCC_4_total))*100

#### plot across both anchors
# Merge anchor1 and anchor2 into one dataframe per replicate
ESR1_pyCC_2_combo <- data.frame(
  Anchor1 = ESR1_pyCC_2_anchor1$Number.of.Overlaps.with.ChIA.PET.Peaks,
  Anchor2 = ESR1_pyCC_2_anchor2$Number.of.Overlaps.with.ChIA.PET.Peaks
)

ESR1_pyCC_4_combo <- data.frame(
  Anchor1 = ESR1_pyCC_4_anchor1$Number.of.Overlaps.with.ChIA.PET.Peaks,
  Anchor2 = ESR1_pyCC_4_anchor2$Number.of.Overlaps.with.ChIA.PET.Peaks
)



# Peak overlaps if Anchor1 > 0 OR Anchor2 > 0
ESR1_pyCC_2_combo$AnyAnchor <- (ESR1_pyCC_2_combo$Anchor1 > 0 | ESR1_pyCC_2_combo$Anchor2 > 0)
ESR1_pyCC_4_combo$AnyAnchor <- (ESR1_pyCC_4_combo$Anchor1 > 0 | ESR1_pyCC_4_combo$Anchor2 > 0)

# Add AnyAnchor column to original peak data
ESR1_pyCC_2_Chiacombined <- ESR1_pyCC_2_anchor1
ESR1_pyCC_2_Chiacombined$AnyAnchor <- ESR1_pyCC_2_combo$AnyAnchor

ESR1_pyCC_4_Chiacombined <- ESR1_pyCC_4_anchor1
ESR1_pyCC_4_Chiacombined$AnyAnchor <- ESR1_pyCC_4_combo$AnyAnchor

ESR1_pyCC_2_Chiacombined$AnyAnchor <- as.numeric(ESR1_pyCC_2_Chiacombined$AnyAnchor)
ESR1_pyCC_4_Chiacombined$AnyAnchor <- as.numeric(ESR1_pyCC_4_Chiacombined$AnyAnchor)

# Convert p-values to predicted probabilities
predicted_prob_2 <- 1 - ESR1_pyCC_2_Chiacombined$pvalue_adj.Reference
predicted_prob_4 <- 1 - ESR1_pyCC_4_Chiacombined$pvalue_adj.Reference

# Create prediction objects
pred2 <- prediction(predicted_prob_2, ESR1_pyCC_2_Chiacombined$AnyAnchor)
pred4 <- prediction(predicted_prob_4, ESR1_pyCC_4_Chiacombined$AnyAnchor)

# Generate ROC curves
ROC2 <- performance(pred2, measure = "tpr", x.measure = "fpr")
ROC4 <- performance(pred4, measure = "tpr", x.measure = "fpr")

# Plot ROC curve for Peakset 2
plot(ROC2, main = "ROC Curve: Peakset 2 vs ChIA-PET Anchors", col = "blue", lwd = 2)
abline(a = 0, b = 1, col = "red", lty = 2)
auc2 <- performance(pred2, measure = "auc")
text(0.7, 0.2, paste("AUC = ", round(auc2@y.values[[1]], 2)), col = "black", cex = 1.5)

# Plot ROC curve for Peakset 4
plot(ROC4, main = "ROC Curve: Peakset 4 vs ChIA-PET Anchors", col = "blue", lwd = 2)
abline(a = 0, b = 1, col = "red", lty = 2)
auc4 <- performance(pred4, measure = "auc")
text(0.7, 0.2, paste("AUC = ", round(auc4@y.values[[1]], 2)), col = "black", cex = 1.5)