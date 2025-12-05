#Install and load packages
#install.packages("ROCR")

library(ROCR)
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

ESR1_pyCC_4_Encode2_q = read.table("results/1kbp_A_peak_data_ER_WT4.bed_B_Encode_ER_rep2_ChIP_q0.05_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)

#Change overlap values >0 to 1
ESR1_pyCC_2_Andy_q <- ESR1_pyCC_2_Andy_q %>%
  mutate(Number.of.Overlaps.with.ChIP.Seq.Peaks = ifelse(Number.of.Overlaps.with.ChIP.Seq.Peaks > 0, 1, Number.of.Overlaps.with.ChIP.Seq.Peaks))

ESR1_pyCC_2_Encode2_q <- ESR1_pyCC_2_Encode2_q %>%
  mutate(Number.of.Overlaps.with.ChIP.Seq.Peaks = ifelse(Number.of.Overlaps.with.ChIP.Seq.Peaks > 0, 1, Number.of.Overlaps.with.ChIP.Seq.Peaks))

ESR1_pyCC_4_Andy_q <- ESR1_pyCC_4_Andy_q %>%
  mutate(Number.of.Overlaps.with.ChIP.Seq.Peaks = ifelse(Number.of.Overlaps.with.ChIP.Seq.Peaks > 0, 1, Number.of.Overlaps.with.ChIP.Seq.Peaks))

ESR1_pyCC_4_Encode2_q <- ESR1_pyCC_4_Encode2_q %>%
  mutate(Number.of.Overlaps.with.ChIP.Seq.Peaks = ifelse(Number.of.Overlaps.with.ChIP.Seq.Peaks > 0, 1, Number.of.Overlaps.with.ChIP.Seq.Peaks))

# View the updated dataframe
print(ESR1_pyCC_2_Andy_q)
print(ESR1_pyCC_2_Encode2_q)
print(ESR1_pyCC_4_Andy_q)
print(ESR1_pyCC_4_Encode2_q)

#Using ROCR for peakset 2
#Andy

#Create Prediciton object
#Cannot use pvalue, must use predicted probablility - e.g. Pvalue of 0.0009 has a predicted probaility of 0.9991 as it should be very likely this is a true value

predicted_probability2_G <- (1- ESR1_pyCC_2_Andy_q$pvalue_adj.Reference)

predicted_probability2_E <- (1- ESR1_pyCC_2_Encode2_q$pvalue_adj.Reference)

pred2_G <- prediction(predicted_probability2_G, ESR1_pyCC_2_Andy_q$Number.of.Overlaps.with.ChIP.Seq.Peaks)

pred2_E <- prediction(predicted_probability2_E, ESR1_pyCC_2_Encode2_q$Number.of.Overlaps.with.ChIP.Seq.Peaks)

#Create ROC curve

ROC2_G <- performance(pred2_G, measure = "tpr", x.measure = "fpr")

ROC2_E <- performance(pred2_E, measure = "tpr", x.measure = "fpr")

#plot GEO
plot(ROC2_G, main = "ROC curve for ESR1-HyPB Peakset 2 overlap with GSE109820", col = "blue", , lwd=2,  xaxs="i", yaxs="i")

abline(a = 0, b = 1, col = "red", lty = 2)

#Calculate and print AUC
ROC2_G_auc <- performance(pred2_G, measure = "auc")
print(paste("AUC:", ROC2_G_auc@y.values[[1]]))
text(0.7, 0.2, paste("AUC = ", round(ROC2_G_auc@y.values[[1]], 2)), col = "black", cex = 1.5)

#plot Encode
plot(ROC2_E, main = "ROC curve for ESR1-HyPB Peakset 2 overlap with ENCFF063JMY", col = "blue", , lwd=2,  xaxs="i", yaxs="i")

abline(a = 0, b = 1, col = "red", lty = 2)

#Calculate and print AUC
ROC2_E_auc <- performance(pred2_E, measure = "auc")
print(paste("AUC:", ROC2_E_auc@y.values[[1]]))
text(0.7, 0.2, paste("AUC = ", round(ROC2_E_auc@y.values[[1]], 2)), col = "black", cex = 1.5)

#Using ROCR for peakset 4
#Andy

#Create Prediciton object
#Cannot use pvalue, must use predicted probablility - e.g. Pvalue of 0.0009 has a predicted probaility of 0.9991 as it should be very likely this is a true value

predicted_probability4_G <- (1- ESR1_pyCC_4_Andy_q$pvalue_adj.Reference)

predicted_probability4_E <- (1- ESR1_pyCC_4_Encode2_q$pvalue_adj.Reference)

pred4_G <- prediction(predicted_probability4_G, ESR1_pyCC_4_Andy_q$Number.of.Overlaps.with.ChIP.Seq.Peaks)

pred4_E <- prediction(predicted_probability4_E, ESR1_pyCC_4_Encode2_q$Number.of.Overlaps.with.ChIP.Seq.Peaks)

#Create ROC curve

ROC4_G <- performance(pred4_G, measure = "tpr", x.measure = "fpr")

ROC4_E <- performance(pred4_E, measure = "tpr", x.measure = "fpr")

#plot GEO
plot(ROC4_G, main = "ROC curve for ESR1-HyPB Peakset 4 overlap with GSE109820", col = "blue", , lwd=2,  xaxs="i", yaxs="i")

abline(a = 0, b = 1, col = "red", lty = 2)

#Calculate and print AUC
ROC4_G_auc <- performance(pred4_G, measure = "auc")
print(paste("AUC:", ROC4_G_auc@y.values[[1]]))
text(0.7, 0.2, paste("AUC = ", round(ROC4_G_auc@y.values[[1]], 2)), col = "black", cex = 1.5)

#plot Encode
plot(ROC4_E, main = "ROC curve for ESR1-HyPB Peakset 4 overlap with ENCFF063JMY", col = "blue", , lwd=2,  xaxs="i", yaxs="i")

abline(a = 0, b = 1, col = "red", lty = 2)

#Calculate and print AUC
ROC4_E_auc <- performance(pred4_E, measure = "auc")
print(paste("AUC:", ROC4_E_auc@y.values[[1]]))
text(0.7, 0.2, paste("AUC = ", round(ROC4_E_auc@y.values[[1]], 2)), col = "black", cex = 1.5)



