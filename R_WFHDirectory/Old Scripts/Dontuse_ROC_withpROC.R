#Install packages
if (! requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("xrobin/pROC@develop")

library(pROC)

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


# View the updated dataframe
print(ESR1_pyCC_2_Andy_q)
print(ESR1_pyCC_2_Encode2_q)
print(ESR1_pyCC_4_Andy_q)

#USING pROC with pvaule

#Set up ROC for Peakset 2 with Andy ChIP

Set2_Andy_ROC <- roc(ESR1_pyCC_2_Andy_q$Number.of.Overlaps.with.ChIP.Seq.Peaks, ESR1_pyCC_2_Andy_q$pvalue_adj.Reference)

#plot
plot(Set2_Andy_ROC, main="Peakset2 ROC curve with Andy ChIP", lwd=2,  xaxs="i", yaxs="i")

#Using pROC with predicted probablility
predicted_probability <- (1- ESR1_pyCC_2_Andy_q$pvalue_adj.Reference)

probabilityROC <- roc(ESR1_pyCC_2_Andy_q$Number.of.Overlaps.with.ChIP.Seq.Peaks, predicted_probability)

plot(probabilityROC, main="Peakset2 ROC curve with Andy ChIP", lwd=2,  xaxs="i", yaxs="i")

