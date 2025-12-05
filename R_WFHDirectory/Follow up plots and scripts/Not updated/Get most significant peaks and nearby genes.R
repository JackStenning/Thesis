#Set working directory
setwd("//userfs/jps558/w2k/Desktop/R_WorkingDirectory") 

#Load data
ESR1_pyCC_2_anno = read.table("Plots and Scripts/peakset_2/Peakset2_annotated.csv", header = TRUE, sep = ",")

ESR1_pyCC_4_anno = read.table("Plots and Scripts/peakset_4/Peakset4_annotated.csv", header = TRUE, sep = ",")

#Re-order on pvalue
ESR1_pyCC_2_anno = ESR1_pyCC_2_anno[order(ESR1_pyCC_2_anno$pvalue_adj.Reference), ]

ESR1_pyCC_4_anno = ESR1_pyCC_4_anno[order(ESR1_pyCC_4_anno$pvalue_adj.Reference), ]

#Write 

write.csv(ESR1_pyCC_2_anno, "Plots and Scripts/peakset_2/Peakset2_annotated_By_Pvalue.csv", row.names = FALSE)

write.csv(ESR1_pyCC_4_anno, "Plots and Scripts/peakset_4/Peakset4_annotated_By_Pvalue.csv", row.names = FALSE)