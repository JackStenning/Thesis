############# LOOK HERE ######################
#add or delete the '#' symbol below on install packages
#install.packages("Rtools")
#install.packages("ggthemes")

#DELTE AFTER WRITING

#Load software
library(ggplot2)
library(ggthemes)

#Set working directory
setwd("//userfs/jps558/w2k/Desktop/R_WorkingDirectory") 


#Set column names
pyCCcol_names = c("Chr",	"Start",	"End",	"Center",	"Experiment Insertions",	"Background insertions",	"Reference Insertions",	"pvalue Reference",	"pvalue Background",	"Fraction Experiment",	"TPH Experiment",	"Fraction background",	"TPH background",	"TPH background subtracted",	"pvalue_adj Reference", "Number of Overlaps with ChIP-Seq Peaks")

##Read data
#Variables named <Window '-a'><Window'-b'>
ESR1_pyCC_2_GSE_4 = read.table("results/1kbp_A_peak_data_ER_WT2.bed_B_Andy_ChIP_1e-4_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_2_GSE_6 = read.table("results/1kbp_A_peak_data_ER_WT2.bed_B_Andy_ChIP_1e-6_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_2_GSE_9 = read.table("results/1kbp_A_peak_data_ER_WT2.bed_B_Andy_ChIP_1e-9_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_2_GSE_q = read.table("results/1kbp_A_peak_data_ER_WT2.bed_B_Andy_ChIP_q0.05_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_2_Encode1_4 = read.table("results/1kbp_A_peak_data_ER_WT2.bed_B_Encode_ER_rep1_ChIP_1e-4_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_2_Encode1_6 = read.table("results/1kbp_A_peak_data_ER_WT2.bed_B_Encode_ER_rep1_ChIP_1e-6_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_2_Encode1_9 = read.table("results/1kbp_A_peak_data_ER_WT2.bed_B_Encode_ER_rep1_ChIP_1e-9_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_2_Encode1_q = read.table("results/1kbp_A_peak_data_ER_WT2.bed_B_Encode_ER_rep1_ChIP_q0.05_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_2_Encode2_4 = read.table("results/1kbp_A_peak_data_ER_WT2.bed_B_Encode_ER_rep2_ChIP_1e-4_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_2_Encode2_6 = read.table("results/1kbp_A_peak_data_ER_WT2.bed_B_Encode_ER_rep2_ChIP_1e-6_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_2_Encode2_9 = read.table("results/1kbp_A_peak_data_ER_WT2.bed_B_Encode_ER_rep2_ChIP_1e-9_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_2_Encode2_q = read.table("results/1kbp_A_peak_data_ER_WT2.bed_B_Encode_ER_rep2_ChIP_q0.05_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)


ESR1_pyCC_4_GSE_4 = read.table("results/1kbp_A_peak_data_ER_WT4.bed_B_Andy_ChIP_1e-4_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_4_GSE_6 = read.table("results/1kbp_A_peak_data_ER_WT4.bed_B_Andy_ChIP_1e-6_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_4_GSE_9 = read.table("results/1kbp_A_peak_data_ER_WT4.bed_B_Andy_ChIP_1e-9_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_4_GSE_q = read.table("results/1kbp_A_peak_data_ER_WT4.bed_B_Andy_ChIP_q0.05_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_4_Encode1_4 = read.table("results/1kbp_A_peak_data_ER_WT4.bed_B_Encode_ER_rep1_ChIP_1e-4_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_4_Encode1_6 = read.table("results/1kbp_A_peak_data_ER_WT4.bed_B_Encode_ER_rep1_ChIP_1e-6_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_4_Encode1_9 = read.table("results/1kbp_A_peak_data_ER_WT4.bed_B_Encode_ER_rep1_ChIP_1e-9_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_4_Encode1_q = read.table("results/1kbp_A_peak_data_ER_WT4.bed_B_Encode_ER_rep1_ChIP_q0.05_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_4_Encode2_4 = read.table("results/1kbp_A_peak_data_ER_WT4.bed_B_Encode_ER_rep2_ChIP_1e-4_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_4_Encode2_6 = read.table("results/1kbp_A_peak_data_ER_WT4.bed_B_Encode_ER_rep2_ChIP_1e-6_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_4_Encode2_9 = read.table("results/1kbp_A_peak_data_ER_WT4.bed_B_Encode_ER_rep2_ChIP_1e-9_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_4_Encode2_q = read.table("results/1kbp_A_peak_data_ER_WT4.bed_B_Encode_ER_rep2_ChIP_q0.05_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)

#Counting the total reads
ESR1_pyCC_2_total=nrow(ESR1_pyCC_2_GSE_4)

ESR1_pyCC_4_total=nrow(ESR1_pyCC_4_GSE_4)

#Count the number of reads with no overlap
ESR1_pyCC_2_GSE_4_overlap = as.numeric(length(which(ESR1_pyCC_2_GSE_4$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_2_GSE_6_overlap = as.numeric(length(which(ESR1_pyCC_2_GSE_6$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_2_GSE_9_overlap = as.numeric(length(which(ESR1_pyCC_2_GSE_9$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_2_GSE_q_overlap = as.numeric(length(which(ESR1_pyCC_2_GSE_q$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_2_Encode1_4_overlap = as.numeric(length(which(ESR1_pyCC_2_Encode1_4$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_2_Encode1_6_overlap = as.numeric(length(which(ESR1_pyCC_2_Encode1_6$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_2_Encode1_9_overlap = as.numeric(length(which(ESR1_pyCC_2_Encode1_9$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_2_Encode1_q_overlap = as.numeric(length(which(ESR1_pyCC_2_Encode1_q$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_2_Encode2_4_overlap = as.numeric(length(which(ESR1_pyCC_2_Encode2_4$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_2_Encode2_6_overlap = as.numeric(length(which(ESR1_pyCC_2_Encode2_6$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_2_Encode2_9_overlap = as.numeric(length(which(ESR1_pyCC_2_Encode2_9$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_2_Encode2_q_overlap = as.numeric(length(which(ESR1_pyCC_2_Encode2_q$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))

ESR1_pyCC_4_GSE_4_overlap = as.numeric(length(which(ESR1_pyCC_4_GSE_4$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_4_GSE_6_overlap = as.numeric(length(which(ESR1_pyCC_4_GSE_6$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_4_GSE_9_overlap = as.numeric(length(which(ESR1_pyCC_4_GSE_9$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_4_GSE_q_overlap = as.numeric(length(which(ESR1_pyCC_4_GSE_q$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_4_Encode1_4_overlap = as.numeric(length(which(ESR1_pyCC_4_Encode1_4$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_4_Encode1_6_overlap = as.numeric(length(which(ESR1_pyCC_4_Encode1_6$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_4_Encode1_9_overlap = as.numeric(length(which(ESR1_pyCC_4_Encode1_9$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_4_Encode1_q_overlap = as.numeric(length(which(ESR1_pyCC_4_Encode1_q$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_4_Encode2_4_overlap = as.numeric(length(which(ESR1_pyCC_4_Encode2_4$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_4_Encode2_6_overlap = as.numeric(length(which(ESR1_pyCC_4_Encode2_6$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_4_Encode2_9_overlap = as.numeric(length(which(ESR1_pyCC_4_Encode2_9$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_4_Encode2_q_overlap = as.numeric(length(which(ESR1_pyCC_4_Encode2_q$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))

#Determine overlap percentage
Percentage_ESR1_pyCC_2_GSE_4 =(1-(ESR1_pyCC_2_GSE_4_overlap/ESR1_pyCC_2_total))*100
Percentage_ESR1_pyCC_2_GSE_6 =(1-(ESR1_pyCC_2_GSE_6_overlap/ESR1_pyCC_2_total))*100
Percentage_ESR1_pyCC_2_GSE_9 =(1-(ESR1_pyCC_2_GSE_9_overlap/ESR1_pyCC_2_total))*100
Percentage_ESR1_pyCC_2_GSE_q =(1-(ESR1_pyCC_2_GSE_q_overlap/ESR1_pyCC_2_total))*100
Percentage_ESR1_pyCC_2_Encode1_4 =(1-(ESR1_pyCC_2_Encode1_4_overlap/ESR1_pyCC_2_total))*100
Percentage_ESR1_pyCC_2_Encode1_6 =(1-(ESR1_pyCC_2_Encode1_6_overlap/ESR1_pyCC_2_total))*100
Percentage_ESR1_pyCC_2_Encode1_9 =(1-(ESR1_pyCC_2_Encode1_9_overlap/ESR1_pyCC_2_total))*100
Percentage_ESR1_pyCC_2_Encode1_q =(1-(ESR1_pyCC_2_Encode1_q_overlap/ESR1_pyCC_2_total))*100
Percentage_ESR1_pyCC_2_Encode2_4 =(1-(ESR1_pyCC_2_Encode2_4_overlap/ESR1_pyCC_2_total))*100
Percentage_ESR1_pyCC_2_Encode2_6 =(1-(ESR1_pyCC_2_Encode2_6_overlap/ESR1_pyCC_2_total))*100
Percentage_ESR1_pyCC_2_Encode2_9 =(1-(ESR1_pyCC_2_Encode2_9_overlap/ESR1_pyCC_2_total))*100
Percentage_ESR1_pyCC_2_Encode2_q =(1-(ESR1_pyCC_2_Encode2_q_overlap/ESR1_pyCC_2_total))*100

Percentage_ESR1_pyCC_4_GSE_4 =(1-(ESR1_pyCC_4_GSE_4_overlap/ESR1_pyCC_4_total))*100
Percentage_ESR1_pyCC_4_GSE_6 =(1-(ESR1_pyCC_4_GSE_6_overlap/ESR1_pyCC_4_total))*100
Percentage_ESR1_pyCC_4_GSE_9 =(1-(ESR1_pyCC_4_GSE_9_overlap/ESR1_pyCC_4_total))*100
Percentage_ESR1_pyCC_4_GSE_q =(1-(ESR1_pyCC_4_GSE_q_overlap/ESR1_pyCC_4_total))*100
Percentage_ESR1_pyCC_4_Encode1_4 =(1-(ESR1_pyCC_4_Encode1_4_overlap/ESR1_pyCC_4_total))*100
Percentage_ESR1_pyCC_4_Encode1_6 =(1-(ESR1_pyCC_4_Encode1_6_overlap/ESR1_pyCC_4_total))*100
Percentage_ESR1_pyCC_4_Encode1_9 =(1-(ESR1_pyCC_4_Encode1_9_overlap/ESR1_pyCC_4_total))*100
Percentage_ESR1_pyCC_4_Encode1_q =(1-(ESR1_pyCC_4_Encode1_q_overlap/ESR1_pyCC_4_total))*100
Percentage_ESR1_pyCC_4_Encode2_4 =(1-(ESR1_pyCC_4_Encode2_4_overlap/ESR1_pyCC_4_total))*100
Percentage_ESR1_pyCC_4_Encode2_6 =(1-(ESR1_pyCC_4_Encode2_6_overlap/ESR1_pyCC_4_total))*100
Percentage_ESR1_pyCC_4_Encode2_9 =(1-(ESR1_pyCC_4_Encode2_9_overlap/ESR1_pyCC_4_total))*100
Percentage_ESR1_pyCC_4_Encode2_q =(1-(ESR1_pyCC_4_Encode2_q_overlap/ESR1_pyCC_4_total))*100

#Create table to plot
Finalpercentage_table = data.frame(Name=c('Peak2_GEO-4', 'Peak2_GEO-6', 'Peak2_GEO-9', 'Peak2_GEO_q0.05', 'Peak2_Encod1-4', 'Peak2_Encod1-6', 'Peak2_Encod1-9', 'Peak2_Encod1_q0.05', 'Peak2_Encod2-4', 'Peak2_Encod2-6', 'Peak2_Encod2-9', 'Peak2_Encod2_q0.05', 'Peak4_GEO-4', 'Peak4_GEO-6', 'Peak4_GEO-9', 'Peak4_GEO_q0.05', 'Peak4_Encod1-4', 'Peak4_Encod1-6', 'Peak4_Encod1-9', 'Peak4_Encod1_q0.05', 'Peak4_Encod2-4', 'Peak4_Encod2-6', 'Peak4_Encod2-9', 'Peak4_Encod2_q0.05'),
                                       Total=rep(c(ESR1_pyCC_2_total, ESR1_pyCC_2_total, ESR1_pyCC_2_total, ESR1_pyCC_2_total, ESR1_pyCC_2_total, ESR1_pyCC_2_total, ESR1_pyCC_2_total, ESR1_pyCC_2_total, ESR1_pyCC_2_total, ESR1_pyCC_2_total, ESR1_pyCC_2_total, ESR1_pyCC_2_total, ESR1_pyCC_4_total, ESR1_pyCC_4_total, ESR1_pyCC_4_total, ESR1_pyCC_4_total, ESR1_pyCC_4_total, ESR1_pyCC_4_total, ESR1_pyCC_4_total, ESR1_pyCC_4_total, ESR1_pyCC_4_total, ESR1_pyCC_4_total, ESR1_pyCC_4_total, ESR1_pyCC_4_total)),
                                       Overlap_withChIP= c(ESR1_pyCC_2_GSE_4_overlap, ESR1_pyCC_2_GSE_6_overlap, ESR1_pyCC_2_GSE_9_overlap, ESR1_pyCC_2_GSE_q_overlap, ESR1_pyCC_2_Encode1_4_overlap, ESR1_pyCC_2_Encode1_6_overlap, ESR1_pyCC_2_Encode1_9_overlap, ESR1_pyCC_2_Encode1_q_overlap, ESR1_pyCC_2_Encode2_4_overlap, ESR1_pyCC_2_Encode2_6_overlap, ESR1_pyCC_2_Encode2_9_overlap, ESR1_pyCC_2_Encode2_q_overlap, ESR1_pyCC_4_GSE_4_overlap, ESR1_pyCC_4_GSE_6_overlap, ESR1_pyCC_4_GSE_9_overlap, ESR1_pyCC_4_GSE_q_overlap, ESR1_pyCC_4_Encode1_4_overlap, ESR1_pyCC_4_Encode1_6_overlap, ESR1_pyCC_4_Encode1_9_overlap, ESR1_pyCC_4_Encode1_q_overlap, ESR1_pyCC_4_Encode2_4_overlap, ESR1_pyCC_4_Encode2_6_overlap, ESR1_pyCC_4_Encode2_9_overlap, ESR1_pyCC_4_Encode2_q_overlap),
                                       PercentageOverlapping= c(Percentage_ESR1_pyCC_2_GSE_4, Percentage_ESR1_pyCC_2_GSE_6, Percentage_ESR1_pyCC_2_GSE_9, Percentage_ESR1_pyCC_2_GSE_q, Percentage_ESR1_pyCC_2_Encode1_4, Percentage_ESR1_pyCC_2_Encode1_6, Percentage_ESR1_pyCC_2_Encode1_9, Percentage_ESR1_pyCC_2_Encode1_q, Percentage_ESR1_pyCC_2_Encode2_4, Percentage_ESR1_pyCC_2_Encode2_6, Percentage_ESR1_pyCC_2_Encode2_9, Percentage_ESR1_pyCC_2_Encode2_q, Percentage_ESR1_pyCC_4_GSE_4, Percentage_ESR1_pyCC_4_GSE_6, Percentage_ESR1_pyCC_4_GSE_9, Percentage_ESR1_pyCC_4_GSE_q, Percentage_ESR1_pyCC_4_Encode1_4, Percentage_ESR1_pyCC_4_Encode1_6, Percentage_ESR1_pyCC_4_Encode1_9, Percentage_ESR1_pyCC_4_Encode1_q, Percentage_ESR1_pyCC_4_Encode2_4, Percentage_ESR1_pyCC_4_Encode2_6, Percentage_ESR1_pyCC_4_Encode2_9, Percentage_ESR1_pyCC_4_Encode2_q))

#Then we plot
graphorder= c('Peak2_GEO-4', 'Peak2_GEO-6', 'Peak2_GEO-9', 'Peak2_GEO_q0.05', 'Peak2_Encod1-4', 'Peak2_Encod1-6', 'Peak2_Encod1-9', 'Peak2_Encod1_q0.05', 'Peak2_Encod2-4', 'Peak2_Encod2-6', 'Peak2_Encod2-9', 'Peak2_Encod2_q0.05', 'Peak4_GEO-4', 'Peak4_GEO-6', 'Peak4_GEO-9', 'Peak4_GEO_q0.05', 'Peak4_Encod1-4', 'Peak4_Encod1-6', 'Peak4_Encod1-9', 'Peak4_Encod1_q0.05', 'Peak4_Encod2-4', 'Peak4_Encod2-6', 'Peak4_Encod2-9', 'Peak4_Encod2_q0.05')
pyCC_Figure <- ggplot(data = Finalpercentage_table, 
                     aes(x= factor(Name, graphorder), 
                         y = PercentageOverlapping,
                         fill = Name, 
                         stat("identity"))) +  geom_col(color = "black",
                                                        show.legend = FALSE) +  ylim(0, 100) + labs(x = "Name of CC peak file and ChIP with pvalue cut off", 
                                                                                                    y = "% of CC peaks overlapping with ChIP-Seq peak",
                                                                                                    title = "pyCalling Cards Peaks and ChIP-Seq Overlap",
                                                                                                    subtitle = "Determining percentage of overlapping peaks from pyCalling Card ESR1 peak files within 1 kbp of either \nGSE109820, ENCFF365BIT or ENCFF063JMY ChIP-Seq datasets at varying p-value cut offs") +  theme_dark() + geom_text(aes(label = round(PercentageOverlapping, 1)), # Add percentage as text
                                                                                                                                                                                                                                                                              vjust = -0.5, # Adjust vertical position above the bars
                                                                                                                                                                                                                                                                              size = 2) + coord_cartesian(expand = FALSE, clip = "off") 


pyCC_Figure + theme(axis.text.x = element_text(angle = 90))
########################################################################
