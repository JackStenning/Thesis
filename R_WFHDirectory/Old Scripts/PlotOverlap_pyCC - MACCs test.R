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
ESR1_pyCC_4_Andy_q = read.table("results/1kbp_A_peak_data_ER_WT4.bed_B_Andy_ChIP_q0.05_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_7_Andy_q = read.table("results/1kbp_A_peak_data_ER_WT7.bed_B_Andy_ChIP_q0.05_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_8_Andy_q = read.table("results/1kbp_A_peak_data_ER_WT8.bed_B_Andy_ChIP_q0.05_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_9_Andy_q = read.table("results/1kbp_A_peak_data_ER_WT9.bed_B_Andy_ChIP_q0.05_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_10_Andy_q = read.table("results/1kbp_A_peak_data_ER_WT10.bed_B_Andy_ChIP_q0.05_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)

#Counting the total reads
ESR1_pyCC_4_total=nrow(ESR1_pyCC_4_Andy_q)
ESR1_pyCC_7_total=nrow(ESR1_pyCC_7_Andy_q)
ESR1_pyCC_8_total=nrow(ESR1_pyCC_8_Andy_q)
ESR1_pyCC_9_total=nrow(ESR1_pyCC_9_Andy_q)
ESR1_pyCC_10_total=nrow(ESR1_pyCC_10_Andy_q)




#Count the number of reads with no overlap
ESR1_pyCC_4_Andy_q_overlap = as.numeric(length(which(ESR1_pyCC_4_Andy_q$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_7_Andy_q_overlap = as.numeric(length(which(ESR1_pyCC_7_Andy_q$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_8_Andy_q_overlap = as.numeric(length(which(ESR1_pyCC_8_Andy_q$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_9_Andy_q_overlap = as.numeric(length(which(ESR1_pyCC_9_Andy_q$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_10_Andy_q_overlap = as.numeric(length(which(ESR1_pyCC_10_Andy_q$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))

#Determine overlap percentage
Percentage_ESR1_pyCC_4_Andy_q =(1-(ESR1_pyCC_4_Andy_q_overlap/ESR1_pyCC_4_total))*100
Percentage_ESR1_pyCC_7_Andy_q =(1-(ESR1_pyCC_7_Andy_q_overlap/ESR1_pyCC_7_total))*100
Percentage_ESR1_pyCC_8_Andy_q =(1-(ESR1_pyCC_8_Andy_q_overlap/ESR1_pyCC_8_total))*100
Percentage_ESR1_pyCC_9_Andy_q =(1-(ESR1_pyCC_9_Andy_q_overlap/ESR1_pyCC_9_total))*100
Percentage_ESR1_pyCC_10_Andy_q =(1-(ESR1_pyCC_10_Andy_q_overlap/ESR1_pyCC_10_total))*100

#Create table to plot
#Create table to plot
pyCC_Finalpercentage_table = data.frame(Name=c('ER_4', 'ER_7', 'ER_8', 'ER_9', 'ER_10'),
                                       Total=rep(c(ESR1_pyCC_4_total, ESR1_pyCC_7_total, ESR1_pyCC_8_total, ESR1_pyCC_9_total, ESR1_pyCC_10_total)),
                                       Overlap_withChIP= c(ESR1_pyCC_4_Andy_q_overlap, ESR1_pyCC_7_Andy_q_overlap, ESR1_pyCC_8_Andy_q_overlap, ESR1_pyCC_9_Andy_q_overlap, ESR1_pyCC_10_Andy_q_overlap),
                                       PercentageOverlapping= c(Percentage_ESR1_pyCC_4_Andy_q, Percentage_ESR1_pyCC_7_Andy_q, Percentage_ESR1_pyCC_8_Andy_q, Percentage_ESR1_pyCC_9_Andy_q, Percentage_ESR1_pyCC_10_Andy_q))

#Then we plot
graphorder= c('ER_4', 'ER_7', 'ER_8', 'ER_9', 'ER_10')
pyCC_Figure <- ggplot(data = pyCC_Finalpercentage_table, 
                     aes(x= factor(Name, graphorder), 
                         y = PercentageOverlapping,
                         fill = Name, 
                         stat("identity"))) +  geom_col(color = "black",
                                                        show.legend = FALSE) +  ylim(0, 100) + labs(x = "Name of CC peak file and ChIP with pvalue cut off", 
                                                                                                    y = "% of CC peaks overlapping with ChIP-Seq peak",
                                                                                                    title = "pyCalling Cards Peaks and ChIP-Seq Overlap",
                                                                                                    subtitle = "Determining percentage of overlapping peaks from pyCalling Card ESR1 peak files and either Andy's or \nEncode ChIP-Seq Dataset") +  theme_dark() + coord_cartesian(expand = FALSE, clip = "off") + geom_hline(yintercept = 83, linetype = "dashed") 


pyCC_Figure + theme(axis.text.x = element_text(angle = 90))
########################################################################
