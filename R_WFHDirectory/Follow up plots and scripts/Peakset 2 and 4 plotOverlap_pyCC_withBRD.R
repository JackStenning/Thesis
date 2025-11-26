############# LOOK HERE ######################
#add or delete the '#' symbol below on install packages
#install.packages("Rtools")
#install.packages("ggthemes")

#DELTE AFTER WRITING

#Load software
library(ggplot2)
library(ggthemes)

#Set working directory
setwd("C:/Users/jps558/OneDrive - University of York/Desktop/R_WorkingDirectory/R_WFHDirectory") 


#Set column names
pyCCcol_names = c("Chr",	"Start",	"End",	"Center",	"Experiment Insertions",	"Background insertions",	"Reference Insertions",	"pvalue Reference",	"pvalue Background",	"Fraction Experiment",	"TPH Experiment",	"Fraction background",	"TPH background",	"TPH background subtracted",	"pvalue_adj Reference", "Number of Overlaps with ER ChIP-Seq Peaks", "Number of Overlaps with BRD4 ChIP-Seq Peaks")

##Read data
#Variables named <Window '-a'><Window'-b'>
#Encode Rep 2 was taken as it has more peaks than replicate 1 - even though there is less overlap with my data
ESR1_pyCC_2_Liu = read.table("overlap/cc_BRD4/count_BRD4_1kbp_A_peak_data_ER_WT2.bed_B_Andy_ChIP_q0.05_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_2_Zhang = read.table("overlap/cc_BRD4/count_BRD4Z_1kbp_A_peak_data_ER_WT2.bed_B_Andy_ChIP_q0.05_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)


ESR1_pyCC_4_Liu = read.table("overlap/cc_BRD4/count_BRD4_1kbp_A_peak_data_ER_WT4.bed_B_Andy_ChIP_q0.05_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_4_Zhang = read.table("overlap/cc_BRD4/count_BRD4Z_1kbp_A_peak_data_ER_WT4.bed_B_Andy_ChIP_q0.05_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)


#Counting the total reads

ESR1_pyCC_2_total=nrow(ESR1_pyCC_2_Liu)

ESR1_pyCC_4_total=nrow(ESR1_pyCC_4_Liu)


#Count the number of reads with no overlap
ESR1_pyCC_2_Liu_overlap = as.numeric(length(which(ESR1_pyCC_2_Liu$Number.of.Overlaps.with.BRD4.ChIP.Seq.Peaks==0)))
ESR1_pyCC_2_Zhang_overlap = as.numeric(length(which(ESR1_pyCC_2_Zhang$Number.of.Overlaps.with.BRD4.ChIP.Seq.Peaks==0)))



ESR1_pyCC_4_Liu_overlap = as.numeric(length(which(ESR1_pyCC_4_Liu$Number.of.Overlaps.with.BRD4.ChIP.Seq.Peaks==0)))
ESR1_pyCC_4_Zhang_overlap = as.numeric(length(which(ESR1_pyCC_4_Zhang$Number.of.Overlaps.with.BRD4.ChIP.Seq.Peaks==0)))


#Determine overlap percentage
Percentage_ESR1_pyCC_2_Liu =(1-(ESR1_pyCC_2_Liu_overlap/ESR1_pyCC_2_total))*100
Percentage_ESR1_pyCC_2_Zhang =(1-(ESR1_pyCC_2_Zhang_overlap/ESR1_pyCC_2_total))*100



Percentage_ESR1_pyCC_4_Liu =(1-(ESR1_pyCC_4_Liu_overlap/ESR1_pyCC_4_total))*100
Percentage_ESR1_pyCC_4_Zhang = (1-(ESR1_pyCC_4_Zhang_overlap/ESR1_pyCC_4_total))*100

#Create table to plot
#Create table to plot
pyCC_Finalpercentage_table = data.frame(Name=c('Peakset_2_Liu', 'Peakset_2_Zhang', 'Peakset_4_Liu', 'Peakset_4_Zhang'),
                                       Total=rep(c(ESR1_pyCC_2_total, ESR1_pyCC_2_total, ESR1_pyCC_4_total, ESR1_pyCC_4_total)),
                                       Overlap_withChIP= c(ESR1_pyCC_2_Liu_overlap, ESR1_pyCC_2_Zhang_overlap, ESR1_pyCC_4_Liu_overlap, ESR1_pyCC_4_Zhang_overlap),
                                       PercentageOverlapping= c(Percentage_ESR1_pyCC_2_Liu, Percentage_ESR1_pyCC_2_Zhang, Percentage_ESR1_pyCC_4_Liu, Percentage_ESR1_pyCC_4_Zhang))

#Then we plot
graphorder= c('Peakset_2_Liu', 'Peakset_2_Zhang', 'Peakset_4_Liu', 'Peakset_4_Zhang')
tiff("pyCalling Cards Peaks high and low and BRD4 ChIP-Seq Overlap.tiff", width = 2481, height = 1749, res = 300)
pyCC_Figure <- ggplot(data = pyCC_Finalpercentage_table, 
                     aes(x= factor(Name, graphorder), 
                         y = PercentageOverlapping,
                         fill = Name, 
                         stat("identity"))) +  geom_col(color = "black",
                                                        show.legend = FALSE) +  ylim(0, 100) + labs(x = "Name of CC peak file and ChIP Dataset", 
                                                                                                    y = "% of CC peaks overlapping with ChIP-Seq peak",
                                                                                                    title = "pyCalling Cards Peaks and ChIP-Seq Overlap",
                                                                                                    subtitle = "Determining percentage of overlapping peaks from pyCalling Card ESR1 peak files with either Liu or \nZhnag BRD4 ChIP-Seq datasets within 1 kbp") +  theme_minimal() + geom_text(aes(label = round(PercentageOverlapping, 1)), # Add percentage as text
                                                                                                                                                                                                                                                                             vjust = -0.5, # Adjust vertical position above the bars
                                                                                                                                                                                                                                                                             size = 4) + coord_cartesian(expand = FALSE, clip = "off") + geom_hline(yintercept = 83, linetype = "dashed") 


pyCC_Figure + theme(axis.text.x = element_text(angle = 90))
dev.off()
########################################################################
