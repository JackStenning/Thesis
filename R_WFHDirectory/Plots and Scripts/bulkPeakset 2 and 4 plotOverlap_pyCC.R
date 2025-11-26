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
#Encode Rep 2 was taken as it has more peaks than replicate 1 - even though there is less overlap with my data
ESR1_pyCC_2_GEO_q = read.table("results/1kbp_A_peak_data_ER_WT2.bed_B_Andy_ChIP_q0.05_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_2_Encode2_q = read.table("results/1kbp_A_peak_data_ER_WT2.bed_B_Encode_ER_rep2_ChIP_q0.05_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)


ESR1_pyCC_4_GEO_q = read.table("results/1kbp_A_peak_data_ER_WT4.bed_B_Andy_ChIP_q0.05_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_4_Encode2_q = read.table("results/1kbp_A_peak_data_ER_WT4.bed_B_Encode_ER_rep2_ChIP_q0.05_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)


#Counting the total reads

ESR1_pyCC_2_total=nrow(ESR1_pyCC_2_GEO_q)

ESR1_pyCC_4_total=nrow(ESR1_pyCC_4_GEO_q)


#Count the number of reads with no overlap
ESR1_pyCC_2_GEO_q_overlap = as.numeric(length(which(ESR1_pyCC_2_GEO_q$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_2_Encode2_q_overlap = as.numeric(length(which(ESR1_pyCC_2_Encode2_q$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))



ESR1_pyCC_4_GEO_q_overlap = as.numeric(length(which(ESR1_pyCC_4_GEO_q$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_4_Encode2_q_overlap = as.numeric(length(which(ESR1_pyCC_4_Encode2_q$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))


#Determine overlap percentage
Percentage_ESR1_pyCC_2_GEO_q =(1-(ESR1_pyCC_2_GEO_q_overlap/ESR1_pyCC_2_total))*100
Percentage_ESR1_pyCC_2_Encode2_q =(1-(ESR1_pyCC_2_Encode2_q_overlap/ESR1_pyCC_2_total))*100



Percentage_ESR1_pyCC_4_GEO_q =(1-(ESR1_pyCC_4_GEO_q_overlap/ESR1_pyCC_4_total))*100
Percentage_ESR1_pyCC_4_Encode2_q = (1-(ESR1_pyCC_4_Encode2_q_overlap/ESR1_pyCC_4_total))*100

#Create table to plot
#Create table to plot
pyCC_Finalpercentage_table = data.frame(Name=c('Peakset_2_GEO', 'Peakset_2_Encode', 'Peakset_4_GEO', 'Peakset_4_Encode'),
                                       Total=rep(c(ESR1_pyCC_2_total, ESR1_pyCC_2_total, ESR1_pyCC_4_total, ESR1_pyCC_4_total)),
                                       Overlap_withChIP= c(ESR1_pyCC_2_GEO_q_overlap, ESR1_pyCC_2_Encode2_q_overlap, ESR1_pyCC_4_GEO_q_overlap, ESR1_pyCC_4_Encode2_q_overlap),
                                       PercentageOverlapping= c(Percentage_ESR1_pyCC_2_GEO_q, Percentage_ESR1_pyCC_2_Encode2_q, Percentage_ESR1_pyCC_4_GEO_q, Percentage_ESR1_pyCC_4_Encode2_q))

#Then we plot
graphorder= c('Peakset_2_GEO', 'Peakset_2_Encode', 'Peakset_4_GEO', 'Peakset_4_Encode')
pyCC_Figure <- ggplot(data = pyCC_Finalpercentage_table, 
                     aes(x= factor(Name, graphorder), 
                         y = PercentageOverlapping,
                         fill = Name, 
                         stat("identity"))) +  geom_col(color = "black",
                                                        show.legend = FALSE) +  ylim(0, 100) + labs(x = "Name of CC peak file and ChIP Dataset", 
                                                                                                    y = "% of CC peaks overlapping with ChIP-Seq peak",
                                                                                                    title = "pyCalling Cards Peaks and ChIP-Seq Overlap",
                                                                                                    subtitle = "Determining percentage of overlapping peaks from pyCalling Card ESR1 peak files with either GSE109820 or \nENCFF063JMY ChIP-Seq datasets within 1 kbp") +  theme_dark() + geom_text(aes(label = round(PercentageOverlapping, 1)), # Add percentage as text
                                                                                                                                                                                                                                                                             vjust = -0.5, # Adjust vertical position above the bars
                                                                                                                                                                                                                                                                             size = 4) + coord_cartesian(expand = FALSE, clip = "off") + geom_hline(yintercept = 83, linetype = "dashed") 


pyCC_Figure + theme(axis.text.x = element_text(angle = 90))
########################################################################
