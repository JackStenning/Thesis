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
pyCCcol_names = c("Chr",	"Start",	"End",	"Name",	"Score",	"Strand",	"Signal Value",	"pvalue",	"qvalue",	"Peak", "Number of Overlaps with BRD4 ChIP-Seq Peaks")

##Read data
#Variables named <Window '-a'><Window'-b'>
#Encode Rep 2 was taken as it has more peaks than replicate 1 - even though there is less overlap with my data
GEO_Liu = read.table("overlap/chip_BRD4/count_BRD4_Andy_ChIP_q0.05_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
GEO_Zhang = read.table("overlap/chip_BRD4/count_BRD4Z_Andy_ChIP_q0.05_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)


Encode_Liu = read.table("overlap/chip_BRD4/count_BRD4_Encode_ER_rep2_ChIP_q0.05_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
Encode_Zhang = read.table("overlap/chip_BRD4/count_BRD4Z_Encode_ER_rep2_ChIP_q0.05_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)


#Counting the total reads

ESR1_pyCC_2_total=nrow(GEO_Liu)

ESR1_pyCC_4_total=nrow(Encode_Liu)


#Count the number of reads with no overlap
GEO_Liu_overlap = as.numeric(length(which(GEO_Liu$Number.of.Overlaps.with.BRD4.ChIP.Seq.Peaks==0)))
GEO_Zhang_overlap = as.numeric(length(which(GEO_Zhang$Number.of.Overlaps.with.BRD4.ChIP.Seq.Peaks==0)))



Encode_Liu_overlap = as.numeric(length(which(Encode_Liu$Number.of.Overlaps.with.BRD4.ChIP.Seq.Peaks==0)))
Encode_Zhang_overlap = as.numeric(length(which(Encode_Zhang$Number.of.Overlaps.with.BRD4.ChIP.Seq.Peaks==0)))


#Determine overlap percentage
Percentage_GEO_Liu =(1-(GEO_Liu_overlap/ESR1_pyCC_2_total))*100
Percentage_GEO_Zhang =(1-(GEO_Zhang_overlap/ESR1_pyCC_2_total))*100



Percentage_Encode_Liu =(1-(Encode_Liu_overlap/ESR1_pyCC_4_total))*100
Percentage_Encode_Zhang = (1-(Encode_Zhang_overlap/ESR1_pyCC_4_total))*100

#Create table to plot
#Create table to plot
pyCC_Finalpercentage_table = data.frame(Name=c('GEO_Liu', 'GEO_Zhang', 'Encode_Liu', 'Encode_Zhang'),
                                       Total=rep(c(ESR1_pyCC_2_total, ESR1_pyCC_2_total, ESR1_pyCC_4_total, ESR1_pyCC_4_total)),
                                       Overlap_withChIP= c(GEO_Liu_overlap, GEO_Zhang_overlap, Encode_Liu_overlap, Encode_Zhang_overlap),
                                       PercentageOverlapping= c(Percentage_GEO_Liu, Percentage_GEO_Zhang, Percentage_Encode_Liu, Percentage_Encode_Zhang))

#Then we plot
graphorder= c('GEO_Liu', 'GEO_Zhang', 'Encode_Liu', 'Encode_Zhang')

tiff("ER and BRD4 ChIP-Seq Overlap.tiff", width = 2481, height = 1749, res = 300) 
pyCC_Figure <- ggplot(data = pyCC_Finalpercentage_table, 
                     aes(x= factor(Name, graphorder), 
                         y = PercentageOverlapping,
                         fill = Name, 
                         stat("identity"))) +  geom_col(color = "black",
                                                        show.legend = FALSE) +  ylim(0, 100) + labs(x = "Name of ChIP Dataset", 
                                                                                                    y = "% of ER ChIP peaks overlapping with BRD4 ChIP-Seq peak",
                                                                                                    title = "ER and BRD4 ChIP-Seq Overlap",
                                                                                                    subtitle = "Determining percentage of overlapping peaks from ER ChIP-Seq peak files with either Liu or \nZhnag BRD4 ChIP-Seq datasets within 1 kbp") +  theme_minimal() + geom_text(aes(label = round(PercentageOverlapping, 1)), # Add percentage as text
                                                                                                                                                                                                                                                                             vjust = -0.5, # Adjust vertical position above the bars
                                                                                                                                                                                                                                                                             size = 4) + coord_cartesian(expand = FALSE, clip = "off") + geom_hline(yintercept = 83, linetype = "dashed") 


pyCC_Figure + theme(axis.text.x = element_text(angle = 90))
dev.off()
########################################################################
