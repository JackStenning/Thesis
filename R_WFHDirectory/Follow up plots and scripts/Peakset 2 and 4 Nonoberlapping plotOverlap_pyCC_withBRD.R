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
ESR1_pyCC_2_encode_Liu = read.table("overlap/cc_BRD4/count_BRD4_Peakset2_Encode_no_chipseq_overlaps.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_2_encode_Zhang = read.table("overlap/cc_BRD4/count_BRD4Z_Peakset2_Encode_no_chipseq_overlaps.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)


ESR1_pyCC_2_GEO_Liu = read.table("overlap/cc_BRD4/count_BRD4_Peakset2_GEO_no_chipseq_overlaps.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_2_GEO_Zhang = read.table("overlap/cc_BRD4/count_BRD4Z_Peakset2_GEO_no_chipseq_overlaps.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)


#Counting the total reads

ESR1_pyCC_2_total=nrow(ESR1_pyCC_2_encode_Liu)

ESR1_pyCC_4_total=nrow(ESR1_pyCC_2_GEO_Liu)


#Count the number of reads with no overlap
ESR1_pyCC_2_encode_Liu_overlap = as.numeric(length(which(ESR1_pyCC_2_encode_Liu$Number.of.Overlaps.with.BRD4.ChIP.Seq.Peaks==0)))
ESR1_pyCC_2_encode_Zhang_overlap = as.numeric(length(which(ESR1_pyCC_2_encode_Zhang$Number.of.Overlaps.with.BRD4.ChIP.Seq.Peaks==0)))



ESR1_pyCC_2_GEO_Liu_overlap = as.numeric(length(which(ESR1_pyCC_2_GEO_Liu$Number.of.Overlaps.with.BRD4.ChIP.Seq.Peaks==0)))
ESR1_pyCC_4_Zhang_overlap = as.numeric(length(which(ESR1_pyCC_2_GEO_Zhang$Number.of.Overlaps.with.BRD4.ChIP.Seq.Peaks==0)))


#Determine overlap percentage
Percentage_ESR1_pyCC_2_encode_Liu =(1-(ESR1_pyCC_2_encode_Liu_overlap/ESR1_pyCC_2_total))*100
Percentage_ESR1_pyCC_2_encode_Zhang =(1-(ESR1_pyCC_2_encode_Zhang_overlap/ESR1_pyCC_2_total))*100



Percentage_ESR1_pyCC_2_GEO_Liu =(1-(ESR1_pyCC_2_GEO_Liu_overlap/ESR1_pyCC_4_total))*100
Percentage_ESR1_pyCC_4_Zhang = (1-(ESR1_pyCC_4_Zhang_overlap/ESR1_pyCC_4_total))*100

#Create table to plot
#Create table to plot
pyCC_Finalpercentage_table = data.frame(Name=c('Peakset_2_Encode_Liu', 'Peakset_2_Encode_Zhang', 'Peakset_2_GEO_Liu', 'Peakset_2_GEO_Zhang'),
                                       Total=rep(c(ESR1_pyCC_2_total, ESR1_pyCC_2_total, ESR1_pyCC_4_total, ESR1_pyCC_4_total)),
                                       Overlap_withChIP= c(ESR1_pyCC_2_encode_Liu_overlap, ESR1_pyCC_2_encode_Zhang_overlap, ESR1_pyCC_2_GEO_Liu_overlap, ESR1_pyCC_4_Zhang_overlap),
                                       PercentageOverlapping= c(Percentage_ESR1_pyCC_2_encode_Liu, Percentage_ESR1_pyCC_2_encode_Zhang, Percentage_ESR1_pyCC_2_GEO_Liu, Percentage_ESR1_pyCC_4_Zhang))

#Then we plot
graphorder= c('Peakset_2_Encode_Liu', 'Peakset_2_Encode_Zhang', 'Peakset_2_GEO_Liu', 'Peakset_2_GEO_Zhang')
tiff("High stringency non-overlapping peaks pyCalling Cards Peaks and ChIP-Seq Overlap.tiff", width = 3510, height = 2481, res = 300)
pyCC_Figure <- ggplot(data = pyCC_Finalpercentage_table, 
                     aes(x= factor(Name, graphorder), 
                         y = PercentageOverlapping,
                         fill = Name, 
                         stat("identity"))) +  geom_col(color = "black",
                                                        show.legend = FALSE) +  ylim(0, 100) + labs(x = "Name of CC peak file and ChIP Dataset", 
                                                                                                    y = "% of CC peaks overlapping with ChIP-Seq peak",
                                                                                                    title = "pyCalling Cards Peaks and ChIP-Seq Overlap",
                                                                                                    subtitle = "Determining percentage of overlapping peaks from pyCalling Card ESR1 high stringency peak files with either Liu or \nZhnag BRD4 ChIP-Seq datasets within 1 kbp") +  theme_minimal() + geom_text(aes(label = round(PercentageOverlapping, 1)), # Add percentage as text
                                                                                                                                                                                                                                                                             vjust = -0.5, # Adjust vertical position above the bars
                                                                                                                                                                                                                                                                             size = 4) + coord_cartesian(expand = FALSE, clip = "off") + geom_hline(yintercept = 83, linetype = "dashed") 


pyCC_Figure + theme(axis.text.x = element_text(angle = 90))
dev.off()
########################################################################

#Set column names
pyCCcol_names = c("Chr",	"Start",	"End",	"Center",	"Experiment Insertions",	"Background insertions",	"Reference Insertions",	"pvalue Reference",	"pvalue Background",	"Fraction Experiment",	"TPH Experiment",	"Fraction background",	"TPH background",	"TPH background subtracted",	"pvalue_adj Reference", "Number of Overlaps with ER ChIP-Seq Peaks", "Number of Overlaps with BRD4 ChIP-Seq Peaks")

##Read data
#Variables named <Window '-a'><Window'-b'>
#Encode Rep 2 was taken as it has more peaks than replicate 1 - even though there is less overlap with my data
ESR1_pyCC_2_encode_Liu = read.table("overlap/cc_BRD4/count_BRD4_Peakset4_Encode_no_chipseq_overlaps.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_2_encode_Zhang = read.table("overlap/cc_BRD4/count_BRD4Z_Peakset4_Encode_no_chipseq_overlaps.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)


ESR1_pyCC_2_GEO_Liu = read.table("overlap/cc_BRD4/count_BRD4_Peakset4_GEO_no_chipseq_overlaps.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_2_GEO_Zhang = read.table("overlap/cc_BRD4/count_BRD4Z_Peakset4_GEO_no_chipseq_overlaps.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)


#Counting the total reads

ESR1_pyCC_2_total=nrow(ESR1_pyCC_2_encode_Liu)

ESR1_pyCC_4_total=nrow(ESR1_pyCC_2_GEO_Liu)


#Count the number of reads with no overlap
ESR1_pyCC_2_encode_Liu_overlap = as.numeric(length(which(ESR1_pyCC_2_encode_Liu$Number.of.Overlaps.with.BRD4.ChIP.Seq.Peaks==0)))
ESR1_pyCC_2_encode_Zhang_overlap = as.numeric(length(which(ESR1_pyCC_2_encode_Zhang$Number.of.Overlaps.with.BRD4.ChIP.Seq.Peaks==0)))



ESR1_pyCC_2_GEO_Liu_overlap = as.numeric(length(which(ESR1_pyCC_2_GEO_Liu$Number.of.Overlaps.with.BRD4.ChIP.Seq.Peaks==0)))
ESR1_pyCC_4_Zhang_overlap = as.numeric(length(which(ESR1_pyCC_2_GEO_Zhang$Number.of.Overlaps.with.BRD4.ChIP.Seq.Peaks==0)))


#Determine overlap percentage
Percentage_ESR1_pyCC_2_encode_Liu =(1-(ESR1_pyCC_2_encode_Liu_overlap/ESR1_pyCC_2_total))*100
Percentage_ESR1_pyCC_2_encode_Zhang =(1-(ESR1_pyCC_2_encode_Zhang_overlap/ESR1_pyCC_2_total))*100



Percentage_ESR1_pyCC_2_GEO_Liu =(1-(ESR1_pyCC_2_GEO_Liu_overlap/ESR1_pyCC_4_total))*100
Percentage_ESR1_pyCC_4_Zhang = (1-(ESR1_pyCC_4_Zhang_overlap/ESR1_pyCC_4_total))*100

#Create table to plot
#Create table to plot
pyCC_Finalpercentage_table = data.frame(Name=c('Peakset_4_Encode_Liu', 'Peakset_4_Encode_Zhang', 'Peakset_4_GEO_Liu', 'Peakset_4_GEO_Zhang'),
                                        Total=rep(c(ESR1_pyCC_2_total, ESR1_pyCC_2_total, ESR1_pyCC_4_total, ESR1_pyCC_4_total)),
                                        Overlap_withChIP= c(ESR1_pyCC_2_encode_Liu_overlap, ESR1_pyCC_2_encode_Zhang_overlap, ESR1_pyCC_2_GEO_Liu_overlap, ESR1_pyCC_4_Zhang_overlap),
                                        PercentageOverlapping= c(Percentage_ESR1_pyCC_2_encode_Liu, Percentage_ESR1_pyCC_2_encode_Zhang, Percentage_ESR1_pyCC_2_GEO_Liu, Percentage_ESR1_pyCC_4_Zhang))

#Then we plot
graphorder= c('Peakset_4_Encode_Liu', 'Peakset_4_Encode_Zhang', 'Peakset_4_GEO_Liu', 'Peakset_4_GEO_Zhang')
tiff("Low stringency non-overlapping peaks pyCalling Cards Peaks and ChIP-Seq Overlap.tiff", width = 3510, height = 2481, res = 300)
pyCC_Figure <- ggplot(data = pyCC_Finalpercentage_table, 
                      aes(x= factor(Name, graphorder), 
                          y = PercentageOverlapping,
                          fill = Name, 
                          stat("identity"))) +  geom_col(color = "black",
                                                         show.legend = FALSE) +  ylim(0, 100) + labs(x = "Name of CC peak file and ChIP Dataset", 
                                                                                                     y = "% of CC peaks overlapping with ChIP-Seq peak",
                                                                                                     title = "pyCalling Cards Peaks and ChIP-Seq Overlap",
                                                                                                     subtitle = "Determining percentage of overlapping peaks from pyCalling Card ESR1 Low stringency peak files with either Liu or \nZhnag BRD4 ChIP-Seq datasets within 1 kbp") +  theme_minimal() + geom_text(aes(label = round(PercentageOverlapping, 1)), # Add percentage as text
                                                                                                                                                                                                                                                                                                                 vjust = -0.5, # Adjust vertical position above the bars
                                                                                                                                                                                                                                                                                                                 size = 4) + coord_cartesian(expand = FALSE, clip = "off") + geom_hline(yintercept = 83, linetype = "dashed") 


pyCC_Figure + theme(axis.text.x = element_text(angle = 90))
dev.off()
########################################################################

