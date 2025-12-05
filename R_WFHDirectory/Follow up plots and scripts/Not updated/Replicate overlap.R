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
ESR1_pyCC_2_ER1 = read.table("rep_overlap/a_peakset2_b_filtered_sortedER1_EKDL230008858.qbed.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_2_ER2 = read.table("rep_overlap/a_peakset2_b_filtered_sortedER2_EKDL230008858.qbed.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_2_ER3 = read.table("rep_overlap/a_peakset2_b_filtered_sortedER3_EKDL230008858.qbed.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_2_ER5 = read.table("rep_overlap/a_peakset2_b_filtered_sortedER5_EKDL230008858.qbed.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_2_ER6 = read.table("rep_overlap/a_peakset2_b_filtered_sortedER6_EKDL230008858.qbed.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)


ESR1_pyCC_4_ER1 = read.table("rep_overlap/a_peakset4_b_filtered_sortedER1_EKDL230008858.qbed.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_4_ER2 = read.table("rep_overlap/a_peakset4_b_filtered_sortedER2_EKDL230008858.qbed.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_4_ER3 = read.table("rep_overlap/a_peakset4_b_filtered_sortedER3_EKDL230008858.qbed.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_4_ER5 = read.table("rep_overlap/a_peakset4_b_filtered_sortedER5_EKDL230008858.qbed.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_4_ER6 = read.table("rep_overlap/a_peakset4_b_filtered_sortedER6_EKDL230008858.qbed.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)

#Counting the total reads

ESR1_pyCC_2_total=nrow(ESR1_pyCC_2_ER1)

ESR1_pyCC_4_total=nrow(ESR1_pyCC_4_ER1)





#Count the number of reads with no non_overlap
ESR1_pyCC_2_ER1_non_overlap = as.numeric(length(which(ESR1_pyCC_2_ER1$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_2_ER2_non_overlap = as.numeric(length(which(ESR1_pyCC_2_ER2$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_2_ER3_non_overlap = as.numeric(length(which(ESR1_pyCC_2_ER3$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_2_ER5_non_overlap = as.numeric(length(which(ESR1_pyCC_2_ER5$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_2_ER6_non_overlap = as.numeric(length(which(ESR1_pyCC_2_ER6$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))


ESR1_pyCC_4_ER1_non_overlap = as.numeric(length(which(ESR1_pyCC_4_ER1$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_4_ER2_non_overlap = as.numeric(length(which(ESR1_pyCC_4_ER2$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_4_ER3_non_overlap = as.numeric(length(which(ESR1_pyCC_4_ER3$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_4_ER5_non_overlap = as.numeric(length(which(ESR1_pyCC_4_ER5$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_4_ER6_non_overlap = as.numeric(length(which(ESR1_pyCC_4_ER6$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))


#Determine overlap percentage
Percentage_ESR1_pyCC_2_ER1 =(1-(ESR1_pyCC_2_ER1_non_overlap/ESR1_pyCC_2_total))*100
Percentage_ESR1_pyCC_2_ER2 =(1-(ESR1_pyCC_2_ER2_non_overlap/ESR1_pyCC_2_total))*100
Percentage_ESR1_pyCC_2_ER3 =(1-(ESR1_pyCC_2_ER3_non_overlap/ESR1_pyCC_2_total))*100
Percentage_ESR1_pyCC_2_ER5 =(1-(ESR1_pyCC_2_ER5_non_overlap/ESR1_pyCC_2_total))*100
Percentage_ESR1_pyCC_2_ER6 =(1-(ESR1_pyCC_2_ER6_non_overlap/ESR1_pyCC_2_total))*100

Percentage_ESR1_pyCC_4_ER1 =(1-(ESR1_pyCC_4_ER1_non_overlap/ESR1_pyCC_4_total))*100
Percentage_ESR1_pyCC_4_ER2 =(1-(ESR1_pyCC_4_ER2_non_overlap/ESR1_pyCC_4_total))*100
Percentage_ESR1_pyCC_4_ER3 =(1-(ESR1_pyCC_4_ER3_non_overlap/ESR1_pyCC_4_total))*100
Percentage_ESR1_pyCC_4_ER5 =(1-(ESR1_pyCC_4_ER5_non_overlap/ESR1_pyCC_4_total))*100
Percentage_ESR1_pyCC_4_ER6 =(1-(ESR1_pyCC_4_ER6_non_overlap/ESR1_pyCC_4_total))*100



#Create table to plot
#Create table to plot
pyCC_Finalpercentage_table = data.frame(Name=c('Peak_2_ER1', 'Peak_2_ER2', 'Peak_2_ER3', 'Peak_2_ER5', 'Peak_2_ER6', 'Peak_4_ER1', 'Peak_4_ER2', 'Peak_4_ER3', 'Peak_4_ER5', 'Peak_4_ER6'),
                                       Total=rep(c(ESR1_pyCC_2_total, ESR1_pyCC_2_total, ESR1_pyCC_2_total, ESR1_pyCC_2_total, ESR1_pyCC_2_total, ESR1_pyCC_4_total, ESR1_pyCC_4_total, ESR1_pyCC_4_total, ESR1_pyCC_4_total, ESR1_pyCC_4_total)),
                                       Overlapping_with_Replicates = c(ESR1_pyCC_2_ER1_overlaping, ESR1_pyCC_2_ER2_overlaping, ESR1_pyCC_2_ER3_overlaping, ESR1_pyCC_2_ER5_overlaping, ESR1_pyCC_2_ER6_overlaping, ESR1_pyCC_4_ER1_overlaping, ESR1_pyCC_4_ER2_overlaping, ESR1_pyCC_4_ER3_overlaping, ESR1_pyCC_4_ER5_overlaping, ESR1_pyCC_2_ER6_overlaping),
                                       Non_overlapping_with_Replicates= c(ESR1_pyCC_2_ER1_non_overlap, ESR1_pyCC_2_ER2_non_overlap, ESR1_pyCC_2_ER3_non_overlap, ESR1_pyCC_2_ER5_non_overlap, ESR1_pyCC_2_ER6_non_overlap, ESR1_pyCC_4_ER1_non_overlap, ESR1_pyCC_4_ER2_non_overlap, ESR1_pyCC_4_ER3_non_overlap, ESR1_pyCC_4_ER5_non_overlap, ESR1_pyCC_4_ER6_non_overlap),
                                       PercentageOverlapping= c(Percentage_ESR1_pyCC_2_ER1, Percentage_ESR1_pyCC_2_ER2, Percentage_ESR1_pyCC_2_ER3, Percentage_ESR1_pyCC_2_ER5, Percentage_ESR1_pyCC_2_ER6, Percentage_ESR1_pyCC_4_ER1, Percentage_ESR1_pyCC_4_ER2, Percentage_ESR1_pyCC_4_ER3, Percentage_ESR1_pyCC_4_ER5, Percentage_ESR1_pyCC_4_ER6))

write.csv(pyCC_Finalpercentage_table, "rep_overlap/Peakset overlap with Replicates.csv", row.names = FALSE)

#Then we plot
graphorder= c('Peak_2_ER1', 'Peak_2_ER2', 'Peak_2_ER3', 'Peak_2_ER5', 'Peak_2_ER6', 'Peak_4_ER1', 'Peak_4_ER2', 'Peak_4_ER3', 'Peak_4_ER5', 'Peak_4_ER6')
pyCC_Figure <- ggplot(data = pyCC_Finalpercentage_table, 
                     aes(x= factor(Name, graphorder), 
                         y = PercentageOverlapping,
                         fill = Name, 
                         stat("identity"))) +  geom_col(color = "black",
                                                        show.legend = FALSE) +  ylim(0, 100) + labs(x = "Name of CC peak file and ER Calling Card replicate", 
                                                                                                    y = "% of CC peaks overlapping with replicate",
                                                                                                    title = "pyCalling Cards Peaks overlap with replicates",
                                                                                                    subtitle = "Determining percentage of overlapping peaks from pyCalling Card \nESR1 peak files and the replicates used to generate them") +  theme_dark() + coord_cartesian(expand = FALSE, clip = "off")


pyCC_Figure + theme(axis.text.x = element_text(angle = 90))


#VennDiagram
if (!requireNamespace("ggVennDiagram", quietly = TRUE)) {
  install.packages("ggVennDiagram")
}

# Load ggVennDiagram
library(ggVennDiagram)
library(ggplot2)

#Create list

ESR1_pyCC_2_ER1_overlaping = seq.int(1, (ESR1_pyCC_2_total-(as.numeric(length(which(ESR1_pyCC_2_ER1$Number.of.Overlaps.with.ChIP.Seq.Peaks==0))))))
ESR1_pyCC_2_ER2_overlaping = seq.int(1, (ESR1_pyCC_2_total-(as.numeric(length(which(ESR1_pyCC_2_ER2$Number.of.Overlaps.with.ChIP.Seq.Peaks==0))))))
ESR1_pyCC_2_ER3_overlaping = seq.int(1, (ESR1_pyCC_2_total-(as.numeric(length(which(ESR1_pyCC_2_ER3$Number.of.Overlaps.with.ChIP.Seq.Peaks==0))))))
ESR1_pyCC_2_ER5_overlaping = seq.int(1, (ESR1_pyCC_2_total-(as.numeric(length(which(ESR1_pyCC_2_ER5$Number.of.Overlaps.with.ChIP.Seq.Peaks==0))))))
ESR1_pyCC_2_ER6_overlaping = seq.int(1, (ESR1_pyCC_2_total-(as.numeric(length(which(ESR1_pyCC_2_ER6$Number.of.Overlaps.with.ChIP.Seq.Peaks==0))))))


ESR1_pyCC_4_ER1_overlaping = seq.int(1, (ESR1_pyCC_4_total-(as.numeric(length(which(ESR1_pyCC_4_ER1$Number.of.Overlaps.with.ChIP.Seq.Peaks==0))))))
ESR1_pyCC_4_ER2_overlaping = seq.int(1, (ESR1_pyCC_4_total-(as.numeric(length(which(ESR1_pyCC_4_ER2$Number.of.Overlaps.with.ChIP.Seq.Peaks==0))))))
ESR1_pyCC_4_ER3_overlaping = seq.int(1, (ESR1_pyCC_4_total-(as.numeric(length(which(ESR1_pyCC_4_ER3$Number.of.Overlaps.with.ChIP.Seq.Peaks==0))))))
ESR1_pyCC_4_ER5_overlaping = seq.int(1, (ESR1_pyCC_4_total-(as.numeric(length(which(ESR1_pyCC_4_ER5$Number.of.Overlaps.with.ChIP.Seq.Peaks==0))))))
ESR1_pyCC_4_ER6_overlaping = seq.int(1, (ESR1_pyCC_4_total-(as.numeric(length(which(ESR1_pyCC_4_ER6$Number.of.Overlaps.with.ChIP.Seq.Peaks==0))))))

#Peakset 2

venn_list <- list(
  Rep1 = ESR1_pyCC_2_ER1_overlaping,
  Rep2 = ESR1_pyCC_2_ER2_overlaping,
  Rep3 = ESR1_pyCC_2_ER3_overlaping,
  Rep5 = ESR1_pyCC_2_ER5_overlaping,
  Rep6 = ESR1_pyCC_2_ER6_overlaping
)

# Plot the Venn diagram
Peak2 <- ggVennDiagram(venn_list, label = "count", label_alpha = 0, label_color = "white")
Peak2 + scale_fill_distiller(palette = "RdBu") + labs(title = "ESR1-HyPB Peakset 2 Overlap with each replicate") + theme(plot.title = element_text(hjust = 0.5))




#Peakset 4

venn_list <- list(
  Rep1 = ESR1_pyCC_4_ER1_overlaping,
  Rep2 = ESR1_pyCC_4_ER2_overlaping,
  Rep3 = ESR1_pyCC_4_ER3_overlaping,
  Rep5 = ESR1_pyCC_4_ER5_overlaping,
  Rep6 = ESR1_pyCC_4_ER6_overlaping
)

# Plot the Venn diagram
Peak4 <- ggVennDiagram(venn_list, label = "count", label_alpha = 0, label_color = "white")
Peak4 + scale_fill_distiller(palette = "RdBu") + labs(title = "ESR1-HyPB Peakset 4 Overlap with each replicate") + theme(plot.title = element_text(hjust = 0.5))


########################################################################
