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
pyCCcol_names = c("Chr",	"Start",	"End",	"Center",	"Experiment Insertions",	"Background insertions",	"Reference Insertions",	"pvalue Reference",	"pvalue Background",	"Fraction Experiment",	"TPH Experiment",	"Fraction background",	"TPH background",	"TPH background subtracted",	"pvalue_adj Reference", "Number of Overlaps with ER ChIP-Seq Peaks", "Number of Overlaps with ChIA-PET Peaks")

##Read data
#Variables named <Window '-a'><Window'-b'>
#Encode Rep 2 was taken as it has more peaks than replicate 1 - even though there is less overlap with my data
ESR1_pyCC_2_encode_Anchor1 = read.table("overlap/cc_Chia/count_Chia_Peakset2_Encode_no_chipseq_overlaps_Chiapet_hg38_combined.bed_anchor1.txt", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_2_encode_Anchor2 = read.table("overlap/cc_Chia/count_Chia_Peakset2_Encode_no_chipseq_overlaps_Chiapet_hg38_combined.bed_anchor2.txt", header = TRUE, sep = "\t", col.names = pyCCcol_names)


ESR1_pyCC_2_GEO_Anchor1 = read.table("overlap/cc_Chia/count_Chia_Peakset2_GEO_no_chipseq_overlaps_Chiapet_hg38_combined.bed_anchor1.txt", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_2_GEO_Anchor2 = read.table("overlap/cc_Chia/count_Chia_Peakset2_GEO_no_chipseq_overlaps_Chiapet_hg38_combined.bed_anchor2.txt", header = TRUE, sep = "\t", col.names = pyCCcol_names)


#Counting the total reads

ESR1_pyCC_2_total=nrow(ESR1_pyCC_2_encode_Anchor1)

ESR1_pyCC_4_total=nrow(ESR1_pyCC_2_GEO_Anchor1)


#Count the number of reads with no overlap
ESR1_pyCC_2_encode_Anchor1_overlap = as.numeric(length(which(ESR1_pyCC_2_encode_Anchor1$Number.of.Overlaps.with.ChIA.PET.Peaks==0)))
ESR1_pyCC_2_encode_Anchor2_overlap = as.numeric(length(which(ESR1_pyCC_2_encode_Anchor2$Number.of.Overlaps.with.ChIA.PET.Peaks==0)))



ESR1_pyCC_2_GEO_Anchor1_overlap = as.numeric(length(which(ESR1_pyCC_2_GEO_Anchor1$Number.of.Overlaps.with.ChIA.PET.Peaks==0)))
ESR1_pyCC_4_Anchor2_overlap = as.numeric(length(which(ESR1_pyCC_2_GEO_Anchor2$Number.of.Overlaps.with.ChIA.PET.Peaks==0)))


#Determine overlap percentage
Percentage_ESR1_pyCC_2_encode_Anchor1 =(1-(ESR1_pyCC_2_encode_Anchor1_overlap/ESR1_pyCC_2_total))*100
Percentage_ESR1_pyCC_2_encode_Anchor2 =(1-(ESR1_pyCC_2_encode_Anchor2_overlap/ESR1_pyCC_2_total))*100



Percentage_ESR1_pyCC_2_GEO_Anchor1 =(1-(ESR1_pyCC_2_GEO_Anchor1_overlap/ESR1_pyCC_4_total))*100
Percentage_ESR1_pyCC_4_Anchor2 = (1-(ESR1_pyCC_4_Anchor2_overlap/ESR1_pyCC_4_total))*100

#Create table to plot
#Create table to plot
pyCC_Finalpercentage_table = data.frame(Name=c('Peakset_2_Encode_Anchor1', 'Peakset_2_Encode_Anchor2', 'Peakset_2_GEO_Anchor1', 'Peakset_2_GEO_Anchor2'),
                                       Total=rep(c(ESR1_pyCC_2_total, ESR1_pyCC_2_total, ESR1_pyCC_4_total, ESR1_pyCC_4_total)),
                                       Overlap_withChIP= c(ESR1_pyCC_2_encode_Anchor1_overlap, ESR1_pyCC_2_encode_Anchor2_overlap, ESR1_pyCC_2_GEO_Anchor1_overlap, ESR1_pyCC_4_Anchor2_overlap),
                                       PercentageOverlapping= c(Percentage_ESR1_pyCC_2_encode_Anchor1, Percentage_ESR1_pyCC_2_encode_Anchor2, Percentage_ESR1_pyCC_2_GEO_Anchor1, Percentage_ESR1_pyCC_4_Anchor2))

#Then we plot
graphorder= c('Peakset_2_Encode_Anchor1', 'Peakset_2_Encode_Anchor2', 'Peakset_2_GEO_Anchor1', 'Peakset_2_GEO_Anchor2')
tiff("High stringency non-overlapping peaks pyCalling Cards Peaks and ChIA-PET Overlap.tiff", width = 3510, height = 2481, res = 300)
pyCC_Figure <- ggplot(data = pyCC_Finalpercentage_table, 
                     aes(x= factor(Name, graphorder), 
                         y = PercentageOverlapping,
                         fill = Name, 
                         stat("identity"))) +  geom_col(color = "black",
                                                        show.legend = FALSE) +  ylim(0, 100) + labs(x = "Name of CC peak file and ChIP Dataset", 
                                                                                                    y = "% of CC peaks overlapping with ChIP-Seq peak",
                                                                                                    title = "pyCalling Cards Peaks and ChIP-Seq Overlap",
                                                                                                    subtitle = "Determining percentage of overlapping peaks from pyCalling Card ESR1 high stringency peak files with either Anchor1 or \nZhnag BRD4 ChIP-Seq datasets within 1 kbp") +  theme_minimal() + geom_text(aes(label = round(PercentageOverlapping, 1)), # Add percentage as text
                                                                                                                                                                                                                                                                             vjust = -0.5, # Adjust vertical position above the bars
                                                                                                                                                                                                                                                                             size = 4) + coord_cartesian(expand = FALSE, clip = "off") + geom_hline(yintercept = 83, linetype = "dashed") 


pyCC_Figure + theme(axis.text.x = element_text(angle = 90))
dev.off()
########################################################################
#### plot across both anchors
# Merge anchor1 and anchor2 into one dataframe per replicate
ESR1_pyCC_2_combined <- data.frame(
  Anchor1 = ESR1_pyCC_2_encode_Anchor1$Number.of.Overlaps.with.ChIA.PET.Peaks,
  Anchor2 = ESR1_pyCC_2_encode_Anchor2$Number.of.Overlaps.with.ChIA.PET.Peaks
)

ESR1_pyCC_4_combined <- data.frame(
  Anchor1 = ESR1_pyCC_2_GEO_Anchor1$Number.of.Overlaps.with.ChIA.PET.Peaks,
  Anchor2 = ESR1_pyCC_2_GEO_Anchor2$Number.of.Overlaps.with.ChIA.PET.Peaks
)

# Peak overlaps if Anchor1 > 0 OR Anchor2 > 0
ESR1_pyCC_2_combined$AnyAnchor <- (ESR1_pyCC_2_combined$Anchor1 > 0 | ESR1_pyCC_2_combined$Anchor2 > 0)
ESR1_pyCC_4_combined$AnyAnchor <- (ESR1_pyCC_4_combined$Anchor1 > 0 | ESR1_pyCC_4_combined$Anchor2 > 0)

# Rep 2
total2 <- nrow(ESR1_pyCC_2_combined)
overlap2 <- sum(ESR1_pyCC_2_combined$AnyAnchor)
perc2 <- (overlap2 / total2) * 100

# Rep 4
total4 <- nrow(ESR1_pyCC_4_combined)
overlap4 <- sum(ESR1_pyCC_4_combined$AnyAnchor)
perc4 <- (overlap4 / total4) * 100

# Build summary table
pyCC_AnyAnchor_table <- data.frame(
  Replicate = c("Peakset_2", "Peakset_4"),
  Total = c(total2, total4),
  Overlap_withAnyAnchor = c(overlap2, overlap4),
  PercentageOverlapping = c(perc2, perc4)
)

tiff("High stringency non-overlapping peaks pyCalling Cards Peaks and ChIA-PET Overlap both anchors.tiff", width = 3510, height = 2481, res = 300)
ggplot(pyCC_AnyAnchor_table, aes(x = Replicate, y = PercentageOverlapping, fill = Replicate)) +
  geom_col(color = "black", show.legend = FALSE) +
  ylim(0, 100) +
  labs(x = "Replicate",
       y = "% of CC peaks overlapping with ChIA-PET (any anchor)",
       title = "pyCalling Cards Peaks and ChIA-PET Overlap",
       subtitle = "Percentage of CC peaks overlapping with either Anchor1 or Anchor2 that do not overlap with ChIP-Seq") +
  theme_minimal() +
  geom_text(aes(label = round(PercentageOverlapping, 1)), vjust = -0.5, size = 4)
dev.off()

###################################################### PLOT FOR LOW STRINGENCY PEAKS
#Set column names
pyCCcol_names = c("Chr",	"Start",	"End",	"Center",	"Experiment Insertions",	"Background insertions",	"Reference Insertions",	"pvalue Reference",	"pvalue Background",	"Fraction Experiment",	"TPH Experiment",	"Fraction background",	"TPH background",	"TPH background subtracted",	"pvalue_adj Reference", "Number of Overlaps with ER ChIP-Seq Peaks", "Number of Overlaps with ChIA-PET Peaks")

##Read data
#Variables named <Window '-a'><Window'-b'>
#Encode Rep 2 was taken as it has more peaks than replicate 1 - even though there is less overlap with my data
ESR1_pyCC_2_encode_Anchor1 = read.table("overlap/cc_Chia/count_Chia_Peakset4_Encode_no_chipseq_overlaps_Chiapet_hg38_combined.bed_anchor1.txt", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_2_encode_Anchor2 = read.table("overlap/cc_Chia/count_Chia_Peakset4_Encode_no_chipseq_overlaps_Chiapet_hg38_combined.bed_anchor1.txt", header = TRUE, sep = "\t", col.names = pyCCcol_names)


ESR1_pyCC_2_GEO_Anchor1 = read.table("overlap/cc_Chia/count_Chia_Peakset4_GEO_no_chipseq_overlaps_Chiapet_hg38_combined.bed_anchor1.txt", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_2_GEO_Anchor2 = read.table("overlap/cc_Chia/count_Chia_Peakset4_GEO_no_chipseq_overlaps_Chiapet_hg38_combined.bed_anchor2.txt", header = TRUE, sep = "\t", col.names = pyCCcol_names)


#Counting the total reads

ESR1_pyCC_2_total=nrow(ESR1_pyCC_2_encode_Anchor1)

ESR1_pyCC_4_total=nrow(ESR1_pyCC_2_GEO_Anchor1)


#Count the number of reads with no overlap
ESR1_pyCC_2_encode_Anchor1_overlap = as.numeric(length(which(ESR1_pyCC_2_encode_Anchor1$Number.of.Overlaps.with.ChIA.PET.Peaks==0)))
ESR1_pyCC_2_encode_Anchor2_overlap = as.numeric(length(which(ESR1_pyCC_2_encode_Anchor2$Number.of.Overlaps.with.ChIA.PET.Peaks==0)))



ESR1_pyCC_2_GEO_Anchor1_overlap = as.numeric(length(which(ESR1_pyCC_2_GEO_Anchor1$Number.of.Overlaps.with.ChIA.PET.Peaks==0)))
ESR1_pyCC_4_Anchor2_overlap = as.numeric(length(which(ESR1_pyCC_2_GEO_Anchor2$Number.of.Overlaps.with.ChIA.PET.Peaks==0)))


#Determine overlap percentage
Percentage_ESR1_pyCC_2_encode_Anchor1 =(1-(ESR1_pyCC_2_encode_Anchor1_overlap/ESR1_pyCC_2_total))*100
Percentage_ESR1_pyCC_2_encode_Anchor2 =(1-(ESR1_pyCC_2_encode_Anchor2_overlap/ESR1_pyCC_2_total))*100



Percentage_ESR1_pyCC_2_GEO_Anchor1 =(1-(ESR1_pyCC_2_GEO_Anchor1_overlap/ESR1_pyCC_4_total))*100
Percentage_ESR1_pyCC_4_Anchor2 = (1-(ESR1_pyCC_4_Anchor2_overlap/ESR1_pyCC_4_total))*100

#Create table to plot
pyCC_Finalpercentage_table = data.frame(Name=c('Peakset_4_Encode_Anchor1', 'Peakset_4_Encode_Anchor2', 'Peakset_4_GEO_Anchor1', 'Peakset_4_GEO_Anchor2'),
                                        Total=rep(c(ESR1_pyCC_2_total, ESR1_pyCC_2_total, ESR1_pyCC_4_total, ESR1_pyCC_4_total)),
                                        Overlap_withChIP= c(ESR1_pyCC_2_encode_Anchor1_overlap, ESR1_pyCC_2_encode_Anchor2_overlap, ESR1_pyCC_2_GEO_Anchor1_overlap, ESR1_pyCC_4_Anchor2_overlap),
                                        PercentageOverlapping= c(Percentage_ESR1_pyCC_2_encode_Anchor1, Percentage_ESR1_pyCC_2_encode_Anchor2, Percentage_ESR1_pyCC_2_GEO_Anchor1, Percentage_ESR1_pyCC_4_Anchor2))

#Then we plot
graphorder= c('Peakset_4_Encode_Anchor1', 'Peakset_4_Encode_Anchor2', 'Peakset_4_GEO_Anchor1', 'Peakset_4_GEO_Anchor2')
tiff("Low stringency non-overlapping peaks pyCalling Cards Peaks and ChIA-PET Overlap.tiff", width = 3510, height = 2481, res = 300)
pyCC_Figure <- ggplot(data = pyCC_Finalpercentage_table, 
                      aes(x= factor(Name, graphorder), 
                          y = PercentageOverlapping,
                          fill = Name, 
                          stat("identity"))) +  geom_col(color = "black",
                                                         show.legend = FALSE) +  ylim(0, 100) + labs(x = "Name of CC peak file and ChIP Dataset", 
                                                                                                     y = "% of CC peaks overlapping with ChIP-Seq peak",
                                                                                                     title = "pyCalling Cards Peaks and ChIP-Seq Overlap",
                                                                                                     subtitle = "Determining percentage of overlapping peaks from pyCalling Card ESR1 Low stringency peak files with either Anchor1 or \nZhnag BRD4 ChIP-Seq datasets within 1 kbp") +  theme_minimal() + geom_text(aes(label = round(PercentageOverlapping, 1)), # Add percentage as text
                                                                                                                                                                                                                                                                                                                 vjust = -0.5, # Adjust vertical position above the bars
                                                                                                                                                                                                                                                                                                                 size = 4) + coord_cartesian(expand = FALSE, clip = "off") + geom_hline(yintercept = 83, linetype = "dashed") 


pyCC_Figure + theme(axis.text.x = element_text(angle = 90))
dev.off()
########################################################################

#### plot across both anchors
# Merge anchor1 and anchor2 into one dataframe per replicate
ESR1_pyCC_2_combined <- data.frame(
  Anchor1 = ESR1_pyCC_2_encode_Anchor1$Number.of.Overlaps.with.ChIA.PET.Peaks,
  Anchor2 = ESR1_pyCC_2_encode_Anchor2$Number.of.Overlaps.with.ChIA.PET.Peaks
)

ESR1_pyCC_4_combined <- data.frame(
  Anchor1 = ESR1_pyCC_2_GEO_Anchor1$Number.of.Overlaps.with.ChIA.PET.Peaks,
  Anchor2 = ESR1_pyCC_2_GEO_Anchor2$Number.of.Overlaps.with.ChIA.PET.Peaks
)

# Peak overlaps if Anchor1 > 0 OR Anchor2 > 0
ESR1_pyCC_2_combined$AnyAnchor <- (ESR1_pyCC_2_combined$Anchor1 > 0 | ESR1_pyCC_2_combined$Anchor2 > 0)
ESR1_pyCC_4_combined$AnyAnchor <- (ESR1_pyCC_4_combined$Anchor1 > 0 | ESR1_pyCC_4_combined$Anchor2 > 0)

# Rep 2
total2 <- nrow(ESR1_pyCC_2_combined)
overlap2 <- sum(ESR1_pyCC_2_combined$AnyAnchor)
perc2 <- (overlap2 / total2) * 100

# Rep 4
total4 <- nrow(ESR1_pyCC_4_combined)
overlap4 <- sum(ESR1_pyCC_4_combined$AnyAnchor)
perc4 <- (overlap4 / total4) * 100

# Build summary table
pyCC_AnyAnchor_table <- data.frame(
  Replicate = c("Peakset_2", "Peakset_4"),
  Total = c(total2, total4),
  Overlap_withAnyAnchor = c(overlap2, overlap4),
  PercentageOverlapping = c(perc2, perc4)
)

tiff("Low stringency non-overlapping peaks pyCalling Cards Peaks and ChIA-PET Overlap both anchors.tiff", width = 3510, height = 2481, res = 300)
ggplot(pyCC_AnyAnchor_table, aes(x = Replicate, y = PercentageOverlapping, fill = Replicate)) +
  geom_col(color = "black", show.legend = FALSE) +
  ylim(0, 100) +
  labs(x = "Replicate",
       y = "% of CC peaks overlapping with ChIA-PET (any anchor)",
       title = "pyCalling Cards Peaks and ChIA-PET Overlap",
       subtitle = "Percentage of CC peaks overlapping with either Anchor1 or Anchor2 that do not overlap with ChIP-Seq") +
  theme_minimal() +
  geom_text(aes(label = round(PercentageOverlapping, 1)), vjust = -0.5, size = 4)
dev.off()