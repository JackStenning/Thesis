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
pyCCcol_names = c("Chr",	"Start",	"End",	"Center",	"Experiment Insertions",	"Background insertions",	"Reference Insertions",	"pvalue Reference",	"pvalue Background",	"Fraction Experiment",	"TPH Experiment",	"Fraction background",	"TPH background",	"TPH background subtracted",	"pvalue_adj Reference", "Number of Overlaps with ChIA-PET Peaks")
ChIP_names = c("Chr",	"Start",	"End",	"Name",	"Score",	"Strand",	"Signal Value",	"pvalue",	"qvalue",	"Peak", "Number of Overlaps with ChIA-PET Peaks")

##Read data
#Variables named <Window '-a'><Window'-b'>
#Encode Rep 2 was taken as it has more peaks than replicate 1 - even though there is less overlap with my data
ESR1_pyCC_2_anchor1 = read.table("overlap/cc_Chia/count_Chia_peak_data_ER_WT2_Chiapet_hg38_combined.bed_anchor1.txt", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_2_anchor2 = read.table("overlap/cc_Chia/count_Chia_peak_data_ER_WT2_Chiapet_hg38_combined.bed_anchor2.txt", header = TRUE, sep = "\t", col.names = pyCCcol_names)


ESR1_pyCC_4_anchor1 = read.table("overlap/cc_Chia/count_Chia_peak_data_ER_WT4_Chiapet_hg38_combined.bed_anchor1.txt", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_4_anchor2 = read.table("overlap/cc_Chia/count_Chia_peak_data_ER_WT4_Chiapet_hg38_combined.bed_anchor2.txt", header = TRUE, sep = "\t", col.names = pyCCcol_names)

GEO_anchor1 = read.table("overlap/chip_Chia/count_Chia_Andy_ChIP_q0.05_peaksNP.bedChiapet_hg38_combined.bed_anchor1.txt", header = TRUE, sep = "\t", col.names = ChIP_names)
GEO_anchor2 = read.table("overlap/chip_Chia/count_Chia_Andy_ChIP_q0.05_peaksNP.bedChiapet_hg38_combined.bed_anchor2.txt", header = TRUE, sep = "\t", col.names = ChIP_names)


Encode_Anchor1 = read.table("overlap/chip_Chia/count_Chia_Encode_ER_rep2_ChIP_q0.05_peaksNP.bedChiapet_hg38_combined.bed_anchor1.txt", header = TRUE, sep = "\t", col.names = ChIP_names)
Encode_Anchor2 = read.table("overlap/chip_Chia/count_Chia_Encode_ER_rep2_ChIP_q0.05_peaksNP.bedChiapet_hg38_combined.bed_anchor2.txt", header = TRUE, sep = "\t", col.names = ChIP_names)


#Counting the total reads

ESR1_pyCC_2_total=nrow(ESR1_pyCC_2_anchor1)

ESR1_pyCC_4_total=nrow(ESR1_pyCC_4_anchor1)

GEO_total=nrow(GEO_anchor1)

Encode_total=nrow(Encode_Anchor1)

#Count the number of reads with no overlap
ESR1_pyCC_2_anchor1_overlap = as.numeric(length(which(ESR1_pyCC_2_anchor1$Number.of.Overlaps.with.ChIA.PET.Peaks==0)))
ESR1_pyCC_2_anchor2_overlap = as.numeric(length(which(ESR1_pyCC_2_anchor2$Number.of.Overlaps.with.ChIA.PET.Peaks==0)))



ESR1_pyCC_4_anchor1_overlap = as.numeric(length(which(ESR1_pyCC_4_anchor1$Number.of.Overlaps.with.ChIA.PET.Peaks==0)))
ESR1_pyCC_4_anchor2_overlap = as.numeric(length(which(ESR1_pyCC_4_anchor2$Number.of.Overlaps.with.ChIA.PET.Peaks==0)))

GEO_anchor1_overlap = as.numeric(length(which(GEO_anchor1$Number.of.Overlaps.with.ChIA.PET.Peaks==0)))
GEO_anchor2_overlap = as.numeric(length(which(GEO_anchor2$Number.of.Overlaps.with.ChIA.PET.Peaks==0)))



Encode_Anchor1_overlap = as.numeric(length(which(Encode_Anchor1$Number.of.Overlaps.with.ChIA.PET.Peaks==0)))
Encode_Anchor2_overlap = as.numeric(length(which(Encode_Anchor2$Number.of.Overlaps.with.ChIA.PET.Peaks==0)))


#Determine overlap percentage
Percentage_ESR1_pyCC_2_anchor1 =(1-(ESR1_pyCC_2_anchor1_overlap/ESR1_pyCC_2_total))*100
Percentage_ESR1_pyCC_2_anchor2 =(1-(ESR1_pyCC_2_anchor2_overlap/ESR1_pyCC_2_total))*100



Percentage_ESR1_pyCC_4_anchor1 =(1-(ESR1_pyCC_4_anchor1_overlap/ESR1_pyCC_4_total))*100
Percentage_ESR1_pyCC_4_anchor2 = (1-(ESR1_pyCC_4_anchor2_overlap/ESR1_pyCC_4_total))*100

Percentage_GEO_anchor1 =(1-(GEO_anchor1_overlap/GEO_total))*100
Percentage_GEO_anchor2 =(1-(GEO_anchor2_overlap/GEO_total))*100



Percentage_Encode_Anchor1 =(1-(Encode_Anchor1_overlap/Encode_total))*100
Percentage_Encode_Anchor2 = (1-(Encode_Anchor2_overlap/Encode_total))*100

#Create table to plot
#Create table to plot
pyCC_Finalpercentage_table = data.frame(Name=c('Peakset_2_Anchor1', 'Peakset_2_Anchor2', 'Peakset_4_Anchor1', 'Peakset_4_Anchor2', 'GEO_anchor1', 'GEO_anchor2', 'Encode_Anchor1', 'Encode_Anchor2'),
                                       Total=rep(c(ESR1_pyCC_2_total, ESR1_pyCC_2_total, ESR1_pyCC_4_total, ESR1_pyCC_4_total, GEO_total, GEO_total, Encode_total, Encode_total)),
                                       Overlap_withChIP= c(ESR1_pyCC_2_anchor1_overlap, ESR1_pyCC_2_anchor2_overlap, ESR1_pyCC_4_anchor1_overlap, ESR1_pyCC_4_anchor2_overlap, GEO_anchor1_overlap, GEO_anchor2_overlap, Encode_Anchor1_overlap, Encode_Anchor2_overlap),
                                       PercentageOverlapping= c(Percentage_ESR1_pyCC_2_anchor1, Percentage_ESR1_pyCC_2_anchor2, Percentage_ESR1_pyCC_4_anchor1, Percentage_ESR1_pyCC_4_anchor2,Percentage_GEO_anchor1, Percentage_GEO_anchor2, Percentage_Encode_Anchor1, Percentage_Encode_Anchor2))

#Then we plot
graphorder= c('Peakset_2_Anchor1', 'Peakset_4_Anchor1', 'GEO_anchor1', 'Encode_Anchor1', 'Peakset_2_Anchor2','Peakset_4_Anchor2', 'GEO_anchor2', 'Encode_Anchor2')
tiff("pyCalling Cards Peaks AND ChIP with ChIA-PET Overlap across all anchors.tiff", width = 2481, height = 2000, res = 300)

pyCC_Figure <- ggplot(data = pyCC_Finalpercentage_table, 
                      aes(x = factor(Name, graphorder), 
                          y = PercentageOverlapping,
                          fill = Name)) +  
  geom_col(color = "black", show.legend = FALSE) +  
  ylim(0, 100) + 
  labs(x = "Name of CC peak file and ChIA-PET Dataset anchor", 
       y = "% of CC peaks overlapping with ChIA-PET anchor",
       title = "pyCalling Cards Peaks and ChIA-PET Overlap",
       subtitle = "Determining percentage of overlapping peaks from pyCalling Card ESR1 peak files with either ChIA-PET \nAnchor within 1 kbp") +  
  theme_minimal(base_size = 15) +  # increased from 14 → 16
  geom_text(aes(label = round(PercentageOverlapping, 1)), 
            vjust = -0.5, 
            size = 4) +
  coord_cartesian(expand = FALSE, clip = "off") +
  theme(axis.text.x = element_text(angle = 90),
        axis.title.y = element_text(size = 15))   # also +2

pyCC_Figure
dev.off()
########################################################################

#### plot across both anchors
# Merge anchor1 and anchor2 into one dataframe per replicate
ESR1_pyCC_2_combo <- data.frame(
  Anchor1 = ESR1_pyCC_2_anchor1$Number.of.Overlaps.with.ChIA.PET.Peaks,
  Anchor2 = ESR1_pyCC_2_anchor2$Number.of.Overlaps.with.ChIA.PET.Peaks
)

ESR1_pyCC_4_combo <- data.frame(
  Anchor1 = ESR1_pyCC_4_anchor1$Number.of.Overlaps.with.ChIA.PET.Peaks,
  Anchor2 = ESR1_pyCC_4_anchor2$Number.of.Overlaps.with.ChIA.PET.Peaks
)

GEO_combo <- data.frame(
  Anchor1 = GEO_anchor1$Number.of.Overlaps.with.ChIA.PET.Peaks,
  Anchor2 = GEO_anchor2$Number.of.Overlaps.with.ChIA.PET.Peaks
)

Encode_combo <- data.frame(
  Anchor1 = Encode_Anchor1$Number.of.Overlaps.with.ChIA.PET.Peaks,
  Anchor2 = Encode_Anchor2$Number.of.Overlaps.with.ChIA.PET.Peaks
)

# Peak overlaps if Anchor1 > 0 OR Anchor2 > 0
ESR1_pyCC_2_combo$AnyAnchor <- (ESR1_pyCC_2_combo$Anchor1 > 0 | ESR1_pyCC_2_combo$Anchor2 > 0)
ESR1_pyCC_4_combo$AnyAnchor <- (ESR1_pyCC_4_combo$Anchor1 > 0 | ESR1_pyCC_4_combo$Anchor2 > 0)
GEO_combo$AnyAnchor <- (GEO_combo$Anchor1 > 0 | GEO_combo$Anchor2 > 0)
Encode_combo$AnyAnchor <- (Encode_combo$Anchor1 > 0 | Encode_combo$Anchor2 > 0)

# Rep 2
total2 <- nrow(ESR1_pyCC_2_combo)
overlap2 <- sum(ESR1_pyCC_2_combo$AnyAnchor)
perc2 <- (overlap2 / total2) * 100

# Rep 4
total4 <- nrow(ESR1_pyCC_4_combo)
overlap4 <- sum(ESR1_pyCC_4_combo$AnyAnchor)
perc4 <- (overlap4 / total4) * 100

# GEO
total2_2 <- nrow(GEO_combo)
overlap2_2 <- sum(GEO_combo$AnyAnchor)
perc2_2 <- (overlap2_2 / total2_2) * 100

# Encode
total4_2 <- nrow(Encode_combo)
overlap4_2 <- sum(Encode_combo$AnyAnchor)
perc4_2 <- (overlap4_2 / total4_2) * 100

# Build summary table
pyCC_AnyAnchor_table <- data.frame(
  Replicate = c("Peakset_2", "Peakset_4", "GEO", "Encode"),
  Total = c(total2, total4, total2_2, total4_2),
  Overlap_withAnyAnchor = c(overlap2, overlap4, overlap2_2, overlap4_2),
  PercentageOverlapping = c(perc2, perc4, perc2_2, perc4_2)
)
graphorder= c("Peakset_2", "Peakset_4", "GEO", "Encode")

tiff("pyCalling Cards Peaks AND ChIP with ChIA-PET Overlap across both anchors1 nd 2.tiff",
     width = 2481, height = 1749, res = 300)

ggplot(pyCC_AnyAnchor_table,
       aes(x = factor(Replicate, graphorder),
           y = PercentageOverlapping,
           fill = Replicate)) +
  geom_col(color = "black", show.legend = FALSE) +
  ylim(0, 100) +
  labs(x = "Replicate",
       y = "% of CC peaks overlapping with ChIA-PET (any anchor)",
       title = "pyCalling Cards Peaks and ChIA-PET Overlap",
       subtitle = "Percentage of CC peaks overlapping with either Anchor1 or Anchor2") +
  theme_minimal(base_size = 14) +   # ← +2 pt text everywhere
  geom_text(aes(label = round(PercentageOverlapping, 1)),
            vjust = -0.5,
            size = 5)                # labels also +2 (was 4)

dev.off()
