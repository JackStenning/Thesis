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
ChIP_names = c("Chr",	"Start",	"End",	"Name",	"Score",	"Strand",	"Signal Value",	"pvalue",	"qvalue",	"Peak", "Number of Overlaps with ChIA-PET Peaks")

##Read data
#Variables named <Window '-a'><Window'-b'>
#Encode Rep 2 was taken as it has more peaks than replicate 1 - even though there is less overlap with my data
GEO_anchor1 = read.table("overlap/chip_Chia/count_Chia_Andy_ChIP_q0.05_peaksNP.bedChiapet_hg38_combined.bed_anchor1.txt", header = TRUE, sep = "\t", col.names = ChIP_names)
GEO_anchor2 = read.table("overlap/chip_Chia/count_Chia_Andy_ChIP_q0.05_peaksNP.bedChiapet_hg38_combined.bed_anchor2.txt", header = TRUE, sep = "\t", col.names = ChIP_names)


Encode_Anchor1 = read.table("overlap/chip_Chia/count_Chia_Encode_ER_rep2_ChIP_q0.05_peaksNP.bedChiapet_hg38_combined.bed_anchor1.txt", header = TRUE, sep = "\t", col.names = ChIP_names)
Encode_Anchor2 = read.table("overlap/chip_Chia/count_Chia_Encode_ER_rep2_ChIP_q0.05_peaksNP.bedChiapet_hg38_combined.bed_anchor2.txt", header = TRUE, sep = "\t", col.names = ChIP_names)


#Counting the total reads

GEO_total=nrow(GEO_anchor1)

Encode_total=nrow(Encode_Anchor1)


#Count the number of reads with no overlap
GEO_anchor1_overlap = as.numeric(length(which(GEO_anchor1$Number.of.Overlaps.with.ChIA.PET.Peaks==0)))
GEO_anchor2_overlap = as.numeric(length(which(GEO_anchor2$Number.of.Overlaps.with.ChIA.PET.Peaks==0)))



Encode_Anchor1_overlap = as.numeric(length(which(Encode_Anchor1$Number.of.Overlaps.with.ChIA.PET.Peaks==0)))
Encode_Anchor2_overlap = as.numeric(length(which(Encode_Anchor2$Number.of.Overlaps.with.ChIA.PET.Peaks==0)))


#Determine overlap percentage
Percentage_GEO_anchor1 =(1-(GEO_anchor1_overlap/GEO_total))*100
Percentage_GEO_anchor2 =(1-(GEO_anchor2_overlap/GEO_total))*100



Percentage_Encode_Anchor1 =(1-(Encode_Anchor1_overlap/Encode_total))*100
Percentage_Encode_Anchor2 = (1-(Encode_Anchor2_overlap/Encode_total))*100

#Create table to plot
#Create table to plot
pyCC_Finalpercentage_table = data.frame(Name=c('GEO_anchor1', 'GEO_anchor2', 'Encode_Anchor1', 'Encode_Anchor2'),
                                       Total=rep(c(GEO_total, GEO_total, Encode_total, Encode_total)),
                                       Overlap_withChIP= c(GEO_anchor1_overlap, GEO_anchor2_overlap, Encode_Anchor1_overlap, Encode_Anchor2_overlap),
                                       PercentageOverlapping= c(Percentage_GEO_anchor1, Percentage_GEO_anchor2, Percentage_Encode_Anchor1, Percentage_Encode_Anchor2))

#Then we plot
graphorder= c('GEO_anchor1', 'GEO_anchor2', 'Encode_Anchor1', 'Encode_Anchor2')

tiff("ChIP-Seq and ChIA-PET Overlap.tiff", width = 3510, height = 2481, res = 300) 
pyCC_Figure <- ggplot(data = pyCC_Finalpercentage_table, 
                     aes(x= factor(Name, graphorder), 
                         y = PercentageOverlapping,
                         fill = Name, 
                         stat("identity"))) +  geom_col(color = "black",
                                                        show.legend = FALSE) +  ylim(0, 100) + labs(x = "Name of ChIP peak file and ChIA-PET Dataset anchor", 
                                                                                                    y = "% of CC peaks overlapping with ChIA-PET anchor",
                                                                                                    title = "ChIP-Seq and ChIA-PET Overlap",
                                                                                                    subtitle = "Determining percentage of overlapping peaks from ER ChIP-Seq peak files with either Anchor1 or \nAnchor2 of ER ChIA-PET within 1 kbp") +  theme_minimal() + geom_text(aes(label = round(PercentageOverlapping, 1)), # Add percentage as text
                                                                                                                                                                                                                                                                             vjust = -0.5, # Adjust vertical position above the bars
                                                                                                                                                                                                                                                                             size = 4) + coord_cartesian(expand = FALSE, clip = "off") + geom_hline(yintercept = 83, linetype = "dashed") 


pyCC_Figure + theme(axis.text.x = element_text(angle = 90))
dev.off()
########################################################################

#### plot across both anchors
# Merge anchor1 and anchor2 into one dataframe per replicate
GEO_combined <- data.frame(
  Anchor1 = GEO_anchor1$Number.of.Overlaps.with.ChIA.PET.Peaks,
  Anchor2 = GEO_anchor2$Number.of.Overlaps.with.ChIA.PET.Peaks
)

Encode_combined <- data.frame(
  Anchor1 = Encode_Anchor1$Number.of.Overlaps.with.ChIA.PET.Peaks,
  Anchor2 = Encode_Anchor2$Number.of.Overlaps.with.ChIA.PET.Peaks
)

# Peak overlaps if Anchor1 > 0 OR Anchor2 > 0
GEO_combined$AnyAnchor <- (GEO_combined$Anchor1 > 0 | GEO_combined$Anchor2 > 0)
Encode_combined$AnyAnchor <- (Encode_combined$Anchor1 > 0 | Encode_combined$Anchor2 > 0)

# GEO
total2_2 <- nrow(GEO_combined)
overlap2_2 <- sum(GEO_combined$AnyAnchor)
perc2_2 <- (overlap2_2 / total2_2) * 100

# Encode
total4_2 <- nrow(Encode_combined)
overlap4_2 <- sum(Encode_combined$AnyAnchor)
perc4_2 <- (overlap4_2 / total4_2) * 100

# Build summary table
pyCC_AnyAnchor_table <- data.frame(
  Replicate = c("GEO", "Encode"),
  Total = c(total2_2, total4_2),
  Overlap_withAnyAnchor = c(overlap2_2, overlap4_2),
  PercentageOverlapping = c(perc2_2, perc4_2)
)

tiff("ChIP-Seq and ChIA-PET Overlap at both anchors.tiff", width = 3510, height = 2481, res = 300) 
ggplot(pyCC_AnyAnchor_table, aes(x = Replicate, y = PercentageOverlapping, fill = Replicate)) +
  geom_col(color = "black", show.legend = FALSE) +
  ylim(0, 100) +
  labs(x = "Replicate",
       y = "% of CC peaks overlapping with ChIA-PET (any anchor)",
       title = "ChIP-Seq and ChIA-PET Overlap at both anchors",
       subtitle = "Percentage of CC peaks overlapping with either Anchor1 or Anchor2") +
  theme_minimal() +
  geom_text(aes(label = round(PercentageOverlapping, 1)), vjust = -0.5, size = 4)
dev.off()