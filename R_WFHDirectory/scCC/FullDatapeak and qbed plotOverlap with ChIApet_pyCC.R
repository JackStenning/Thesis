############# LOOK HERE ######################
#add or delete the '#' symbol below on install packages
#install.packages("Rtools")
#install.packages("ggthemes")

#DELTE AFTER WRITING

#Load software
library(ggplot2)
library(ggthemes)
library(dplyr)

#Set working directory
setwd("C:/Users/jps558/OneDrive - University of York/Desktop/R_WorkingDirectory/R_WFHDirectory/scCC") 


#Set column names
pyCCcol_names = c("Chr",	"Start",	"End",	"Experiment Insertions", "Reference Inserions",	"Expected Insertions",	"pvalue",	"Fraction Experiment",	"TPH Experiment",	"TPH background subtracted",	"pvalue_adj Reference", "Number of Overlaps with ChIA-PET Peaks")
pyCCcol_q_names = c("Chr",	"Start",	"End",	"Experiment Insertions", "Strand", "barcode", "Number of Overlaps with ChIA-PET Peaks")
##Read data
#Variables named <Window '-a'><Window'-b'>
#Encode Rep 2 was taken as it has more peaks than replicate 1 - even though there is less overlap with my data
ESR1_pyCC_full_Peak_anchor1 = read.table("FinalCCFresult/1kbp_A_FullResultpeak_data_ER_test4.bed_B_chiaanchor1.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_full_Peak_anchor2 = read.table("FinalCCFresult/1kbp_A_FullResultpeak_data_ER_test4.bed_B_chiaanchor2.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)


ESR1_pyCC_full_UMI_anchor1 = read.table("FinalCCFresult/1kbp_A_FullUMIpeak_data_ER_test4.bed_B_chiaanchor1.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_full_UMI_anchor2 = read.table("FinalCCFresult/1kbp_A_FullUMIpeak_data_ER_test4.bed_B_chiaanchor2.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)

ESR1_pyCC_full_Peak_anchor1_qbed = read.table("FinalCCFresult/1kbp_A_Fullresult_scCC_final.sorted.ccf_B_chiaanchor1.bed", header = TRUE, sep = "\t", col.names = pyCCcol_q_names)
ESR1_pyCC_full_Peak_anchor2_qbed = read.table("FinalCCFresult/1kbp_A_Fullresult_scCC_final.sorted.ccf_B_chiaanchor2.bed", header = TRUE, sep = "\t", col.names = pyCCcol_q_names)


ESR1_pyCC_full_UMI_anchor1_qbed = read.table("FinalCCFresult/1kbp_A_Fullresult_UMIFilt_scCC_final.ccf_B_chiaanchor1.bed", header = TRUE, sep = "\t", col.names = pyCCcol_q_names)
ESR1_pyCC_full_UMI_anchor2_qbed = read.table("FinalCCFresult/1kbp_A_Fullresult_UMIFilt_scCC_final.ccf_B_chiaanchor2.bed", header = TRUE, sep = "\t", col.names = pyCCcol_q_names)

#Counting the total reads

ESR1_pyCC_2_total=nrow(ESR1_pyCC_full_Peak_anchor1)

ESR1_pyCC_4_total=nrow(ESR1_pyCC_full_UMI_anchor1)

ESR1_pyCC_2_total_q=nrow(ESR1_pyCC_full_Peak_anchor1_qbed)

ESR1_pyCC_4_total_q=nrow(ESR1_pyCC_full_UMI_anchor1_qbed)


#Count the number of reads with no overlap
ESR1_pyCC_full_Peak_anchor1_overlap = as.numeric(length(which(ESR1_pyCC_full_Peak_anchor1$Number.of.Overlaps.with.ChIA.PET.Peaks==0)))
ESR1_pyCC_full_Peak_anchor2_overlap = as.numeric(length(which(ESR1_pyCC_full_Peak_anchor2$Number.of.Overlaps.with.ChIA.PET.Peaks==0)))



ESR1_pyCC_full_UMI_anchor1_overlap = as.numeric(length(which(ESR1_pyCC_full_UMI_anchor1$Number.of.Overlaps.with.ChIA.PET.Peaks==0)))
ESR1_pyCC_full_UMI_anchor2_overlap = as.numeric(length(which(ESR1_pyCC_full_UMI_anchor2$Number.of.Overlaps.with.ChIA.PET.Peaks==0)))

ESR1_pyCC_full_Peak_anchor1_qbed_overlap_q = as.numeric(length(which(ESR1_pyCC_full_Peak_anchor1_qbed$Number.of.Overlaps.with.ChIA.PET.Peaks==0)))
ESR1_pyCC_full_Encode_q_qbed_overlap_q = as.numeric(length(which(ESR1_pyCC_full_Peak_anchor2_qbed$Number.of.Overlaps.with.ChIA.PET.Peaks==0)))



ESR1_pyCC_full_UMI_anchor1_qbed_overlap_q = as.numeric(length(which(ESR1_pyCC_full_UMI_anchor1_qbed$Number.of.Overlaps.with.ChIA.PET.Peaks==0)))
ESR1_pyCC_full_UMI_anchor2_qbed_overlap_q = as.numeric(length(which(ESR1_pyCC_full_UMI_anchor2_qbed$Number.of.Overlaps.with.ChIA.PET.Peaks==0)))

#Determine overlap percentage
Percentage_ESR1_pyCC_full_Peak_anchor1 =(1-(ESR1_pyCC_full_Peak_anchor1_overlap/ESR1_pyCC_2_total))*100
Percentage_ESR1_pyCC_full_Peak_anchor2 =(1-(ESR1_pyCC_full_Peak_anchor2_overlap/ESR1_pyCC_2_total))*100



Percentage_ESR1_pyCC_full_UMI_anchor1 =(1-(ESR1_pyCC_full_UMI_anchor1_overlap/ESR1_pyCC_4_total))*100
Percentage_ESR1_pyCC_full_UMI_anchor2 = (1-(ESR1_pyCC_full_UMI_anchor2_overlap/ESR1_pyCC_4_total))*100

Percentage_ESR1_pyCC_full_Peak_anchor1_qbed_q =(1-(ESR1_pyCC_full_Peak_anchor1_qbed_overlap_q/ESR1_pyCC_2_total_q))*100
Percentage_ESR1_pyCC_full_Peak_anchor2bed_q =(1-(ESR1_pyCC_full_Encode_q_qbed_overlap_q/ESR1_pyCC_2_total_q))*100



Percentage_ESR1_pyCC_full_UMI_anchor1_qbed_q =(1-(ESR1_pyCC_full_UMI_anchor1_qbed_overlap_q/ESR1_pyCC_4_total_q))*100
Percentage_ESR1_pyCC_full_UMI_anchor2_qbed_q = (1-(ESR1_pyCC_full_UMI_anchor2_qbed_overlap_q/ESR1_pyCC_4_total_q))*100

#Create table to plot
#Create table to plot
pyCC_Finalpercentage_table <- data.frame(
  Name = c(
    'Full_Result_Peaks_Anchor1', 'Full_Result_Peaks_Anchor2', 
    'Full_Result_Peaks_UMI_Anchor1', 'Full_Result_Peaks_UMI_Anchor2',
    'Full_Result_qbed_Anchor1', 'Full_Result_qbed_Anchor2', 
    'Full_Result_qbed_UMI_Anchor1', 'Full_Result_qbed_UMI_Anchor2'
  ),
  Total = c(
    ESR1_pyCC_2_total, ESR1_pyCC_4_total, 
    ESR1_pyCC_2_total_q, ESR1_pyCC_4_total_q,
    ESR1_pyCC_2_total, ESR1_pyCC_4_total, 
    ESR1_pyCC_2_total_q, ESR1_pyCC_4_total_q
  ),
  Overlap_withChIA = c(
    ESR1_pyCC_full_Peak_anchor1_overlap, ESR1_pyCC_full_UMI_anchor1_overlap, 
    ESR1_pyCC_full_Peak_anchor1_qbed_overlap_q, ESR1_pyCC_full_UMI_anchor1_qbed_overlap_q,
    ESR1_pyCC_full_Peak_anchor2_overlap, ESR1_pyCC_full_UMI_anchor2_overlap, 
    ESR1_pyCC_full_Encode_q_qbed_overlap_q, ESR1_pyCC_full_UMI_anchor2_qbed_overlap_q
  ),
  PercentageOverlapping = c(
    Percentage_ESR1_pyCC_full_Peak_anchor1, Percentage_ESR1_pyCC_full_UMI_anchor1, 
    Percentage_ESR1_pyCC_full_Peak_anchor1_qbed_q, Percentage_ESR1_pyCC_full_UMI_anchor1_qbed_q,
    Percentage_ESR1_pyCC_full_Peak_anchor2, Percentage_ESR1_pyCC_full_UMI_anchor2, 
    Percentage_ESR1_pyCC_full_Peak_anchor2bed_q, Percentage_ESR1_pyCC_full_UMI_anchor2_qbed_q
  )
)

#Then we plot
graphorder= c(
  'Full_Result_Peaks_Anchor1', 'Full_Result_Peaks_Anchor2', 
  'Full_Result_Peaks_UMI_Anchor1', 'Full_Result_Peaks_UMI_Anchor2',
  'Full_Result_qbed_Anchor1', 'Full_Result_qbed_Anchor2', 
  'Full_Result_qbed_UMI_Anchor1', 'Full_Result_qbed_UMI_Anchor2'
)
pyCC_Figure <- ggplot(data = pyCC_Finalpercentage_table, 
                      aes(x= factor(Name, graphorder), 
                          y = PercentageOverlapping,
                          fill = Name, 
                          stat("identity"))) +  geom_col(color = "black",
                                                         show.legend = FALSE) +  ylim(0, 100) + labs(x = "Name of CC peak file and ChIA Dataset", 
                                                                                                     y = "% of CC peaks overlapping with ChIA-PET peak",
                                                                                                     title = "Full result - single cell pyCalling Cards Peaks and ChIA-PET Overlap",
                                                                                                     subtitle = "Determining percentage of overlapping peaks from pyCalling Card single-cell ESR1 peak files with ChIA-PET datasets within 1 kbp") +  theme_minimal() + geom_text(aes(label = round(PercentageOverlapping, 1)), # Add percentage as text
                                                                                                                                                                                                                                                                                                                 vjust = -0.5, # Adjust vertical position above the bars
                                                                                                                                                                                                                                                                                                                 size = 4) + coord_cartesian(expand = FALSE, clip = "off") + geom_hline(yintercept = 83, linetype = "dashed") 


pyCC_Figure + theme(axis.text.x = element_text(angle = 90))
########################################################################

#### plot across both anchors
# Merge anchor1 and anchor2 into one dataframe per replicate
FullPeaks <- data.frame(
  Anchor1 = ESR1_pyCC_full_Peak_anchor1$Number.of.Overlaps.with.ChIA.PET.Peaks,
  Anchor2 = ESR1_pyCC_full_Peak_anchor2$Number.of.Overlaps.with.ChIA.PET.Peaks
)

FullPeaks_UMI <- data.frame(
  Anchor1 = ESR1_pyCC_full_UMI_anchor1$Number.of.Overlaps.with.ChIA.PET.Peaks,
  Anchor2 = ESR1_pyCC_full_UMI_anchor2$Number.of.Overlaps.with.ChIA.PET.Peaks
)

FullPeaks_qbed <- data.frame(
  Anchor1 = ESR1_pyCC_full_Peak_anchor1_qbed$Number.of.Overlaps.with.ChIA.PET.Peaks,
  Anchor2 = ESR1_pyCC_full_Peak_anchor2_qbed$Number.of.Overlaps.with.ChIA.PET.Peaks
)

FullPeaks_qbed_UMI <- data.frame(
  Anchor1 = ESR1_pyCC_full_UMI_anchor1_qbed$Number.of.Overlaps.with.ChIA.PET.Peaks,
  Anchor2 = ESR1_pyCC_full_UMI_anchor2_qbed$Number.of.Overlaps.with.ChIA.PET.Peaks
)

# Peak overlaps if Anchor1 > 0 OR Anchor2 > 0
FullPeaks$AnyAnchor <- (FullPeaks$Anchor1 > 0 | FullPeaks$Anchor2 > 0)
FullPeaks_UMI$AnyAnchor <- (FullPeaks_UMI$Anchor1 > 0 | FullPeaks_UMI$Anchor2 > 0)
FullPeaks_qbed$AnyAnchor <- (FullPeaks_qbed$Anchor1 > 0 | FullPeaks_qbed$Anchor2 > 0)
FullPeaks_qbed_UMI$AnyAnchor <- (FullPeaks_qbed_UMI$Anchor1 > 0 | FullPeaks_qbed_UMI$Anchor2 > 0)

# Rep 2
total2 <- nrow(FullPeaks)
overlap2 <- sum(FullPeaks$AnyAnchor)
perc2 <- (overlap2 / total2) * 100

# Rep 4
total4 <- nrow(FullPeaks_UMI)
overlap4 <- sum(FullPeaks_UMI$AnyAnchor)
perc4 <- (overlap4 / total4) * 100

# GEO
total2_2 <- nrow(FullPeaks_qbed)
overlap2_2 <- sum(FullPeaks_qbed$AnyAnchor)
perc2_2 <- (overlap2_2 / total2_2) * 100

# Encode
total4_2 <- nrow(FullPeaks_qbed_UMI)
overlap4_2 <- sum(FullPeaks_qbed_UMI$AnyAnchor)
perc4_2 <- (overlap4_2 / total4_2) * 100

# Build summary table
pyCC_AnyAnchor_table <- data.frame(
  Replicate = c("Full_result_Peaks", "Full_result_Peaks_UMI", "Full_result_qbed", "Full_result_qbed_UMI"),
  Total = c(total2, total4, total2_2, total4_2),
  Overlap_withAnyAnchor = c(overlap2, overlap4, overlap2_2, overlap4_2),
  PercentageOverlapping = c(perc2, perc4, perc2_2, perc4_2)
)

graphorder= c("Full_result_Peaks", "Full_result_Peaks_UMI", "Full_result_qbed", "Full_result_qbed_UMI")

ggplot(pyCC_AnyAnchor_table, aes(x = factor(Replicate, graphorder), y = PercentageOverlapping, fill = Replicate)) +
  geom_col(color = "black", show.legend = FALSE) +
  ylim(0, 100) +
  labs(x = "Replicate",
       y = "% of CC peaks overlapping with ChIA-PET (any anchor)",
       title = "pyCalling Cards Peaks and ChIA-PET Overlap",
       subtitle = "Percentage of CC peaks overlapping with either Anchor1 or Anchor2") +
  theme_minimal() +
  geom_text(aes(label = round(PercentageOverlapping, 1)), vjust = -0.5, size = 4)



# Remove multiple entries in qbed
# Function to collapse qBED by genomic coordinates
collapse_qbed <- function(df) {
  df %>%
    group_by(Chr, Start, End) %>%
    summarise(
      `Experiment Insertions` = sum(`Experiment.Insertions`),  # sum insertions
      Strand = first(Strand),                                   # keep strand if needed
      barcode = paste(unique(barcode), collapse = ";"),        # optional: combine barcodes
      `Number of Overlaps with ChIA-PET Peaks` = sum(`Number.of.Overlaps.with.ChIA.PET.Peaks`)  # sum overlaps
    ) %>%
    ungroup()
}

# Collapse all qBED dfs
Peak_anchor1_qbed_coll <- collapse_qbed(ESR1_pyCC_full_Peak_anchor1_qbed)
Peak_anchor2_qbed_coll <- collapse_qbed(ESR1_pyCC_full_Peak_anchor2_qbed)
UMI_anchor1_qbed_coll <- collapse_qbed(ESR1_pyCC_full_UMI_anchor1_qbed)
UMI_anchor2_qbed_coll <- collapse_qbed(ESR1_pyCC_full_UMI_anchor2_qbed)

# Recalculate the total number of loci and overlaps
calc_percentage <- function(df) {
  total <- nrow(df)
  perc <- (1 - sum(df$`Number of Overlaps with ChIA-PET Peaks` == 0)/total) * 100
  return(perc)
}

Percentage_ESR1_pyCC_full_Peak_anchor1_qbed_q <- calc_percentage(Peak_anchor1_qbed_coll)
Percentage_ESR1_pyCC_full_Peak_anchor2_qbed_q <- calc_percentage(Peak_anchor2_qbed_coll)
Percentage_ESR1_pyCC_full_UMI_anchor1_qbed_q <- calc_percentage(UMI_anchor1_qbed_coll)
Percentage_ESR1_pyCC_full_UMI_anchor2_qbed_q <- calc_percentage(UMI_anchor2_qbed_coll)

pyCC_Finalpercentage_table_coll <- data.frame(
  Name = c(
    'Full_Result_Peaks_Anchor1', 'Full_Result_Peaks_Anchor2', 
    'Full_Result_Peaks_UMI_Anchor1', 'Full_Result_Peaks_UMI_Anchor2',
    'Full_Result_qbed_Anchor1', 'Full_Result_qbed_Anchor2', 
    'Full_Result_qbed_UMI_Anchor1', 'Full_Result_qbed_UMI_Anchor2'
  ),
  PercentageOverlapping = c(
    Percentage_ESR1_pyCC_full_Peak_anchor1, Percentage_ESR1_pyCC_full_Peak_anchor2,
    Percentage_ESR1_pyCC_full_UMI_anchor1, Percentage_ESR1_pyCC_full_UMI_anchor2,
    Percentage_ESR1_pyCC_full_Peak_anchor1_qbed_q, Percentage_ESR1_pyCC_full_Peak_anchor2_qbed_q,
    Percentage_ESR1_pyCC_full_UMI_anchor1_qbed_q, Percentage_ESR1_pyCC_full_UMI_anchor2_qbed_q
  )
)

graphorder = c(
  'Full_Result_Peaks_Anchor1', 'Full_Result_Peaks_Anchor2', 
  'Full_Result_Peaks_UMI_Anchor1', 'Full_Result_Peaks_UMI_Anchor2',
  'Full_Result_qbed_Anchor1', 'Full_Result_qbed_Anchor2', 
  'Full_Result_qbed_UMI_Anchor1', 'Full_Result_qbed_UMI_Anchor2'
)

ggplot(pyCC_Finalpercentage_table_coll, aes(x = factor(Name, graphorder), y = PercentageOverlapping, fill = Name)) +
  geom_col(color = "black", show.legend = FALSE) +
  ylim(0, 100) +
  labs(x = "Name of CC peak file and ChIA Dataset", 
       y = "% of CC peaks overlapping with ChIA-PET peak",
       title = "Collapsed qBED + Full Peaks: Single-cell pyCalling Cards Overlap") +
  theme_minimal() +
  geom_text(aes(label = round(PercentageOverlapping, 1)), vjust = -0.5, size = 4) +
  theme(axis.text.x = element_text(angle = 90))



# Combine anchors: AnyAnchor
FullPeaks_qbed_coll <- data.frame(
  Anchor1 = Peak_anchor1_qbed_coll$`Number of Overlaps with ChIA-PET Peaks`,
  Anchor2 = Peak_anchor2_qbed_coll$`Number of Overlaps with ChIA-PET Peaks`
)
FullPeaks_qbed_coll$AnyAnchor <- (FullPeaks_qbed_coll$Anchor1 > 0 | FullPeaks_qbed_coll$Anchor2 > 0)
perc_qbed_coll <- (sum(FullPeaks_qbed_coll$AnyAnchor)/nrow(FullPeaks_qbed_coll))*100

FullPeaks_qbed_UMI_coll <- data.frame(
  Anchor1 = UMI_anchor1_qbed_coll$`Number of Overlaps with ChIA-PET Peaks`,
  Anchor2 = UMI_anchor2_qbed_coll$`Number of Overlaps with ChIA-PET Peaks`
)
FullPeaks_qbed_UMI_coll$AnyAnchor <- (FullPeaks_qbed_UMI_coll$Anchor1 > 0 | FullPeaks_qbed_UMI_coll$Anchor2 > 0)
perc_qbed_UMI_coll <- (sum(FullPeaks_qbed_UMI_coll$AnyAnchor)/nrow(FullPeaks_qbed_UMI_coll))*100

pyCC_AnyAnchor_table_collapsed <- data.frame(
  Replicate = c("Full_result_Peaks", "Full_result_Peaks_UMI", "Full_result_qbed", "Full_result_qbed_UMI"),
  PercentageOverlapping = c(perc2, perc4, perc_qbed_coll, perc_qbed_UMI_coll)
)

graphorder= c("Full_result_Peaks", "Full_result_Peaks_UMI", "Full_result_qbed", "Full_result_qbed_UMI")

ggplot(pyCC_AnyAnchor_table_collapsed, aes(x = factor(Replicate, graphorder), y = PercentageOverlapping, fill = Replicate)) +
  geom_col(color = "black", show.legend = FALSE) +
  ylim(0, 100) +
  labs(x = "Replicate",
       y = "% of CC peaks overlapping with ChIA-PET (any anchor)",
       title = "Collapsed qBED + Full Peaks: Any Anchor Overlap") +
  theme_minimal() +
  geom_text(aes(label = round(PercentageOverlapping, 1)), vjust = -0.5, size = 4)