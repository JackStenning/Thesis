############# LOOK HERE ######################
#add or delete the '#' symbol below on install packages
#install.packages("Rtools")
#install.packages("ggthemes")


############# SETUP ######################
# Load packages
library(ggplot2)
library(ggthemes)
library(dplyr)

# Set working directory
setwd("C:/Users/jps558/OneDrive - University of York/Desktop/R_WorkingDirectory/R_WFHDirectory/scCC") 

# Column names
pyCCcol_q_names = c("Chr", "Start", "End", "Experiment Insertions", "Strand", "barcode", "Number of Overlaps with ChIP-Seq Peaks")

############# READ DATA ##################
ESR1_pyCC_full_Highcc = read.table("result_trim_minimum/1kbp_A_JS_ES1_allpass_adapter_trimmed.tagged.bam_final.sorted.ccf_B_highCC.bed", 
                                   header = TRUE, sep = "\t", col.names = pyCCcol_q_names)
ESR1_pyCC_full_Lowcc_q = read.table("result_trim_minimum/1kbp_A_JS_ES1_allpass_adapter_trimmed.tagged.bam_final.sorted.ccf_B_lowCC.bed", 
                                    header = TRUE, sep = "\t", col.names = pyCCcol_q_names)
ESR1_pyCC_full_UMI_Highcc_q = read.table("result_trim_minimum/1kbp_A_JS_ES1_allpass_adapter_trimmed.tagged.bam_UMIFilt_final.ccf_B_highCC.bed", 
                                         header = TRUE, sep = "\t", col.names = pyCCcol_q_names)
ESR1_pyCC_full_UMI_Lowcc_q = read.table("result_trim_minimum/1kbp_A_JS_ES1_allpass_adapter_trimmed.tagged.bam_UMIFilt_final.ccf_B_lowCC.bed", 
                                        header = TRUE, sep = "\t", col.names = pyCCcol_q_names)

############# TOTAL COUNTS ##################
ESR1_pyCC_2_total_q = nrow(ESR1_pyCC_full_Highcc)
ESR1_pyCC_3_total_q = nrow(ESR1_pyCC_full_Lowcc_q)
ESR1_pyCC_4_total_q = nrow(ESR1_pyCC_full_UMI_Highcc_q)
ESR1_pyCC_5_total_q = nrow(ESR1_pyCC_full_UMI_Lowcc_q)

############# OVERLAP COUNTS ##################
ESR1_pyCC_full_Highcc_overlap = sum(ESR1_pyCC_full_Highcc$Number.of.Overlaps.with.ChIP.Seq.Peaks > 0, na.rm = TRUE)
ESR1_pyCC_full_Lowcc_q_overlap = sum(ESR1_pyCC_full_Lowcc_q$Number.of.Overlaps.with.ChIP.Seq.Peaks > 0, na.rm = TRUE)
ESR1_pyCC_full_UMI_Highcc_q_overlap = sum(ESR1_pyCC_full_UMI_Highcc_q$Number.of.Overlaps.with.ChIP.Seq.Peaks > 0, na.rm = TRUE)
ESR1_pyCC_full_UMI_Lowcc_q_overlap = sum(ESR1_pyCC_full_UMI_Lowcc_q$Number.of.Overlaps.with.ChIP.Seq.Peaks > 0, na.rm = TRUE)

############# PERCENTAGE OVERLAP ##################
Percentage_ESR1_pyCC_full_Highcc = (ESR1_pyCC_full_Highcc_overlap / ESR1_pyCC_2_total_q) * 100
Percentage_ESR1_pyCC_full_Lowcc_q = (ESR1_pyCC_full_Lowcc_q_overlap / ESR1_pyCC_3_total_q) * 100
Percentage_ESR1_pyCC_full_UMI_Highcc_q = (ESR1_pyCC_full_UMI_Highcc_q_overlap / ESR1_pyCC_4_total_q) * 100
Percentage_ESR1_pyCC_full_UMI_Lowcc_q = (ESR1_pyCC_full_UMI_Lowcc_q_overlap / ESR1_pyCC_5_total_q) * 100

############# CREATE FINAL TABLE ##################
pyCC_Finalpercentage_table <- data.frame(
  Name = c('Full_Result_Peaks_High', 'Full_UMI_Peaks_High', 
           'Full_Result_Peaks_Low', 'Full_UMI_Peaks_Low'),
  Total = c(ESR1_pyCC_2_total_q, ESR1_pyCC_4_total_q, 
            ESR1_pyCC_3_total_q, ESR1_pyCC_5_total_q),
  Overlap_withChIP = c(ESR1_pyCC_full_Highcc_overlap, ESR1_pyCC_full_UMI_Highcc_q_overlap,
                       ESR1_pyCC_full_Lowcc_q_overlap, ESR1_pyCC_full_UMI_Lowcc_q_overlap),
  PercentageOverlapping = c(Percentage_ESR1_pyCC_full_Highcc, Percentage_ESR1_pyCC_full_UMI_Highcc_q,
                            Percentage_ESR1_pyCC_full_Lowcc_q, Percentage_ESR1_pyCC_full_UMI_Lowcc_q)
)

############# PLOT FULL PEAKS ##################
graphorder = c('Full_Result_Peaks_High', 'Full_UMI_Peaks_High', 
               'Full_Result_Peaks_Low', 'Full_UMI_Peaks_Low')

pyCC_Figure <- ggplot(data = pyCC_Finalpercentage_table, 
                      aes(x = factor(Name, graphorder), 
                          y = PercentageOverlapping,
                          fill = Name)) +  
  geom_col(color = "black", show.legend = FALSE) +  
  ylim(0, 100) + 
  labs(x = "Name of CC peak file and single-cell Dataset", 
       y = "% of single-cell CC peaks overlapping with bulk CC peaks files",
       title = "Single-cell pyCalling Cards Peaks and Bulk CC Overlap") +
  theme_minimal() + 
  geom_text(aes(label = round(PercentageOverlapping, 3)), vjust = -0.5, size = 4) + 
  coord_cartesian(expand = FALSE, clip = "off") + 
  geom_hline(yintercept = 83, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 90))

pyCC_Figure


############# COLLAPSE QBED ##################
collapse_qbed <- function(df) {
  df %>%
    group_by(Chr, Start, End) %>%
    summarise(
      `Experiment Insertions` = sum(`Experiment.Insertions`),
      Strand = first(Strand),
      barcode = paste(unique(barcode), collapse = ";"),
      `Number of Overlaps with ChIP-Seq Peaks` = sum(`Number.of.Overlaps.with.ChIP.Seq.Peaks`)
    ) %>%
    ungroup()
}

Full_High_qbed_coll <- collapse_qbed(ESR1_pyCC_full_Highcc)
UMI_High_qbed_coll <- collapse_qbed(ESR1_pyCC_full_UMI_Highcc_q)
Full_Low_qbed_coll <- collapse_qbed(ESR1_pyCC_full_Lowcc_q)
UMI_Low_qbed_coll <- collapse_qbed(ESR1_pyCC_full_UMI_Lowcc_q)

# Calc percentages for collapsed
calc_percentage <- function(df) {
  total <- nrow(df)
  perc <- (sum(df$`Number of Overlaps with ChIP-Seq Peaks` != 0) / total) * 100
  overlap_count <- sum(df$`Number of Overlaps with ChIP-Seq Peaks` != 0)
  return(list(total = total, perc = perc, overlap = overlap_count))
}

perc_Full_High_qbed <- calc_percentage(Full_High_qbed_coll)$perc
perc_UMI_High_qbed <- calc_percentage(UMI_High_qbed_coll)$perc
perc_Full_Low_qbed <- calc_percentage(Full_Low_qbed_coll)$perc
perc_UMI_Low_qbed <- calc_percentage(UMI_Low_qbed_coll)$perc

overlap_Full_High_qbed <- calc_percentage(Full_High_qbed_coll)$overlap
overlap_UMI_High_qbed <- calc_percentage(UMI_High_qbed_coll)$overlap
overlap_Full_Low_qbed <- calc_percentage(Full_Low_qbed_coll)$overlap
overlap_UMI_Low_qbed <- calc_percentage(UMI_Low_qbed_coll)$overlap

total_Full_High_qbed <- calc_percentage(Full_High_qbed_coll)$total
total_UMI_High_qbed <- calc_percentage(UMI_High_qbed_coll)$total
total_Full_Low_qbed <- calc_percentage(Full_Low_qbed_coll)$total
total_UMI_Low_qbed <- calc_percentage(UMI_Low_qbed_coll)$total

############# CREATE TABLE FOR COLLAPSED QBED ##################
pyCC_Finalpercentage_table_coll <- data.frame(
  Name = c('Full_Result_qBed_High', 'Full_UMI_qBed_High', 'Full_Result_qBed_Low', 'Full_UMI_qBed_Low'),
  Total = c(total_Full_High_qbed, total_UMI_High_qbed, total_Full_Low_qbed, total_UMI_Low_qbed),
  Overlap_withChIP = c(overlap_Full_High_qbed, overlap_UMI_High_qbed, overlap_Full_Low_qbed, overlap_UMI_Low_qbed),
  PercentageOverlapping = c(perc_Full_High_qbed, perc_UMI_High_qbed, perc_Full_Low_qbed, perc_UMI_Low_qbed)
)

############# PLOT COLLAPSED QBED ##################
graphorder = c('Full_Result_qBed_High', 'Full_UMI_qBed_High', 'Full_Result_qBed_Low', 'Full_UMI_qBed_Low')
pyCC_Figure_coll <- ggplot(data = pyCC_Finalpercentage_table_coll, 
                           aes(x = factor(Name, graphorder), 
                               y = PercentageOverlapping,
                               fill = Name)) +  
  geom_col(color = "black", show.legend = FALSE) +  
  ylim(0, 100) + 
  labs(x = "Collapsed qBED Dataset", 
       y = "% of loci overlapping with bulk CC peaks",
       title = "Collapsed qBED: Single-cell pyCalling Cards and Bulk CC Overlap") +
  theme_minimal() + 
  geom_text(aes(label = round(PercentageOverlapping, 3)), vjust = -0.5, size = 4) + 
  coord_cartesian(expand = FALSE, clip = "off") + 
  geom_hline(yintercept = 83, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 90))

pyCC_Figure_coll
