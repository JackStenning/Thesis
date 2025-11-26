############# LOOK HERE ######################
#add or delete the '#' symbol below on install packages
#install.packages("Rtools")
#install.packages("ggthemes")

#DELTE AFTER WRITING

#Load software
library(ggplot2)
library(ggthemes)

#Set working directory
setwd("C:/Users/jps558/OneDrive - University of York/Desktop/R_WorkingDirectory/R_WFHDirectory/scCC") 


#Set column names
pyCCcol_names = c("Chr",	"Start",	"End",	"Experiment Insertions", "Reference Inserions",	"Expected Insertions",	"pvalue",	"Fraction Experiment",	"TPH Experiment",	"TPH background subtracted",	"pvalue_adj Reference", "Number of Overlaps with ChIP-Seq Peaks")
pyCCcol_q_names = c("Chr",	"Start",	"End",	"Experiment Insertions", "Strand", "barcode", "Number of Overlaps with ChIP-Seq Peaks")
##Read data
#Variables named <Window '-a'><Window'-b'>
#Encode Rep 2 was taken as it has more peaks than replicate 1 - even though there is less overlap with my data
ESR1_pyCC_full_Highcc = read.table("Fullresult/1kbp_A_peak_data_ER_test4.bed_B_highCC.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_full_Lowcc_q = read.table("Fullresult/1kbp_A_peak_data_ER_test4.bed_B_lowCC.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)


ESR1_pyCC_full_UMI_Highcc_q = read.table("FullresultUMI/1kbp_A_peak_data_ER_test4.bed_B_highCC.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_full_UMI_Lowcc_q = read.table("FullresultUMI/1kbp_A_peak_data_ER_test4.bed_B_lowCC.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)

ESR1_pyCC_full_Highcc_qbed = read.table("result_full/1kbp_A_JS_ES1_allpass_adapter_trimmed.tagged.bam_final.sorted.ccf_B_highCC.bed", header = TRUE, sep = "\t", col.names = pyCCcol_q_names)
ESR1_pyCC_full_Lowcc_qbed = read.table("result_full/1kbp_A_JS_ES1_allpass_adapter_trimmed.tagged.bam_final.sorted.ccf_B_lowCC.bed", header = TRUE, sep = "\t", col.names = pyCCcol_q_names)


ESR1_pyCC_full_UMI_Highcc_q_qbed = read.table("result_full/1kbp_A_JS_ES1_allpass_adapter_trimmed.tagged.bam_UMIFilt_final.ccf_B_highCC.bed", header = TRUE, sep = "\t", col.names = pyCCcol_q_names)
ESR1_pyCC_full_UMI_Lowcc_q_qbed = read.table("result_full/1kbp_A_JS_ES1_allpass_adapter_trimmed.tagged.bam_UMIFilt_final.ccf_B_lowCC.bed", header = TRUE, sep = "\t", col.names = pyCCcol_q_names)

#Counting the total reads

ESR1_pyCC_2_total=nrow(ESR1_pyCC_full_Highcc)

ESR1_pyCC_4_total=nrow(ESR1_pyCC_full_UMI_Highcc_q)

ESR1_pyCC_2_total_q=nrow(ESR1_pyCC_full_Highcc_qbed)

ESR1_pyCC_4_total_q=nrow(ESR1_pyCC_full_UMI_Highcc_q_qbed)


#Count the number of reads with no overlap
ESR1_pyCC_full_Highcc_overlap = as.numeric(length(which(ESR1_pyCC_full_Highcc$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_full_Lowcc_q_overlap = as.numeric(length(which(ESR1_pyCC_full_Lowcc_q$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))



ESR1_pyCC_full_UMI_Highcc_q_overlap = as.numeric(length(which(ESR1_pyCC_full_UMI_Highcc_q$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_full_UMI_Lowcc_q_overlap = as.numeric(length(which(ESR1_pyCC_full_UMI_Lowcc_q$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))

ESR1_pyCC_full_Highcc_qbed_overlap_q = as.numeric(length(which(ESR1_pyCC_full_Highcc_qbed$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_full_Encode_q_qbed_overlap_q = as.numeric(length(which(ESR1_pyCC_full_Lowcc_qbed$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))



ESR1_pyCC_full_UMI_Highcc_q_qbed_overlap_q = as.numeric(length(which(ESR1_pyCC_full_UMI_Highcc_q_qbed$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_full_UMI_Lowcc_q_qbed_overlap_q = as.numeric(length(which(ESR1_pyCC_full_UMI_Lowcc_q_qbed$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))

#Determine overlap percentage
Percentage_ESR1_pyCC_full_Highcc =(1-(ESR1_pyCC_full_Highcc_overlap/ESR1_pyCC_2_total))*100
Percentage_ESR1_pyCC_full_Lowcc_q =(1-(ESR1_pyCC_full_Lowcc_q_overlap/ESR1_pyCC_2_total))*100



Percentage_ESR1_pyCC_full_UMI_Highcc_q =(1-(ESR1_pyCC_full_UMI_Highcc_q_overlap/ESR1_pyCC_4_total))*100
Percentage_ESR1_pyCC_full_UMI_Lowcc_q = (1-(ESR1_pyCC_full_UMI_Lowcc_q_overlap/ESR1_pyCC_4_total))*100

Percentage_ESR1_pyCC_full_Highcc_qbed_q =(1-(ESR1_pyCC_full_Highcc_qbed_overlap_q/ESR1_pyCC_2_total_q))*100
Percentage_ESR1_pyCC_full_Lowcc_qbed_q =(1-(ESR1_pyCC_full_Encode_q_qbed_overlap_q/ESR1_pyCC_2_total_q))*100



Percentage_ESR1_pyCC_full_UMI_Highcc_q_qbed_q =(1-(ESR1_pyCC_full_UMI_Highcc_q_qbed_overlap_q/ESR1_pyCC_4_total_q))*100
Percentage_ESR1_pyCC_full_UMI_Lowcc_q_qbed_q = (1-(ESR1_pyCC_full_UMI_Lowcc_q_qbed_overlap_q/ESR1_pyCC_4_total_q))*100

#Create table to plot
#Create table to plot
pyCC_Finalpercentage_table <- data.frame(
  Name = c(
    'Full_Result_Peaks_High', 'Full_UMI_Peaks_High', 
    'Full_Result_qBed_High', 'Full_UMI_qBed_High',
    'Full_Result_Peaks_Low', 'Full_UMI_Peaks_Low', 
    'Full_Result_qBed_Low', 'Full_UMI_qBed_Low'
  ),
  Total = c(
    ESR1_pyCC_2_total, ESR1_pyCC_4_total, 
    ESR1_pyCC_2_total_q, ESR1_pyCC_4_total_q,
    ESR1_pyCC_2_total, ESR1_pyCC_4_total, 
    ESR1_pyCC_2_total_q, ESR1_pyCC_4_total_q
  ),
  Overlap_withChIP = c(
    ESR1_pyCC_full_Highcc_overlap, ESR1_pyCC_full_UMI_Highcc_q_overlap, 
    ESR1_pyCC_full_Highcc_qbed_overlap_q, ESR1_pyCC_full_UMI_Highcc_q_qbed_overlap_q,
    ESR1_pyCC_full_Lowcc_q_overlap, ESR1_pyCC_full_UMI_Lowcc_q_overlap, 
    ESR1_pyCC_full_Encode_q_qbed_overlap_q, ESR1_pyCC_full_UMI_Lowcc_q_qbed_overlap_q
  ),
  PercentageOverlapping = c(
    Percentage_ESR1_pyCC_full_Highcc, Percentage_ESR1_pyCC_full_UMI_Highcc_q, 
    Percentage_ESR1_pyCC_full_Highcc_qbed_q, Percentage_ESR1_pyCC_full_UMI_Highcc_q_qbed_q,
    Percentage_ESR1_pyCC_full_Lowcc_q, Percentage_ESR1_pyCC_full_UMI_Lowcc_q, 
    Percentage_ESR1_pyCC_full_Lowcc_qbed_q, Percentage_ESR1_pyCC_full_UMI_Lowcc_q_qbed_q
  )
)

#Then we plot
graphorder= c(
  'Full_Result_Peaks_High', 'Full_UMI_Peaks_High', 
  'Full_Result_qBed_High', 'Full_UMI_qBed_High',
  'Full_Result_Peaks_Low', 'Full_UMI_Peaks_Low', 
  'Full_Result_qBed_Low', 'Full_UMI_qBed_Low'
)
pyCC_Figure <- ggplot(data = pyCC_Finalpercentage_table, 
                      aes(x= factor(Name, graphorder), 
                          y = PercentageOverlapping,
                          fill = Name, 
                          stat("identity"))) +  geom_col(color = "black",
                                                         show.legend = FALSE) +  ylim(0, 100) + labs(x = "Name of CC peak file and single-cell Dataset", 
                                                                                                     y = "% of single-cell CC peaks and qbed overlapping with bulk CC peaks files",
                                                                                                     title = "Full result - single cell pyCalling Cards Peaks and Bulk CC Overlap",
                                                                                                     subtitle = "Determining percentage of overlapping peaks from pyCalling Card single-cell ESR1 peak files with bulk CC datasets within 1 kbp") +  theme_minimal() + geom_text(aes(label = round(PercentageOverlapping, 1)), # Add percentage as text
                                                                                                                                                                                                                                                                                                                 vjust = -0.5, # Adjust vertical position above the bars
                                                                                                                                                                                                                                                                                                                 size = 4) + coord_cartesian(expand = FALSE, clip = "off") + geom_hline(yintercept = 83, linetype = "dashed") 


pyCC_Figure + theme(axis.text.x = element_text(angle = 90))
########################################################################


#Without multiple qbed entries
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

Full_GEO_qbed_coll <- collapse_qbed(ESR1_pyCC_full_Highcc_qbed)
UMI_GEO_qbed_coll  <- collapse_qbed(ESR1_pyCC_full_UMI_Highcc_q_qbed)
Full_Encode_qbed_coll <- collapse_qbed(ESR1_pyCC_full_Lowcc_qbed)
UMI_Encode_qbed_coll  <- collapse_qbed(ESR1_pyCC_full_UMI_Lowcc_q_qbed)

#Calc percentages
calc_percentage <- function(df) {
  total <- nrow(df)
  perc <- (1 - sum(df$`Number of Overlaps with ChIP-Seq Peaks` == 0)/total) * 100
  return(list(total = total, perc = perc))
}

perc_Full_GEO_qbed <- calc_percentage(Full_GEO_qbed_coll)$perc
perc_UMI_GEO_qbed  <- calc_percentage(UMI_GEO_qbed_coll)$perc
perc_Full_Encode_qbed <- calc_percentage(Full_Encode_qbed_coll)$perc
perc_UMI_Encode_qbed  <- calc_percentage(UMI_Encode_qbed_coll)$perc

total_Full_GEO_qbed <- calc_percentage(Full_GEO_qbed_coll)$total
total_UMI_GEO_qbed  <- calc_percentage(UMI_GEO_qbed_coll)$total
total_Full_Encode_qbed <- calc_percentage(Full_Encode_qbed_coll)$total
total_UMI_Encode_qbed  <- calc_percentage(UMI_Encode_qbed_coll)$total

#Build table
pyCC_Finalpercentage_table_coll <- data.frame(
  Name = c(
    'Full_Result_Peaks_High', 'Full_UMI_Peaks_High', 
    'Full_Result_qBed_High', 'Full_UMI_qBed_High',
    'Full_Result_Peaks_Low', 'Full_UMI_Peaks_Low', 
    'Full_Result_qBed_Low', 'Full_UMI_qBed_Low'
  ),
  Total = c(
    ESR1_pyCC_2_total, ESR1_pyCC_4_total,
    total_Full_GEO_qbed, total_UMI_GEO_qbed,
    ESR1_pyCC_2_total, ESR1_pyCC_4_total,
    total_Full_Encode_qbed, total_UMI_Encode_qbed
  ),
  Overlap_withChIP = c(
    ESR1_pyCC_full_Highcc_overlap, ESR1_pyCC_full_UMI_Highcc_q_overlap,
    sum(Full_GEO_qbed_coll$`Number of Overlaps with ChIP-Seq Peaks` > 0),
    sum(UMI_GEO_qbed_coll$`Number of Overlaps with ChIP-Seq Peaks` > 0),
    ESR1_pyCC_full_Lowcc_q_overlap, ESR1_pyCC_full_UMI_Lowcc_q_overlap,
    sum(Full_Encode_qbed_coll$`Number of Overlaps with ChIP-Seq Peaks` > 0),
    sum(UMI_Encode_qbed_coll$`Number of Overlaps with ChIP-Seq Peaks` > 0)
  ),
  PercentageOverlapping = c(
    Percentage_ESR1_pyCC_full_Highcc, Percentage_ESR1_pyCC_full_UMI_Highcc_q,
    perc_Full_GEO_qbed, perc_UMI_GEO_qbed,
    Percentage_ESR1_pyCC_full_Lowcc_q, Percentage_ESR1_pyCC_full_UMI_Lowcc_q,
    perc_Full_Encode_qbed, perc_UMI_Encode_qbed
  )
)

#Plot
graphorder = c(
  'Full_Result_Peaks_High', 'Full_UMI_Peaks_High', 
  'Full_Result_qBed_High', 'Full_UMI_qBed_High',
  'Full_Result_Peaks_Low', 'Full_UMI_Peaks_Low', 
  'Full_Result_qBed_Low', 'Full_UMI_qBed_Low'
)

pyCC_Figure <- ggplot(data = pyCC_Finalpercentage_table_coll, 
                      aes(x = factor(Name, graphorder), 
                          y = PercentageOverlapping,
                          fill = Name)) +  
  geom_col(color = "black", show.legend = FALSE) +  
  labs(x = "Name of CC peak file and single-cell Dataset", 
       y = "% of single-cell CC insertions overlapping with bulk CC peaks",
       title = "Collapsed qBED + Full Peaks: Single-cell pyCalling Cards and Bulk CC Overlap",
       subtitle = "Determining percentage of overlapping peaks from pyCalling Card single-cell ESR1 peak files with bulk CC datasets within 1 kbp") +
  theme_minimal() + 
  geom_text(aes(label = round(PercentageOverlapping, 1)), vjust = -0.5, size = 4) + 
  coord_cartesian(expand = FALSE, clip = "off") + 
  geom_hline(yintercept = 83, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 90))

pyCC_Figure
