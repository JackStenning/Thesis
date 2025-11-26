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
pyCCcol_names = c("Chr",	"Start",	"End",	"Experiment Insertions", "Reference Inserions",	"Expected Insertions",	"pvalue",	"Fraction Experiment",	"TPH Experiment",	"TPH background subtracted",	"pvalue_adj Reference", "Number of Overlaps with ChIP-Seq Peaks")
pyCCcol_q_names = c("Chr",	"Start",	"End",	"Experiment Insertions", "Strand", "barcode", "Number of Overlaps with ChIP-Seq Peaks")
##Read data
#Variables named <Window '-a'><Window'-b'>
#Encode Rep 2 was taken as it has more peaks than replicate 1 - even though there is less overlap with my data
ESR1_pyCC_full_GEO_q = read.table("Fullresult/1kbp_A_peak_data_ER_test4.bed_B_Andy.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_full_Encode2_q = read.table("Fullresult/1kbp_A_peak_data_ER_test4.bed_B_BRD4.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)


ESR1_pyCC_UMI_GEO_q = read.table("FullresultUMI/1kbp_A_peak_data_ER_test4.bed_B_Andy.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_UMI_Encode2_q = read.table("FullresultUMI/1kbp_A_peak_data_ER_test4.bed_B_BRD4.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)

ESR1_pyCC_full_GEO_q_qbed = read.table("result_full/1kbp_A_JS_ES1_allpass_adapter_trimmed.tagged.bam_final.sorted.ccf_B_Andy.bed", header = TRUE, sep = "\t", col.names = pyCCcol_q_names)
ESR1_pyCC_full_Encode2_qbed = read.table("result_full/1kbp_A_JS_ES1_allpass_adapter_trimmed.tagged.bam_final.sorted.ccf_B_BRD4.bed", header = TRUE, sep = "\t", col.names = pyCCcol_q_names)


ESR1_pyCC_UMI_GEO_q_qbed = read.table("result_full/1kbp_A_JS_ES1_allpass_adapter_trimmed.tagged.bam_UMIFilt_final.ccf_B_Andy.bed", header = TRUE, sep = "\t", col.names = pyCCcol_q_names)
ESR1_pyCC_UMI_Encode2_q_qbed = read.table("result_full/1kbp_A_JS_ES1_allpass_adapter_trimmed.tagged.bam_UMIFilt_final.ccf_B_BRD4.bed", header = TRUE, sep = "\t", col.names = pyCCcol_q_names)

#Counting the total reads

ESR1_pyCC_2_total=nrow(ESR1_pyCC_full_GEO_q)

ESR1_pyCC_4_total=nrow(ESR1_pyCC_UMI_GEO_q)

ESR1_pyCC_2_total_q=nrow(ESR1_pyCC_full_GEO_q_qbed)

ESR1_pyCC_4_total_q=nrow(ESR1_pyCC_UMI_GEO_q_qbed)


#Count the number of reads with no overlap
ESR1_pyCC_full_GEO_q_overlap = as.numeric(length(which(ESR1_pyCC_full_GEO_q$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_full_Encode2_q_overlap = as.numeric(length(which(ESR1_pyCC_full_Encode2_q$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))



ESR1_pyCC_UMI_GEO_q_overlap = as.numeric(length(which(ESR1_pyCC_UMI_GEO_q$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_UMI_Encode2_q_overlap = as.numeric(length(which(ESR1_pyCC_UMI_Encode2_q$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))

ESR1_pyCC_full_GEO_q_qbed_overlap_q = as.numeric(length(which(ESR1_pyCC_full_GEO_q_qbed$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_full_Encode_q_qbed_overlap_q = as.numeric(length(which(ESR1_pyCC_full_Encode2_qbed$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))



ESR1_pyCC_UMI_GEO_q_qbed_overlap_q = as.numeric(length(which(ESR1_pyCC_UMI_GEO_q_qbed$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_UMI_Encode2_q_qbed_overlap_q = as.numeric(length(which(ESR1_pyCC_UMI_Encode2_q_qbed$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))

#Determine overlap percentage
Percentage_ESR1_pyCC_full_GEO_q =(1-(ESR1_pyCC_full_GEO_q_overlap/ESR1_pyCC_2_total))*100
Percentage_ESR1_pyCC_full_Encode2_q =(1-(ESR1_pyCC_full_Encode2_q_overlap/ESR1_pyCC_2_total))*100



Percentage_ESR1_pyCC_UMI_GEO_q =(1-(ESR1_pyCC_UMI_GEO_q_overlap/ESR1_pyCC_4_total))*100
Percentage_ESR1_pyCC_UMI_Encode2_q = (1-(ESR1_pyCC_UMI_Encode2_q_overlap/ESR1_pyCC_4_total))*100

Percentage_ESR1_pyCC_full_GEO_q_qbed_q =(1-(ESR1_pyCC_full_GEO_q_qbed_overlap_q/ESR1_pyCC_2_total_q))*100
Percentage_ESR1_pyCC_full_Encode2_qbed_q =(1-(ESR1_pyCC_full_Encode_q_qbed_overlap_q/ESR1_pyCC_2_total_q))*100



Percentage_ESR1_pyCC_UMI_GEO_q_qbed_q =(1-(ESR1_pyCC_UMI_GEO_q_qbed_overlap_q/ESR1_pyCC_4_total_q))*100
Percentage_ESR1_pyCC_UMI_Encode2_q_qbed_q = (1-(ESR1_pyCC_UMI_Encode2_q_qbed_overlap_q/ESR1_pyCC_4_total_q))*100

#Create table to plot
#Create table to plot
pyCC_Finalpercentage_table <- data.frame(
  Name = c(
    'Full_Result_Peaks_GEO', 'Full_UMI_Peaks_GEO', 
    'Full_Result_qBed_GEO', 'Full_UMI_qBed_GEO',
    'Full_Result_Peaks_Encode', 'Full_UMI_Peaks_Encode', 
    'Full_Result_qBed_Encode', 'Full_UMI_qBed_Encode'
  ),
  Total = c(
    ESR1_pyCC_2_total, ESR1_pyCC_4_total, 
    ESR1_pyCC_2_total_q, ESR1_pyCC_4_total_q,
    ESR1_pyCC_2_total, ESR1_pyCC_4_total, 
    ESR1_pyCC_2_total_q, ESR1_pyCC_4_total_q
  ),
  Overlap_withChIP = c(
    ESR1_pyCC_full_GEO_q_overlap, ESR1_pyCC_UMI_GEO_q_overlap, 
    ESR1_pyCC_full_GEO_q_qbed_overlap_q, ESR1_pyCC_UMI_GEO_q_qbed_overlap_q,
    ESR1_pyCC_full_Encode2_q_overlap, ESR1_pyCC_UMI_Encode2_q_overlap, 
    ESR1_pyCC_full_Encode_q_qbed_overlap_q, ESR1_pyCC_UMI_Encode2_q_qbed_overlap_q
  ),
  PercentageOverlapping = c(
    Percentage_ESR1_pyCC_full_GEO_q, Percentage_ESR1_pyCC_UMI_GEO_q, 
    Percentage_ESR1_pyCC_full_GEO_q_qbed_q, Percentage_ESR1_pyCC_UMI_GEO_q_qbed_q,
    Percentage_ESR1_pyCC_full_Encode2_q, Percentage_ESR1_pyCC_UMI_Encode2_q, 
    Percentage_ESR1_pyCC_full_Encode2_qbed_q, Percentage_ESR1_pyCC_UMI_Encode2_q_qbed_q
  )
)

#Then we plot
# Landscape A4 TIFF export
tiff("Full_pyCC_vs_ChIP_A4_Landscape.tiff",
     width = 3508, height = 2480, res = 300)

graphorder = c(
  'Full_Result_Peaks_GEO', 'Full_UMI_Peaks_GEO', 
  'Full_Result_qBed_GEO', 'Full_UMI_qBed_GEO',
  'Full_Result_Peaks_Encode', 'Full_UMI_Peaks_Encode', 
  'Full_Result_qBed_Encode', 'Full_UMI_qBed_Encode'
)

pyCC_Figure <- ggplot(data = pyCC_Finalpercentage_table, 
                      aes(x = factor(Name, graphorder), 
                          y = PercentageOverlapping,
                          fill = Name)) +
  geom_col(color = "black", show.legend = FALSE) +
  ylim(0, 100) +
  labs(
    x = "Name of CC peak file and ChIP Dataset",
    y = "% of CC peaks overlapping with ChIP-Seq peak",
    title = "Full result - single cell pyCalling Cards Peaks and ChIP-Seq Overlap",
    subtitle = "Determining percentage of overlapping peaks from pyCalling Card single-cell ESR1 peak files with either GSE109820 or \nENCFF063JMY ChIP-Seq datasets within 1 kbp"
  ) +
  
  # Increase bar-label text +2 pt (4 → 6)
  geom_text(aes(label = round(PercentageOverlapping, 1)),
            vjust = -0.5,
            size = 6) +
  
  coord_cartesian(expand = FALSE, clip = "off") +
  geom_hline(yintercept = 83, linetype = "dashed") +
  
  # +2 pt to all theme text except title/subtitle
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, size = 14),
    axis.text.y = element_text(size = 17),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 18),
    
    # Leave title/subtitle unchanged
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

print(pyCC_Figure)

dev.off()

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

Full_GEO_qbed_coll <- collapse_qbed(ESR1_pyCC_full_GEO_q_qbed)
UMI_GEO_qbed_coll  <- collapse_qbed(ESR1_pyCC_UMI_GEO_q_qbed)
Full_Encode_qbed_coll <- collapse_qbed(ESR1_pyCC_full_Encode2_qbed)
UMI_Encode_qbed_coll  <- collapse_qbed(ESR1_pyCC_UMI_Encode2_q_qbed)

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
    'Full_Result_Peaks_GEO', 'Full_UMI_Peaks_GEO', 
    'Full_Result_qBed_GEO', 'Full_UMI_qBed_GEO',
    'Full_Result_Peaks_Encode', 'Full_UMI_Peaks_Encode', 
    'Full_Result_qBed_Encode', 'Full_UMI_qBed_Encode'
  ),
  Total = c(
    ESR1_pyCC_2_total, ESR1_pyCC_4_total,
    total_Full_GEO_qbed, total_UMI_GEO_qbed,
    ESR1_pyCC_2_total, ESR1_pyCC_4_total,
    total_Full_Encode_qbed, total_UMI_Encode_qbed
  ),
  Overlap_withChIP = c(
    ESR1_pyCC_full_GEO_q_overlap, ESR1_pyCC_UMI_GEO_q_overlap,
    sum(Full_GEO_qbed_coll$`Number of Overlaps with ChIP-Seq Peaks` > 0),
    sum(UMI_GEO_qbed_coll$`Number of Overlaps with ChIP-Seq Peaks` > 0),
    ESR1_pyCC_full_Encode2_q_overlap, ESR1_pyCC_UMI_Encode2_q_overlap,
    sum(Full_Encode_qbed_coll$`Number of Overlaps with ChIP-Seq Peaks` > 0),
    sum(UMI_Encode_qbed_coll$`Number of Overlaps with ChIP-Seq Peaks` > 0)
  ),
  PercentageOverlapping = c(
    Percentage_ESR1_pyCC_full_GEO_q, Percentage_ESR1_pyCC_UMI_GEO_q,
    perc_Full_GEO_qbed, perc_UMI_GEO_qbed,
    Percentage_ESR1_pyCC_full_Encode2_q, Percentage_ESR1_pyCC_UMI_Encode2_q,
    perc_Full_Encode_qbed, perc_UMI_Encode_qbed
  )
)

#Plot
# Export to A4 landscape TIFF
# Export to A4 landscape TIFF
tiff("Collapsed_qBED_Full_Peaks_vs_ChIP_A4_Landscape.tiff",
     width = 3508, height = 2480, res = 300)

graphorder = c(
  'Full_Result_Peaks_GEO', 'Full_UMI_Peaks_GEO', 
  'Full_Result_qBed_GEO', 'Full_UMI_qBed_GEO',
  'Full_Result_Peaks_Encode', 'Full_UMI_Peaks_Encode', 
  'Full_Result_qBed_Encode', 'Full_UMI_qBed_Encode'
)

pyCC_Figure <- ggplot(data = pyCC_Finalpercentage_table_coll, 
                      aes(x = factor(Name, graphorder), 
                          y = PercentageOverlapping,
                          fill = Name)) +  
  geom_col(color = "black", show.legend = FALSE) +  
  ylim(0, 100) + 
  labs(
    x = "Name of CC peak file and ChIP Dataset", 
    y = "% of CC peaks overlapping with ChIP-Seq peak",
    title = "Collapsed qBED + Full Peaks: Single-cell pyCalling Cards and ChIP-Seq Overlap",
    subtitle = "Percentage of overlapping peaks from single-cell ESR1 pyCalling Cards with GEO or Encode ChIP-Seq"
  ) +  
  theme_minimal(base_size = 14) +   # base size increased by +2
  geom_text(aes(label = round(PercentageOverlapping, 1)),
            vjust = -0.5,
            size = 6) +     # label size increased (4 → 6)
  coord_cartesian(expand = FALSE, clip = "off") +
  theme(
    axis.text.x = element_text(angle = 90, size = 14),
    axis.text.y = element_text(size = 17),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 18),
    
    # Titles stay unchanged
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

print(pyCC_Figure)

dev.off()
