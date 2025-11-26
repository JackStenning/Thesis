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
pyCCcol_names = c("Chr",	"Start",	"End",	"Experiment Insertions", "Strand", "barcode", "Number of Overlaps with ChIP-Seq Peaks")

##Read data
#Variables named <Window '-a'><Window'-b'>
#Encode Rep 2 was taken as it has more peaks than replicate 1 - even though there is less overlap with my data
ESR1_pyCC_full_GEO_q = read.table("result_trim_minimum/1kbp_A_JS_ES1_allpass_adapter_trimmed.tagged.bam_final.sorted.ccf_B_Andy.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_full_Encode2_q = read.table("result_trim_minimum/1kbp_A_JS_ES1_allpass_adapter_trimmed.tagged.bam_final.sorted.ccf_B_BRD4.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)


ESR1_pyCC_UMI_GEO_q = read.table("result_trim_minimum/1kbp_A_JS_ES1_allpass_adapter_trimmed.tagged.bam_UMIFilt_final.ccf_B_Andy.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_UMI_Encode2_q = read.table("result_trim_minimum/1kbp_A_JS_ES1_allpass_adapter_trimmed.tagged.bam_UMIFilt_final.ccf_B_BRD4.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)


#Counting the total reads

ESR1_pyCC_2_total=nrow(ESR1_pyCC_full_GEO_q)

ESR1_pyCC_4_total=nrow(ESR1_pyCC_UMI_GEO_q)


#Count the number of reads with no overlap
ESR1_pyCC_full_GEO_q_overlap = as.numeric(length(which(ESR1_pyCC_full_GEO_q$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_full_Encode2_q_overlap = as.numeric(length(which(ESR1_pyCC_full_Encode2_q$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))



ESR1_pyCC_UMI_GEO_q_overlap = as.numeric(length(which(ESR1_pyCC_UMI_GEO_q$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_UMI_Encode2_q_overlap = as.numeric(length(which(ESR1_pyCC_UMI_Encode2_q$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))


#Determine overlap percentage
Percentage_ESR1_pyCC_full_GEO_q =(1-(ESR1_pyCC_full_GEO_q_overlap/ESR1_pyCC_2_total))*100
Percentage_ESR1_pyCC_full_Encode2_q =(1-(ESR1_pyCC_full_Encode2_q_overlap/ESR1_pyCC_2_total))*100



Percentage_ESR1_pyCC_UMI_GEO_q =(1-(ESR1_pyCC_UMI_GEO_q_overlap/ESR1_pyCC_4_total))*100
Percentage_ESR1_pyCC_UMI_Encode2_q = (1-(ESR1_pyCC_UMI_Encode2_q_overlap/ESR1_pyCC_4_total))*100

#Create table to plot
#Create table to plot
pyCC_Finalpercentage_table = data.frame(Name=c('Pre-Trimmed_qbed_GEO', 'Pre-Trimmed_qbed_UMI_GEO', 'Pre-Trimmed_qbed_Encode', 'Pre-Trimmed_qbed_UMI_Encode'),
                                       Total=rep(c(ESR1_pyCC_2_total,  ESR1_pyCC_4_total, ESR1_pyCC_2_total, ESR1_pyCC_4_total)),
                                       Overlap_withChIP= c(ESR1_pyCC_full_GEO_q_overlap, ESR1_pyCC_UMI_GEO_q_overlap, ESR1_pyCC_full_Encode2_q_overlap, ESR1_pyCC_UMI_Encode2_q_overlap),
                                       PercentageOverlapping= c(Percentage_ESR1_pyCC_full_GEO_q, Percentage_ESR1_pyCC_UMI_GEO_q, Percentage_ESR1_pyCC_full_Encode2_q, Percentage_ESR1_pyCC_UMI_Encode2_q))

#Then we plot
graphorder = c('Pre-Trimmed_qbed_GEO', 'Pre-Trimmed_qbed_UMI_GEO',
               'Pre-Trimmed_qbed_Encode', 'Pre-Trimmed_qbed_UMI_Encode')

graphorder = c('Pre-Trimmed_qbed_GEO', 'Pre-Trimmed_qbed_UMI_GEO',
               'Pre-Trimmed_qbed_Encode', 'Pre-Trimmed_qbed_UMI_Encode')

tiff("PyCCtrimmed_qBed_vs_ERChIP.tiff",
     width = 3508, height = 2480, res = 300)
pyCC_Figure <- ggplot(data = pyCC_Finalpercentage_table, 
                      aes(x = factor(Name, graphorder),
                          y = PercentageOverlapping,
                          fill = Name)) +
  
  geom_col(color = "black", show.legend = FALSE) +
  ylim(0, 100) +
  labs(x = "Name of CC peak file and ChIP Dataset", 
       y = "% of CC qBed overlapping with ChIP-Seq peak",
       title = "Pre-trimmed result - single cell pyCalling Cards qBed vs ChIP-Seq",
       subtitle = "Overlap within 1kbp") +
  
  # <<< +2 font size here
  theme_minimal(base_size = 14) +
  
  geom_text(aes(label = round(PercentageOverlapping, 1)),
            vjust = -0.5,
            size = 6) +  # was 4 → now +2
  
  coord_cartesian(expand = FALSE, clip = "off") +
  
  # Axis text also larger
  theme(axis.text.x = element_text(angle = 90, size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 14))

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

Full_GEO_qbed_coll <- collapse_qbed(ESR1_pyCC_full_GEO_q)
UMI_GEO_qbed_coll  <- collapse_qbed(ESR1_pyCC_UMI_GEO_q)
Full_Encode_qbed_coll <- collapse_qbed(ESR1_pyCC_full_Encode2_q)
UMI_Encode_qbed_coll  <- collapse_qbed(ESR1_pyCC_UMI_Encode2_q)

# Calculate overlap percentages for collapsed qBED
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

# Build table for plotting collapsed qBED
pyCC_Finalpercentage_table_coll <- data.frame(
  Name = c(
    'Collapsed_qBed_GEO', 'Collapsed_qBed_UMI_GEO', 
    'Collapsed_qBed_Encode', 'Collapsed_qBed_UMI_Encode'
  ),
  Total = c(
    total_Full_GEO_qbed, total_UMI_GEO_qbed,
    total_Full_Encode_qbed, total_UMI_Encode_qbed
  ),
  Overlap_withChIP = c(
    sum(Full_GEO_qbed_coll$`Number of Overlaps with ChIP-Seq Peaks` > 0),
    sum(UMI_GEO_qbed_coll$`Number of Overlaps with ChIP-Seq Peaks` > 0),
    sum(Full_Encode_qbed_coll$`Number of Overlaps with ChIP-Seq Peaks` > 0),
    sum(UMI_Encode_qbed_coll$`Number of Overlaps with ChIP-Seq Peaks` > 0)
  ),
  PercentageOverlapping = c(
    perc_Full_GEO_qbed, perc_UMI_GEO_qbed,
    perc_Full_Encode_qbed, perc_UMI_Encode_qbed
  )
)

# Plot collapsed qBED overlap
graphorder = c(
  'Collapsed_qBed_GEO', 'Collapsed_qBed_UMI_GEO', 
  'Collapsed_qBed_Encode', 'Collapsed_qBed_UMI_Encode'
)

pyCC_Figure_coll <- ggplot(data = pyCC_Finalpercentage_table_coll, 
                           aes(x = factor(Name, graphorder), 
                               y = PercentageOverlapping,
                               fill = Name)) +  
  geom_col(color = "black", show.legend = FALSE) +  
  ylim(0, 100) + 
  labs(x = "Name of CC peak file and ChIP Dataset", 
       y = "% of CC qBed overlapping with ChIP-Seq peak",
       title = "Collapsed pre-trimmed qBED: Single-cell pyCalling Cards and ChIP-Seq Overlap",
       subtitle = "Percentage of overlapping peaks from pre-trimmed single-cell ESR1 pyCalling Cards qBED with GEO or Encode ChIP-Seq") +  
  theme_minimal() + 
  geom_text(aes(label = round(PercentageOverlapping, 1)), vjust = -0.5, size = 4) + 
  coord_cartesian(expand = FALSE, clip = "off") + 
  geom_hline(yintercept = 83, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 90))

pyCC_Figure_coll
dev.off()
