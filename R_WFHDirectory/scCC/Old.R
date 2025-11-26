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
pyCCcol_names = c("Chr",	"Start",	"End",	"Experiment Insertions", "Strand", "barcode", "Number of Overlaps with ChIP-Seq Peaks")

##Read data
#Variables named <Window '-a'><Window'-b'>
#Encode Rep 2 was taken as it has more peaks than replicate 1 - even though there is less overlap with my data
ESR1_pyCC_full_GEO_q_qbed = read.table("result_full/1kbp_A_JS_ES1_allpass_adapter_trimmed.tagged.bam_final.sorted.ccf_B_Andy.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_full_Encode2_qbed = read.table("result_full/1kbp_A_JS_ES1_allpass_adapter_trimmed.tagged.bam_final.sorted.ccf_B_BRD4.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)


ESR1_pyCC_UMI_GEO_q_qbed = read.table("result_full/1kbp_A_JS_ES1_allpass_adapter_trimmed.tagged.bam_UMIFilt_final.ccf_B_Andy.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_UMI_Encode2_q_qbed = read.table("result_full/1kbp_A_JS_ES1_allpass_adapter_trimmed.tagged.bam_UMIFilt_final.ccf_B_BRD4.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)


#Counting the total reads

ESR1_pyCC_2_total_q=nrow(ESR1_pyCC_full_GEO_q_qbed)

ESR1_pyCC_4_total_q=nrow(ESR1_pyCC_UMI_GEO_q_qbed)


#Count the number of reads with no overlap
ESR1_pyCC_full_GEO_q_qbed_overlap_q = as.numeric(length(which(ESR1_pyCC_full_GEO_q_qbed$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_full_Encode_q_qbed_overlap_q = as.numeric(length(which(ESR1_pyCC_full_Encode2_qbed$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))



ESR1_pyCC_UMI_GEO_q_qbed_overlap_q = as.numeric(length(which(ESR1_pyCC_UMI_GEO_q_qbed$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))
ESR1_pyCC_UMI_Encode2_q_qbed_overlap_q = as.numeric(length(which(ESR1_pyCC_UMI_Encode2_q_qbed$Number.of.Overlaps.with.ChIP.Seq.Peaks==0)))


#Determine overlap percentage
Percentage_ESR1_pyCC_full_GEO_q_qbed_q =(1-(ESR1_pyCC_full_GEO_q_qbed_overlap_q/ESR1_pyCC_2_total_q))*100
Percentage_ESR1_pyCC_full_Encode2_qbed_q =(1-(ESR1_pyCC_full_Encode_q_qbed_overlap_q/ESR1_pyCC_2_total_q))*100



Percentage_ESR1_pyCC_UMI_GEO_q_qbed_q =(1-(ESR1_pyCC_UMI_GEO_q_qbed_overlap_q/ESR1_pyCC_4_total_q))*100
Percentage_ESR1_pyCC_UMI_Encode2_q_qbed_q = (1-(ESR1_pyCC_UMI_Encode2_q_qbed_overlap_q/ESR1_pyCC_4_total_q))*100

#Create table to plot
#Create table to plot
pyCC_Finalpercentage_table <- data.frame(
  Name = c('Full_Result_qBed_GEO', 'Full_UMI_qBed_GEO', 'Full_Result_qBed_Encode', 'Full_UMI_qBed_Encode'),
  Total = c(ESR1_pyCC_2_total_q, ESR1_pyCC_4_total_q, ESR1_pyCC_2_total_q, ESR1_pyCC_4_total_q),
  Overlap_withChIP = c(ESR1_pyCC_full_GEO_q_qbed_overlap_q, ESR1_pyCC_UMI_GEO_q_qbed_overlap_q, ESR1_pyCC_full_Encode_q_qbed_overlap_q, ESR1_pyCC_UMI_Encode2_q_qbed_overlap_q),
  PercentageOverlapping = c(Percentage_ESR1_pyCC_full_GEO_q_qbed_q, Percentage_ESR1_pyCC_UMI_GEO_q_qbed_q, Percentage_ESR1_pyCC_full_Encode2_qbed_q, Percentage_ESR1_pyCC_UMI_Encode2_q_qbed_q)
)

#Then we plot
graphorder= c('Full_Result_qBed_GEO', 'Full_UMI_qBed_GEO', 'Full_Result_qBed_Encode', 'Full_UMI_qBed_Encode')
pyCC_Figure <- ggplot(data = pyCC_Finalpercentage_table, 
                     aes(x= factor(Name, graphorder), 
                         y = PercentageOverlapping,
                         fill = Name, 
                         stat("identity"))) +  geom_col(color = "black",
                                                        show.legend = FALSE) +  ylim(0, 100) + labs(x = "Name of CC peak file and ChIP Dataset", 
                                                                                                    y = "% of CC qBed overlapping with ChIP-Seq peak",
                                                                                                    title = "Full result - single cell pyCalling Cards qBed and ChIP-Seq Overlap",
                                                                                                    subtitle = "Determining percentage of overlapping qBed from pyCalling Card ESR1 peak files with either GSE109820 or \nENCFF063JMY ChIP-Seq datasets within 1 kbp") +  theme_dark() + geom_text(aes(label = round(PercentageOverlapping, 1)), # Add percentage as text
                                                                                                                                                                                                                                                                             vjust = -0.5, # Adjust vertical position above the bars
                                                                                                                                                                                                                                                                             size = 4) + coord_cartesian(expand = FALSE, clip = "off") + geom_hline(yintercept = 83, linetype = "dashed") 


pyCC_Figure + theme(axis.text.x = element_text(angle = 90))
########################################################################
