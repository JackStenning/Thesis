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
pyCCcol_names = c("Chr",	"Start",	"End",	"Center",	"Experiment Insertions",	"Background insertions",	"Reference Insertions",	"pvalue Reference",	"pvalue Background",	"Fraction Experiment",	"TPH Experiment",	"Fraction background",	"TPH background",	"TPH background subtracted",	"pvalue_adj Reference", "Number of Overlaps with ChIP-Seq Peaks")

##Read data
#Variables named <Window '-a'><Window'-b'>
#Encode Rep 2 was taken as it has more peaks than replicate 1 - even though there is less overlap with my data
ESR1_pyCC_2_GEO_q = read.table("results/1kbp_A_peak_data_ER_WT2.bed_B_Andy_ChIP_q0.05_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_2_Encode2_q = read.table("results/1kbp_A_peak_data_ER_WT2.bed_B_Encode_ER_rep2_ChIP_q0.05_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)


ESR1_pyCC_4_GEO_q = read.table("results/1kbp_A_peak_data_ER_WT4.bed_B_Andy_ChIP_q0.05_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)
ESR1_pyCC_4_Encode2_q = read.table("results/1kbp_A_peak_data_ER_WT4.bed_B_Encode_ER_rep2_ChIP_q0.05_peaksNP.bed", header = TRUE, sep = "\t", col.names = pyCCcol_names)

# Load your BED file (if not already loaded)
# Example: df <- read.table("input.bed", header = TRUE, sep = "\t")

# Filter rows where overlaps == 0
filtered_2_GEO_q <- ESR1_pyCC_2_GEO_q[ESR1_pyCC_2_GEO_q$`Number.of.Overlaps.with.ChIP.Seq.Peaks` == 0, ]

# Write to a new BED file
write.table(
  filtered_2_GEO_q,
  file = "Peakset2_GEO_no_chipseq_overlaps.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE# set to FALSE if you don’t want the header
)

# Filter rows where overlaps == 0
filtered_2_Encode_q <- ESR1_pyCC_2_Encode2_q[ESR1_pyCC_2_Encode2_q$`Number.of.Overlaps.with.ChIP.Seq.Peaks` == 0, ]

# Write to a new BED file
write.table(
  filtered_2_Encode_q,
  file = "Peakset2_Encode_no_chipseq_overlaps.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE  # set to FALSE if you don’t want the header
)

# Filter rows where overlaps == 0
filtered_4_GEO_q <- ESR1_pyCC_4_GEO_q[ESR1_pyCC_4_GEO_q$`Number.of.Overlaps.with.ChIP.Seq.Peaks` == 0, ]

# Write to a new BED file
write.table(
  filtered_4_GEO_q,
  file = "Peakset4_GEO_no_chipseq_overlaps.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE  # set to FALSE if you don’t want the header
)

# Filter rows where overlaps == 0
filtered_4_ENCODE_q <- ESR1_pyCC_4_Encode2_q[ESR1_pyCC_4_Encode2_q$`Number.of.Overlaps.with.ChIP.Seq.Peaks` == 0, ]

# Write to a new BED file
write.table(
  filtered_4_ENCODE_q,
  file = "Peakset4_Encode_no_chipseq_overlaps.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE  # set to FALSE if you don’t want the header
)
