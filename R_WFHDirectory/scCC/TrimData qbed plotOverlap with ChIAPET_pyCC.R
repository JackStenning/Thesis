############# LOOK HERE ######################
#add or delete the '#' symbol below on install packages
#install.packages("Rtools")
#install.packages("ggthemes")

#DELTE AFTER WRITING

# Load packages
library(ggplot2)
library(ggthemes)
library(dplyr)

# Set working directory
setwd("C:/Users/jps558/OneDrive - University of York/Desktop/R_WorkingDirectory/R_WFHDirectory/scCC") 

# Column names
pyCCcol_names = c(
  "Chr","Start","End","Experiment Insertions","Reference Inserions",
  "Expected Insertions","pvalue","Fraction Experiment","TPH Experiment",
  "TPH background subtracted","pvalue_adj Reference","Number of Overlaps with ChIA-PET Peaks"
)
pyCCcol_q_names = c("Chr","Start","End","Experiment Insertions","Strand","barcode","Number of Overlaps with ChIA-PET Peaks")

# Read Trimmed qBED data
ESR1_pyCC_Trim_Anchor1 <- read.table(
  "result_trim_minimum/1kbp_A_JS_ES1_allpass_adapter_trimmed.tagged.bam_final.sorted.ccf_B_chiaanchor1.bed",
  header = TRUE, sep = "\t", col.names = pyCCcol_q_names
)
ESR1_pyCC_Trim_Anchor2 <- read.table(
  "result_trim_minimum/1kbp_A_JS_ES1_allpass_adapter_trimmed.tagged.bam_final.sorted.ccf_B_chiaanchor2.bed",
  header = TRUE, sep = "\t", col.names = pyCCcol_q_names
)
ESR1_pyCC_Trim_UMI_Anchor1 <- read.table(
  "result_trim_minimum/1kbp_A_JS_ES1_allpass_adapter_trimmed.tagged.bam_UMIFilt_final.ccf_B_chiaanchor1.bed",
  header = TRUE, sep = "\t", col.names = pyCCcol_q_names
)
ESR1_pyCC_Trim_UMI_Anchor2 <- read.table(
  "result_trim_minimum/1kbp_A_JS_ES1_allpass_adapter_trimmed.tagged.bam_UMIFilt_final.ccf_B_chiaanchor2.bed",
  header = TRUE, sep = "\t", col.names = pyCCcol_q_names
)

# Count totals and calculate overlaps

total_trim <- nrow(ESR1_pyCC_Trim_Anchor1)
total_trim_UMI <- nrow(ESR1_pyCC_Trim_UMI_Anchor1)

overlap_trim_anchor1 <- sum(ESR1_pyCC_Trim_Anchor1$Number.of.Overlaps.with.ChIA.PET.Peaks > 0)
overlap_trim_anchor2 <- sum(ESR1_pyCC_Trim_Anchor2$Number.of.Overlaps.with.ChIA.PET.Peaks > 0)
overlap_trim_UMI_anchor1 <- sum(ESR1_pyCC_Trim_UMI_Anchor1$Number.of.Overlaps.with.ChIA.PET.Peaks > 0)
overlap_trim_UMI_anchor2 <- sum(ESR1_pyCC_Trim_UMI_Anchor2$Number.of.Overlaps.with.ChIA.PET.Peaks > 0)

perc_trim_anchor1 <- (overlap_trim_anchor1 / total_trim) * 100
perc_trim_anchor2 <- (overlap_trim_anchor2 / total_trim) * 100
perc_trim_UMI_anchor1 <- (overlap_trim_UMI_anchor1 / total_trim_UMI) * 100
perc_trim_UMI_anchor2 <- (overlap_trim_UMI_anchor2 / total_trim_UMI) * 100

# Build per-anchor summary table
trim_table <- data.frame(
  Name = c(
    "Trim_Anchor1", "Trim_Anchor2", 
    "Trim_UMI_Anchor1", "Trim_UMI_Anchor2"
  ),
  PercentageOverlapping = c(
    perc_trim_anchor1, perc_trim_anchor2,
    perc_trim_UMI_anchor1, perc_trim_UMI_anchor2
  )
)

# Plot per-anchor overlap
graphorder = c("Trim_Anchor1", "Trim_Anchor2", "Trim_UMI_Anchor1", "Trim_UMI_Anchor2")

ggplot(trim_table, aes(x = factor(Name, graphorder), y = PercentageOverlapping, fill = Name)) +
  geom_col(color = "black", show.legend = FALSE) +
  ylim(0, 100) +
  labs(x = "Dataset", y = "% Overlap with ChIA-PET Peak",
       title = "Trimmed pyCalling Cards Overlap: Per Anchor") +
  geom_text(aes(label = round(PercentageOverlapping, 1)), vjust = -0.5, size = 4) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

# Combine anchors: AnyAnchor

Trim_df <- data.frame(
  A1 = ESR1_pyCC_Trim_Anchor1$Number.of.Overlaps.with.ChIA.PET.Peaks,
  A2 = ESR1_pyCC_Trim_Anchor2$Number.of.Overlaps.with.ChIA.PET.Peaks
)
Trim_df$Any <- Trim_df$A1 > 0 | Trim_df$A2 > 0
perc_trim_any <- (sum(Trim_df$Any) / nrow(Trim_df)) * 100

TrimUMI_df <- data.frame(
  A1 = ESR1_pyCC_Trim_UMI_Anchor1$Number.of.Overlaps.with.ChIA.PET.Peaks,
  A2 = ESR1_pyCC_Trim_UMI_Anchor2$Number.of.Overlaps.with.ChIA.PET.Peaks
)
TrimUMI_df$Any <- TrimUMI_df$A1 > 0 | TrimUMI_df$A2 > 0
perc_trimUMI_any <- (sum(TrimUMI_df$Any) / nrow(TrimUMI_df)) * 100

any_table <- data.frame(
  Dataset = c("Trimmed", "Trimmed_UMI"),
  PercentageOverlapping = c(perc_trim_any, perc_trimUMI_any)
)

ggplot(any_table, aes(x = Dataset, y = PercentageOverlapping, fill = Dataset)) +
  geom_col(color = "black", show.legend = FALSE) +
  ylim(0, 100) +
  labs(x = "Dataset", y = "% Overlap with ChIA-PET (Any Anchor)",
       title = "Trimmed pyCalling Cards Overlap: Any Anchor") +
  geom_text(aes(label = round(PercentageOverlapping, 1)), vjust = -0.5, size = 4) +
  theme_minimal()

# Collapse qBED (if you want to deduplicate loci)

collapse_qbed <- function(df) {
  df %>%
    group_by(Chr, Start, End) %>%
    summarise(
      Experiment_Insertions = sum(`Experiment.Insertions`),
      Strand = first(Strand),
      barcode = paste(unique(barcode), collapse = ";"),
      Overlaps = sum(`Number.of.Overlaps.with.ChIA.PET.Peaks`)
    ) %>%
    ungroup()
}

# Collapse trimmed datasets
Trim_anchor1_coll <- collapse_qbed(ESR1_pyCC_Trim_Anchor1)
Trim_anchor2_coll <- collapse_qbed(ESR1_pyCC_Trim_Anchor2)
Trim_UMI_anchor1_coll <- collapse_qbed(ESR1_pyCC_Trim_UMI_Anchor1)
Trim_UMI_anchor2_coll <- collapse_qbed(ESR1_pyCC_Trim_UMI_Anchor2)

# Calculate percentage overlaps for collapsed datasets
calc_percentage <- function(df) {
  total <- nrow(df)
  perc <- (1 - sum(df$Overlaps == 0) / total) * 100
  return(perc)
}

perc_trim_coll_a1 <- calc_percentage(Trim_anchor1_coll)
perc_trim_coll_a2 <- calc_percentage(Trim_anchor2_coll)
perc_trimUMI_coll_a1 <- calc_percentage(Trim_UMI_anchor1_coll)
perc_trimUMI_coll_a2 <- calc_percentage(Trim_UMI_anchor2_coll)

coll_table <- data.frame(
  Name = c("Trim_Coll_Anchor1", "Trim_Coll_Anchor2", "Trim_UMI_Coll_Anchor1", "Trim_UMI_Coll_Anchor2"),
  PercentageOverlapping = c(perc_trim_coll_a1, perc_trim_coll_a2, perc_trimUMI_coll_a1, perc_trimUMI_coll_a2)
)

ggplot(coll_table, aes(x = factor(Name, levels = Name), y = PercentageOverlapping, fill = Name)) +
  geom_col(color = "black", show.legend = FALSE) +
  ylim(0, 100) +
  labs(x = "Dataset", y = "% Overlap with ChIA-PET",
       title = "Trimmed pyCalling Cards (Collapsed qBED)") +
  geom_text(aes(label = round(PercentageOverlapping, 1)), vjust = -0.5, size = 4) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

#Any anchor calculation
collapsed_any_df <- data.frame(
  A1 = Trim_anchor1_coll$Overlaps,
  A2 = Trim_anchor2_coll$Overlaps
)
collapsed_any_df$Any <- collapsed_any_df$A1 > 0 | collapsed_any_df$A2 > 0
perc_trim_coll_any <- (sum(collapsed_any_df$Any) / nrow(collapsed_any_df)) * 100

collapsed_any_UMI_df <- data.frame(
  A1 = Trim_UMI_anchor1_coll$Overlaps,
  A2 = Trim_UMI_anchor2_coll$Overlaps
)
collapsed_any_UMI_df$Any <- collapsed_any_UMI_df$A1 > 0 | collapsed_any_UMI_df$A2 > 0
perc_trimUMI_coll_any <- (sum(collapsed_any_UMI_df$Any) / nrow(collapsed_any_UMI_df)) * 100

# Build collapsed any-anchor table
collapsed_any_table <- data.frame(
  Dataset = c("Trimmed_Collapsed", "Trimmed_UMI_Collapsed"),
  PercentageOverlapping = c(perc_trim_coll_any, perc_trimUMI_coll_any)
)

# Plot collapsed any-anchor overlaps
# Save horizontal A4 TIFF (landscape)
tiff("PyCCCollapsed_qBED_AnyAnchor_A4_Landscape.tiff",
     width = 3508, height = 2480, res = 300)

ggplot(collapsed_any_table,
       aes(x = Dataset, y = PercentageOverlapping, fill = Dataset)) +
  geom_col(color = "black", show.legend = FALSE) +
  ylim(0, 100) +
  
  labs(x = "Dataset",
       y = "% Overlap with ChIA-PET (Any Anchor)",
       title = "Trimmed pyCalling Cards Overlap (Collapsed qBED, Any Anchor)",
       subtitle = NULL) +
  
  # +2 pt bar labels
  geom_text(aes(label = round(PercentageOverlapping, 1)),
            vjust = -0.5,
            size = 6) +     # was 4 → now +2
  
  # Increase all font sizes EXCEPT title/subtitle
  theme_minimal(base_size = 14) +   # was ~12 → now +2
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x  = element_text(size = 14),
    axis.text.y  = element_text(size = 14),
    plot.title   = element_text(size = 14, face = "bold"),    # unchanged
    plot.subtitle = element_text(size = 12)                   # unchanged
  )

dev.off()
