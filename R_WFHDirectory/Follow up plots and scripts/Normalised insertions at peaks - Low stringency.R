###############################################################
### Full Calling Card IPM Reproducibility Script (All Pairs) ###
### Low Stringency (Peakset 4) + TIFF Output Added
###############################################################

library(ggplot2)
library(dplyr)
library(IRanges)

setwd("C:/Users/jps558/OneDrive - University of York/Desktop/R_WorkingDirectory/R_WFHDirectory") 

###############################################################
### Create output directory for TIFF files
###############################################################

output_dir <- "Peakset4_LowStringency_Plots"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

###############################################################
### Load Peak Data
###############################################################

qbedcol_names = c("Chr","Start","End","Number.of.Insertions","Strand","Barcode")

pyCCcol_names = c(
  "Chr","Start","End","Center","Experiment.Insertions","Background.insertions",
  "Reference.Insertions","pvalue.Reference","pvalue.Background","Fraction.Experiment",
  "TPH.Experiment","Fraction.background","TPH.background","TPH.background.subtracted",
  "pvalue_adj.Reference"
)

ER_peak2 = read.table("Normalising replicates/Peaks/peak_data_ER_WT2.bed",
                      header = TRUE, sep = "\t", col.names = pyCCcol_names)

ER_peak4 = read.table("Normalising replicates/Peaks/peak_data_ER_WT4.bed",
                      header = TRUE, sep = "\t", col.names = pyCCcol_names)

ER_peak4_subset <- ER_peak4[, c("Chr","Start","End","TPH.Experiment")]


###############################################################
### Load Replicate Data (Peakset 4)
###############################################################

ER_Rep1_4 = read.table("Normalising replicates/Peaks/a_filtered_sortedER1_EKDL230008858.qbed_b_peakset4_.bed",
                       header = TRUE, sep = "\t", col.names = qbedcol_names)

ER_Rep2_4 = read.table("Normalising replicates/Peaks/a_filtered_sortedER2_EKDL230008858.qbed_b_peakset4_.bed",
                       header = TRUE, sep = "\t", col.names = qbedcol_names)

ER_Rep3_4 = read.table("Normalising replicates/Peaks/a_filtered_sortedER3_EKDL230008858.qbed_b_peakset4_.bed",
                       header = TRUE, sep = "\t", col.names = qbedcol_names)

ER_Rep5_4 = read.table("Normalising replicates/Peaks/a_filtered_sortedER5_EKDL230008858.qbed_b_peakset4_.bed",
                       header = TRUE, sep = "\t", col.names = qbedcol_names)

ER_Rep6_4 = read.table("Normalising replicates/Peaks/a_filtered_sortedER6_EKDL230008858.qbed_b_peakset4_.bed",
                       header = TRUE, sep = "\t", col.names = qbedcol_names)

# subset first 4 columns
ER_Rep1_4_subset = ER_Rep1_4[, 1:4]
ER_Rep2_4_subset = ER_Rep2_4[, 1:4]
ER_Rep3_4_subset = ER_Rep3_4[, 1:4]
ER_Rep5_4_subset = ER_Rep5_4[, 1:4]
ER_Rep6_4_subset = ER_Rep6_4[, 1:4]


###############################################################
### Combine Replicate Insertions
###############################################################

combine_rep <- function(df) {
  df %>%
    group_by(Chr, Start, End) %>%
    summarise(Number.of.Insertions = sum(Number.of.Insertions, na.rm = TRUE),
              .groups = "drop")
}

replicate_list <- list(
  Rep1 = combine_rep(ER_Rep1_4_subset),
  Rep2 = combine_rep(ER_Rep2_4_subset),
  Rep3 = combine_rep(ER_Rep3_4_subset),
  Rep5 = combine_rep(ER_Rep5_4_subset),
  Rep6 = combine_rep(ER_Rep6_4_subset)
)


###############################################################
### Map Insertions to Peaks
###############################################################

ir_peak <- IRanges(start = ER_peak4_subset$Start, end = ER_peak4_subset$End)

map_to_peaks <- function(rep_df) {
  ir_rep <- IRanges(start = rep_df$Start, end = rep_df$End)
  ov <- findOverlaps(ir_peak, ir_rep, maxgap = 1000)
  
  insert_vec <- numeric(length(ir_peak))
  insert_vec[unique(queryHits(ov))] <- tapply(
    rep_df$Number.of.Insertions[subjectHits(ov)],
    queryHits(ov),
    sum
  )
  
  return(insert_vec)
}

mapped_insertions <- lapply(replicate_list, map_to_peaks)


###############################################################
### Function: Compare two replicates, calculate IPM, save TIFFs
###############################################################

plot_replicate_pair <- function(repA, repB, nameA, nameB, output_dir) {
  
  df <- data.frame(
    Chr = ER_peak4_subset$Chr,
    Start = ER_peak4_subset$Start,
    End = ER_peak4_subset$End,
    A = repA,
    B = repB
  ) %>%
    filter(A != 0 & B != 0)
  
  # IPM
  df$IPM_A <- (df$A / sum(df$A)) * 1e6
  df$IPM_B <- (df$B / sum(df$B)) * 1e6
  
  df$log10_A <- log10(df$IPM_A)
  df$log10_B <- log10(df$IPM_B)
  
  corr <- cor(df$IPM_A, df$IPM_B)
  corr_log <- cor(df$log10_A, df$log10_B)
  
  ### ---- IPM Plot ---- ###
  p1 <- ggplot(df, aes(IPM_A, IPM_B)) +
    geom_point(color="blue", size=2) +
    geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
    annotate("text", x=max(df$IPM_A)*0.7, y=max(df$IPM_B)*0.9,
             label=paste("Pearson r =", round(corr,2)), size=5) +
    labs(title=paste("IPM Reproducibility:", nameA, "vs", nameB),
         x=paste(nameA, "(IPM)"), y=paste(nameB, "(IPM)")) +
    theme_minimal()
  
  ### ---- log10(IPM) Plot ---- ###
  p2 <- ggplot(df, aes(log10_A, log10_B)) +
    geom_point(color="blue", size=2) +
    geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
    annotate("text", x=max(df$log10_A)*0.7, y=max(df$log10_B)*0.9,
             label=paste("Pearson r =", round(corr_log,2)), size=5) +
    labs(title=paste("log10(IPM) Reproducibility:", nameA, "vs", nameB),
         x=paste(nameA, "(log10 IPM)"), y=paste(nameB, "(log10 IPM)")) +
    theme_light()
  
  ### ---- Save TIFFs ---- ###
  tiff(
    filename = file.path(output_dir, paste0("IPM_", nameA, "_vs_", nameB, ".tiff")),
    width = 7, height = 7, units = "in", res = 300, compression = "lzw"
  )
  print(p1)
  dev.off()
  
  tiff(
    filename = file.path(output_dir, paste0("log10IPM_", nameA, "_vs_", nameB, ".tiff")),
    width = 7, height = 7, units = "in", res = 300, compression = "lzw"
  )
  print(p2)
  dev.off()
  
  return(list(data=df, plot_ipm=p1, plot_log=p2))
}

###############################################################
### All Pairwise Replicate Comparisons
###############################################################

replicate_pairs <- combn(names(mapped_insertions), 2, simplify = FALSE)

results <- list()

for (pair in replicate_pairs) {
  A <- pair[1]
  B <- pair[2]
  
  results[[paste(A, B, sep="_vs_")]] <- plot_replicate_pair(
    mapped_insertions[[A]],
    mapped_insertions[[B]],
    A,
    B,
    output_dir
  )
}

###############################################################
### Print plots to screen (optional)
###############################################################

for (nm in names(results)) {
  print(results[[nm]]$plot_ipm)
  print(results[[nm]]$plot_log)
}

facet4 <- build_facet_plots(results, "Peakset 4 (Low Stringency)")
print(facet4$IPM_facet)
print(facet4$log10_facet)


# Peakset 4
facet4 <- build_facet_plots_with_corr(
  results,
  "Peakset 4 (Low Stringency)",
  tiff_path_ipm = "Peakset4_IPM_facet.tiff",
  tiff_path_log = "Peakset4_log10IPM_facet.tiff"
)


# Print to screen if desired
print(facet4$IPM_facet)
print(facet4$log10_facet)



###############################################################
### END OF SCRIPT
###############################################################
