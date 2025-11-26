###############################################################
### Full Calling Card IPM Reproducibility Script (Peakset 2)
### Generates All Pairwise Replicate Plots + Saves TIFF Files
###############################################################

library(ggplot2)
library(dplyr)
library(IRanges)

setwd("C:/Users/jps558/OneDrive - University of York/Desktop/R_WorkingDirectory/R_WFHDirectory") 

###############################################################
### Create output directory for TIFFs
###############################################################
output_dir <- "Peakset2_Plot"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

###############################################################
### Load column names
###############################################################

qbedcol_names = c("Chr","Start","End","Number.of.Insertions","Strand","Barcode")

pyCCcol_names = c(
  "Chr","Start","End","Center","Experiment.Insertions","Background.insertions",
  "Reference.Insertions","pvalue.Reference","pvalue.Background","Fraction.Experiment",
  "TPH.Experiment","Fraction.background","TPH.background","TPH.background.subtracted",
  "pvalue_adj.Reference"
)

###############################################################
### Load Peak Data (Peakset 2)
###############################################################

ER_peak2 = read.table(
  "Normalising replicates/Peaks/peak_data_ER_WT2.bed",
  header = TRUE, sep = "\t", col.names = pyCCcol_names
)

ER_peak2_subset <- ER_peak2[, c("Chr","Start","End","TPH.Experiment")]

# For compatibility with the rest of the script:
ER_peak4_subset <- ER_peak2_subset

###############################################################
### Load Replicate Data for Peakset 2
###############################################################

ER_Rep1_2 = read.table(
  "Normalising replicates/Peaks/a_filtered_sortedER1_EKDL230008858.qbed_b_peakset2_.bed",
  header = TRUE, sep = "\t", col.names = qbedcol_names
)

ER_Rep2_2 = read.table(
  "Normalising replicates/Peaks/a_filtered_sortedER2_EKDL230008858.qbed_b_peakset2_.bed",
  header = TRUE, sep = "\t", col.names = qbedcol_names
)

ER_Rep3_2 = read.table(
  "Normalising replicates/Peaks/a_filtered_sortedER3_EKDL230008858.qbed_b_peakset2_.bed",
  header = TRUE, sep = "\t", col.names = qbedcol_names
)

ER_Rep5_2 = read.table(
  "Normalising replicates/Peaks/a_filtered_sortedER5_EKDL230008858.qbed_b_peakset2_.bed",
  header = TRUE, sep = "\t", col.names = qbedcol_names
)

ER_Rep6_2 = read.table(
  "Normalising replicates/Peaks/a_filtered_sortedER6_EKDL230008858.qbed_b_peakset2_.bed",
  header = TRUE, sep = "\t", col.names = qbedcol_names
)

###############################################################
### Subset replicate columns (Chr, Start, End, Insertions)
###############################################################

ER_Rep1_2_subset = ER_Rep1_2[, 1:4]
ER_Rep2_2_subset = ER_Rep2_2[, 1:4]
ER_Rep3_2_subset = ER_Rep3_2[, 1:4]
ER_Rep5_2_subset = ER_Rep5_2[, 1:4]
ER_Rep6_2_subset = ER_Rep6_2[, 1:4]


###############################################################
### Combine replicate insertion rows by Chr/Start/End
###############################################################

combine_rep <- function(df) {
  df %>%
    group_by(Chr, Start, End) %>%
    summarise(Number.of.Insertions = sum(Number.of.Insertions, na.rm = TRUE),
              .groups = "drop")
}

replicate_list <- list(
  Rep1 = combine_rep(ER_Rep1_2_subset),
  Rep2 = combine_rep(ER_Rep2_2_subset),
  Rep3 = combine_rep(ER_Rep3_2_subset),
  Rep5 = combine_rep(ER_Rep5_2_subset),
  Rep6 = combine_rep(ER_Rep6_2_subset)
)


###############################################################
### Map replicate insertions to peaks (IRanges)
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
### Function: Compare replicate pair + produce plots + save TIFFs
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
  
  ### IPM plot
  p1 <- ggplot(df, aes(IPM_A, IPM_B)) +
    geom_point(color="blue", size=2) +
    geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
    annotate("text", x=max(df$IPM_A)*0.7, y=max(df$IPM_B)*0.9,
             label=paste("Pearson r =", round(corr,2)), size=5) +
    labs(title=paste("IPM Reproducibility:", nameA, "vs", nameB),
         x=paste(nameA, "(IPM)"), y=paste(nameB, "(IPM)")) +
    theme_minimal()
  
  ### log10(IPM) plot
  p2 <- ggplot(df, aes(log10_A, log10_B)) +
    geom_point(color="blue", size=2) +
    geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
    annotate("text", x=max(df$log10_A)*0.7, y=max(df$log10_B)*0.9,
             label=paste("Pearson r =", round(corr_log,2)), size=5) +
    labs(title=paste("log10(IPM) Reproducibility:", nameA, "vs", nameB),
         x=paste(nameA, "(log10 IPM)"), y=paste(nameB, "(log10 IPM)")) +
    theme_light()
  
  ### Save as TIFF
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
### Perform all pairwise comparisons
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
### FACET PLOT BUILDER FOR ALL PAIRS
###############################################################

build_facet_plots <- function(results_list, peakset_name) {
  
  # Collect all comparison data frames into one tibble
  all_df <- bind_rows(
    lapply(names(results_list), function(nm) {
      if (is.null(results_list[[nm]])) return(NULL)
      df <- results_list[[nm]]$data
      df$Comparison <- nm
      df
    }),
    .id = NULL
  )
  
  # -- IPM facet plot --
  p_ipm <- ggplot(all_df, aes(IPM_A, IPM_B)) +
    geom_point(size = 1.2, alpha = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    facet_wrap(~ Comparison, scales = "free") +
    labs(title = paste("IPM Reproducibility:", peakset_name),
         x = "Rep A (IPM)", y = "Rep B (IPM)") +
    theme_bw() +
    theme(strip.text = element_text(size = 10))
  
  # -- log10(IPM) facet plot --
  p_log <- ggplot(all_df, aes(log10_A, log10_B)) +
    geom_point(size = 1.2, alpha = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    facet_wrap(~ Comparison, scales = "free") +
    labs(title = paste("log10(IPM) Reproducibility:", peakset_name),
         x = "Rep A (log10 IPM)", y = "Rep B (log10 IPM)") +
    theme_bw() +
    theme(strip.text = element_text(size = 10))
  
  return(list(IPM_facet = p_ipm, log10_facet = p_log))
}

#Facet
facet4 <- build_facet_plots(results, "Peakset 2 (High Stringency)")
print(facet4$IPM_facet)
print(facet4$log10_facet)

build_facet_plots_with_corr <- function(results_list, peakset_name, tiff_path_ipm=NULL, tiff_path_log=NULL) {
  
  # Collect all comparison data frames and compute correlations
  all_df <- bind_rows(
    lapply(names(results_list), function(nm) {
      if (is.null(results_list[[nm]])) return(NULL)
      df <- results_list[[nm]]$data
      # Compute Pearson correlations
      corr <- cor(df$IPM_A, df$IPM_B)
      corr_log <- cor(df$log10_A, df$log10_B)
      df$Comparison <- paste0(nm, " (r=", round(corr, 2), ")")
      df$Comparison_log <- paste0(nm, " (r=", round(corr_log, 2), ")")
      df
    }),
    .id = NULL
  )
  
  # -- IPM facet plot --
  p_ipm <- ggplot(all_df, aes(IPM_A, IPM_B)) +
    geom_point(size = 1.5, alpha = 0.7, color="blue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    facet_wrap(~ Comparison, scales = "free") +
    labs(title = paste("IPM Reproducibility:", peakset_name),
         x = "Rep A (IPM)", y = "Rep B (IPM)") +
    theme_bw() +
    theme(strip.text = element_text(size = 10))
  
  # -- log10(IPM) facet plot --
  p_log <- ggplot(all_df, aes(log10_A, log10_B)) +
    geom_point(size = 1.5, alpha = 0.7, color="blue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    facet_wrap(~ Comparison_log, scales = "free") +
    labs(title = paste("log10(IPM) Reproducibility:", peakset_name),
         x = "Rep A (log10 IPM)", y = "Rep B (log10 IPM)") +
    theme_bw() +
    theme(strip.text = element_text(size = 10))
  
  # -- Save to TIFF if paths provided --
  if (!is.null(tiff_path_ipm)) {
    tiff(tiff_path_ipm, width = 12, height = 10, units = "in", res = 300, compression = "lzw")
    print(p_ipm)
    dev.off()
  }
  if (!is.null(tiff_path_log)) {
    tiff(tiff_path_log, width = 12, height = 10, units = "in", res = 300, compression = "lzw")
    print(p_log)
    dev.off()
  }
  
  return(list(IPM_facet = p_ipm, log10_facet = p_log))
}

# Peakset 2
facet2 <- build_facet_plots_with_corr(
  results,
  "Peakset 2 (High Stringency)",
  tiff_path_ipm = "Peakset2_IPM_facet.tiff",
  tiff_path_log = "Peakset2_log10IPM_facet.tiff"
)


# Print to screen if desired
print(facet2$IPM_facet)
print(facet2$log10_facet)
###############################################################
### END OF SCRIPT
###############################################################
