################################################################################
# BRD4: TTAA-based hypergeometric pipeline + violin + bar plots (saved)
#   - Uses BRD4 calling-card samples (WT2_Liu, WT2_Zhang, WT4_Liu, WT4_Zhang)
#   - Uses BRD4 ChIP union as target
#   - Outputs Excel/CSV and plots (violin of TTAA overlap; bar of enrichment & -log10p)
################################################################################

# ---- Packages ----
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("GenomicRanges", quietly = TRUE)) BiocManager::install("GenomicRanges")
if (!requireNamespace("IRanges", quietly = TRUE)) BiocManager::install("IRanges")
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("ggpubr", quietly = TRUE)) install.packages("ggpubr")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(GenomicRanges)
library(IRanges)
library(openxlsx)
library(ggplot2)
library(ggpubr)
library(dplyr)

setwd("C:/Users/jps558/OneDrive - University of York/Desktop/R_WorkingDirectory/R_WFHDirectory") 

# ---- Helper: load 3-col BED -> GRanges (0-based BED -> convert to 1-based) ----
load_bed3 <- function(path, zero_based = TRUE) {
  if (!file.exists(path)) stop("File not found: ", path)
  df <- read.table(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE, comment.char = "")
  if (ncol(df) < 3) stop("BED file has fewer than 3 columns: ", path)
  starts <- as.integer(df[[2]]); ends <- as.integer(df[[3]])
  if (zero_based) starts <- starts + 1L
  GRanges(seqnames = df[[1]], ranges = IRanges(start = starts, end = ends))
}

# ---- Paths & sample files (from your BRD4 script) ----
ttaa_file <- "C:/Users/jps558/OneDrive - University of York/Desktop/hg38_TTAA_canon.bed"

# Calling-card peak BEDs (3-column)
cc_wt2_liu_bed   <- "overlap/cc_BRD4/count_BRD4_1kbp_A_peak_data_ER_WT2.bed_B_Andy_ChIP_q0.05_peaksNP.bed"
cc_wt2_zhang_bed <- "overlap/cc_BRD4/count_BRD4Z_1kbp_A_peak_data_ER_WT2.bed_B_Andy_ChIP_q0.05_peaksNP.bed"
cc_wt4_liu_bed   <- "overlap/cc_BRD4/count_BRD4_1kbp_A_peak_data_ER_WT4.bed_B_Andy_ChIP_q0.05_peaksNP.bed"
cc_wt4_zhang_bed <- "overlap/cc_BRD4/count_BRD4Z_1kbp_A_peak_data_ER_WT4.bed_B_Andy_ChIP_q0.05_peaksNP.bed"

# BRD4 ChIP peak BEDs (3-column)
brd4_liu_bed   <- "overlap/BRD4_control_q0.05peaks.bed"
brd4_zheng_bed <- "overlap/BRD4_control_Z_q0.05peaks.bed"

# Output files
out_xlsx <- "BRD4_Hypergeom_results.xlsx"
out_csv  <- "BRD4_Hypergeom_results.csv"
plot_dir <- "BRD4_plots"

# ---- Quick file existence checks ----
needed <- c(ttaa_file,
            cc_wt2_liu_bed, cc_wt2_zhang_bed, cc_wt4_liu_bed, cc_wt4_zhang_bed,
            brd4_liu_bed, brd4_zheng_bed)
missing <- needed[!file.exists(needed)]
if (length(missing) > 0) {
  stop("Missing files (fix paths) before running:\n", paste(missing, collapse = "\n"))
}

# ---- Load TTAA universe and samples ----
ttaaSites <- load_bed3(ttaa_file)
N <- length(ttaaSites)
cat("TTAA sites (N):", N, "\n")

CC_WT2_Liu   <- load_bed3(cc_wt2_liu_bed)
CC_WT2_Zhang <- load_bed3(cc_wt2_zhang_bed)
CC_WT4_Liu   <- load_bed3(cc_wt4_liu_bed)
CC_WT4_Zhang <- load_bed3(cc_wt4_zhang_bed)

cc_list <- list(
  WT2_Liu   = CC_WT2_Liu,
  WT2_Zhang = CC_WT2_Zhang,
  WT4_Liu   = CC_WT4_Liu,
  WT4_Zhang = CC_WT4_Zhang
)

cat("Calling-card samples loaded (#regions):\n")
print(sapply(cc_list, length))

# ---- Load BRD4 ChIP and union ----
BRD4_Liu   <- load_bed3(brd4_liu_bed)
BRD4_Zheng <- load_bed3(brd4_zheng_bed)
cat("BRD4 peaks: Liu =", length(BRD4_Liu), ", Zheng =", length(BRD4_Zheng), "\n")
BRD4_union <- reduce(c(BRD4_Liu, BRD4_Zheng))
cat("BRD4 union regions:", length(BRD4_union), "\n")

# Normalize seqlevels style (choose BRD4_union as reference)
seqlevelsStyle(ttaaSites) <- seqlevelsStyle(BRD4_union)
for (nm in names(cc_list)) seqlevelsStyle(cc_list[[nm]]) <- seqlevelsStyle(BRD4_union)

# ---- Hypergeometric test function (returns one-row DF) ----
hyper_brd4 <- function(sample_name, peakGR, ttaaSites, brd4GR) {
  Nloc <- length(ttaaSites)
  ov_brd4 <- findOverlaps(ttaaSites, brd4GR)
  K <- length(unique(queryHits(ov_brd4)))
  ov_sample <- findOverlaps(ttaaSites, peakGR)
  n <- length(unique(queryHits(ov_sample)))
  k <- length(intersect(unique(queryHits(ov_brd4)), unique(queryHits(ov_sample))))
  # numeric stability: compute log.p and convert
  logp <- phyper(k - 1L, K, Nloc - K, n, lower.tail = FALSE, log.p = TRUE)
  neglog10p <- if (is.infinite(logp) && logp < 0) Inf else -as.numeric(logp) / log(10)
  pval <- if (is.infinite(logp) && logp < 0) 0 else exp(logp)
  fisher_res <- tryCatch(fisher.test(matrix(c(k, n - k, K - k, Nloc - K - (n - k)), nrow = 2, byrow = TRUE), alternative = "greater"), error = function(e) list(p.value = NA, estimate = NA))
  enrichment <- if (n == 0 || K == 0) NA else (k / n) / (K / Nloc)
  data.frame(
    Sample = sample_name,
    N = Nloc, K = K, n = n, k = k,
    Percent_CC_TTAA_in_BRD4 = ifelse(n == 0, NA, round(100 * k / n, 4)),
    Enrichment = enrichment,
    Hypergeom_p = pval,
    NegLog10_Hyper = neglog10p,
    Fisher_p = fisher_res$p.value,
    Fisher_OR = if (!is.null(fisher_res$estimate)) as.numeric(fisher_res$estimate) else NA,
    stringsAsFactors = FALSE
  )
}

# ---- Run CC x BRD4_union tests ----
res_list <- lapply(names(cc_list), function(nm) hyper_brd4(nm, cc_list[[nm]], ttaaSites, BRD4_union))
results <- do.call(rbind, res_list)
rownames(results) <- NULL
print(results)

# ---- Save results ----
write.xlsx(results, file = out_xlsx, overwrite = TRUE)
write.csv(results, file = out_csv, row.names = FALSE)
cat("Results written to:", out_xlsx, "and", out_csv, "\n")

# ---- Build per-TTAA table for violin plotting (0/1 overlap per TTAA per sample) ----
ttaa_ids <- seq_len(length(ttaaSites))
ttaa_flags <- lapply(names(cc_list), function(nm) {
  ov <- findOverlaps(ttaaSites, cc_list[[nm]])
  hits <- unique(queryHits(ov))
  data.frame(
    TTAA = ttaa_ids,
    Overlap = as.integer(ttaa_ids %in% hits),
    Sample = nm,
    stringsAsFactors = FALSE
  )
})
ttaa_plot_df <- bind_rows(ttaa_flags)

# ---- Plot saving helper ----
save_plot <- function(plot_obj, filename_base, outdir = plot_dir, width = 7, height = 5, dpi = 300) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  pdf_file <- file.path(outdir, paste0(filename_base, ".pdf"))
  png_file <- file.path(outdir, paste0(filename_base, ".png"))
  ggsave(pdf_file, plot = plot_obj, width = width, height = height)
  ggsave(png_file, plot = plot_obj, width = width, height = height, dpi = dpi)
  message("Saved plots: ", pdf_file, " , ", png_file)
}

# ---- Violin plot: TTAA overlap (0/1) per sample ----
# Use proportions (since values are 0 or 1, violin shows density of overlaps)
p_violin <- ggplot(ttaa_plot_df, aes(x = Sample, y = Overlap, fill = Sample)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.7) +
  geom_boxplot(width = 0.09, outlier.size = 0.5, fill = "white") +
  stat_compare_means(comparisons = list(c("WT2_Liu", "WT2_Zhang"), c("WT4_Liu", "WT4_Zhang"),
                                        c("WT2_Liu", "WT4_Liu"), c("WT2_Zhang", "WT4_Zhang")),
                     method = "wilcox.test", label = "p.signif") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "TTAA overlap distribution (calling-card peaks)", y = "Overlap (0/1)", x = "")

save_plot(p_violin, "TTAA_overlap_violin")

# ---- Bar plots: Enrichment and -log10(p) ----
p_enrich <- ggplot(results, aes(x = Sample, y = Enrichment, fill = Sample)) +
  geom_col(color = "black") +
  theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Fold enrichment (observed/expected)", y = "Fold enrichment", x = "")

save_plot(p_enrich, "Fold_enrichment")

p_logp <- ggplot(results, aes(x = Sample, y = NegLog10_Hyper, fill = Sample)) +
  geom_col(color = "black") +
  theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "-log10 Hypergeometric p-value", y = "-log10(p)", x = "")

save_plot(p_logp, "NegLog10_Hyper")

cat("\nAll done — results and plots saved in current directory (plots in '", plot_dir, "').\n", sep = "")
