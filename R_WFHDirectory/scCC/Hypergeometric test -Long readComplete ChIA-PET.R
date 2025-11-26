###############################
# Long-read calling-card vs ChIA-PET
# Hypergeometric pipeline (TTAA-universe)
# Compares: Full / trimmed / UMI / full peak sets

#### !!!!! FOR SOME REASON I CALLED THE ENCCODE FILE BRD4 INSTEAD- THIS IS ENCODE OVERLAP WHICH CAN BE VERIFIED IN V2_075_qbed_Overlap.SH !!!!!!!! #########


###############################

## -------- Dependencies --------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("GenomicRanges", quietly = TRUE)) BiocManager::install("GenomicRanges")
if (!requireNamespace("IRanges", quietly = TRUE)) BiocManager::install("IRanges")
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")

library(GenomicRanges)
library(IRanges)
library(openxlsx)

#Set working directory
setwd("C:/Users/jps558/OneDrive - University of York/Desktop/R_WorkingDirectory/R_WFHDirectory/scCC") 

## -------- Helpers: robust BED/table loaders --------
# Load a 3-column BED-like file (no header assumed)
load_bed3 <- function(path, zero_based = TRUE) {
  if (!file.exists(path)) stop("File not found: ", path)
  df <- read.table(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE, comment.char = "")
  if (ncol(df) < 3) stop("BED-like file must have >= 3 columns: ", path)
  starts <- as.integer(df[[2]]); ends <- as.integer(df[[3]])
  if (zero_based) starts <- starts + 1L
  GRanges(seqnames = df[[1]], ranges = IRanges(start = starts, end = ends))
}

# Load a tabular peak file that may have header and many columns.
# We assume the first 3 columns are chr,start,end if header=TRUE.
load_table_to_gr <- function(path, header = TRUE, zero_based = TRUE) {
  if (!file.exists(path)) stop("File not found: ", path)
  df <- read.table(path, header = header, sep = "\t", stringsAsFactors = FALSE, comment.char = "")
  # if there are fewer than 3 columns, error
  if (ncol(df) < 3) stop("Peak table must have >=3 columns (chr,start,end): ", path)
  # when header=TRUE and the file contains named columns, we still grab first 3
  chr <- df[[1]]; start <- as.integer(df[[2]]); end <- as.integer(df[[3]])
  if (zero_based) start <- start + 1L
  GRanges(seqnames = chr, ranges = IRanges(start = start, end = end))
}

# Combine two anchor files (anchor1 + anchor2) into a single GRanges (union)
combine_anchor_pair_gr <- function(anchor1_path, anchor2_path) {
  gr1 <- load_table_to_gr(anchor1_path, header = TRUE, zero_based = TRUE)
  gr2 <- load_table_to_gr(anchor2_path, header = TRUE, zero_based = TRUE)
  reduce(c(gr1, gr2))
}

## -------- Paths you provided (adjust if necessary) --------
# TTAA universe
ttaa_file <- "C:/Users/jps558/OneDrive - University of York/Desktop/hg38_TTAA_canon.bed"

# Long-read calling cards - GEO and ENCODE overlap samples (Untrimmed)
p_Untrimmed_GEO <- "Fullresult/1kbp_A_peak_data_ER_test4.bed_B_Andy.bed"
p_Untrimmed_ENCODE <- "Fullresult/1kbp_A_peak_data_ER_test4.bed_B_BRD4.bed"

p_Untrimmed_UMI_GEO <- "FullresultUMI/1kbp_A_peak_data_ER_test4.bed_B_Andy.bed"
p_Untrimmed_UMI_ENCODE <- "FullresultUMI/1kbp_A_peak_data_ER_test4.bed_B_BRD4.bed"

# Trimmed (peak tables)
p_trim_GEO        <- "result_trim_minimum/1kbp_A_JS_ES1_allpass_adapter_trimmed.tagged.bam_final.sorted.ccf_B_Andy.bed"
p_trim_ENCODE       <- "result_trim_minimum/1kbp_A_JS_ES1_allpass_adapter_trimmed.tagged.bam_final.sorted.ccf_B_BRD4.bed"

p_trim_UMI_GEO    <- "result_trim_minimum/1kbp_A_JS_ES1_allpass_adapter_trimmed.tagged.bam_UMIFilt_final.ccf_B_Andy.bed"
p_trim_UMI_ENCODE   <- "result_trim_minimum/1kbp_A_JS_ES1_allpass_adapter_trimmed.tagged.bam_UMIFilt_final.ccf_B_BRD4.bed"

# Untrimmed / anchor-style (Untrimmed anchors and UMI anchors)
Untrimmed_anchor1_path       <- "FinalCCFresult/1kbp_A_FullResultpeak_data_ER_test4.bed_B_chiaanchor1.bed"
Untrimmed_anchor2_path       <- "FinalCCFresult/1kbp_A_FullResultpeak_data_ER_test4.bed_B_chiaanchor2.bed"
Untrimmed_umi_anchor1_path   <- "FinalCCFresult/1kbp_A_FullUMIpeak_data_ER_test4.bed_B_chiaanchor1.bed"
Untrimmed_umi_anchor2_path   <- "FinalCCFresult/1kbp_A_FullUMIpeak_data_ER_test4.bed_B_chiaanchor2.bed"

# qbed anchor-style peak files (scCC final qbed-style)
Untrimmed_peak_qbed_anchor1  <- "FinalCCFresult/1kbp_A_Fullresult_scCC_final.sorted.ccf_B_chiaanchor1.bed"
Untrimmed_peak_qbed_anchor2  <- "FinalCCFresult/1kbp_A_Fullresult_scCC_final.sorted.ccf_B_chiaanchor2.bed"
Untrimmed_umi_qbed_anchor1   <- "FinalCCFresult/1kbp_A_Fullresult_UMIFilt_scCC_final.ccf_B_chiaanchor1.bed"
Untrimmed_umi_qbed_anchor2   <- "FinalCCFresult/1kbp_A_Fullresult_UMIFilt_scCC_final.ccf_B_chiaanchor2.bed"

# ChIA-PET anchor set (same as before)
chia1_file <- "overlap/BEDPE/Chiapet_hg38_combined.bed_anchor1.bed"
chia2_file <- "overlap/BEDPE/Chiapet_hg38_combined.bed_anchor2.bed"

# Output
out_xlsx <- "LongRead_CC_vs_ChIA_hypergeom.xlsx"

## -------- Load TTAA and ChIA-PET combined anchors --------
cat("Loading TTAA universe...\n")
ttaaSites <- load_bed3(ttaa_file, zero_based = TRUE)
N <- length(ttaaSites)
cat("N (TTAA sites):", N, "\n\n")

cat("Loading ChIA-PET anchors and combining...\n")
chiaGR1 <- load_bed3(chia1_file, zero_based = TRUE)
chiaGR2 <- load_bed3(chia2_file, zero_based = TRUE)
combinedChia <- reduce(c(chiaGR1, chiaGR2))
cat("ChIA-PET combined regions:", length(combinedChia), "\n\n")

## -------- Build sample GRanges list from the files you provided --------
samples_gr <- list()

# 1) Untrimmed
samples_gr$Untrimmed_GEO    <- load_table_to_gr(p_Untrimmed_GEO, header = TRUE, zero_based = TRUE)
samples_gr$Untrimmed_ENCODE    <- load_table_to_gr(p_Untrimmed_ENCODE, header = TRUE, zero_based = TRUE)

# 2) Untrimmed UMI (if files are same structure)
samples_gr$Untrimmed_UMI_GEO <- load_table_to_gr(p_Untrimmed_UMI_GEO, header = TRUE, zero_based = TRUE)
samples_gr$Untrimmed_UMI_ENCODE <- load_table_to_gr(p_Untrimmed_UMI_ENCODE, header = TRUE, zero_based = TRUE)

# 3) trimmed (final.ccf style)
samples_gr$trimmed_GEO      <- load_table_to_gr(p_trim_GEO, header = TRUE, zero_based = TRUE)
samples_gr$trimmed_ENCODE     <- load_table_to_gr(p_trim_ENCODE, header = TRUE, zero_based = TRUE)

samples_gr$trimmed_UMI_GEO  <- load_table_to_gr(p_trim_UMI_GEO, header = TRUE, zero_based = TRUE)
samples_gr$trimmed_UMI_ENCODE <- load_table_to_gr(p_trim_UMI_ENCODE, header = TRUE, zero_based = TRUE)

# 4) Untrimmed anchor pairs -> union per sample (Untrimmed)
samples_gr$Untrimmed_peak_union <- combine_anchor_pair_gr(Untrimmed_anchor1_path, Untrimmed_anchor2_path)

# 5) Untrimmed UMI anchor pairs -> union
samples_gr$Untrimmed_UMI_peak_union <- combine_anchor_pair_gr(Untrimmed_umi_anchor1_path, Untrimmed_umi_anchor2_path)

# 6) qbed-style anchor pairs (scCC final) -> union
samples_gr$Untrimmed_peak_qbed_union <- combine_anchor_pair_gr(Untrimmed_peak_qbed_anchor1, Untrimmed_peak_qbed_anchor2)

# 7) Untrimmed UMI qbed anchor union
samples_gr$Untrimmed_UMI_qbed_union <- combine_anchor_pair_gr(Untrimmed_umi_qbed_anchor1, Untrimmed_umi_qbed_anchor2)

# Print counts for sanity check
cat("Sample ranges (counts):\n")
print(sapply(samples_gr, length))
cat("\n")

## -------- Normalize seqlevels to combinedChia style to avoid mismatches --------
for (nm in names(samples_gr)) seqlevelsStyle(samples_gr[[nm]]) <- seqlevelsStyle(combinedChia)
seqlevelsStyle(ttaaSites) <- seqlevelsStyle(combinedChia)

## -------- hypergeometric test function (TTAA-universe) --------
hyper_test_sample <- function(sample_name, sampleGR, ttaaSites, targetChiaGR) {
  Nloc <- length(ttaaSites)
  
  # K = TTAA in ChIA-PET anchors
  ov_chia <- findOverlaps(ttaaSites, targetChiaGR)
  K <- length(unique(queryHits(ov_chia)))
  
  # n = TTAA in sample peaks
  ov_sample <- findOverlaps(ttaaSites, sampleGR)
  n <- length(unique(queryHits(ov_sample)))
  
  # k = TTAA in both
  k <- length(intersect(unique(queryHits(ov_chia)), unique(queryHits(ov_sample))))
  
  # log p-value
  logp <- phyper(k - 1L, K, Nloc - K, n, lower.tail = FALSE, log.p = TRUE)
  neglog10p <- -as.numeric(logp) / log(10)
  
  # Fisher exact (one-sided)
  cont_tab <- matrix(c(k, n - k, K - k, Nloc - K - (n - k)), nrow = 2, byrow = TRUE)
  fisher_res <- fisher.test(cont_tab, alternative = "greater")
  
  enrichment <- if (n == 0) NA else (k / n) / (K / Nloc)
  percent_overlap <- if (n == 0) NA else 100 * k / n
  
  data.frame(
    Sample = sample_name,
    N_TTAA = Nloc,
    K_in_ChIA = K,
    n_in_sample = n,
    k_overlap = k,
    Percent_overlap = round(percent_overlap, 4),
    Fold_enrichment = round(enrichment, 5),
    NegLog10_hyper = ifelse(is.infinite(neglog10p), Inf, round(neglog10p, 3)),
    Fisher_p = fisher_res$p.value,
    Fisher_OR = if (!is.null(fisher_res$estimate)) as.numeric(fisher_res$estimate) else NA,
    stringsAsFactors = FALSE
  )
}

## -------- Run test for all samples --------
res_list <- lapply(names(samples_gr), function(nm) {
  cat("Testing sample:", nm, "...\n")
  hyper_test_sample(nm, samples_gr[[nm]], ttaaSites, combinedChia)
})
results_df <- do.call(rbind, res_list)

## -------- Save output --------
write.xlsx(results_df, file = out_xlsx, overwrite = TRUE)
cat("\n=== Done. Results written to:", out_xlsx, "===\n")
print(results_df)
