### ============================================================
###  MASTER HYPERGEOMETRIC PIPELINE FOR MULTIPLE PEAK DATASETS
### ============================================================

## ---------------- Required packages ----------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("GenomicRanges", quietly = TRUE)) BiocManager::install("GenomicRanges")
if (!requireNamespace("IRanges", quietly = TRUE)) BiocManager::install("IRanges")
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")

library(GenomicRanges)
library(IRanges)
library(openxlsx)

## ---------------- Helper: load 3-column BED ----------------
load_bed_to_gr <- function(path, chrom.col = 1, start.col = 2, end.col = 3, zero_based_bed = TRUE) {
  df <- read.table(path, header = FALSE, stringsAsFactors = FALSE, sep = "\t")
  if (zero_based_bed) df[[start.col]] <- df[[start.col]] + 1L
  GRanges(seqnames = df[[chrom.col]],
          ranges = IRanges(start = df[[start.col]], end = df[[end.col]]))
}

## ============================================================
##  Load TTAA universe
## ============================================================
ttaa_file <- "C:/Users/jps558/OneDrive - University of York/Desktop/hg38_TTAA_canon.bed"
ttaaSites <- load_bed_to_gr(ttaa_file)

## ============================================================
##  Load ChIA-PET anchors and merge
## ============================================================
chia1_file <- "overlap/BEDPE/Chiapet_hg38_combined.bed_anchor1.bed"
chia2_file <- "overlap/BEDPE/Chiapet_hg38_combined.bed_anchor2.bed"

chiaGR1 <- load_bed_to_gr(chia1_file)
chiaGR2 <- load_bed_to_gr(chia2_file)
combinedChia <- reduce(c(chiaGR1, chiaGR2))

## ============================================================
##  Function to combine Anchor1+Anchor2 and return GRanges
## ============================================================

combine_anchor_pair <- function(anchor1_file, anchor2_file, colnames_vec) {
  A1 <- read.table(anchor1_file, header = TRUE, sep = "\t", col.names = colnames_vec)
  A2 <- read.table(anchor2_file, header = TRUE, sep = "\t", col.names = colnames_vec)
  
  stopifnot(length(A1[[length(colnames_vec)]]) == length(A2[[length(colnames_vec)]]))
  
  OverEither <- ifelse(A1[[length(colnames_vec)]] > 0 | A2[[length(colnames_vec)]] > 0, 1, 0)
  
  df3 <- data.frame(A1[,1:3], OverEither)
  
  GRanges(
    seqnames = df3[,1],
    ranges = IRanges(start = df3[,2] + 1, end = df3[,3]),
    OverlapEither = df3$OverEither
  )
}

## ============================================================
##  Load all datasets
## ============================================================

### ---- Calling cards: Low stringency ----
pyCC_names <- c("Chr","Start","End","Center","ExpIns","BackIns","RefIns",
                "pRef","pBack","FracExp","TPHexp","FracBack","TPHback",
                "TPHbackSub","padjRef","Overlaps")

Low_ccGR <- combine_anchor_pair(
  "overlap/cc_Chia/count_Chia_peak_data_ER_WT4_Chiapet_hg38_combined.bed_anchor1.txt",
  "overlap/cc_Chia/count_Chia_peak_data_ER_WT4_Chiapet_hg38_combined.bed_anchor2.txt",
  pyCC_names
)

### ---- Calling cards: High stringency ----
High_ccGR <- combine_anchor_pair(
  "overlap/cc_Chia/count_Chia_peak_data_ER_WT2_Chiapet_hg38_combined.bed_anchor1.txt",
  "overlap/cc_Chia/count_Chia_peak_data_ER_WT2_Chiapet_hg38_combined.bed_anchor2.txt",
  pyCC_names
)

### ---- ChIP-seq: ENCODE ----
ChIP_names <- c("Chr","Start","End","Name","Score","Strand","Signal","p","q",
                "Peak","Overlaps")

Encode_GR <- combine_anchor_pair(
  "overlap/chip_Chia/count_Chia_Encode_ER_rep2_ChIP_q0.05_peaksNP.bedChiapet_hg38_combined.bed_anchor1.txt",
  "overlap/chip_Chia/count_Chia_Encode_ER_rep2_ChIP_q0.05_peaksNP.bedChiapet_hg38_combined.bed_anchor2.txt",
  ChIP_names
)

### ---- ChIP-seq: GEO ----
GEO_GR <- combine_anchor_pair(
  "overlap/chip_Chia/count_Chia_Andy_ChIP_q0.05_peaksNP.bedChiapet_hg38_combined.bed_anchor1.txt",
  "overlap/chip_Chia/count_Chia_Andy_ChIP_q0.05_peaksNP.bedChiapet_hg38_combined.bed_anchor2.txt",
  ChIP_names
)

## Normalize seqlevels style
seqlevelsStyle(ttaaSites) <- seqlevelsStyle(combinedChia)
seqlevelsStyle(Low_ccGR) <- seqlevelsStyle(combinedChia)
seqlevelsStyle(High_ccGR) <- seqlevelsStyle(combinedChia)
seqlevelsStyle(Encode_GR) <- seqlevelsStyle(combinedChia)
seqlevelsStyle(GEO_GR) <- seqlevelsStyle(combinedChia)


## ============================================================
##  Function: hypergeometric analysis for one dataset
## ============================================================

hypertest <- function(sample_name, peakGR, ttaaSites, chiaGR) {
  N <- length(ttaaSites)
  
  ov_chia <- findOverlaps(ttaaSites, chiaGR)
  K <- length(unique(queryHits(ov_chia)))
  
  ov_peaks <- findOverlaps(ttaaSites, peakGR)
  n <- length(unique(queryHits(ov_peaks)))
  
  k <- length(intersect(unique(queryHits(ov_chia)),
                        unique(queryHits(ov_peaks))))
  
  logp <- phyper(k - 1, K, N - K, n, lower.tail = FALSE, log.p = TRUE)
  log10p <- -logp / log(10)
  fisher_p <- fisher.test(matrix(c(k, n - k, K - k, N - K - (n - k)),
                                 nrow = 2, byrow = TRUE),
                          alternative = "greater")$p.value
  
  enrichment <- (k / n) / (K / N)
  
  data.frame(
    Sample = sample_name,
    N_TTAA = N,
    K_TTAA_in_ChIA = K,
    n_TTAA_in_sample = n,
    k_overlap = k,
    Fold_Enrichment = enrichment,
    log10p_hyper = log10p,
    FisherP = fisher_p
  )
}

## ============================================================
##  Run for all datasets
## ============================================================

results <- rbind(
  hypertest("Low_CC",   Low_ccGR,   ttaaSites, combinedChia),
  hypertest("High_CC",  High_ccGR,  ttaaSites, combinedChia),
  hypertest("ENCODE",   Encode_GR,  ttaaSites, combinedChia),
  hypertest("GEO",      GEO_GR,     ttaaSites, combinedChia)
)

## ============================================================
##  Save Excel output
## ============================================================

write.xlsx(results, file = "Hypergeom_summary.xlsx", overwrite = TRUE)

cat("\n--- DONE! Results written to: Hypergeom_summary.xlsx ---\n")
print(results)
