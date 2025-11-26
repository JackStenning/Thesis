###############################################################################
# Calling-card enrichment pipeline (integrated with pyCC-style input)
# - Reads pyCC peaks (many columns) using user-provided column names
# - Expands peaks +/- 1000 bp
# - Generates 3 reproducible matched shuffled peak sets
# - Extracts mean ChIP signal from BigWig (weighted mean across overlaps)
# - Wilcoxon tests + empirical p-values against 3 shuffles
# - Optional dummy ChIP control (score permutation)
###############################################################################

# libraries
library(GenomicRanges)
library(rtracklayer)
library(GenomicAlignments) # not strictly required here but often useful
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

#Set working directory

setwd("C:/Users/jps558/OneDrive - University of York/Desktop/R_WorkingDirectory/R_WFHDirectory")
# -----------------------------
# 0. Parameters & file paths
# -----------------------------
pyCCcol_names = c("Chr",	"Start",	"End",	"Experiment Insertions", "Reference Inserions",	"Expected Insertions",	"pvalue",	"Fraction Experiment",	"TPH Experiment",	"TPH background subtracted",	"pvalue_adj Reference")

full_path <- "scCC/Fullresult/peak_data_ER_test4.bed"     # Full calling-card peaks 
UMI_path  <- "scCC/FullresultUMI/peak_data_ER_test4.bed"     # UMI calling-card peaks 

bw_geo    <- "peaks/SRR6652020_A90.bw"        
bw_encode <- "peaks/ENCFF063JMY.bw"           

expand_flank <- 1000    # ±1kb expansion
shuffle_seeds_Full <- c(101, 102, 103)
shuffle_seeds_UMI  <- c(201, 202, 203)

# -----------------------------
# 1. Read pyCC peak files as tables (preserve float columns)
# -----------------------------
message("Reading pyCC peak files as tables...")
Full_df <- read.table(full_path, header = TRUE, sep = "\t", col.names = pyCCcol_names, comment.char = "", stringsAsFactors = FALSE)
UMI_df  <- read.table(UMI_path,  header = TRUE, sep = "\t", col.names = pyCCcol_names, comment.char = "", stringsAsFactors = FALSE)

# Basic sanity checks
if(!all(c("Chr","Start","End") %in% colnames(Full_df))) stop("Full peaks missing required columns")
if(!all(c("Chr","Start","End") %in% colnames(UMI_df)))  stop("UMI peaks missing required columns")

# Convert to GRanges (coerce start/end to integer to be safe)
Full_peaks <- GRanges(seqnames = as.character(Full_df$Chr),
                      ranges   = IRanges(start = as.integer(round(Full_df$Start)),
                                         end   = as.integer(round(Full_df$End))),
                      mcols = Full_df[ , setdiff(colnames(Full_df), c("Chr","Start","End")), drop = FALSE])

UMI_peaks <- GRanges(seqnames = as.character(UMI_df$Chr),
                     ranges   = IRanges(start = as.integer(round(UMI_df$Start)),
                                        end   = as.integer(round(UMI_df$End))),
                     mcols = UMI_df[ , setdiff(colnames(UMI_df), c("Chr","Start","End")), drop = FALSE])

message("Number of peaks: Full = ", length(Full_peaks), ", UMI = ", length(UMI_peaks))

# -----------------------------
# 2. Expand peaks ±1kb
# -----------------------------
expand_peaks <- function(gr, flank = 1000) {
  # ensure boundaries are >=1 after expansion
  gr_exp <- resize(gr, width = width(gr) + 2*flank, fix = "center")
  start(gr_exp)[start(gr_exp) < 1] <- 1L
  gr_exp
}

Full_peaks_exp <- expand_peaks(Full_peaks, expand_flank)
UMI_peaks_exp  <- expand_peaks(UMI_peaks,  expand_flank)

# -----------------------------
# 3. Get genome seqlengths from one of the bigWigs
# -----------------------------
message("Reading BigWig header to obtain genome seqlengths (this does not import all signal)...")
seqlens <- seqlengths(import(bw_geo, which=GRanges())) # import(...) with empty which reads header
# fallback if that returns NA or empty, try import without which (may be heavy)
if(is.null(seqlens) || length(seqlens) == 0) {
  seqlens <- seqlengths(import(bw_geo))
}
# filter to remove NA lengths
seqlens <- seqlens[!is.na(seqlens)]
message("Chromosomes available from bigWig: ", paste(names(seqlens), collapse=", "))

# -----------------------------
# 4. Function to create matched shuffled peaks (match chr counts & widths)
#    Note: shuffle is done BEFORE expansion; we then expand each shuffled set
# -----------------------------
shuffle_matched <- function(peaks, seqlens, seed=1) {
  set.seed(seed)
  # match chromosome distribution
  chr_counts <- as.numeric(table(seqnames(peaks)))
  chr_names  <- names(table(seqnames(peaks)))
  # build vector of chromosomes to sample from (preserve counts)
  seqs <- unlist(mapply(function(ch, n) rep(ch, n), chr_names, chr_counts, SIMPLIFY = FALSE), use.names = FALSE)
  # use widths of input peaks (unexpanded widths)
  widths <- width(peaks)
  if(length(widths) != length(seqs)) {
    # if for some reason counts mismatch, fall back to sampling chromosomes proportional to seqlens
    seqs <- sample(names(seqlens), length(widths), replace = TRUE)
  }
  # generate starts respecting chromosome length
  starts <- integer(length(widths))
  for(i in seq_along(widths)) {
    ch <- seqs[i]
    w  <- widths[i]
    max_start <- as.integer(seqlens[ch] - w)
    if(is.na(max_start) || max_start < 1) {
      stop("Chromosome ", ch, " is too small or missing in seqlens.")
    }
    starts[i] <- sample.int(max_start, 1)
  }
  gr <- GRanges(seqnames = seqs, ranges = IRanges(start = starts, width = widths))
  # remove any ranges that accidentally overlap end-of-chromosome (defensive)
  valid <- start(gr) >= 1 & end(gr) <= seqlens[as.character(seqnames(gr))]
  gr[valid]
}

# create 3 shuffles for Full and UMI (before expansion), then expand
Full_shufs <- lapply(shuffle_seeds_Full, function(s) expand_peaks(shuffle_matched(Full_peaks, seqlens, seed=s), expand_flank))
UMI_shufs  <- lapply(shuffle_seeds_UMI,  function(s) expand_peaks(shuffle_matched(UMI_peaks,  seqlens, seed=s), expand_flank))

# ensure the length (#regions) match original expanded sets; if some were dropped, warn
if(length(Full_shufs[[1]]) != length(Full_peaks_exp)) warning("Some Full shuffled peaks were dropped due to chr length constraints; counts may differ.")
if(length(UMI_shufs[[1]])  != length(UMI_peaks_exp))  warning("Some UMI shuffled peaks were dropped due to chr length constraints; counts may differ.")

# -----------------------------
# 5. Function: extract mean ChIP signal per peak from a BigWig
#    This imports only the regions overlapping peaks (fastish) and computes a
#    length-weighted mean score per peak (so wide overlaps handled properly).
# -----------------------------
extract_mean_signal_from_bigwig <- function(bw_file, peaks) {
  # import the bigWig only for the ranges of interest (this produces a GRanges of segments with 'score')
  message("Importing BigWig for ", length(peaks), " peaks from ", bw_file, " (this may still take some time)...")
  bw_gr <- import(bw_file, format="BigWig", which = reduce(peaks))
  if(length(bw_gr) == 0) {
    warning("No signal imported from bigWig; returning zeros.")
    return(rep(0, length(peaks)))
  }
  # find overlaps between bw segments and peaks
  ol <- findOverlaps(peaks, bw_gr)
  # initialize result
  res <- numeric(length(peaks))
  # for each peak, sum score * overlap_width and divide by peak width (weighted mean)
  if(length(ol) > 0) {
    # compute overlap widths
    qhits <- queryHits(ol)
    shits <- subjectHits(ol)
    ov_widths <- pmin(end(peaks)[qhits], end(bw_gr)[shits]) - pmax(start(peaks)[qhits], start(bw_gr)[shits]) + 1
    # accumulate weighted sums per peak
    weighted_sum <- tapply(ov_widths * mcols(bw_gr)$score[shits], qhits, sum)
    tot_width <- tapply(ov_widths, qhits, sum)
    # assign weighted means
    indices <- as.integer(names(weighted_sum))
    res[indices] <- as.numeric(weighted_sum / tot_width)
  }
  # For peaks with no overlapping bigWig segments -> 0 (or NA if you prefer)
  # Optionally: set NA where no overlap: res[res == 0] <- NA
  return(res)
}

# -----------------------------
# 6. Extract signals
# -----------------------------
message("Extracting signals for Full peaks (GEO & ENCODE)...")
Full_geo_signal   <- extract_mean_signal_from_bigwig(bw_geo, Full_peaks_exp)
Full_encode_signal<- extract_mean_signal_from_bigwig(bw_encode, Full_peaks_exp)

message("Extracting signals for UMI peaks (GEO & ENCODE)...")
UMI_geo_signal    <- extract_mean_signal_from_bigwig(bw_geo, UMI_peaks_exp)
UMI_encode_signal <- extract_mean_signal_from_bigwig(bw_encode, UMI_peaks_exp)

# shuffled signals (concatenate three shuffles for background distribution)
Full_geo_sh_signals <- unlist(lapply(Full_shufs, function(sh) extract_mean_signal_from_bigwig(bw_geo, sh)))
UMI_geo_sh_signals  <- unlist(lapply(UMI_shufs,  function(sh) extract_mean_signal_from_bigwig(bw_geo, sh)))

Full_encode_sh_signals <- unlist(lapply(Full_shufs, function(sh) extract_mean_signal_from_bigwig(bw_encode, sh)))
UMI_encode_sh_signals  <- unlist(lapply(UMI_shufs,  function(sh) extract_mean_signal_from_bigwig(bw_encode, sh)))

# -----------------------------
# 7. Statistical tests
#    Wilcoxon rank-sum (Mann-Whitney) and simple empirical p-value vs 3 shuffles
# -----------------------------
message("Performing statistical tests...")

# function to compute empirical p-value given observed vector and list/vector of shuffled signals
empirical_pval <- function(obs_vals, shuf_vals) {
  # here shuf_vals is a vector of combined shuffled per-peak means
  obs_mean <- mean(obs_vals, na.rm = TRUE)
  # compute how many shuffled means >= obs_mean on a per-shuffle aggregated basis:
  # but we only have many per-peak shuffled values; we will compute per-shuffle means if needed.
  # Simpler: compute combined distribution of per-peak shuffled means and compare:
  p_emp <- (sum(shuf_vals >= obs_mean, na.rm=TRUE) + 1) / (length(shuf_vals) + 1) # +1 correction
  list(obs_mean = obs_mean, p_emp = p_emp)
}

# Wilcoxon (non-parametric) comparing per-peak distributions
Full_wilcox_geo <- wilcox.test(Full_geo_signal, Full_geo_sh_signals, alternative = "greater")
UMI_wilcox_geo  <- wilcox.test(UMI_geo_signal,  UMI_geo_sh_signals,  alternative = "greater")

Full_wilcox_encode <- wilcox.test(Full_encode_signal, Full_encode_sh_signals, alternative = "greater")
UMI_wilcox_encode  <- wilcox.test(UMI_encode_signal,  UMI_encode_sh_signals,  alternative = "greater")

# empirical p-values
Full_emp_geo <- empirical_pval(Full_geo_signal, Full_geo_sh_signals)
UMI_emp_geo  <- empirical_pval(UMI_geo_signal,  UMI_geo_sh_signals)

Full_emp_encode <- empirical_pval(Full_encode_signal, Full_encode_sh_signals)
UMI_emp_encode  <- empirical_pval(UMI_encode_signal,  UMI_encode_sh_signals)

# report
message("---- Results (GEO bigWig) ----")
message("Full peaks: Wilcoxon p = ", signif(Full_wilcox_geo$p.value, 4),
        "; empirical p (vs 3 shuffles combined) = ", signif(Full_emp_geo$p_emp, 4),
        "; obs mean = ", signif(Full_emp_geo$obs_mean,4))
message("UMI peaks:  Wilcoxon p = ", signif(UMI_wilcox_geo$p.value, 4),
        "; empirical p (vs 3 shuffles combined) = ", signif(UMI_emp_geo$p_emp, 4),
        "; obs mean = ", signif(UMI_emp_geo$obs_mean,4))

message("---- Results (ENCODE bigWig) ----")
message("Full peaks: Wilcoxon p = ", signif(Full_wilcox_encode$p.value, 4),
        "; empirical p = ", signif(Full_emp_encode$p_emp, 4),
        "; obs mean = ", signif(Full_emp_encode$obs_mean,4))
message("UMI peaks:  Wilcoxon p = ", signif(UMI_wilcox_encode$p.value, 4),
        "; empirical p = ", signif(UMI_emp_encode$p_emp, 4),
        "; obs mean = ", signif(UMI_emp_encode$obs_mean,4))

# -----------------------------
# 9. Plots
# -----------------------------
# -------------------------------
# 1. Prepare long-format dataframe
# -------------------------------
plot_df <- bind_rows(
  data.frame(type="Full_obs_GEO",       value=Full_geo_signal),
  data.frame(type="Full_shuf_GEO",      value=Full_geo_sh_signals),
  data.frame(type="UMI_obs_GEO",        value=UMI_geo_signal),
  data.frame(type="UMI_shuf_GEO",       value=UMI_geo_sh_signals),
  data.frame(type="Full_obs_ENCODE",    value=Full_encode_signal),
  data.frame(type="Full_shuf_ENCODE",   value=Full_encode_sh_signals),
  data.frame(type="UMI_obs_ENCODE",     value=UMI_encode_signal),
  data.frame(type="UMI_shuf_ENCODE",    value=UMI_encode_sh_signals)
)

# -------------------------------
# 2. Define relevant comparisons
# -------------------------------
obs_vs_shuf <- list(
  c("Full_obs_GEO", "Full_shuf_GEO"),
  c("UMI_obs_GEO",  "UMI_shuf_GEO"),
  c("Full_obs_ENCODE", "Full_shuf_ENCODE"),
  c("UMI_obs_ENCODE",  "UMI_shuf_ENCODE")
)

Full_vs_UMI <- list(
  c("Full_obs_GEO", "UMI_obs_GEO"),
  c("Full_obs_ENCODE", "UMI_obs_ENCODE")
)

my_comparisons <- c(obs_vs_shuf, Full_vs_UMI)

# -------------------------------
# 3a. Plot all data (observed + shuffled)
# -------------------------------
p_all <- ggplot(plot_df, aes(x = type, y = value, fill = type)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, fill = "white") +
  geom_jitter(width = 0.1, alpha = 0.4, size = 0.8) +
  stat_compare_means(comparisons = my_comparisons, 
                     method = "wilcox.test",
                     exact = FALSE,       
                     label = "p.signif") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(x="Condition", y="Signal enrichment", title="ChIP-seq enrichment (All samples)") +
  scale_fill_brewer(palette = "Set2")

# -------------------------------
# 3b. Plot only observed calling card data
# -------------------------------
obs_df <- plot_df %>% 
  filter(grepl("_obs_", type))  # only observed

# Only compare Full vs UMI within observed
obs_comparisons <- list(
  c("Full_obs_GEO", "UMI_obs_GEO"),
  c("Full_obs_ENCODE", "UMI_obs_ENCODE")
)

p_obs <- ggplot(obs_df, aes(x = type, y = value, fill = type)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, fill = "white") +
  geom_jitter(width = 0.1, alpha = 0.4, size = 0.8) +
  stat_compare_means(comparisons = obs_comparisons, 
                     method = "wilcox.test",
                     exact = FALSE,       # <<< HERE
                     label = "p.signif") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(x="Condition", y="Signal enrichment", title="ChIP-seq enrichment (Observed only)") +
  scale_fill_brewer(palette = "Set2")

# -------------------------------
# 4. Print plots
# -------------------------------
print(p_all)
print(p_obs)

###############################################################################
# Simple function to save plots only (from Script 1, adapted for Script 2)
###############################################################################
save_plots_tiff <- function(outdir = "scCC_ChIP_enrichment_plots_tiff", plots = list()){
  if(!dir.exists(outdir)){
    dir.create(outdir, recursive = TRUE)
  }
  for(nm in names(plots)){
    ggsave(
      filename = file.path(outdir, paste0(nm, ".tiff")),
      plot = plots[[nm]],
      width = 10, height = 6,
      dpi = 300,            # high-quality publication resolution
      compression = "lzw",  # standard compression for TIFF
      device = "tiff"
    )
  }
  message("TIFF plots saved to: ", outdir)
}

# Call the saving function using the Script 2 plot objects
save_plots_tiff(
  outdir = "scCC_ChIP_enrichment_plots_tiff",
  plots = list(
    all_data = p_all,
    observed_only = p_obs
  )
)