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

# -----------------------------
# 0. Parameters & file paths
# -----------------------------
pyCCcol_names = c("Chr",	"Start",	"End",	"Experiment Insertions", "Reference Inserions",	"Expected Insertions",	"pvalue",	"Fraction Experiment",	"TPH Experiment",	"TPH background subtracted",	"pvalue_adj Reference")

full_path <- "scCC/Fullresult/peak_data_ER_test4.bed"     # Full calling-card peaks 
UMI_path  <- "scCC/FullresultUMI/peak_data_ER_test4.bed"     # UMI calling-card peaks 

bw_1    <- "peaks/wgEncodeGisChiaPetMcf7EraaSigRep1.bigWig"        
bw_2 <- "peaks/wgEncodeGisChiaPetMcf7EraaSigRep2.bigWig"       
bw_3    <- "peaks/wgEncodeGisChiaPetMcf7EraaSigRep3.bigWig"  


expand_flank <- 1000    # ±1kb expansion
shuffle_seeds_Full <- c(101, 102, 103)
shuffle_seeds_UMI  <- c(201, 202, 203)

# -----------------------------
# 1. Read pyCC peaks
# -----------------------------
Full_df <- read.table(full_path, header = TRUE, sep="\t", col.names=pyCCcol_names,
                      comment.char="", stringsAsFactors=FALSE)
UMI_df <- read.table(UMI_path, header = TRUE, sep="\t", col.names=pyCCcol_names,
                     comment.char="", stringsAsFactors=FALSE)

# Convert to GRanges
Full_peaks <- GRanges(seqnames = as.character(Full_df$Chr),
                      ranges = IRanges(start=as.integer(round(Full_df$Start)),
                                       end=as.integer(round(Full_df$End))),
                      mcols=Full_df[, setdiff(colnames(Full_df), c("Chr","Start","End")), drop=FALSE])

UMI_peaks <- GRanges(seqnames = as.character(UMI_df$Chr),
                     ranges = IRanges(start=as.integer(round(UMI_df$Start)),
                                      end=as.integer(round(UMI_df$End))),
                     mcols=UMI_df[, setdiff(colnames(UMI_df), c("Chr","Start","End")), drop=FALSE])

# -----------------------------
# Restrict peaks to chromosomes present in all BigWigs (safe)
# -----------------------------
valid_chrs_list <- lapply(list(bw_1, bw_2, bw_3), function(bw) seqlevels(import(bw, which=GRanges())))
# Only keep chromosomes that are present in both the BigWigs and in the peaks
valid_chrs <- Reduce(intersect, valid_chrs_list)

# Compute intersection with the actual peaks seqlevels to avoid invalid seqlevels error
Full_valid_chrs <- intersect(seqlevels(Full_peaks), valid_chrs)
UMI_valid_chrs  <- intersect(seqlevels(UMI_peaks),  valid_chrs)

Full_peaks <- keepSeqlevels(Full_peaks, Full_valid_chrs, pruning.mode="coarse")
UMI_peaks  <- keepSeqlevels(UMI_peaks,  UMI_valid_chrs, pruning.mode="coarse")

# -----------------------------
# 2. Safe chromosome lengths (hg38)
# -----------------------------
seqlens <- c(
  chr1=248956422, chr2=242193529, chr3=198295559, chr4=190214555,
  chr5=181538259, chr6=170805979, chr7=159345973, chr8=145138636,
  chr9=138394717, chr10=133797422, chr11=135086622, chr12=133275309,
  chr13=114364328, chr14=107043718, chr15=101991189, chr16=90338345,
  chr17=83257441, chr18=80373285, chr19=58617616, chr20=64444167,
  chr21=46709983, chr22=50818468, chrX=156040895, chrY=57227415
)

# -----------------------------
# 3. Expand peaks safely
# -----------------------------
expand_peaks <- function(gr, flank=1000) {
  gr_exp <- resize(gr, width=width(gr)+2*flank, fix="center")
  start(gr_exp)[start(gr_exp)<1] <- 1L
  # Clip to chromosome length
  chr_end <- seqlens[as.character(seqnames(gr_exp))]
  end(gr_exp)[end(gr_exp) > chr_end] <- chr_end[end(gr_exp) > chr_end]
  gr_exp
}

Full_peaks_exp <- expand_peaks(Full_peaks, expand_flank)
UMI_peaks_exp  <- expand_peaks(UMI_peaks, expand_flank)

# -----------------------------
# 4. Shuffle peaks safely
# -----------------------------
shuffle_matched <- function(peaks, seqlens, seed=1) {
  set.seed(seed)
  chr_counts <- as.numeric(table(seqnames(peaks)))
  chr_names  <- names(table(seqnames(peaks)))
  seqs <- unlist(mapply(function(ch,n) rep(ch,n), chr_names, chr_counts, SIMPLIFY=FALSE), use.names=FALSE)
  widths <- width(peaks)
  
  # Clip widths to chromosome if necessary
  for(i in seq_along(widths)) {
    ch <- seqs[i]
    w <- widths[i]
    max_start <- as.integer(seqlens[ch] - w)
    if(is.na(max_start) || max_start < 1) {
      w <- min(w, seqlens[ch]-1)
      widths[i] <- w
      max_start <- as.integer(seqlens[ch] - w)
    }
    widths[i] <- w
  }
  
  starts <- integer(length(widths))
  for(i in seq_along(widths)) {
    ch <- seqs[i]
    w <- widths[i]
    max_start <- as.integer(seqlens[ch]-w)
    starts[i] <- sample.int(max_start,1)
  }
  gr <- GRanges(seqnames=seqs, ranges=IRanges(start=starts, width=widths))
  valid <- start(gr)>=1 & end(gr) <= seqlens[as.character(seqnames(gr))]
  gr[valid]
}

Full_shufs <- lapply(shuffle_seeds_Full, function(s) expand_peaks(shuffle_matched(Full_peaks, seqlens, seed=s), expand_flank))
UMI_shufs  <- lapply(shuffle_seeds_UMI, function(s) expand_peaks(shuffle_matched(UMI_peaks, seqlens, seed=s), expand_flank))

# -----------------------------
# 5. Extract mean signal from BigWig
# -----------------------------
extract_mean_signal_from_bigwig <- function(bw_file, peaks) {
  message("Importing ", bw_file)
  bw_gr <- import(bw_file, format="BigWig", which=reduce(peaks))
  if(length(bw_gr)==0) return(rep(0, length(peaks)))
  ol <- findOverlaps(peaks, bw_gr)
  res <- numeric(length(peaks))
  if(length(ol)>0){
    qhits <- queryHits(ol)
    shits <- subjectHits(ol)
    ov_widths <- pmin(end(peaks)[qhits], end(bw_gr)[shits]) - pmax(start(peaks)[qhits], start(bw_gr)[shits]) + 1
    weighted_sum <- tapply(ov_widths * mcols(bw_gr)$score[shits], qhits, sum)
    tot_width <- tapply(ov_widths, qhits, sum)
    indices <- as.integer(names(weighted_sum))
    res[indices] <- as.numeric(weighted_sum/tot_width)
  }
  return(res)
}

# -----------------------------
# 6. Extract signals for 3 BigWigs
# -----------------------------
Full_signal_1 <- extract_mean_signal_from_bigwig(bw_1, Full_peaks_exp)
Full_signal_2 <- extract_mean_signal_from_bigwig(bw_2, Full_peaks_exp)
Full_signal_3 <- extract_mean_signal_from_bigwig(bw_3, Full_peaks_exp)

UMI_signal_1 <- extract_mean_signal_from_bigwig(bw_1, UMI_peaks_exp)
UMI_signal_2 <- extract_mean_signal_from_bigwig(bw_2, UMI_peaks_exp)
UMI_signal_3 <- extract_mean_signal_from_bigwig(bw_3, UMI_peaks_exp)

Full_shuf_1 <- unlist(lapply(Full_shufs, function(sh) extract_mean_signal_from_bigwig(bw_1, sh)))
Full_shuf_2 <- unlist(lapply(Full_shufs, function(sh) extract_mean_signal_from_bigwig(bw_2, sh)))
Full_shuf_3 <- unlist(lapply(Full_shufs, function(sh) extract_mean_signal_from_bigwig(bw_3, sh)))

UMI_shuf_1 <- unlist(lapply(UMI_shufs, function(sh) extract_mean_signal_from_bigwig(bw_1, sh)))
UMI_shuf_2 <- unlist(lapply(UMI_shufs, function(sh) extract_mean_signal_from_bigwig(bw_2, sh)))
UMI_shuf_3 <- unlist(lapply(UMI_shufs, function(sh) extract_mean_signal_from_bigwig(bw_3, sh)))

# -----------------------------
# 7. Statistical tests
# -----------------------------
empirical_pval <- function(obs, shuf){
  obs_mean <- mean(obs, na.rm=TRUE)
  p_emp <- (sum(shuf >= obs_mean, na.rm=TRUE)+1)/(length(shuf)+1)
  list(obs_mean=obs_mean, p_emp=p_emp)
}

# Wilcoxon
Full_wilcox_1 <- wilcox.test(Full_signal_1, Full_shuf_1, alternative="greater")
Full_wilcox_2 <- wilcox.test(Full_signal_2, Full_shuf_2, alternative="greater")
Full_wilcox_3 <- wilcox.test(Full_signal_3, Full_shuf_3, alternative="greater")

UMI_wilcox_1 <- wilcox.test(UMI_signal_1, UMI_shuf_1, alternative="greater")
UMI_wilcox_2 <- wilcox.test(UMI_signal_2, UMI_shuf_2, alternative="greater")
UMI_wilcox_3 <- wilcox.test(UMI_signal_3, UMI_shuf_3, alternative="greater")

# Empirical
Full_emp_1 <- empirical_pval(Full_signal_1, Full_shuf_1)
Full_emp_2 <- empirical_pval(Full_signal_2, Full_shuf_2)
Full_emp_3 <- empirical_pval(Full_signal_3, Full_shuf_3)
UMI_emp_1 <- empirical_pval(UMI_signal_1, UMI_shuf_1)
UMI_emp_2 <- empirical_pval(UMI_signal_2, UMI_shuf_2)
UMI_emp_3 <- empirical_pval(UMI_signal_3, UMI_shuf_3)

# -----------------------------
# 8. Prepare plotting dataframe
# -----------------------------
plot_df <- bind_rows(
  data.frame(type = paste0("Full_obs_BW1"),   value = Full_signal_1),
  data.frame(type = paste0("Full_obs_BW2"),   value = Full_signal_2),
  data.frame(type = paste0("Full_obs_BW3"),   value = Full_signal_3),
  data.frame(type = paste0("Full_shuf_BW1"),  value = Full_shuf_1),
  data.frame(type = paste0("Full_shuf_BW2"),  value = Full_shuf_2),
  data.frame(type = paste0("Full_shuf_BW3"),  value = Full_shuf_3),
  data.frame(type = paste0("UMI_obs_BW1"),    value = UMI_signal_1),
  data.frame(type = paste0("UMI_obs_BW2"),    value = UMI_signal_2),
  data.frame(type = paste0("UMI_obs_BW3"),    value = UMI_signal_3),
  data.frame(type = paste0("UMI_shuf_BW1"),   value = UMI_shuf_1),
  data.frame(type = paste0("UMI_shuf_BW2"),   value = UMI_shuf_2),
  data.frame(type = paste0("UMI_shuf_BW3"),   value = UMI_shuf_3)
)

obs_vs_shuf <- list(
  c("Full_obs_BW1", "Full_shuf_BW1"),
  c("Full_obs_BW2", "Full_shuf_BW2"),
  c("Full_obs_BW3", "Full_shuf_BW3"),
  c("UMI_obs_BW1",  "UMI_shuf_BW1"),
  c("UMI_obs_BW2",  "UMI_shuf_BW2"),
  c("UMI_obs_BW3",  "UMI_shuf_BW3")
)

Full_vs_UMI <- list(
  c("Full_obs_BW1", "UMI_obs_BW1"),
  c("Full_obs_BW2", "UMI_obs_BW2"),
  c("Full_obs_BW3", "UMI_obs_BW3")
)

my_comparisons <- c(obs_vs_shuf, Full_vs_UMI)
# -------------------------------
# 3a. Plot all data (observed + shuffled)
# -------------------------------
p_all <- ggplot(plot_df, aes(x = type, y = value, fill = type)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, fill = "white") +
  geom_jitter(width = 0.1, alpha = 0.4, size = 0.8) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Condition", y = "Signal enrichment", title = "ChIP-seq enrichment (All samples)") +
  scale_fill_brewer(palette = "Paired")  # supports up to 12 colors

# -------------------------------
# 3b. Plot only observed calling card data
# -------------------------------
obs_df <- plot_df %>% filter(grepl("_obs_", type))

obs_comparisons <- list(
  c("Full_obs_BW1", "UMI_obs_BW1"),
  c("Full_obs_BW2", "UMI_obs_BW2"),
  c("Full_obs_BW3", "UMI_obs_BW3")
)

p_obs <- ggplot(obs_df, aes(x = type, y = value, fill = type)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, fill = "white") +
  geom_jitter(width = 0.1, alpha = 0.4, size = 0.8) +
  stat_compare_means(comparisons = obs_comparisons, method = "wilcox.test", label = "p.signif") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Condition", y = "Signal enrichment", title = "ChIP-seq enrichment (Observed only)") +
  scale_fill_brewer(palette = "Paired")

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