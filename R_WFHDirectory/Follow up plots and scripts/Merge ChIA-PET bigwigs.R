#!/usr/bin/env Rscript

# ============================================================
# Merge three bigWig files into one (mean signal)
# ============================================================

#Set working directory
setwd("C:/Users/jps558/OneDrive - University of York/Desktop/R_WorkingDirectory/R_WFHDirectory") 


suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(IRanges)
})

# -----------------------------
# Input files (edit these)
# -----------------------------
bw1 <- "peaks/wgEncodeGisChiaPetMcf7EraaSigRep1.bigWig"
bw2 <- "peaks/wgEncodeGisChiaPetMcf7EraaSigRep2.bigWig"
bw3 <- "peaks/wgEncodeGisChiaPetMcf7EraaSigRep3.bigWig"

output_bw <- "merged_output.bw"
message("1) Importing BigWig files as GRanges (for intervals + seqlengths)...")

gr1 <- import(bw1, as = "GRanges")
gr2 <- import(bw2, as = "GRanges")
gr3 <- import(bw3, as = "GRanges")

# Get seqlengths from first file (needed for export)
sl <- seqlengths(gr1)

message("2) Creating unified interval set (reduce over all GRanges)...")

all_gr <- c(gr1, gr2, gr3)
common_gr <- reduce(all_gr)

# Prepare vector for results
merged_scores <- numeric(length(common_gr))
merged_scores[] <- NA_real_

message("3) Importing BigWig coverage as RleList...")
rle1 <- import(bw1, as = "RleList")
rle2 <- import(bw2, as = "RleList")
rle3 <- import(bw3, as = "RleList")

# Seqnames present in the unified intervals
seqs_in_common <- unique(as.character(seqnames(common_gr)))

message("4) Computing interval means chromosome-by-chromosome...")

get_means_for_seq <- function(rle_list, seqname, idx, gr) {
  if (!seqname %in% names(rle_list)) 
    return(rep(NA_real_, length(idx)))
  
  rle <- rle_list[[seqname]]
  irs <- ranges(gr[idx])
  v <- Views(rle, start(irs), end(irs))
  return(as.numeric(viewMeans(v)))
}

for (seqname in seqs_in_common) {
  idx <- which(as.character(seqnames(common_gr)) == seqname)
  if (length(idx) == 0) next
  
  m1 <- get_means_for_seq(rle1, seqname, idx, common_gr)
  m2 <- get_means_for_seq(rle2, seqname, idx, common_gr)
  m3 <- get_means_for_seq(rle3, seqname, idx, common_gr)
  
  merged_scores[idx] <- rowMeans(cbind(m1, m2, m3), na.rm = TRUE)
}

message("5) Building merged GRanges object...")

out_gr <- common_gr
out_gr$score <- merged_scores

# Harmonize seqlevels with available seqlengths
common_seqs <- intersect(seqlevels(out_gr), names(sl))

# prune ranges for dropped seqlevels
out_gr <- keepSeqlevels(out_gr, common_seqs, pruning.mode = "coarse")

# assign seqlengths
seqlengths(out_gr) <- sl[common_seqs]

message("6) Exporting merged BigWig...")

export(out_gr, output_bw, format = "bigWig")

message("Done! Merged file written to: ", output_bw)