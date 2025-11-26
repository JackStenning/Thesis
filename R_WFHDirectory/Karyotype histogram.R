

## Make Histogram

# Filter to include only standard chromosomes (optional) 
chrom_sizes <- chrom_sizes[names(chrom_sizes) %in% paste0("chr", c(1:22, "X", "Y"))]

# Create bins of desired size, e.g., 1Mb
bin_size <- 1e5 # 1 million bases = 1 Mb

bins <- tileGenome(chrom_sizes,
                   tilewidth=bin_size,
                   cut.last.tile.in.chrom=TRUE)

# Count insertions per bin
counts <- countOverlaps(bins, gr)

# Add counts as metadata
mcols(bins)$score <- counts

# Plot
kp <- plotKaryotype(genome = "hg38")
kpBars(kp, data = bins, y0 = 0, y1 = mcols(bins)$score / max(counts), col = "#377EB8")
