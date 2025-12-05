# Install packages if needed
if (!requireNamespace("GenomicRanges", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("IRanges", quietly = TRUE)) BiocManager::install("IRanges")
if (!requireNamespace("GenomicRanges", quietly = TRUE)) BiocManager::install("GenomicRanges")

library(GenomicRanges)


# Set working directory
setwd("C:/Users/jps558/OneDrive - University of York/Desktop/R_WorkingDirectory/R_WFHDirectory") 


# ===  Load BED files ===
load_bed <- function(file) {
  bed <- read.table(file, stringsAsFactors = FALSE, header = FALSE)
  # BED: chrom, start, end
  GRanges(seqnames = bed[,1], ranges = IRanges(start = bed[,2]+1, end = bed[,3]))
}

#Set column names
ChIP_names = c("Chr",	"Start",	"End",	"Name",	"Score",	"Strand",	"Signal Value",	"pvalue",	"qvalue",	"Peak", "Number of Overlaps with ChIA-PET Peaks")
Encode_Anchor1 = read.table("overlap/chip_Chia/count_Chia_Encode_ER_rep2_ChIP_q0.05_peaksNP.bedChiapet_hg38_combined.bed_anchor1.txt", header = TRUE, sep = "\t", col.names = ChIP_names)
Encode_Anchor2 = read.table("overlap/chip_Chia/count_Chia_Encode_ER_rep2_ChIP_q0.05_peaksNP.bedChiapet_hg38_combined.bed_anchor2.txt", header = TRUE, sep = "\t", col.names = ChIP_names)

# Make sure the two overlap columns are the same length as the number of CC peaks
length(Encode_Anchor1$Number.of.Overlaps.with.ChIA.PET.Peaks)
length(Encode_Anchor2$Number.of.Overlaps.with.ChIA.PET.Peaks)

# If they are the same, create the combined column
OverlapEither <- ifelse((Encode_Anchor1$Number.of.Overlaps.with.ChIA.PET.Peaks > 0) |
                          (Encode_Anchor2$Number.of.Overlaps.with.ChIA.PET.Peaks > 0), 1, 0)

# Build new data frame with first 3 columns of CC peaks and combined column
ChIP_with_overlap_E <- data.frame(Encode_Anchor1[,1:3], OverlapEither)

# View result
head(ChIP_with_overlap_E)

ChIP_Peaks_E <- ChIP_with_overlap_E

ChIP_Peaks_EGR <- GRanges(
  seqnames = ChIP_with_overlap_E$Chr,
  ranges = IRanges(start = ChIP_with_overlap_E$Start + 1, end = ChIP_with_overlap_E$End),
  OverlapEither = ChIP_with_overlap_E$OverlapEither  # optional metadata column
)

# Check
ChIP_Peaks_EGR

## Load ChIA
chiaPeaks1 <- read.table("overlap/BEDPE/Chiapet_hg38_combined.bed_anchor1.bed", sep = "\t")
chiaPeaks2 <- read.table("overlap/BEDPE/Chiapet_hg38_combined.bed_anchor2.bed", sep = "\t")

chiaGR1 <- GRanges(seqnames = chiaPeaks1$V1,
                   ranges = IRanges(start = chiaPeaks1$V2 + 1, end = chiaPeaks1$V3))

chiaGR2 <- GRanges(seqnames = chiaPeaks2$V1,
                   ranges = IRanges(start = chiaPeaks2$V2 + 1, end = chiaPeaks2$V3))

combinedChia <- c(chiaGR1, chiaGR2)


combinedChia <- reduce(combinedChia)


length(combinedChia)


## Load TTAA

ttaaSites <- load_bed("C:/Users/jps558/OneDrive - University of York/Desktop/hg38_TTAA_canon.bed")

# TTAA sites in ChIA-PET peaks
ov_ChIA <- findOverlaps(ttaaSites, combinedChia)
K <- length(unique(queryHits(ov_ChIA)))
cat("TTAA sites in ChIA-PET peaks (K):", K, "\n")

# TTAA sites in calling-card peaks
ov_cc <- findOverlaps(ttaaSites, ChIP_Peaks_EGR)
n <- length(unique(queryHits(ov_cc)))
cat("TTAA sites in calling-card peaks (n):", n, "\n")

# TTAA sites in both
k <- length(intersect(unique(queryHits(ov_cc)), unique(queryHits(ov_ChIA))))
cat("TTAA sites in both (k):", k, "\n")

# Hypergeometric test
N <- length(ttaaSites)
pval <- phyper(k-1, K, N-K, n, lower.tail = FALSE)
cat("Hypergeometric enrichment p-value:", pval, "\n")


# Optional: Fisher exact test 
cont_table <- matrix(c(k, n-k, K-k, N-K-(n-k)),
                     nrow = 2, byrow = TRUE)
rownames(cont_table) <- c("Calling-cards", "Other_TTAA")
colnames(cont_table) <- c("ChIA-PET", "Not_ChIA-PET")
fisher_res <- fisher.test(cont_table, alternative = "greater")
cat("Fisher exact p-value (enrichment):", fisher_res$p.value, "\n")

# Enrichment fold 
enrichment <- (k/n) / (K/N)
cat("Enrichment (observed / expected):", enrichment, "\n")
