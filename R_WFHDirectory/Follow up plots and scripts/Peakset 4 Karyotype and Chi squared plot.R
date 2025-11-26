#THIS CODE WAS RUN USING R 4.4

#Load Libraries
library(karyoploteR)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)


#Set working directory
setwd("C:/Users/jps558/OneDrive - University of York/Desktop/R_WorkingDirectory/R_WFHDirectory") 

# Example BED file with insertions
insertions <- read.table("peaks/peak_data_ER_WT4.bed", header=FALSE)
colnames(insertions) <- c("Chr",	"Start",	"End",	"Center",	"Experiment Insertions",	"Background insertions",	"Reference Insertions",	"pvalue Reference",	"pvalue Background",	"Fraction Experiment",	"TPH Experiment",	"Fraction background",	"TPH background",	"TPH background subtracted",	"pvalue_adj Reference")

# Convert to GRanges
gr <- GRanges(
  seqnames = insertions$Chr,
  ranges = IRanges(start = insertions$Start, end = insertions$End),
  ExperimentInsertions = insertions$`Experiment Insertions`
)

## Make line plot
tiff("Lowstringency_karyoplot.tiff", width = 2481, height = 1749, res = 300)
# Create genome ideogram (e.g. human hg38)
karyoplot1 <- plotKaryotype(genome = "hg38")

# Add points for insertions
kpSegments(karyoplot1,
           chr = insertions$Chr,
           x0 = insertions$Start,
           x1 = insertions$Start,
           y0 = 0,
           y1 = 0.5,
           col = "blue")

dev.off()
#Calculate Chi-Squared
#Get Chrom lengths
genome <- BSgenome.Hsapiens.UCSC.hg38
chrom_sizes <- seqlengths(genome)

# Count observed insertions per chromosome
obs <- table(insertions$Chr) 

# Keep only matching chromosome sizes
chroms <- names(obs)
chrom_sizes_filtered <- chrom_sizes[chroms]

# Calculate expected counts proportional to chromosome length
expected <- sum(obs) * (chrom_sizes_filtered / sum(chrom_sizes_filtered))

# Run chi-squared test
chisq.test(x = as.numeric(obs), p = chrom_sizes_filtered / sum(chrom_sizes_filtered))

#Identify residuals
test <- chisq.test(x = as.numeric(obs), p = chrom_sizes_filtered / sum(chrom_sizes_filtered))
residuals <- test$stdres # standardized residuals
names(residuals) <- names(obs)
residuals

#Create order for plot
# Define desired chromosome order
desired_order <- c(paste0("chr", 1:22), "chrX", "chrY")

# Match to that order
chr_order <- match(desired_order, names(obs))

# Drop NAs (in case some chromosomes are missing)
chr_order <- chr_order[!is.na(chr_order)]

# Reorder and plot
obs <- obs[chr_order]
expected <- expected[chr_order]

tiff("Lowstringency_Chisquared_peaks.tiff", width = 2481, height = 1749, res = 300)

barplot(rbind(obs, expected), beside = TRUE,
        col = c("skyblue", "orange"),
        legend.text = c("Observed", "Expected"),
        main = "Low stringency peakset χ² across all chromosomes",
        ylab = "Number of Calling Card Peaks",
        ylim = c(0, 100),
        names.arg = names(obs), las = 2, cex.names = 0.8)

dev.off()

#Save data across scripts
# Your tag for this script's data
tag <- "Lowstringency_Peak_Chisquared"  # replace with a unique label for each script

# Create this script's dataframe (replace obs/expected/chromosomes with your actual vectors)
df_new <- data.frame(
  chromosome = names(obs),
  obs = as.numeric(obs),
  expected = as.numeric(expected),
  tag = tag
)

# Load existing data if present
if (file.exists(data_file)) {
  df_all <- read.csv(data_file, stringsAsFactors = FALSE)
  df_all <- rbind(df_all, df_new)
} else {
  df_all <- df_new
}

# Save back to CSV
write.csv(df_all, data_file, row.names = FALSE)
