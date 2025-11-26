#THIS CODE WAS RUN USING R 4.4
# Install BiocManager if you don't have it yet
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#Reset repositories
options(repos = BiocManager::repositories())

# Install BiocParallel first — this is critical
BiocManager::install("BiocParallel", ask = FALSE, update = FALSE)

# Now install the rest of your packages
BiocManager::install(c(
  "rtracklayer",
  "GenomicFeatures",
  "BSgenome.Hsapiens.UCSC.hg38",
  "karyoploteR"
), force = TRUE, ask = FALSE, update = FALSE)

#Load Libraries
library(karyoploteR)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(openxlsx)


#Set working directory
setwd("C:/Users/jps558/OneDrive - University of York/Desktop/R_WorkingDirectory/R_WFHDirectory") 

# Example BED file with insertions
insertions <- read.table("peaks/Andy_ChIP_q0.05_peaksNP.bed", header=FALSE)
colnames(insertions) <- c("Chr",	"Start",	"End",	"Name",	"Score",	"Strand",	"Signal Value",	"pvalue",	"qvalue",	"Peak")

# Convert to GRanges
gr <- GRanges(
  seqnames = insertions$Chr,
  ranges = IRanges(start = insertions$Start, end = insertions$End),
  ExperimentInsertions = insertions$`Experiment Insertions`
)

## Make line plot
tiff("GEO_Karyoplot.tiff", width = 2481, height = 1749, res = 300)
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

tiff("GEO_Peak_Chisquared.tiff", width = 2481, height = 1749, res = 300)

barplot(rbind(obs, expected), beside = TRUE,
        col = c("skyblue", "orange"),
        legend.text = c("Observed", "Expected"),
        ylab = "Number of ChIP-Seq Peaks",
        main = "GEO peakset χ² across all chromosomes",
        ylim = c(0, 2000),
        names.arg = names(obs), las = 2, cex.names = 0.8)

dev.off()

#Save data across scripts
# File for storing the shared data
data_file <- "shared_chisq.csv"

# Your tag for this script's data
tag <- "GEO_Peak_Chisquared"  # replace with a unique label for each script

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

#Encode
# Example BED file with insertions
insertions <- read.table("peaks/!Encode_ER_rep2_ChIP_q0.05_peaksNP.bed", header=FALSE)
colnames(insertions) <- c("Chr",	"Start",	"End",	"Name",	"Score",	"Strand",	"Signal Value",	"pvalue",	"qvalue",	"Peak")

# Convert to GRanges
gr <- GRanges(
  seqnames = insertions$Chr,
  ranges = IRanges(start = insertions$Start, end = insertions$End),
  ExperimentInsertions = insertions$`Experiment Insertions`
)

## Make line plot
tiff("Encode_Karyoplot.tiff", width = 2481, height = 1749, res = 300)
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

tiff("Encode_Peak_Chisquared.tiff", width = 2481, height = 1749, res = 300)

barplot(rbind(obs, expected), beside = TRUE,
        col = c("skyblue", "orange"),
        legend.text = c("Observed", "Expected"),
        main = "Encode peakset χ² across all chromosomes",
        ylab = "Number of ChIP-Seq Peaks",
        ylim = c(0, 2500),
        names.arg = names(obs), las = 2, cex.names = 0.8)

dev.off()

#Save data across scripts
# Your tag for this script's data
tag <- "Encode_Peak_Chisquared"  # replace with a unique label for each script

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



# Align expected and test outputs to observed chromosomes
expected <- expected[names(obs)]
test$stdres <- test$stdres[names(obs)]

# Calculate per-chromosome p-values from standardized residuals
pvals_per_chr <- 2 * (1 - pnorm(abs(test$stdres)))
names(pvals_per_chr) <- names(obs)

# Build data frame of per-chromosome results
results_df <- data.frame(
  Chromosome = names(obs),
  Observed = as.numeric(obs),
  Expected = as.numeric(expected),
  StdResidual = test$stdres,
  Pvalue_per_chr = pvals_per_chr
)

# Add significance asterisks column
results_df$Significance <- cut(
  results_df$Pvalue_per_chr,
  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
  labels = c("***", "**", "*", "")
)


# Add overall chi-square p-value summary row with significance
summary_row <- data.frame(
  Chromosome = "Overall_ChiSq",
  Observed = NA,
  Expected = NA,
  StdResidual = NA,
  Pvalue_per_chr = test$p.value,
  Significance = ifelse(test$p.value < 0.001, "***",
                        ifelse(test$p.value < 0.01, "**",
                               ifelse(test$p.value < 0.05, "*", "")))
)

# Combine
results_all <- rbind(results_df, summary_row)

# Use dataset tag if available
output_file <- if (exists("tag")) {
  paste0(tag, "_ChiSquare_Results.xlsx")
} else {
  "ChiSquare_Results.xlsx"
}

write.xlsx(results_all, file = output_file, rowNames = FALSE)

cat("✅ Chi-square test results saved for dataset:", output_file, "\n")


# Set unique tag for this dataset
tag <- "GEO_Peak_Chisquared"  # or "Encode_Peak_Chisquared", "High_Stringency", etc.

# Align expected and residuals to observed chromosomes
expected <- expected[names(obs)]
test$stdres <- test$stdres[names(obs)]

# Per-chromosome p-values
pvals_per_chr <- 2 * (1 - pnorm(abs(test$stdres)))
names(pvals_per_chr) <- names(obs)

# Build dataframe
results_df <- data.frame(
  Chromosome = names(obs),
  Observed = as.numeric(obs),
  Expected = as.numeric(expected),
  StdResidual = test$stdres,
  Pvalue_per_chr = pvals_per_chr
)

# Add significance asterisks column
results_df$Significance <- cut(
  results_df$Pvalue_per_chr,
  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
  labels = c("***", "**", "*", "")
)


# Add overall chi-square p-value summary row with significance
summary_row <- data.frame(
  Chromosome = "Overall_ChiSq",
  Observed = NA,
  Expected = NA,
  StdResidual = NA,
  Pvalue_per_chr = test$p.value,
  Significance = ifelse(test$p.value < 0.001, "***",
                        ifelse(test$p.value < 0.01, "**",
                               ifelse(test$p.value < 0.05, "*", "")))
)

# Combine
results_all <- rbind(results_df, summary_row)

# Export Excel
library(openxlsx)
output_file <- paste0(tag, "_ChiSquare_Results.xlsx")
write.xlsx(results_all, file = output_file, rowNames = FALSE)
cat("✅ Chi-square results saved for dataset:", output_file, "\n")
