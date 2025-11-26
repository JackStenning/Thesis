#THIS CODE WAS RUN USING R 4.4

#Load Libraries
library(karyoploteR)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)


#Set working directory
setwd("C:/Users/jps558/OneDrive - University of York/Desktop/R_WorkingDirectory/R_WFHDirectory") 

# Example BED file with insertions
insertions <- read.table("peaks/peak_data_ER_WT2.bed", header=FALSE)
colnames(insertions) <- c("Chr",	"Start",	"End",	"Center",	"Experiment Insertions",	"Background insertions",	"Reference Insertions",	"pvalue Reference",	"pvalue Background",	"Fraction Experiment",	"TPH Experiment",	"Fraction background",	"TPH background",	"TPH background subtracted",	"pvalue_adj Reference")

# Convert to GRanges
gr <- GRanges(
  seqnames = insertions$Chr,
  ranges = IRanges(start = insertions$Start, end = insertions$End),
  ExperimentInsertions = insertions$`Experiment Insertions`
)

## Make line plot
tiff("Highstringency_Karyoplot.tiff", width = 2481, height = 1749, res = 300)
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

tiff("Highstringency_Chisquared_peaks.tiff", width = 2481, height = 1749, res = 300)

barplot(rbind(obs, expected), beside = TRUE,
        col = c("skyblue", "orange"),
        legend.text = c("Observed", "Expected"),
        main = "High stringency peakset χ² across all chromosomes",
        ylab = "Number of Calling Card Peaks",
        ylim = c(0, 50),
        names.arg = names(obs), las = 2, cex.names = 0.8)


dev.off()

# Function to save chi-square results with significance
save_chisq_results <- function(obs, expected, test, tag) {
  expected <- expected[names(obs)]
  test$stdres <- test$stdres[names(obs)]
  
  pvals_per_chr <- 2 * (1 - pnorm(abs(test$stdres)))
  names(pvals_per_chr) <- names(obs)
  
  results_df <- data.frame(
    Chromosome = names(obs),
    Observed = as.numeric(obs),
    Expected = as.numeric(expected),
    StdResidual = test$stdres,
    Pvalue_per_chr = pvals_per_chr
  )
  
  results_df$Significance <- cut(
    results_df$Pvalue_per_chr,
    breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
    labels = c("***", "**", "*", "")
  )
  
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
  
  results_all <- rbind(results_df, summary_row)
  
  output_file <- paste0(tag, "_ChiSquare_Results.xlsx")
  write.xlsx(results_all, file = output_file, rowNames = FALSE)
  cat("✅ Chi-square results saved for dataset:", output_file, "\n")
}

# -------------------------------
# Highstringency dataset
insertions_high <- read.table("peaks/peak_data_ER_WT2.bed", header = FALSE)
colnames(insertions_high) <- c("Chr","Start","End","Center","Experiment Insertions","Background insertions",
                               "Reference Insertions","pvalue Reference","pvalue Background",
                               "Fraction Experiment","TPH Experiment","Fraction background",
                               "TPH background","TPH background subtracted","pvalue_adj Reference")

obs_high <- table(insertions_high$Chr)
chrom_sizes <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)
chrom_sizes_filtered <- chrom_sizes[names(obs_high)]
expected_high <- sum(obs_high) * (chrom_sizes_filtered / sum(chrom_sizes_filtered))
test_high <- chisq.test(x = as.numeric(obs_high), p = chrom_sizes_filtered / sum(chrom_sizes_filtered))

save_chisq_results(obs_high, expected_high, test_high, tag = "Highstringency_Peak_Chisquared")

# -------------------------------
# Lowstringency dataset
insertions_low <- read.table("peaks/peak_data_ER_WT4.bed", header = FALSE)
colnames(insertions_low) <- colnames(insertions_high)  # same columns as Highstringency

obs_low <- table(insertions_low$Chr)
chrom_sizes_filtered <- chrom_sizes[names(obs_low)]
expected_low <- sum(obs_low) * (chrom_sizes_filtered / sum(chrom_sizes_filtered))
test_low <- chisq.test(x = as.numeric(obs_low), p = chrom_sizes_filtered / sum(chrom_sizes_filtered))

save_chisq_results(obs_low, expected_low, test_low, tag = "Lowstringency_Peak_Chisquared")