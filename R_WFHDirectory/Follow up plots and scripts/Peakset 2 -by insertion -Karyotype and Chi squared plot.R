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

# If your column name has spaces, convert it to a safe name
colnames(insertions) <- make.names(colnames(insertions))

# Sum experimental insertions per chromosome (observed total insertions per chr)
obs_insertions <- tapply(insertions$Experiment.Insertions, insertions$Chr, sum)
obs_insertions <- as.numeric(obs_insertions)
names(obs_insertions) <- names(tapply(insertions$Experiment.Insertions, insertions$Chr, sum))

# Get chromosome sizes from BSgenome (you already have chrom_sizes_filtered)
chrom_sizes <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)
chroms <- names(obs_insertions)
chrom_sizes_filtered <- chrom_sizes[chroms]

# Expected counts proportional to chromosome length
total_insertions <- sum(obs_insertions)
expected_probs <- chrom_sizes_filtered / sum(chrom_sizes_filtered)
expected_counts <- total_insertions * expected_probs

# Check Chi-square assumptions: expected_counts should generally be >= 5
data.frame(chr = names(obs_insertions),
           observed = obs_insertions,
           expected = expected_counts)

# Run chi-squared test (use the numeric vectors in the same chromosome order)
chisq_result <- chisq.test(x = as.numeric(obs_insertions), p = expected_probs)
chisq_result


# Create desired order and subset to chromosomes present
desired_order <- c(paste0("chr", 1:22), "chrX", "chrY")
ord <- match(desired_order, names(obs_insertions))
ord <- ord[!is.na(ord)]

obs_insertions_ord <- obs_insertions[ord]
expected_counts_ord <- expected_counts[ord]

tiff("Highstringency_Chisquared_insertions.tiff", width = 3510, height = 2481, res = 300)

barplot(rbind(obs_insertions_ord, expected_counts_ord), beside = TRUE,
        col = c("skyblue", "orange"),
        legend.text = c("Observed (insertion counts)", "Expected (length-weighted)"),
        main = "High stringency peakset χ² across all chromosomes",
        names.arg = names(obs_insertions_ord), las = 2, cex.names = 0.8,
        ylab = "Number of insertions",
        ylim = c(0, 2000))

dev.off()

#Save data across scripts
# Your tag for this script's data
tag <- "Highstringency_insertions_Chisquared"  # replace with a unique label for each script

# Create this script's dataframe (replace obs/expected/chromosomes with your actual vectors)
df_new <- data.frame(
  chromosome = names(obs_insertions_ord),
  obs = as.numeric(obs_insertions_ord),
  expected = as.numeric(expected_counts_ord),
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