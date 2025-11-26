if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene", force = TRUE)
BiocManager::install("ChIPseeker")
BiocManager::install("EnsDb.Hsapiens.v86")
BiocManager::install("clusterProfiler")
BiocManager::install("AnnotationDbi")
BiocManager::install("org.Hs.eg.db")


# Load libraries
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(EnsDb.Hsapiens.v86)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(dplyr)
library(tidyr)
library(ggplot2)

#Set working directory

setwd("C:/Users/jps558/OneDrive - University of York/Desktop/R_WorkingDirectory/R_WFHDirectory") 


#Load data
samplefiles <- list.files("peaks/", pattern='.bed', full.names = T)
samplefiles <- as.list(samplefiles)
names(samplefiles) <- c("Encode ChIP", "GEO ChIP","High Stringency Calling Cards", "Low Stringency Calling Cards")

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

#Get Annotations
peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb, 
                       tssRegion=c(-1000, 1000), verbose=FALSE)
peakAnnoList

#Plot data barchart
# Generate the annotation bar plot
p <- plotAnnoBar(peakAnnoList, title = "Feature distribution of ER ChIP-seq and HyPB-ESR1 calling cards")

# Modify y-axis to include additional ticks
p + scale_y_continuous(breaks = seq(0, 100, by = 10), labels = function(x) paste0(x, "%"))


#Distance to TSS
plotDistToTSS(peakAnnoList, title = "Distribution of transcription factor-binding loci \n relative to TSS")


#Save percentages
# Extract annotation data and compute percentages
feature_percentages <- lapply(peakAnnoList, function(anno) {
  anno_df <- as.data.frame(anno)  # Convert annotation to dataframe
  feature_counts <- table(anno_df$annotation)  # Count occurrences of each feature
  feature_percent <- prop.table(feature_counts) * 100  # Convert to percentages
  return(feature_percent)
})

# Print feature percentages
print(feature_percentages)

# Convert list to a tidy dataframe
feature_df <- bind_rows(feature_percentages, .id = "Sample") %>%
  pivot_longer(cols = -Sample, names_to = "Feature", values_to = "Percentage")

#Simplify catagories
feature_df_simplified <- feature_df %>%
  mutate(
    Feature = case_when(
      grepl("Promoter", Feature, ignore.case = TRUE) ~ "Promoter",
      grepl("5' UTR", Feature, ignore.case = TRUE) ~ "5' UTR",
      grepl("3' UTR", Feature, ignore.case = TRUE) ~ "3' UTR",
      
      # 1st Exon first
      grepl("Exon\\s*1", Feature, ignore.case = TRUE) ~ "1st Exon",
      # Other Exons: number after exon is NOT 1
      grepl("Exon\\s*[2-9]", Feature, ignore.case = TRUE) ~ "Other Exon",
      
      # 1st Intron first
      grepl("Intron\\s*1", Feature, ignore.case = TRUE) ~ "1st Intron",
      # Other Introns: number after intron is NOT 1
      grepl("Intron\\s*[2-9]", Feature, ignore.case = TRUE) ~ "Other Intron",
      
      grepl("Downstream", Feature, ignore.case = TRUE) ~ "Downstream (<=300bp)",
      grepl("Intergenic", Feature, ignore.case = TRUE) ~ "Distal Intergenic",
      TRUE ~ Feature
    )
  ) %>%
  group_by(Sample, Feature) %>%
  summarise(Percentage = sum(Percentage, na.rm = TRUE), .groups = "drop")
# Print simplified table
print(feature_df_simplified)

graphorder = c("3' UTR", "Promoter", "1st Exon", "Other Exon", "1st Intron", "Other Intron", "Distal Intergenic", "Downstream (<=300bp)", "5' UTR")

write.csv(feature_df_simplified, "Percentage of peak representation in genomic features.csv", row.names = TRUE)

# Define chromosome order
chrom_order <- c(as.character(1:22), "X", "Y")

# Convert Feature into an ordered factor
feature_df_simplified$Feature <- factor(feature_df_simplified$Feature, levels = chrom_order, ordered = TRUE)

#Plot
tiff("Feature Distribution Across Samples.tiff", width = 2481, height = 1749, res = 300)
ggplot(feature_df_simplified, aes(x = factor(Feature, graphorder), y = Percentage, fill = Sample)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Feature Distribution Across Samples", x = "Genomic Feature", y = "Percentage (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_discrete(name = "Feature")
dev.off()