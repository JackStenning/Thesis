if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene", force = TRUE)
BiocManager::install("ChIPseeker", force = TRUE)
BiocManager::install("EnsDb.Hsapiens.v86", force = TRUE)
BiocManager::install("clusterProfiler", force = TRUE)
BiocManager::install("AnnotationDbi", force = TRUE)
BiocManager::install("org.Hs.eg.db", force = TRUE)

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

# Set working directory
setwd("C:/Users/jps558/OneDrive - University of York/Desktop/R_WorkingDirectory/R_WFHDirectory") 

# Load data
samplefiles <- list.files("peaks/", pattern='.bed', full.names = TRUE)
samplefiles <- as.list(samplefiles)
names(samplefiles) <- c("Encode ChIP", "GEO ChIP", "High Stringency Calling Cards", "Low Stringency Calling Cards")

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Get Annotations
peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb = txdb, tssRegion = c(-1000, 1000), verbose = FALSE)

# Plot annotation bar chart
p <- plotAnnoBar(peakAnnoList, title = "Feature distribution of ER ChIP-seq and HyPB-ESR1 calling cards")
p + scale_y_continuous(breaks = seq(0, 100, by = 10), labels = function(x) paste0(x, "%"))

# Distance to TSS
plotDistToTSS(peakAnnoList, title = "Distribution of transcription factor-binding loci \n relative to TSS")

# Extract annotation percentages
feature_percentages <- lapply(peakAnnoList, function(anno) {
  anno_df <- as.data.frame(anno)
  feature_counts <- table(anno_df$annotation)
  prop.table(feature_counts) * 100
})

# Convert list to tidy dataframe
feature_df <- bind_rows(feature_percentages, .id = "Sample") %>%
  pivot_longer(cols = -Sample, names_to = "Feature", values_to = "Percentage") %>%
  mutate(Feature = as.character(Feature))  # ensure character

# Simplify features
feature_df_simplified <- feature_df %>%
  mutate(
    Feature = case_when(
      grepl("Promoter", Feature, ignore.case = TRUE) ~ "Promoter",
      grepl("5'?\\s?UTR", Feature, ignore.case = TRUE) ~ "5' UTR",
      grepl("3'?\\s?UTR", Feature, ignore.case = TRUE) ~ "3' UTR",
      grepl("Exon\\s*1", Feature, ignore.case = TRUE) ~ "1st Exon",
      grepl("Exon\\s*[2-9]", Feature, ignore.case = TRUE) ~ "Other Exon",
      grepl("intron\\s*1", Feature, ignore.case = TRUE) ~ "1st Intron",
      grepl("intron\\s*[2-9]", Feature, ignore.case = TRUE) ~ "Other Intron",
      grepl("Downstream", Feature, ignore.case = TRUE) ~ "Downstream (<=300bp)",
      grepl("Intergenic", Feature, ignore.case = TRUE) ~ "Distal Intergenic",
      TRUE ~ Feature
    )
  ) %>%
  group_by(Sample, Feature) %>%
  summarise(Percentage = sum(Percentage, na.rm = TRUE), .groups = "drop")

# Define plotting order for features
graphorder <- c("3' UTR", "Promoter", "1st Exon", "Other Exon", 
                "1st Intron", "Other Intron", "Distal Intergenic", 
                "Downstream (<=300bp)", "5' UTR")
feature_df_simplified$Feature <- factor(feature_df_simplified$Feature, levels = graphorder, ordered = TRUE)

# Save simplified percentages
write.csv(feature_df_simplified, "Percentage of peak representation in genomic features.csv", row.names = FALSE)

# Plot feature distribution
tiff("Feature Distribution Across Samples.tiff", width = 2481, height = 1749, res = 300)
ggplot(feature_df_simplified, aes(x = Feature, y = Percentage, fill = Sample)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Feature Distribution Across Samples", x = "Genomic Feature", y = "Percentage (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_discrete(name = "Sample")
dev.off()

library(ggplot2)
library(dplyr)

# 1. Assign Group
feature_df_simplified <- feature_df_simplified %>%
  mutate(Group = ifelse(grepl("ChIP", Sample), "ChIP", "CallingCards"))

# 2. Safe statistical test per Feature
safe_test <- function(df) {
  groups <- unique(df$Group)
  if(length(groups) < 2) {
    # Feature present in only 1 group → cannot test
    return(data.frame(
      Feature = unique(df$Feature),
      t_test_p = NA,
      wilcox_p = NA
    ))
  } else {
    return(data.frame(
      Feature = unique(df$Feature),
      t_test_p = t.test(Percentage ~ Group, data = df)$p.value,
      wilcox_p = wilcox.test(Percentage ~ Group, data = df)$p.value
    ))
  }
}

stat_results <- feature_df_simplified %>%
  group_by(Feature) %>%
  summarise(
    t_test_p = tryCatch(t.test(Percentage ~ Group)$p.value, error = function(e) NA),
    wilcox_p = tryCatch(wilcox.test(Percentage ~ Group)$p.value, error = function(e) NA),
    .groups = "drop"
  ) %>%
  mutate(
    t_test_adj = p.adjust(t_test_p, method = "BH"),
    wilcox_adj = p.adjust(wilcox_p, method = "BH"),
    sig = case_when(
      wilcox_adj < 0.001 ~ "***",
      wilcox_adj < 0.01 ~ "**",
      wilcox_adj < 0.05 ~ "*",
      TRUE ~ ""
    )
  )


# 3. Prepare plotting order
graphorder <- c("3' UTR", "Promoter", "1st Exon", "Other Exon", 
                "1st Intron", "Other Intron", "Distal Intergenic", 
                "Downstream (<=300bp)", "5' UTR")
feature_df_simplified$Feature <- factor(feature_df_simplified$Feature, levels = graphorder, ordered = TRUE)
stat_results$Feature <- factor(stat_results$Feature, levels = graphorder, ordered = TRUE)

# 4. Compute y-position for significance stars
y_max <- feature_df_simplified %>%
  group_by(Feature) %>%
  summarise(y_pos = max(Percentage, na.rm = TRUE) * 1.05)  # 5% above max for stars

stat_results <- left_join(stat_results, y_max, by = "Feature")

# 5. Plot with significance stars
tiff("Feature Distribution with Significance.tiff", width = 2481, height = 1749, res = 300)
ggplot(feature_df_simplified, aes(x = Feature, y = Percentage, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  theme_minimal() +
  labs(title = "Feature Distribution Across Samples",
       x = "Genomic Feature",
       y = "Percentage (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("CallingCards" = "#1f77b4", "ChIP" = "#ff7f0e")) +
  # Add significance stars
  geom_text(data = stat_results, aes(x = Feature, y = y_pos, label = sig),
            inherit.aes = FALSE, vjust = 0, size = 6)
dev.off()
