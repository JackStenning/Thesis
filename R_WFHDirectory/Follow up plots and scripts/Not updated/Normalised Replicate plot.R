# Load necessary library
library(ggplot2)

#Set working directory
setwd("//userfs/jps558/w2k/Desktop/R_WorkingDirectory") 


#Set column names
qbedcol_names = c("Chr",	"Start",	"End",	"Number of Insertions",	"Strand",	"Barcode")

pyCCcol_names = c("Chr",	"Start",	"End",	"Center",	"Experiment Insertions",	"Background insertions",	"Reference Insertions",	"pvalue Reference",	"pvalue Background",	"Fraction Experiment",	"TPH Experiment",	"Fraction background",	"TPH background",	"TPH background subtracted",	"pvalue_adj Reference")


#Load data - replicate .qbed files changed to .bed

ER_Rep1 = read.table("Normalising replicates/filtered_sortedER1_EKDL230008858.bed", header = TRUE, sep = '\t', col.names = qbedcol_names)
ER_Rep2 = read.table("Normalising replicates/filtered_sortedER2_EKDL230008858.bed", header = TRUE, sep = '\t', col.names = qbedcol_names)
ER_Rep3 = read.table("Normalising replicates/filtered_sortedER3_EKDL230008858.bed", header = TRUE, sep = '\t', col.names = qbedcol_names)
ER_Rep5 = read.table("Normalising replicates/filtered_sortedER5_EKDL230008858.bed", header = TRUE, sep = '\t', col.names = qbedcol_names)
ER_Rep6 = read.table("Normalising replicates/filtered_sortedER6_EKDL230008858.bed", header = TRUE, sep = '\t', col.names = qbedcol_names)

ER_Rep1_3 = read.table("Normalising replicates/filtered_sortedER1-3_EKDL230008858.bed", header = TRUE, sep = '\t', col.names = qbedcol_names)
ER_Rep5_6 = read.table("Normalising replicates/filtered_sortedER5-6_EKDL230008858.bed", header = TRUE, sep = '\t', col.names = qbedcol_names)


#Total insertions
ER1_total = sum(ER_Rep1$Number.of.Insertions)
ER2_total = sum(ER_Rep2$Number.of.Insertions)
ER3_total = sum(ER_Rep3$Number.of.Insertions)
ER5_total = sum(ER_Rep5$Number.of.Insertions)
ER6_total = sum(ER_Rep6$Number.of.Insertions)

ER_Rep1_3_total = sum(ER_Rep1_3$Number.of.Insertions) 
ER_Rep5_6_total = sum(ER_Rep5_6$Number.of.Insertions) 


#Calculate IPM
ER_Rep1$IPM = (ER_Rep1$Number.of.Insertions/ER1_total)*1e6
ER_Rep2$IPM = (ER_Rep2$Number.of.Insertions/ER2_total)*1e6
ER_Rep3$IPM = (ER_Rep3$Number.of.Insertions/ER3_total)*1e6
ER_Rep5$IPM = (ER_Rep5$Number.of.Insertions/ER5_total)*1e6
ER_Rep6$IPM = (ER_Rep6$Number.of.Insertions/ER6_total)*1e6


ER_Rep1_3$IPM = (ER_Rep1_3$Number.of.Insertions/ER1_total)*1e6
ER_Rep5_6$IPM = (ER_Rep5_6$Number.of.Insertions/ER1_total)*1e6

ER_peak2$IPM = (ER_peak2$Experiment.Insertions/ER_peak2_total)*1e6
ER_peak4$IPM = (ER_peak4$Experiment.Insertions/ER_peak4_total)*1e6

#Convert to log10
ER_Rep1$log10_IPM <- log10(ER_Rep1$IPM)
ER_Rep2$log10_IPM <- log10(ER_Rep2$IPM)
ER_Rep3$log10_IPM <- log10(ER_Rep3$IPM)
ER_Rep5$log10_IPM <- log10(ER_Rep5$IPM)
ER_Rep6$log10_IPM <- log10(ER_Rep6$IPM)

ER_Rep1_3$log10_IPM = log10(ER_Rep1_3$IPM)
ER_Rep5_6$log10_IPM = log10(ER_Rep5_6$IPM)


ER_peak2$log10_IPM = log10(ER_peak2$IPM)
ER_peak4$log10_IPM = log10(ER_peak4$IPM)

#### Correlation of IPM Rep1 and 2 to the number of peaks in peakset 2
#Subset data
# Create a new dataframe with the first four columns
ER_peak2_subset <- ER_peak2[, c(1, 2, 3, 5)]
ER_Rep1_subset = ER_Rep1[, 1:4]
ER_Rep2_subset = ER_Rep2[, 1:4]
ER_Rep3_subset = ER_Rep3[, 1:4]
ER_Rep5_subset = ER_Rep5[, 1:4]
ER_Rep6_subset = ER_Rep6[, 1:4]

# View the first few rows of the new dataframe
head(ER_peak2_subset)

#Align data
# Merge the first two dataframes
aligned_data_temp <- merge(ER_peak2_subset, ER_Rep1_subset, by = "Start", suffixes = c("_Peak2", "_Rep1"))

# Merge the result with the third dataframe
aligned_dataPeak2_12 <- merge(aligned_data_temp, ER_Rep2_subset, by = "Start", suffixes = c("", "_Rep2"))

# Check the dimensions of the merged data
print(dim(aligned_dataPeak2_12))

















##############BELOW IS REPRODUCIBILITY OF ONLY REPLICATES



#### Correlation of IPM Rep1 and 2
#Align data

aligned_data1_2 <- merge(ER_Rep1, ER_Rep2, by = "Start", suffixes = c("_Rep1", "_Rep2"))

# Check the dimensions of the merged data
print(dim(aligned_data1_2))

#Calculate correlation and plot

correlation1_2 <- cor(aligned_data1_2$log10_IPM_Rep1, aligned_data1_2$log10_IPM_Rep2, method = "pearson")

ggplot(aligned_data1_2, aes(x = log10_IPM_Rep1, y = log10_IPM_Rep2)) +
  geom_point(color = "blue", size = 3) +  # Scatter plot points
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +  # Line of unity
  annotate("text", x = max(aligned_data1_2$log10_IPM_Rep1) * 0.05, 
           y = max(aligned_data1_2$log10_IPM_Rep2) * 0.9, 
           label = paste("Pearson r =", round(correlation1_2, 2)), 
           color = "black", size = 5) +  # Add correlation text
  labs(
    title = "Reproducibility of Insertions from Replicates 1 and 2 (log10_IPM)",
    x = "Replicate 1 (log10 IPM)",
    y = "Replicate 2 (log10 IPM)"
  ) +
  theme_minimal()  # Minimal theme

#### Correlation of IPM Rep1 and 6
#Align data

aligned_data1_6 <- merge(ER_Rep1, ER_Rep6, by = "Start", suffixes = c("_Rep1", "_Rep6"))

# Check the dimensions of the merged data
print(dim(aligned_data1_6))

#Calculate correlation and plot

correlation1_6 <- cor(aligned_data1_6$log10_IPM_Rep1, aligned_data1_6$log10_IPM_Rep6, method = "pearson")

ggplot(aligned_data1_6, aes(x = log10_IPM_Rep1, y = log10_IPM_Rep6)) +
  geom_point(color = "blue", size = 3) +  # Scatter plot points
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +  # Line of unity
  annotate("text", x = max(aligned_data1_6$log10_IPM_Rep1) * 0.05, 
           y = max(aligned_data1_6$log10_IPM_Rep6) * 0.9, 
           label = paste("Pearson r =", round(correlation1_6, 2)), 
           color = "black", size = 5) +  # Add correlation text
  labs(
    title = "Reproducibility of Insertions from Replicates 1 and 6 (log10_IPM)",
    x = "Replicate 1 (log10 IPM)",
    y = "Replicate 6 (log10 IPM)"
  ) +
  theme_minimal()  # Minimal theme


#### Correlation of IPM Rep2 and 6
#Align data

aligned_data2_6 <- merge(ER_Rep2, ER_Rep6, by = "Start", suffixes = c("_Rep2", "_Rep6"))

# Check the dimensions of the merged data
print(dim(aligned_data2_6))

#Calculate correlation and plot

correlation2_6 <- cor(aligned_data2_6$log10_IPM_Rep2, aligned_data2_6$log10_IPM_Rep6, method = "pearson")

ggplot(aligned_data2_6, aes(x = log10_IPM_Rep2, y = log10_IPM_Rep6)) +
  geom_point(color = "blue", size = 3) +  # Scatter plot points
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +  # Line of unity
  annotate("text", x = max(aligned_data2_6$log10_IPM_Rep2) * 0.05, 
           y = max(aligned_data2_6$log10_IPM_Rep6) * 0.9, 
           label = paste("Pearson r =", round(correlation2_6, 2)), 
           color = "black", size = 5) +  # Add correlation text
  labs(
    title = "Reproducibility of Insertions from Replicates 2 and 6 (log10_IPM)",
    x = "Replicate 2 (log10 IPM)",
    y = "Replicate 6 (log10 IPM)"
  ) +
  theme_minimal()  # Minimal theme

#### Correlation of IPM Rep2 and 5
#Align data

aligned_data2_5 <- merge(ER_Rep2, ER_Rep5, by = "Start", suffixes = c("_Rep2", "_Rep5"))

# Check the dimensions of the merged data
print(dim(aligned_data2_5))

#Calculate correlation and plot

correlation2_5 <- cor(aligned_data2_5$log10_IPM_Rep2, aligned_data2_5$log10_IPM_Rep5, method = "pearson")

ggplot(aligned_data2_5, aes(x = log10_IPM_Rep2, y = log10_IPM_Rep5)) +
  geom_point(color = "blue", size = 3) +  # Scatter plot points
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +  # Line of unity
  annotate("text", x = max(aligned_data2_5$log10_IPM_Rep2) * 0.05, 
           y = max(aligned_data2_5$log10_IPM_Rep5) * 0.9, 
           label = paste("Pearson r =", round(correlation2_5, 2)), 
           color = "black", size = 5) +  # Add correlation text
  labs(
    title = "Reproducibility of Insertions from Replicates 2 and 5 (log10_IPM)",
    x = "Replicate 2 (log10 IPM)",
    y = "Replicate 5 (log10 IPM)"
  ) +
  theme_minimal()  # Minimal theme


#### Correlation of IPM Rep3 and 5
#Align data

aligned_data3_5 <- merge(ER_Rep3, ER_Rep5, by = "Start", suffixes = c("_Rep3", "_Rep5"))

# Check the dimensions of the merged data
print(dim(aligned_data3_5))

#Calculate correlation and plot

correlation3_5 <- cor(aligned_data3_5$log10_IPM_Rep3, aligned_data3_5$log10_IPM_Rep5, method = "pearson")

ggplot(aligned_data3_5, aes(x = log10_IPM_Rep3, y = log10_IPM_Rep5)) +
  geom_point(color = "blue", size = 3) +  # Scatter plot points
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +  # Line of unity
    labs(
    title = "Reproducibility of Insertions from Replicates 3 and 5 (log10_IPM)",
    x = "Replicate 3 (log10 IPM)",
    y = "Replicate 5 (log10 IPM)"
  ) +
  theme_minimal()  # Minimal theme

#### Correlation of IPM Rep1-3 and 5-6
#Align data

aligned_dataCombi <- merge(ER_Rep1_3, ER_Rep5_6, by = "Start", suffixes = c("_Rep1_3", "_Rep5_6"))

# Check the dimensions of the merged data
print(dim(aligned_dataCombi))



#Calculate correlation and plot

correlationCombi <- cor(aligned_dataCombi$log10_IPM_Rep1_3, aligned_dataCombi$log10_IPM_Rep5_6, method = "pearson")

ggplot(aligned_dataCombi, aes(x = log10_IPM_Rep1_3, y = log10_IPM_Rep5_6)) +
  geom_point(color = "blue", size = 3) +  # Scatter plot points
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +  # Line of unity
  annotate("text", x = max(aligned_dataCombi$log10_IPM_Rep1_3) * 0.05, 
           y = max(aligned_dataCombi$log10_IPM_Rep5_6) * 0.9, 
           label = paste("Pearson r =", round(correlationCombi, 2)), 
           color = "black", size = 5) +  # Add correlation text
  labs(
    title = "Reproducibility of Insertions from grouped Replicates 1-3 and 5-6 (log10_IPM)",
    x = "Replicate 3 (log10 IPM)",
    y = "Replicate 5 (log10 IPM)"
  ) +
  theme_minimal()  # Minimal theme
