#Load data

Bulk_Seq <- data.frame(
  Sample = c("HyPB-ESR1 C Replicate 1",
             "HyPB-ESR1 C Replicate 2",
             "HyPB-ESR1 C Replicate 3",
             "HyPB-ESR1 C Replicate 4",
             "HyPB-ESR1 C Replicate 5",
             "HyPB-ESR1 C Replicate 6",
             "WT HyPB Replicate 1",
             "WT HyPB Replicate 2",
             "WT HyPB Replicate 3",
             "WT HyPB Replicate 4",
             "WT HyPB Replicate 5",
             "WT HyPB Replicate 6",
             "ESR1-HyPB N Replicate 1",
             "ESR1-HyPB N Replicate 2",
             "ESR1-HyPB N Replicate 3",
             "ESR1-HyPB N Replicate 4",
             "ESR1-HyPB N Replicate 5",
             "ESR1-HyPB N Replicate 6"),
  Group = c("HyPB-ESR1 C-terminal fusion",
            "HyPB-ESR1 C-terminal fusion",
            "HyPB-ESR1 C-terminal fusion",
            "HyPB-ESR1 C-terminal fusion",
            "HyPB-ESR1 C-terminal fusion",
            "HyPB-ESR1 C-terminal fusion",
            "WT HyPB",
            "WT HyPB",
            "WT HyPB",
            "WT HyPB",
            "WT HyPB",
            "WT HyPB",
            "ESR1-HyPB N-terminal fusion",
            "ESR1-HyPB N-terminal fusion",
            "ESR1-HyPB N-terminal fusion",
            "ESR1-HyPB N-terminal fusion",
            "ESR1-HyPB N-terminal fusion",
            "ESR1-HyPB N-terminal fusion"),
  Concentration = c(190,
                    108,
                    172,
                    8,
                    92,
                    78,
                    136,
                    101,
                    163,
                    159,
                    248,
                    154,
                    5,
                    7,
                    3,
                    5,
                    15,
                    5)
)

#Make plot
boxplot(Concentration ~ Group, data=Bulk_Seq, 
        main="RNA Concentration by Group",
        ylab="Concentration (ng/µL)",
        col=c("lightblue", "pink"),
        notch=FALSE)

#Plot using ggplot
library(ggplot2)

library(dplyr)
# Compute box plot stats
#box_stats <- boxplot.stats(Bulk_Seq$Concentration)
#outliers <- Bulk_Seq$Concentration %in% box_stats$out  # Identify outliers

# Create ggplot with labeled outliers

# Original plot
tiff("RNA Concentration by Sample Type.tiff", width = 2481, height = 1749, res = 300)
ggplot(Bulk_Seq, aes(x = reorder(Group, Concentration, FUN = median), y = Concentration)) +
  geom_boxplot(fill = "lightblue", color = "darkblue", outlier.shape = NA) +  
  stat_summary(fun.min = min, fun.max = max, geom = "errorbar", width = 0.3, color = "black") +
  geom_jitter(width = 0.2, alpha = 0.5, color = "black") +
  labs(title = "RNA Concentration by Sample Type", 
       x = "Sample Type", 
       y = "Concentration (ng/µL)") +
  theme_minimal(base_size = 14) +  # increase all text by 2 points
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# Filtered plot (without replicate 4)
Bulk_Seq_filtered <- Bulk_Seq %>%
  filter(Sample != "HyPB-ESR1 C Replicate 4")  # Replace with actual sample name

tiff("RNA Concentration by Sample Type W_o 4.tiff", width = 2481, height = 1749, res = 300)
ggplot(Bulk_Seq_filtered, aes(x = reorder(Group, Concentration, FUN = median), y = Concentration)) +
  geom_boxplot(fill = "lightblue", color = "darkblue", outlier.shape = NA) +  
  stat_summary(fun.min = min, fun.max = max, geom = "errorbar", width = 0.3, color = "black") +
  geom_jitter(width = 0.2, alpha = 0.5, color = "black") +
  labs(title = "RNA Concentration by Sample Type (Without HyPB-ESR1 replicate 4)", 
       x = "Sample Type", 
       y = "Concentration (ng/µL)") +
  theme_minimal(base_size = 14) +  # increase all text by 2 points
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()