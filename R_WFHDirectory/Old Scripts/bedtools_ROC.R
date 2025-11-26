# Install necessary libraries
#install.packages("pROC")
#install.packages("ggplot2")
library(ggplot2)
library(pROC)

# Example data (replace with actual calculated values)
thresholds <- c(1,	3,	5,	6,	7,	8,	9,	10,	11,	12,	13,	0)
sensitivity <- c(0,	0,	0.002287239,	0.005102771,	0.007349081,	0.00911738,	0.010169006,	0.011029938, 0.012082716,	0.012704781,	0.012704781,	0.016979147)  # True Positive Rate
fpr <- c(0,	0,	0.001331811,	0.002849544,	0.003985009,	0.004882211,	0.005259916,	0.005495807,	0.00577843,	0.006013827,	0.006013827,	0.007190161)         # False Positive Rate

data <- data.frame(thresholds, sensitivity, fpr)

# Plot the ROC curve
ggplot(data, aes(x = fpr, y = sensitivity)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(size = 2, color = "red") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  labs(
    title = "ROC Curve",
    x = "False Positive Rate (FPR)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  theme_minimal()

# Calculate AUC using pROC
#roc_obj <- roc(sensitivity = sensitivity, specificity = 1 - fpr)
#auc_value <- auc(roc_obj)
#print(paste("AUC:", auc_value))
