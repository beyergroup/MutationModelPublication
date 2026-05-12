# Load libraries
library(tidyverse)
library(reshape2)

# 1: load your data (replace with your file path)
 load("/cellfile/cellnet/MutationModel/data/Modeling/exomeTrainData/RF/rf_imps_generalModel.RData")

head(rf_imps)

rf_imps <- rf_imps %>%
  mutate(tissue = str_to_title(tissue)) %>%
  filter(tissue != "Generalmodel")


# Pivot to wide format (predictors x tissues)
wide_df <- rf_imps %>%
  select(tissue, predictor, gini_scaled) %>%
  pivot_wider(names_from = tissue, values_from = gini_scaled)

#wide_df <- rf_imps %>%
 # select(tissue, predictor, gini) %>%
  #pivot_wider(names_from = tissue, values_from = gini)

# Compute correlation matrix across tissues
cor_mat <- wide_df %>%
  select(-predictor) %>%
  cor(use = "pairwise.complete.obs", method = "pearson") ##using pearson because Many features have near-zero Gini importance, 
#which creates problems with ranking-based approache (Spearman).

# Convert to long format for ggplot2
cor_df <- melt(cor_mat, varnames = c("Tissue1", "Tissue2"), value.name = "Correlation")

#  Plot
ggplot(cor_df, aes(x = Tissue1, y = Tissue2, fill = Correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "#B2182B", name = "Correlation") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    axis.title = element_blank()
  ) +
  ggtitle("Correlation of Feature Importance Across Tissues")

ggsave(paste0("/cellfile/cellnet/MutationModel/fig/modelEvaluation/20251111_RFginiScaled_featureImportanceCorr.png"),
               height=8, width=8)


library(pheatmap)
library(viridis)

png(file = "/cellfile/cellnet/MutationModel/fig/modelEvaluation/20252511_RFginiScaled_featureImpCorrDendrogram.png",
    width = 3000, height = 3000, res = 350)

pheatmap(cor_mat,
         color = rev(viridis::viridis(100)),
             breaks = seq(0, 1, length.out = 101),
         legend_breaks = seq(0, 1, by = 0.2),
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         # overall font size
fontsize_row = 18,     # row label font size
fontsize_col = 18,
angle_col = 45)


dev.off()