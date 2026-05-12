# =============================================================================
# this script Exports rf_gini data to CSV for Python plotting
#step 1: run this script
# step 2: python3 01b_modelEvaluationWEX_RF_chrWise_TissueWiseGI_forRebuttal.py --data-dir fig/modelEvaluation/WES_predictor_cors_csv/ --output fig/rf_gini_heatmap.png
# =============================================================================

source("scripts/05_analysis/00_NamesAndColors.R")

tissues <- c("brain", "breast", "colon", "esophagus", "kidney", 
             "liver", "lung", "ovary", "prostate", "skin")

# Create output directory
#dir.create("data/rf_gini_csv", showWarnings = FALSE)

# Export each tissue's rf_gini matrix
for (tissue in tissues) {
  load(paste0("data/Modeling/exomeTrainData/RF/", tissue, "_finalModel.RData"))
  # load("data/Modeling/exomeTrainData/RF/RF_imps.RData") --> should it actually be this?
  # rf_gini is a list, get the tissue's data
  if (exists("rf_gini") && tissue %in% names(rf_gini)) {
    gini_df <- rf_gini[[tissue]]
    
    # Apply pretty names to rownames
    if (exists("p2P")) {
      rownames(gini_df) <- p2P[rownames(gini_df)]
    }
    
    write.csv(gini_df, paste0("fig/modelEvaluation/WES_predictor_cors_csv/rf_gini_", tissue, ".csv"))
    message(paste("Exported:", tissue))
  }
}

# Also create a combined long-format CSV
combined_df <- do.call(rbind, lapply(tissues, function(tissue) {
  load(paste0("data/Modeling/exomeTrainData/RF/", tissue, "_finalModel.RData"))
  
  if (exists("rf_gini") && tissue %in% names(rf_gini)) {
    gini_df <- rf_gini[[tissue]]
    
    # Apply pretty names
    if (exists("p2P")) {
      rownames(gini_df) <- p2P[rownames(gini_df)]
    }
    
    # Convert to long format
    gini_long <- reshape2::melt(as.matrix(gini_df))
    colnames(gini_long) <- c("feature", "chromosome", "gini")
    gini_long$tissue <- tissue
    
    return(gini_long)
  }
  return(NULL)
}))

write.csv(combined_df, "fig/modelEvaluation/WES_predictor_cors_csv/rf_gini_combined.csv", row.names = FALSE)
message("Exported combined CSV")