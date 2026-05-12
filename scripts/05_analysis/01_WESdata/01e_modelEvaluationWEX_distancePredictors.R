library(readxl)
source("scripts/05_analysis/00_NamesAndColors.R")
dir.create("fig/modelEvaluation/WES_predictor_distancePredictors", showWarnings = F)
# get predictor order
tab = read_xlsx("data/rawdata/dataMappingAlltissues_distancePredictors.xlsx", 
                sheet="allTissues", col_names=T)
ranges = c("1Mb" = 500000, "100kb" = 50000,"10kb" = 5000, 
           "1kb" = 500, "100bp" = 50, "10bp" = 5, "1bp" = 0)
tab$NA. = NULL
tab[tab == "NA"] = NA
# for predictors where we want multiple ranges, expand table
tab = apply(tab,1,function(x){
  if(is.na(x["range"])){
    return(x)
  } else{
    rangeWindow = strsplit(x["range"],";")[[1]]
    rangeInds = which(names(ranges) == rangeWindow[2]):which(names(ranges) == rangeWindow[1])
    subRanges = ranges[rangeInds]
    t(sapply(names(subRanges), function(r,y){
      y["range"] = subRanges[r]
      if(length(subRanges)>1){
        y["Name"] = paste0(y["Name"]," ",r)
        y["abbreviation"] = paste0(y["abbreviation"],"_",r)
      }
      return(y)
    },y=x))
  }
})
tab = do.call(rbind, tab)
tab = as.data.frame(tab)
predictorOrder = tab[,1:3] # Group, Name, abbreviation
p2P = setNames(predictorOrder$Name, predictorOrder$abbreviation)
predictorGroups = split(predictorOrder$abbreviation, predictorOrder$Group)
p2G = setNames(predictorOrder$Group, predictorOrder$Name)

# get predictor importances #####
list_of_dfs = sapply(tissues, function(tissue){
  load(paste0("data/Modeling/exomeTrainData/RF/", tissue, "_", 
              "finalModel.RData"))
  old_rf = importance[names(p2P)]
  # with distance predictors
  load(paste0("data/Modeling/exomeTrainData_distancePredictors/RF/", tissue, "_", 
              "finalModel.RData"))
  gini_df = data.frame("original rf" = old_rf, 
                       "with distance predictors" = importance[names(p2P)], 
                       row.names = p2P)
}, simplify =F) 
#####

library(ComplexHeatmap)
library(circlize)

# list_of_dfs: named list (names = tissues)
# each df: rows = predictors, cols = c("v1","v2") (or similar)
# p2G: named vector: names = predictors, values = group
SCALE_ORDER <- c('1bp','10bp','100bp','1kb','10kb','100kb','1Mb')

extract_scale <- function(x) {
  for (s in rev(SCALE_ORDER)) {
    if (grepl(s, x, fixed = TRUE)) return(s)
  }
  return(NA)
}
SCALE_COLORS <- c(
  '1bp'   = '#fde725',
  '10bp'  = '#a0da39',
  '100bp' = '#5ec962',
  '1kb'   = '#21918c',
  '10kb'  = '#3b528b',
  '100kb' = '#482878',
  '1Mb'   = '#440154'
)

NO_SCALE_COLOR <- "#d3d3d3"

create_heatmap <- function(list_of_dfs, p2G) {
  
  tissues <- names(list_of_dfs)
  
  # -----------------------------
  # 1. Ensure consistent predictors
  # -----------------------------
  predictors <- Reduce(intersect, lapply(list_of_dfs, rownames))
  
  # reorder everything consistently
  list_of_dfs <- lapply(list_of_dfs, function(df) {
    df[predictors, , drop = FALSE]
  })
  
  # -----------------------------
  # 2. Build combined matrix
  # -----------------------------
  mats <- lapply(tissues, function(t) {
    as.matrix(list_of_dfs[[t]])
  })
  
  combined <- do.call(cbind, mats)
  
  # column labels
  versions <- colnames(list_of_dfs[[1]])
  colnames(combined) <- as.vector(
    outer(tissues, versions, paste, sep = "_")
  )
  
  # -----------------------------
  # 3. Column splitting
  # -----------------------------
  column_split <- rep(tissues, each = length(versions))
  
  # -----------------------------
  # Row grouping
  # -----------------------------
  groups <- p2G[predictors]
  groups[is.na(groups)] <- "Other"
  groups <- factor(groups)
  
  # -----------------------------
  # Resolution annotation
  # -----------------------------
  scales <- sapply(predictors, extract_scale)
  
  # replace NA with "None" for legend
  scales_factor <- factor(
    ifelse(is.na(scales), "None", scales),
    levels = c(SCALE_ORDER, "None")
  )
  
  scale_colors_full <- c(SCALE_COLORS, "None" = NO_SCALE_COLOR)
  
  # -----------------------------
  # Group colors
  # -----------------------------
  group_levels <- levels(groups)
  group_colors <- structure(
    circlize::rand_color(length(group_levels)),
    names = group_levels
  )
  # detect distance predictors
  is_distance <- grepl("distance", predictors)
  
  distance_factor <- factor(
    ifelse(is_distance, "Distance", "Other"),
    levels = c("Distance", "Other")
  )
  
  distance_colors <- c(
    "Distance" = "#000000",   # black highlight
    "Other" = "#f0f0f0"       # light grey
  )
  # -----------------------------
  # Combined row annotation (TWO BARS)
  # -----------------------------
  row_ha <- rowAnnotation(
    Resolution = scales_factor,
    Group = groups,
    Distance = distance_factor,
    
    col = list(
      Resolution = scale_colors_full,
      Group = group_colors,
      Distance = distance_colors
    ),
    
    show_annotation_name = TRUE,
    annotation_name_gp = gpar(fontsize = 10)
  )
  # -----------------------------
  # 5. Scaling (optional but recommended)
  # -----------------------------
  # same logic as your Python script
  combined_scaled <- combined
  sd_val <- sd(combined, na.rm = TRUE)
  if (sd_val > 0) {
    combined_scaled <- combined / sd_val
  }
  
  # -----------------------------
  # 6. Colors
  # -----------------------------
  col_fun <- colorRamp2(
    c(0, max(combined_scaled, na.rm = TRUE)),
    c("#e5e5e5", "#ff0000")
  )
  
  
  # -----------------------------
  # 8. Heatmap
  # -----------------------------
  ht <- Heatmap(
    combined_scaled,
    name = "Gini",
    
    col = col_fun,
    
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    
    show_row_names = FALSE,
    show_column_names = FALSE,
    
    row_split = groups,
    column_split = column_split,
    
    column_title = NULL,
    top_annotation = HeatmapAnnotation(
      Version = rep(versions, length(tissues))),
    left_annotation = row_ha,
    row_title_rot = 0
  )
  
  draw(ht,
       heatmap_legend_side = "right",
       annotation_legend_side = "right")
}
png("fig/comparisonPredictorImportance_distancePredictors.png", width = 1200, height = 1200)
create_heatmap(list_of_dfs, p2G)
dev.off()





# compare performance #####
load("data/Modeling/exomeTrainData/RF/ROC_PR_RF_concat.RData") #ROC_PR_RF_concat
library(ROCR)
predPerTissueRF2 = sapply(tissues, function(tissue){
  load(paste0("data/Modeling/exomeTrainData_distancePredictors/RF/", tissue,
              "_predictions.RData"))
  return(predictions)
}, simplify=F)
ROC_PR_RF_concat_distanceP = sapply(names(predPerTissueRF2), function(tissue){
  predConcat = do.call(rbind,predPerTissueRF2[[tissue]])
  ROC_rf = performance(prediction(predConcat$pred, 
                                  predConcat$label), 
                       "tpr", "fpr")
  AUC_rf = performance(prediction(predConcat$pred, 
                                  predConcat$label), 
                       "auc")
  PR_rf = performance(prediction(predConcat$pred, 
                                 predConcat$label), 
                      "prec", "rec")
  return(list(roc = ROC_rf, pr = PR_rf, auc = AUC_rf))
}, simplify=F)

AUCs = sapply(tissues, function(tissue){
  c(orig = ROC_PR_RF_concat[[tissue]]$auc@y.values[[1]],
    distance = ROC_PR_RF_concat_distanceP[[tissue]]$auc@y.values[[1]])
})
png("fig/comparisonAUC_distancePredictors.png", width = 1300, height = 600, pointsize = 20)
par(mar = c(2,3.5,2,1))
temp = barplot(AUCs, beside = T, legend = F, las =1, names.arg = t2T, ylab = "AUC", ylim = c(0,0.7), mgp = c(2.5,0.8,0))
legend(x=0,y=0.75, legend = c("Original", "With distance predictors"), fill = gray.colors(2), xpd = NA)
dev.off()
#####