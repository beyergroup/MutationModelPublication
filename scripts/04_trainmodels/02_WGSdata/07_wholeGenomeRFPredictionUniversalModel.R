library(ranger)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
part_num <- as.integer(args[1])
part_num_padded <- sprintf("%03d", part_num)

tissues <- c("liver")

#tissues <- c("brain","breast","esophagus", "kidney", "liver", "ovary", "prostate", "skin")
nThreads <- 56

for (tissue in tissues) {
  model_file <- paste0("data/Modeling/WholeGenomeData/RF/TissueCombination_finalModel.RData")
  if (!file.exists(model_file)) {
    cat("SKIP: Model file not found:", model_file, "\n")
    next
  }
  cat("Loading model from:", model_file, "\n")
  load(model_file)
  
  for (chr_num in 1:22) {
    cat("\n==== Processing chromosome", chr_num, "part", part_num, "for tissue:", tissue, "====\n")
    
    input_file <- paste0("data/MutTables/WholeGenomeData/partial/", tissue, "/",
                         tissue, "_MutsResult_chr", chr_num, "_part", part_num_padded, ".RData")
    output_file <- paste0("data/Modeling/WholeGenomeData/RFnew/", tissue, "/",
                          tissue, "_MutsResult_chr", chr_num, "_part", part_num_padded, "_RFpredictionsUniversalModel.RData")
    
    if (!file.exists(input_file)) {
      cat("SKIP: Input file not found:", input_file, "\n")
      next
    }
    if (file.exists(output_file)) {
      cat("SKIP: Output already exists:", output_file, "\n")
      next
    }
    
    cat("Loading data from:", input_file, "\n")
    load(input_file)
    
    cat("Making predictions...\n")
    testDat <- partial_data$pred
    yhat <- predict(rf, testDat, num.threads = nThreads)
    predictions <- data.frame(partial_data$muts, prediction = yhat$predictions[, 2])
    
    save(predictions, file = output_file)
    cat("SUCCESS: Processed chromosome", chr_num, "part", part_num, "for tissue:", tissue, "\n")
  }
}
