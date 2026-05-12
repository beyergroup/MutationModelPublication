library(data.table)

# --- Parse command-line arguments ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript script.R <chr_number> <tissue>")
}
i      <- as.integer(args[1])
tissue <- args[2]

cat("=== Starting: chromosome", i, "| tissue:", tissue, "===\n")

setwd('/cellfile/cellnet/MutationModel/')

# Load full chromosome Muts
cat("Loading Muts for chr", i, "\n")
load(paste0("data/MutTables/WholeGenomeData/Muts_WithContext_chr", i, ".RData"))

# Save memory
Muts$alt     <- NULL
Muts$context <- NULL
Muts$mutated <- NULL
gc()

cat("Total Muts positions:", nrow(Muts), "\n")

# Initialize combined predictions
all_predictions <- data.frame()

# Loop through all 100 parts
for (part in 1:100) {
  part_padded <- sprintf("%03d", part)
  
  pred_file <- paste0("data/Modeling/WholeGenomeData/RFnew/", tissue, "/",
                      tissue, "_MutsResult_chr", i,
                      "_part", part_padded, "_RFpredictions.RData")
  
 
  
  
  if (file.exists(pred_file)) {
    load(pred_file)  # loads 'predictions'
    all_predictions <- rbind(all_predictions, predictions)
    cat("Loaded part", part, "- Total rows so far:", nrow(all_predictions), "\n")
  } else {
    cat("WARNING: Missing file for part", part, "\n")
  }
}

cat("\n=== SUMMARY ===\n")
cat("Total Muts positions:", nrow(Muts), "\n")
cat("Total predictions:   ", nrow(all_predictions), "\n")
cat("Difference:          ", nrow(Muts) - nrow(all_predictions), "\n")

# Check for missing positions
cat("\nChecking for missing positions...\n")
missing         <- !(Muts$pos %in% all_predictions$pos)
missing_indices <- which(missing)
cat("Number of missing positions:", length(missing_indices), "\n")

if (length(missing_indices) > 0) {
  cat("\nFirst 10 missing positions:\n")
  print(Muts[missing_indices[1:min(10, length(missing_indices))], ])
  
  cat("\nMissing position ranges:\n")
  cat("Min missing pos:", min(Muts$pos[missing_indices]), "\n")
  cat("Max missing pos:", max(Muts$pos[missing_indices]), "\n")
} else {
  cat("No missing positions!\n")
}

# Check for duplicates
duplicates <- duplicated(all_predictions$pos)
if (any(duplicates)) {
  cat("\nWARNING:", sum(duplicates), "duplicate positions found!\n")
  cat("First 10 duplicates:\n")
  print(all_predictions[duplicates, ][1:10, ])
}

# Save combined predictions
out_file <- paste0("data/Modeling/WholeGenomeData/RFnew/", tissue, "/",
                   tissue, "_MutsResult_chr", i, "_all_parts_RFcombined.RData")


save(all_predictions, file = out_file)
cat("\nSaved combined predictions to:", out_file, "\n")
cat("=== Done: chromosome", i, "| tissue:", tissue, "===\n\n")