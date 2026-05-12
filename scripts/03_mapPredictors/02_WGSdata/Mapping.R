library(readxl)
library(dplyr)
library(parallel)
library(data.table)
library(readr)
library(vcfR)
library(rtracklayer)

args <- commandArgs(trailingOnly = TRUE)
chr_num <- as.integer(args[1])
part_range <- as.character(args[2])  # "1-3" or "51-100"

chr <- paste0("chr", chr_num)
#tissues <- c("brain","breast","esophagus", "kidney", "liver", "ovary", "prostate", "skin")
tissues <- c("liver")

setwd('/cellfile/cellnet/MutationModel/')
source("lib/dataMapping.R")
#.libPaths(new = "/data/public/cschmalo/R-4.1.2/")
# Create partial directory if it doesn't exist (ADD THIS)
partial_dir <- "/data/MutTables/WholeGenomeData/partial/"
dir.create(partial_dir, showWarnings = FALSE, recursive = TRUE)

ranges <- c("1Mb" = 500000, "100kb" = 50000, "10kb" = 5000, 
            "1kb" = 500, "100bp" = 50, "10bp" = 5, "1bp" = 0)

base_dir <- "data/MutTables/WholeGenomeData/temp2"

cat("Processing parts", part_range, "for chromosome", chr, "\n")

# Load mutation data
muts_file <- paste0("data/MutTables/WholeGenomeData/Muts_WithContext_", chr, ".RData")
load(muts_file)

# Load common metadata 
tab_orig <- read_xlsx("data/rawdata/dataMappingAlltissues_WGS.xlsx", 
                      sheet = "allTissues", col_names = TRUE)
tab_orig$NA. <- NULL
tab_orig[tab_orig == "NA"] <- NA
tab_orig <- as.data.frame(tab_orig)
tab_orig <- tab_orig %>% select(-colon, -lung)

# Function to process one tissue for specific part range
process_tissue_partial <- function(tissue, start_part, end_part) {
  cat("Processing tissue:", tissue, "for parts", start_part, "to", end_part, "\n")
  
  # Select and clean up tissue-specific predictor table
  tab <- tab_orig[, c(colnames(tab_orig)[1:9], tissue)]
  tab <- tab[!is.na(tab[[tissue]]), ]
  
  # Expand multi-range predictors
  tab_exp <- apply(tab, 1, function(x) {
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
  tab_exp <- do.call(rbind, tab_exp)
  tab_exp <- as.data.frame(tab_exp)
  
  # Function to process one part file
  process_part <- function(part_num) {
    part_num_padded <- sprintf("%03d", part_num)  
    part_file <- file.path(base_dir, paste0("genomeMuts_", chr, "_part", part_num_padded, ".bed"))
    
    if (!file.exists(part_file)) {
      cat("Missing file:", part_file, "\n")
      return(NULL)
    }
    
    result <- tryCatch({
      cat("Processing part", part_num, "for", tissue, "\n")
      
      # Read the bed file to get positions for this part
      bed_data <- read.table(part_file, header = FALSE, stringsAsFactors = FALSE)
      colnames(bed_data) <- c("chr", "start", "end", "combined")
      bed_positions <- nrow(bed_data)
      cat("BED file has", bed_positions, "positions\n")
      
      # Map predictors
      pred <- mapPredictors(x = tab_exp, posFile = part_file)
      
      if (is.null(pred) || nrow(pred) == 0) {
        cat("WARNING: Empty result for part", part_num, "in", tissue, "\n")
        return(NULL)
      }
      
      cat("Predictor mapping returned", nrow(pred), "rows\n")
      
      # Check 1: Verify pred and bed file have same number of rows
      if (nrow(pred) != bed_positions) {
        cat("WARNING: Row count mismatch! BED file:", bed_positions, 
            "vs pred:", nrow(pred), "\n")
      } else {
        cat("CHECK PASSED: pred rows match BED file positions\n")
      }
      
      # Filter Muts to only include positions in this part's bed file
      # BED format uses 0-based start, but your Muts uses 1-based positions
      # Typically BED start corresponds to the actual position
      
      # Method 1: Filter by chromosome and position range
      min_pos <- min(bed_data$end)
      max_pos <- max(bed_data$end)
      
      Muts_filtered <- Muts[Muts$chr == chr & 
                              Muts$pos >= min_pos & 
                              Muts$pos <= max_pos, ]
      
      cat("Filtered Muts from", nrow(Muts), "to", nrow(Muts_filtered), "rows\n")
      cat("Position range in BED:", min_pos, "-", max_pos, "\n")
      cat("Position range in filtered Muts:", 
          min(Muts_filtered$pos), "-", max(Muts_filtered$pos), "\n")
      
      # Check 2: Verify consistency between pred and filtered Muts
      if (nrow(pred) != nrow(Muts_filtered)) {
        cat("WARNING: Row count mismatch! pred:", nrow(pred), 
            "vs Muts_filtered:", nrow(Muts_filtered), "\n")
      } else {
        cat("CHECK PASSED: pred and Muts_filtered have consistent rows\n")
      }
      
      # Additional check: Verify positions match exactly
      if (nrow(pred) == nrow(Muts_filtered) && nrow(pred) == bed_positions) {
        # Check if positions align (assuming pred or Muts_filtered have position info)
        if ("pos" %in% colnames(pred)) {
          if (all(pred$pos == Muts_filtered$pos)) {
            cat("CHECK PASSED: Positions match exactly between pred and Muts_filtered\n")
          } else {
            cat("WARNING: Position values don't match between pred and Muts_filtered\n")
          }
        }
      }
      
      # Check 3: Final consistency check
      cat("\n=== FINAL CONSISTENCY CHECK ===\n")
      cat("  BED positions:", bed_positions, "\n")
      cat("  pred rows:", nrow(pred), "\n")
      cat("  Muts_filtered rows:", nrow(Muts_filtered), "\n")
      
      if (bed_positions == nrow(pred) && nrow(pred) == nrow(Muts_filtered)) {
        cat("  ✓ ALL CHECKS PASSED: All datasets have consistent dimensions\n")
      } else {
        cat("  ✗ WARNING: Dimension mismatch detected!\n")
      }
      cat("================================\n\n")
      
      # Save individual part with filtered Muts
      output_file <- paste0("data/MutTables/WholeGenomeData/partial/",
                            tissue, "_MutsResult_chr", chr_num, "_part", part_num_padded, ".RData")
      
      partial_data <- list(meta = tab_exp, pred = pred, muts = Muts_filtered)
      save(partial_data, file = output_file)
      
      cat("Saved part", part_num, "for", tissue, "with", nrow(Muts_filtered), "mutations\n")
      cat(paste(rep("=", 60), collapse = ""), "\n\n")
      
      return(TRUE)
      
    }, error = function(e) {
      cat("ERROR in part", part_num, "for", tissue, ":", conditionMessage(e), "\n")
      traceback()
      return(NULL)
    })
    
    return(result)
  }
  
  # Process specific part range
  part_numbers <- start_part:end_part
  cat("Processing parts", start_part, "to", end_part, "for", tissue, "\n")
  
  # Process each part (saves individually inside process_part)
  results <- mclapply(part_numbers, process_part, mc.cores = 28)
  
  # Count successful saves
  successful <- sum(sapply(results, function(x) !is.null(x)))
  cat("\n=== SUMMARY ===\n")
  cat("Successfully saved", successful, "out of", length(part_numbers), "parts for", tissue, "\n")
  cat("===============\n\n")
  
  return(TRUE)
}

# Determine part range based on argument
if (part_range == "41-60") {
  start_part <- 41
  end_part <- 60
} else if (part_range == "51-100") {
  start_part <- 51
  end_part <- 100
} else {
  stop("Invalid part range. Use '41-60' or '51-100'")
}

# Process all tissues for the specified part range
for (tissue in tissues) {
  process_tissue_partial(tissue, start_part, end_part)
}

cat("Finished processing parts", part_range, "for chromosome", chr, "\n")