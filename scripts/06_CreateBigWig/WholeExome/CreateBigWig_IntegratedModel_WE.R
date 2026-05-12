library(dplyr)
library(data.table)
library(rtracklayer) # For BigWig export

# Define tissue and chromosomes for creating whole exome integrated bw file
tissues <- c("liver")

# Path parameters
input_dir <- "data/Modeling/WholeExomeData/RF/"
output_dir <- "data/CreateBigWig/01_WES/"
hg19_chrom_sizes <- "data/rawdata/hg19.chrom.sizes"

# Function to process one tissue (all parts)
process_tissue <- function(tissue) {
  message("Processing tissue: ", tissue)
  all_parts <- list()
  
  # Process all 50 parts
  for (part in 1:50) {
    file_path <- paste0(input_dir, "exomeMuts_part", part, "_", tissue, "_RFpredictionsExomeGeneric.RData")
    
    if (file.exists(file_path)) {
      load(file_path)
      
      # Format as BedGraph with explicit column names
      bg <- predictions %>%
        select(chr, pos, prediction) %>%
        mutate(start = pos - 1,
               end = pos,
               score = prediction) %>%
        select(seqnames = chr, start, end, score)
      
      all_parts[[part]] <- bg
      message("  Processed part ", part)
    }
  }
  
  # Combine all parts
  combined <- rbindlist(all_parts) %>% 
    arrange(seqnames, start)  # Sort by genomic coordinates
  
  # Convert to GRanges
  gr <- makeGRangesFromDataFrame(combined,
                                 keep.extra.columns = TRUE,
                                 starts.in.df.are.0based = TRUE,
                                 seqinfo = Seqinfo(genome = "hg19"))
  
  # Add score metadata
  mcols(gr)$score <- combined$score
  
  # Save as single BigWig per tissue
  bw_file <- paste0(output_dir, tissue, "_exome_mutation_probIntegratedModel_WE.bw")
  export.bw(gr, con = bw_file)
  message("Saved complete tissue BigWig to: ", bw_file)
}
# Process all tissues
for (tissue in tissues) {
  process_tissue(tissue)
}
