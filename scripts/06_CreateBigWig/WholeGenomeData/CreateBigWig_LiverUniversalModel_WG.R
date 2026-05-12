library(dplyr)
library(data.table)
library(rtracklayer)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript CreateBigWig_PerTissue_WG.R <tissue> <chr>")
}
tissue <- args[1]
chr    <- args[2]   # e.g. "chr18"
message("Processing tissue: ", tissue, ", chromosome: ", chr)

# Define chromosomes
chromosomes <- paste0("chr", 1:22)

# Paths
input_dir  <- paste0("data/Modeling/WholeGenomeData/RFnew/", tissue, "/")
output_dir <- "data/CreateBigWig/02_WGS/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# hg19 seqinfo for BigWig export
hg19_seqinfo <- Seqinfo(genome = "hg19")

# chrom.sizes for bedGraphToBigWig
chrom_sizes_file <- tempfile(pattern = paste0(tissue, "_"), fileext = ".chrom.sizes")
data.table(chr = seqnames(hg19_seqinfo),
           len = seqlengths(hg19_seqinfo))[chr %in% chromosomes] |>
  fwrite(chrom_sizes_file, sep = "\t", col.names = FALSE)


# Temp BedGraph file — one chromosome appended at a time
#   tmp_bg_file <- tempfile(pattern = paste0(tissue, "_"), fileext = ".bedGraph")

positions_total <- 0
tmp_dir <- tempdir()
tmp_bw_files <- character(0)

#for (chr in chromosomes) {
  message("  Processing ", chr)
  
  parts <- sort(list.files(
    input_dir,
    pattern = paste0("^", tissue, "_MutsResult_", chr,
                     "_part\\d{3}_RFpredictionsUniversalModel\\.RData$"),
    full.names = TRUE
  ))
  
  if (length(parts) == 0) {
    message("    Warning: No parts found for ", chr)
    next
  }
  
  tmp_bg <- tempfile(pattern = paste0(tissue, "_", chr, "_"), fileext = ".bedGraph")
  chr_positions <- 0L
  
  for (p in parts) {
    env <- new.env()
    load(p, envir = env)
    predictions <- env$predictions
    rm(env)
    
    bg <- data.table(
      seqnames = predictions$chr,
      start    = predictions$pos - 1L,
      end      = predictions$pos,
      score    = predictions$prediction
    )
    
    fwrite(bg, file = tmp_bg, append = TRUE,
           sep = "\t", col.names = FALSE, row.names = FALSE)
    chr_positions <- chr_positions + nrow(bg)
    
    rm(predictions, bg); gc()
  }
  
  positions_total <- positions_total + chr_positions
  message("    ", chr, ": ", chr_positions, " positions across ",
          length(parts), " parts")
  
  out_bw <- paste0(output_dir, "TissueIntegrated", "_", chr, "_wholeGenome.bw")
  status <- system2("bedGraphToBigWig",
                    args = c(shQuote(tmp_bg),
                             shQuote(chrom_sizes_file),
                             shQuote(out_bw)))
  if (status != 0) stop("bedGraphToBigWig failed for ", chr)
  
  unlink(tmp_bg)
  tmp_bw_files <- c(tmp_bw_files, out_bw)
#}

if (length(tmp_bw_files) == 0) {
  message("No data found for tissue: ", tissue)
  quit(status = 1)
}
message("Total positions written: ", positions_total)

# Read combined BedGraph and convert to BigWig
message("Reading combined BedGraph...")


message("Completed tissue: ", tissue)
