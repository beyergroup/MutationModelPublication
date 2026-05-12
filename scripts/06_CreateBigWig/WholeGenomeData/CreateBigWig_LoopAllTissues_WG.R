library(dplyr)
library(data.table)
library(rtracklayer)

# Get tissue name from command-line argument
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript CreateBigWig_PerTissue_WG.R <tissue_name>")
}
tissue <- args[1]
message("Processing tissue: ", tissue)

# Define chromosomes
chromosomes <- paste0("chr", 1:22)

# Paths
input_dir  <- paste0("data/Modeling/WholeGenomeData/RFnew/", tissue, "/")
output_dir <- "data/CreateBigWig/02_WGS/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# hg19 seqinfo for BigWig export
hg19_seqinfo <- Seqinfo(genome = "hg19")

# Temp BedGraph file — one chromosome appended at a time
#   tmp_bg_file <- tempfile(pattern = paste0(tissue, "_"), fileext = ".bedGraph")

positions_total <- 0
tmp_dir <- tempdir()
tmp_bw_files <- character(0)

for (chr in chromosomes) {
  message("  Processing ", chr)
  
  file_path <- paste0(input_dir, tissue, "_MutsResult_", chr,
                      "_all_parts_RFcombined.RData")
  
  if (!file.exists(file_path)) {
    message("    Warning: File not found - ", file_path)
    next
  }
  
  #load(file_path)
  env <- new.env()
  load(file_path, envir = env)
  predictions <- get(ls(env)[1], envir = env)
  rm(env)
  
  bg <- predictions %>%
    select(chr, pos, prediction) %>%
    mutate(start = pos - 1L,
           end   = pos,
           score = prediction) %>%
    select(seqnames = chr, start, end, score) %>%
    arrange(start)
  n_rows <- nrow(bg)
  positions_total <- positions_total + n_rows
  message("    ", chr, ": ", n_rows, " positions")
  gr <- makeGRangesFromDataFrame(bg,
                                 keep.extra.columns = TRUE,
                                 starts.in.df.are.0based = TRUE,
                                 seqinfo = hg19_seqinfo)
  
 # tmp_bw <- paste0(tmp_dir, "/", tissue, "_", chr, "wholeGenome.bw")
  tmp_bw <- paste0(output_dir, tissue, "_", chr, "_wholeGenome.bw")
  export.bw(gr, con = tmp_bw)
  tmp_bw_files <- c(tmp_bw_files, tmp_bw)
  
  rm(predictions, bg, gr)
  gc()
  
  
 # n_rows <- nrow(bg)
 # positions_total <- positions_total + n_rows
  #message("    ", chr, ": ", n_rows, " positions")
  
  #fwrite(bg, file = tmp_bg_file, append = TRUE,
         #sep = "\t", col.names = FALSE, row.names = FALSE)
  
  #rm(predictions, bg)
 # gc()
}

if (length(tmp_bw_files) == 0) {
  message("No data found for tissue: ", tissue)
  quit(status = 1)
}
message("Total positions written: ", positions_total)

# Read combined BedGraph and convert to BigWig
message("Reading combined BedGraph...")


message("Completed tissue: ", tissue)
