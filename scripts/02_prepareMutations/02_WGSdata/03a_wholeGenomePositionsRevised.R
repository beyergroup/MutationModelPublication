#!/usr/bin/env Rscript
library(GenomicRanges)
library(data.table)
options(scipen = 999)# Faster I/O

# Get chromosome from command line
args <- commandArgs(trailingOnly = TRUE)
task_id <- as.integer(args[1])  # SLURM_ARRAY_TASK_ID (1-22)
cr <- paste0("chr", task_id)    # Maps to chr1-chr22

# Load non-excludable regions (only once)
load("/cellfile/datapublic/ypaul1/Mutations/results/GenomePredictions/autosomes_non_excludable_regions.RData")

# Chromosome lengths (hg19)
chr_lengths <- c(
  chr1 = 249250621, chr2 = 243199373, chr3 = 198022430,
  chr4 = 191154276, chr5 = 180915260, chr6 = 171115067,
  chr7 = 159138663, chr8 = 146364022, chr9 = 141213431,
  chr10 = 135534747, chr11 = 135006516, chr12 = 133851895,
  chr13 = 115169878, chr14 = 107349540, chr15 = 102531392,
  chr16 = 90354753, chr17 = 81195210, chr18 = 78077248,
  chr19 = 59128983, chr20 = 63025520, chr21 = 48129895,
  chr22 = 51304566
)

# Process current chromosome
print(paste("Processing", cr))
chrexons <- result_df[result_df$chr == cr, ]
pos <- unique(unlist(apply(chrexons, 1, function(x) x["start"]:x["end"])))
Muts = data.frame(chr = cr, pos = pos, ref = NA, alt = NA) ##this is saved below

bed = data.frame(Muts$chr, as.integer(pmax(pos - 3L, 1L)), as.integer(pmin(pos + 2L, chr_lengths[cr]))) ##this is saved below

# Write output bed (one file per chromosome)
output_dir <- "data/MutTables/WholeGenomeData/temp/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
output_file <- file.path(output_dir, sprintf("TPs_bed_chr%d.bed", as.integer(args[1])))
fwrite(bed, file = output_file, sep = "\t", col.names=F)


# Write output Muts (one file per chromosome)
output_dir <- "data/MutTables/WholeGenomeData/temp/"
#dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
output_file <- file.path(output_dir, sprintf("Muts1_bed_chr%d.bed", as.integer(args[1])))
fwrite(Muts, file = output_file, sep = "\t", col.names = FALSE)

print(paste("Finished", cr))
