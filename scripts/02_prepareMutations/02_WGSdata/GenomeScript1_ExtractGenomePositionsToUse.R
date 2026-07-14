# Load required libraries
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19) 
library(rtracklayer)

# 1. Load the genome assembly (using BSgenome for hg19/GRCh37)


# 2. Subset to autosomes (chromosomes 1-22)
autosomes <- paste0("chr", 1:22)
genome_autosomes <- getSeq(BSgenome.Hsapiens.UCSC.hg19, autosomes)


# 3. Load the excludable regions 
load("data/predictors/wgEncodeDacMapabilityConsensusExcludable.RData")
##the loaded object is called gr
# Create a GRanges object covering the entire autosomes
autosome_ranges <- GRanges(seqnames = autosomes,
                           ranges = IRanges(start = 1, 
                                            end = seqlengths(genome_autosomes)[autosomes]))

# Find the non-excludable regions (set difference)
non_excludable <- setdiff(autosome_ranges, gr) ## retaining the regions we need

# Convert to data frame format as specified
result_df <- as.data.frame(non_excludable)
result_df <- result_df[, c("seqnames", "start", "end", "width", "strand")]
colnames(result_df)[1] <- "chr"

# Save as RData
save(result_df, file = "data/processedData/autosomes_non_excludable_regions.RData")

# Optional: View the first few rows
head(result_df)