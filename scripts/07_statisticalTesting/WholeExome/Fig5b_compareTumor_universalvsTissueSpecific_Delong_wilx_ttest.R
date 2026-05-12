#!/usr/bin/env Rscript
# =============================================================================
# Fig 5B statistical tests: universal (all-tissue) RF model vs. tissue-specific
# RF models on whole-exome cancer data.
#
# Step 1: Compute per-tissue per-chromosome AUCs from the universal model by
#         applying each per-fold universal RF to each tissue's CWCV test
#         positions for that chromosome. (220 AUCs: 10 tissues x 22 chroms)
# Step 2: Paired Wilcoxon and paired t-test (two-sided), per tissue, comparing
#         the 22 universal-model AUCs vs. the 22 tissue-specific AUCs already
#         saved in ROC_PR_RF_perChr. BH-FDR across the 10 tissues.
#
# Output:
#   Fig5B_universal_perTissue_perChr_AUC.csv         (220 rows)
#   Fig5B_universal_vs_tissuespecific_exome.csv      (10 rows)
# =============================================================================

suppressPackageStartupMessages({
  .libPaths(new = "/data/public/cschmalo/R-4.1.2/")
  library(ranger)
  library(pROC)
})

# ---- paths ------------------------------------------------------------------
mut_dir   <- "data/MutTables/exomeTrainData"
rf_dir    <- "data/Modeling/exomeTrainData/RF"
out_dir   <- "data/Modeling/WholeExomeData/statisticalTestingOfPredictions"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

f_rf_auc  <- file.path(rf_dir, "ROC_PR_RF_perChr.RData")    # tissue-specific per-chr AUCs

universal_model_path <- function(cr) {
  file.path(rf_dir, paste0("generalModel__", cr, "_forPrediction.RData"))
}
tissue_data_path <- function(tis) {
  file.path(mut_dir, paste0(tis, "_Muts_mapped_processed.RData"))
}

n_threads <- 14   # for ranger predict

tissues <- c("brain", "breast", "colon", "esophagus", "kidney", "liver",
             "lung", "ovary", "prostate", "skin")

# =============================================================================
# Step 1: per-tissue per-chromosome AUCs from the universal model
# =============================================================================

# ---- pre-load all tissue data once ------------------------------------------
cat("=== Loading tissue data ===\n")
tissue_dat     <- list()
tissue_chroms  <- list()
for (tis in tissues) {
  cat(sprintf("  %s...\n", tis))
  load(tissue_data_path(tis))   # loads `dat` and `datchroms`
  tissue_dat[[tis]]    <- dat
  tissue_chroms[[tis]] <- datchroms
  rm(dat, datchroms)
}

# Sanity: chromosome sets should be identical across tissues
chrom_sets <- lapply(tissue_chroms, function(x) sort(unique(as.character(x))))
stopifnot(length(unique(lapply(chrom_sets, paste, collapse=","))) == 1)
chroms <- chrom_sets[[1]]
cat(sprintf("Chromosomes: %d (%s)\n", length(chroms), paste(chroms, collapse=",")))

# ---- iterate: outer loop over chromosomes -----------------------------------
cat("\n=== Computing universal-model AUCs per (tissue, chromosome) ===\n")
auc_records <- list()

for (cr in chroms) {
  cat(sprintf("  %s ", cr))
  # Load universal RF for this fold (once per chromosome)
  load(universal_model_path(cr))   # loads `rf`
  univ_features <- rf$forest$independent.variable.names
  
  for (tis in tissues) {
    test_idx  <- which(tissue_chroms[[tis]] == cr)
    if (length(test_idx) == 0) next
    test_data <- tissue_dat[[tis]][test_idx, , drop = FALSE]
    # Restrict to features used by universal model (universal was trained on
    # the intersection of features across tissues; some tissues may have extras)
    missing_feats <- setdiff(univ_features, colnames(test_data))
    if (length(missing_feats) > 0) {
      stop(sprintf("Tissue %s missing universal-model features: %s",
                   tis, paste(missing_feats, collapse=", ")))
    }
    test_data <- test_data[, c(univ_features, "mutated"), drop = FALSE]
    
    p <- predict(rf, data = test_data, num.threads = n_threads, verbose = FALSE)
    pred <- p$predictions[, 2]
    label <- as.numeric(as.character(test_data$mutated))
    
    auc_val <- as.numeric(pROC::auc(label, pred, direction = "<", quiet = TRUE))
    auc_records[[length(auc_records) + 1]] <- data.frame(
      tissue     = tis,
      chromosome = cr,
      n_test     = length(label),
      AUC        = auc_val,
      stringsAsFactors = FALSE
    )
  }
  rm(rf); gc(verbose = FALSE)
}
cat("\n")

univ_auc_df <- do.call(rbind, auc_records)
out_univ_aucs <- file.path(out_dir, "Fig5B_universal_perTissue_perChr_AUC.csv")
write.csv(univ_auc_df, out_univ_aucs, row.names = FALSE)
cat(sprintf("Wrote universal per-chr AUCs: %s\n", out_univ_aucs))

# ---- Spot-check: our per-chr AUCs vs existing pooled crossTissuePerformance -
# crossTissuePerformance was computed on positions concatenated across all 22
# chromosomes (one AUC per tissue). It should be CLOSE to (not identical to)
# the mean of our 22 per-chromosome AUCs. Large divergence => something wrong.
cat("\n=== Spot check: our per-chr universal AUCs vs existing pooled AUCs ===\n")
cross_perf_path <- "/cellfile/cellnet/MutationModel/data/Modeling/exomeTrainData/CrossTissue/generalModel_crossTissuePerformance.RData"
if (file.exists(cross_perf_path)) {
  load(cross_perf_path)   # crossTissuePerformance: named numeric vector
  spot <- data.frame(
    tissue            = tissues,
    mean_perChr_AUC   = vapply(tissues, function(t) {
      mean(univ_auc_df$AUC[univ_auc_df$tissue == t])
    }, numeric(1)),
    pooled_AUC_stored = as.numeric(crossTissuePerformance[tissues]),
    stringsAsFactors  = FALSE
  )
  spot$abs_diff <- abs(spot$mean_perChr_AUC - spot$pooled_AUC_stored)
  print(spot, digits = 4, row.names = FALSE)
  if (any(spot$abs_diff > 0.02, na.rm = TRUE)) {
    warning("Some tissues show mean(per-chr AUC) and pooled AUC differing by >0.02 ",
            "- worth a manual look before trusting Step 2 results.")
  } else {
    cat("All tissues within 0.02 of stored pooled AUC. Sanity check passed.\n")
  }
} else {
  cat(sprintf("Skipped: %s not found.\n", cross_perf_path))
}

# =============================================================================
# Step 2: paired Wilcoxon + t-test, per tissue
# =============================================================================

cat("\n=== Loading tissue-specific per-chromosome AUCs ===\n")
load(f_rf_auc)   # ROC_PR_RF_perChr

extract_chr_aucs <- function(auc_obj, tis) {
  chr_list <- auc_obj[[tis]]
  vapply(chr_list, function(x) as.numeric(x$auc), numeric(1))
}

cat("\n=== Paired tests: universal vs tissue-specific (per tissue) ===\n")
test_records <- lapply(tissues, function(tis) {
  ts_aucs   <- extract_chr_aucs(ROC_PR_RF_perChr, tis)
  univ_sub  <- univ_auc_df[univ_auc_df$tissue == tis, ]
  univ_aucs <- setNames(univ_sub$AUC, univ_sub$chromosome)
  
  # Align by chromosome name (paired)
  common <- intersect(names(ts_aucs), names(univ_aucs))
  stopifnot(length(common) == length(ts_aucs))
  ts_aucs   <- ts_aucs[common]
  univ_aucs <- univ_aucs[common]
  
  diffs <- ts_aucs - univ_aucs   # positive => tissue-specific better
  
  wt <- wilcox.test(ts_aucs, univ_aucs, paired = TRUE,
                    alternative = "two.sided", exact = FALSE)
  tt <- t.test(ts_aucs, univ_aucs, paired = TRUE,
               alternative = "two.sided")
  
  data.frame(
    tissue         = tis,
    n_chroms       = length(common),
    mean_AUC_TS    = mean(ts_aucs),
    mean_AUC_UNIV  = mean(univ_aucs),
    mean_delta     = mean(diffs),       # TS - Universal
    median_delta   = median(diffs),
    sd_delta       = sd(diffs),
    wilcox_stat    = as.numeric(wt$statistic),
    wilcox_p       = wt$p.value,
    t_stat         = as.numeric(tt$statistic),
    t_p            = tt$p.value,
    stringsAsFactors = FALSE
  )
})
results <- do.call(rbind, test_records)
results$wilcox_p_adj_BH <- p.adjust(results$wilcox_p, method = "BH")
results$t_p_adj_BH      <- p.adjust(results$t_p,      method = "BH")

# ---- readability helpers ----------------------------------------------------
format_p <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p >= 1e-4,
                formatC(p, format = "f", digits = 4),
                formatC(p, format = "e", digits = 2)))
}
signif_stars <- function(p_adj) {
  as.character(cut(p_adj,
                   breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                   labels = c("***", "**", "*", "ns"),
                   right  = TRUE))
}

results$wilcox_p_fmt        <- format_p(results$wilcox_p)
results$wilcox_p_adj_BH_fmt <- format_p(results$wilcox_p_adj_BH)
results$wilcox_signif       <- signif_stars(results$wilcox_p_adj_BH)
results$t_p_fmt             <- format_p(results$t_p)
results$t_p_adj_BH_fmt      <- format_p(results$t_p_adj_BH)
results$t_signif            <- signif_stars(results$t_p_adj_BH)

# ---- write output -----------------------------------------------------------
out_results <- file.path(out_dir, "Fig5B_universal_vs_tissuespecific_exome.csv")
write.csv(results, out_results, row.names = FALSE)
cat(sprintf("Wrote test results: %s\n", out_results))

# ---- console summary --------------------------------------------------------
cat("\n=== Fig 5B summary ===\n")
cat("Significance: *** p<0.001, ** p<0.01, * p<0.05, ns >=0.05 (BH-adjusted)\n")
cat("Sign of mean_delta: positive = tissue-specific > universal\n\n")
print(results[, c("tissue", "mean_AUC_TS", "mean_AUC_UNIV", "mean_delta",
                  "wilcox_p_fmt", "wilcox_p_adj_BH_fmt", "wilcox_signif",
                  "t_p_fmt", "t_p_adj_BH_fmt", "t_signif")],
      digits = 4, row.names = FALSE)

cat("\nDone.\n")