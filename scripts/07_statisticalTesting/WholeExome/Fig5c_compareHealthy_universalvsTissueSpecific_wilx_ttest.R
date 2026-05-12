#!/usr/bin/env Rscript
# =============================================================================
# Fig 5C statistical tests: healthy tissue AUCs.
#
# Test 1 (paired, n = 9 tissues): paired Wilcoxon + paired t-test, two-sided,
#   universal vs tissue-specific, on the 9 tissues where both AUCs exist.
# Test 2 (one-sample, n = N healthy tissues): one-sample Wilcoxon + t-test,
#   two-sided, universal-model AUC vs 0.5, across all healthy tissues with a
#   universal-model AUC available.
#
# Output:
#   Fig5C_universal_vs_tissuespecific_healthy_paired.csv   (1 row, n=9 tissues)
#   Fig5C_universal_AUC_vs_random_healthy.csv              (per-tissue AUCs +
#                                                           one across-tissue
#                                                           summary row)
# to run script: cd /cellfile/cellnet/MutationModel    ##also start the conda R env
#  Rscript scripts/07_statisticalTesting/WholeExome/Fig5c_compareHealthy_universalvsTissueSpecific_wilx_ttest.R
# =============================================================================

suppressPackageStartupMessages({
  .libPaths(new = "/data/public/cschmalo/R-4.1.2/")
})

# ---- paths ------------------------------------------------------------------
healthy_dir <- "data/Modeling/healthyTissues"
out_dir     <- "data/Modeling/WholeExomeData/statisticalTestingOfPredictions"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

f_univ <- file.path(healthy_dir, "perfsAllTissueModel.RData")
f_ts   <- file.path(healthy_dir, "perfsTissueSpecific.RData")

# ---- load -------------------------------------------------------------------
cat("Loading healthy tissue performance objects...\n")
load(f_univ)   # perfsAllTissueModel
load(f_ts)     # perfsTissueSpecific

# ---- extract AUC scalars from S4 ROCR objects -------------------------------
extract_auc <- function(perf_list) {
  vapply(perf_list, function(x) x$auc@y.values[[1]], numeric(1))
}
auc_universal       <- extract_auc(perfsAllTissueModel)   # ~25 tissues
auc_tissuespecific  <- extract_auc(perfsTissueSpecific)   # ~9 tissues

cat(sprintf("Universal-model AUCs available: %d tissues\n",  length(auc_universal)))
cat(sprintf("Tissue-specific AUCs available: %d tissues\n",  length(auc_tissuespecific)))

# ---- formatting helpers -----------------------------------------------------
format_p <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p >= 1e-4,
                formatC(p, format = "f", digits = 4),
                formatC(p, format = "e", digits = 2)))
}
signif_stars <- function(p) {
  as.character(cut(p,
                   breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                   labels = c("***", "**", "*", "ns"),
                   right  = TRUE))
}

# =============================================================================
# Test 1: paired test, n = 9 tissues
# =============================================================================

cat("\n=== Test 1: paired test (universal vs tissue-specific, n=9 tissues) ===\n")

# Match by tolower()
univ_lc <- setNames(auc_universal, tolower(names(auc_universal)))
ts_lc   <- setNames(auc_tissuespecific, tolower(names(auc_tissuespecific)))

paired_tissues <- intersect(names(univ_lc), names(ts_lc))
stopifnot(length(paired_tissues) > 0)
cat(sprintf("Paired tissues (n=%d): %s\n",
            length(paired_tissues), paste(paired_tissues, collapse = ", ")))

univ_paired <- univ_lc[paired_tissues]
ts_paired   <- ts_lc[paired_tissues]
diffs <- ts_paired - univ_paired   # positive => tissue-specific better

wt <- wilcox.test(ts_paired, univ_paired, paired = TRUE,
                  alternative = "two.sided", exact = FALSE)
tt <- t.test(ts_paired, univ_paired, paired = TRUE, alternative = "two.sided")

paired_results <- data.frame(
  test_description    = "Paired test across healthy tissues: tissue-specific vs universal model AUC",
  n_tissues           = length(paired_tissues),
  tissues             = paste(paired_tissues, collapse = ";"),
  mean_AUC_TS         = mean(ts_paired),
  mean_AUC_UNIV       = mean(univ_paired),
  mean_delta          = mean(diffs),
  median_delta        = median(diffs),
  sd_delta            = sd(diffs),
  wilcox_stat         = as.numeric(wt$statistic),
  wilcox_p            = wt$p.value,
  t_stat              = as.numeric(tt$statistic),
  t_p                 = tt$p.value,
  stringsAsFactors    = FALSE
)
paired_results$wilcox_p_fmt <- format_p(paired_results$wilcox_p)
paired_results$wilcox_signif <- signif_stars(paired_results$wilcox_p)
paired_results$t_p_fmt       <- format_p(paired_results$t_p)
paired_results$t_signif      <- signif_stars(paired_results$t_p)

out_paired <- file.path(out_dir, "Fig5C_universal_vs_tissuespecific_healthy_paired.csv")
write.csv(paired_results, out_paired, row.names = FALSE)
cat(sprintf("Wrote: %s\n", out_paired))

# =============================================================================
# Test 2: one-sample test of universal-model AUC vs 0.5 across healthy tissues
# =============================================================================

cat("\n=== Test 2: one-sample test (universal-model AUC vs 0.5, across healthy tissues) ===\n")

aucs <- auc_universal
n <- length(aucs)
cat(sprintf("N healthy tissues with universal-model AUC: %d\n", n))

wt2 <- wilcox.test(aucs, mu = 0.5, alternative = "two.sided", exact = FALSE)
tt2 <- t.test(aucs, mu = 0.5, alternative = "two.sided")

# Per-tissue rows for transparency, plus one summary row
per_tissue <- data.frame(
  tissue   = names(aucs),
  AUC      = as.numeric(aucs),
  stringsAsFactors = FALSE
)

summary_row <- data.frame(
  tissue = sprintf("[ACROSS-TISSUE TEST n=%d]", n),
  AUC    = NA_real_,
  stringsAsFactors = FALSE
)

across_tissue_stats <- data.frame(
  test_description    = sprintf("One-sample test across %d healthy tissues: universal-model AUC vs 0.5", n),
  n_tissues           = n,
  mean_AUC            = mean(aucs),
  median_AUC          = median(aucs),
  sd_AUC              = sd(aucs),
  min_AUC             = min(aucs),
  max_AUC             = max(aucs),
  wilcox_stat         = as.numeric(wt2$statistic),
  wilcox_p            = wt2$p.value,
  t_stat              = as.numeric(tt2$statistic),
  t_p                 = tt2$p.value,
  stringsAsFactors    = FALSE
)
across_tissue_stats$wilcox_p_fmt  <- format_p(across_tissue_stats$wilcox_p)
across_tissue_stats$wilcox_signif <- signif_stars(across_tissue_stats$wilcox_p)
across_tissue_stats$t_p_fmt       <- format_p(across_tissue_stats$t_p)
across_tissue_stats$t_signif      <- signif_stars(across_tissue_stats$t_p)

# Write per-tissue table and across-tissue summary as separate CSVs (cleaner)
out_pertissue <- file.path(out_dir, "Fig5C_universal_AUC_vs_random_healthy.csv")
out_summary   <- file.path(out_dir, "Fig5C_universal_AUC_vs_random_healthy_summary.csv")
write.csv(per_tissue,           out_pertissue, row.names = FALSE)
write.csv(across_tissue_stats,  out_summary,   row.names = FALSE)
cat(sprintf("Wrote per-tissue AUCs:           %s\n", out_pertissue))
cat(sprintf("Wrote across-tissue test result: %s\n", out_summary))

# ---- console summary --------------------------------------------------------
cat("\n=== Fig 5C summary ===\n")
cat("Significance: *** p<0.001, ** p<0.01, * p<0.05, ns >=0.05\n\n")

cat("--- Test 1 (paired, n=9 tissues) ---\n")
cat("Sign of mean_delta: positive = tissue-specific > universal\n")
print(paired_results[, c("n_tissues","mean_AUC_TS","mean_AUC_UNIV","mean_delta",
                         "wilcox_p_fmt","wilcox_signif","t_p_fmt","t_signif")],
      digits = 4, row.names = FALSE)

cat(sprintf("\n--- Test 2 (one-sample vs 0.5, n=%d tissues) ---\n", n))
print(across_tissue_stats[, c("n_tissues","mean_AUC","sd_AUC","min_AUC","max_AUC",
                              "wilcox_p_fmt","wilcox_signif","t_p_fmt","t_signif")],
      digits = 4, row.names = FALSE)

cat("\nDone.\n")


##later added for 5C additionally for the claim:
#The predictions were better than random (AUC>0.5) for all tissues that
#were used for the training. The only exception were mutations in the kidney, 
#where the kidney-specific model reached suboptimal performance (AUC < 0.5).

auc_trained <- c(brain=0.581, breast=0.570, colon=0.603, esophagus=0.531,
                 kidney=0.480, liver=0.559, lung=0.538, prostate=0.583, skin=0.562)

# Two-sided
t.test(auc_trained, mu = 0.5, alternative = "two.sided")

# One-sided (greater)
t.test(auc_trained, mu = 0.5, alternative = "greater")

