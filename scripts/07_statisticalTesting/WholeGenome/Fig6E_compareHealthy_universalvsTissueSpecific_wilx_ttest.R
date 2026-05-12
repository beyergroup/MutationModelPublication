#!/usr/bin/env Rscript
# =============================================================================
# Fig 6E statistical tests: whole-genome healthy tissue AUCs.
#
# Test 1 (paired, n = 7 tissues): paired Wilcoxon + paired t-test, two-sided,
#   universal vs tissue-specific WG model, on the 7 tissues with both AUCs.
# Test 2 (one-sample, n = 28 healthy tissues): one-sample Wilcoxon + t-test,
#   two-sided, universal-WG-model AUC vs 0.5.
#
# Output:
#   Fig6E_universal_vs_tissuespecific_healthy_paired_WG.csv
#   Fig6E_universal_AUC_vs_random_healthy_WG.csv
#   Fig6E_universal_AUC_vs_random_healthy_WG_summary.csv
# to run script: cd /cellfile/cellnet/MutationModel    ##also start the conda R env
#  Rscript scripts/07_statisticalTesting/WholeGenome/Fig6E_compareHealthy_universalvsTissueSpecific_wilx_ttest.R
# =============================================================================

suppressPackageStartupMessages({
  .libPaths(new = "/data/public/cschmalo/R-4.1.2/")
})

# ---- paths ------------------------------------------------------------------
healthy_dir <- "data/Modeling/healthyTissues"
out_dir     <- "data/Modeling/WholeGenomeData/statisticalTestingOfPredictions"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

f_univ <- file.path(healthy_dir, "perfsAllTissueModel_WGS.RData")
f_ts   <- file.path(healthy_dir, "perfsTissueSpecific_WGS.RData")

# ---- load -------------------------------------------------------------------
cat("Loading WG healthy tissue performance objects...\n")
load(f_univ)   # perfsAllTissueModel  (also brings perfsTissueSpecific from same file)
load(f_ts)     # perfsTissueSpecific  (overrides if present, ensures canonical version)

# ---- extract AUC scalars from S4 ROCR objects -------------------------------
extract_auc <- function(perf_list) {
  vapply(perf_list, function(x) x$auc@y.values[[1]], numeric(1))
}
auc_universal      <- extract_auc(perfsAllTissueModel)   # 28 tissues
auc_tissuespecific <- extract_auc(perfsTissueSpecific)   # 7 tissues

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
# Test 1: paired test, n = 7 tissues
# =============================================================================

cat("\n=== Test 1: paired test (universal vs tissue-specific WG, n=7 tissues) ===\n")

# Match by tolower()
univ_lc <- setNames(auc_universal,      tolower(names(auc_universal)))
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
  test_description = "Paired test across healthy tissues (WG): tissue-specific vs universal model AUC",
  n_tissues        = length(paired_tissues),
  tissues          = paste(paired_tissues, collapse = ";"),
  mean_AUC_TS      = mean(ts_paired),
  mean_AUC_UNIV    = mean(univ_paired),
  mean_delta       = mean(diffs),
  median_delta     = median(diffs),
  sd_delta         = sd(diffs),
  wilcox_stat      = as.numeric(wt$statistic),
  wilcox_p         = wt$p.value,
  t_stat           = as.numeric(tt$statistic),
  t_p              = tt$p.value,
  stringsAsFactors = FALSE
)
paired_results$wilcox_p_fmt  <- format_p(paired_results$wilcox_p)
paired_results$wilcox_signif <- signif_stars(paired_results$wilcox_p)
paired_results$t_p_fmt       <- format_p(paired_results$t_p)
paired_results$t_signif      <- signif_stars(paired_results$t_p)

out_paired <- file.path(out_dir, "Fig6E_universal_vs_tissuespecific_healthy_paired_WG.csv")
write.csv(paired_results, out_paired, row.names = FALSE)
cat(sprintf("Wrote: %s\n", out_paired))

# =============================================================================
# Test 2: one-sample, universal AUC vs 0.5 across all healthy tissues
# =============================================================================

cat("\n=== Test 2: one-sample test (universal-WG-model AUC vs 0.5) ===\n")

aucs <- auc_universal
n <- length(aucs)
cat(sprintf("N healthy tissues with universal-WG-model AUC: %d\n", n))

wt2 <- wilcox.test(aucs, mu = 0.5, alternative = "two.sided", exact = FALSE)
tt2 <- t.test(aucs, mu = 0.5, alternative = "two.sided")

per_tissue <- data.frame(
  tissue = names(aucs),
  AUC    = as.numeric(aucs),
  stringsAsFactors = FALSE
)

across_tissue_stats <- data.frame(
  test_description = sprintf("One-sample test across %d healthy tissues (WG): universal-model AUC vs 0.5", n),
  n_tissues        = n,
  mean_AUC         = mean(aucs),
  median_AUC       = median(aucs),
  sd_AUC           = sd(aucs),
  min_AUC          = min(aucs),
  max_AUC          = max(aucs),
  wilcox_stat      = as.numeric(wt2$statistic),
  wilcox_p         = wt2$p.value,
  t_stat           = as.numeric(tt2$statistic),
  t_p              = tt2$p.value,
  stringsAsFactors = FALSE
)
across_tissue_stats$wilcox_p_fmt  <- format_p(across_tissue_stats$wilcox_p)
across_tissue_stats$wilcox_signif <- signif_stars(across_tissue_stats$wilcox_p)
across_tissue_stats$t_p_fmt       <- format_p(across_tissue_stats$t_p)
across_tissue_stats$t_signif      <- signif_stars(across_tissue_stats$t_p)

out_pertissue <- file.path(out_dir, "Fig6E_universal_AUC_vs_random_healthy_WG.csv")
out_summary   <- file.path(out_dir, "Fig6E_universal_AUC_vs_random_healthy_WG_summary.csv")
write.csv(per_tissue,          out_pertissue, row.names = FALSE)
write.csv(across_tissue_stats, out_summary,   row.names = FALSE)
cat(sprintf("Wrote per-tissue AUCs:           %s\n", out_pertissue))
cat(sprintf("Wrote across-tissue test result: %s\n", out_summary))

# ---- console summary --------------------------------------------------------
cat("\n=== Fig 6E summary ===\n")
cat("Significance: *** p<0.001, ** p<0.01, * p<0.05, ns >=0.05\n\n")

cat("--- Test 1 (paired, n=7 tissues) ---\n")
cat("Sign of mean_delta: positive = tissue-specific > universal\n")
print(paired_results[, c("n_tissues","mean_AUC_TS","mean_AUC_UNIV","mean_delta",
                         "wilcox_p_fmt","wilcox_signif","t_p_fmt","t_signif")],
      digits = 4, row.names = FALSE)

cat(sprintf("\n--- Test 2 (one-sample vs 0.5, n=%d tissues) ---\n", n))
print(across_tissue_stats[, c("n_tissues","mean_AUC","sd_AUC","min_AUC","max_AUC",
                              "wilcox_p_fmt","wilcox_signif","t_p_fmt","t_signif")],
      digits = 4, row.names = FALSE)

cat("\nDone.\n")