#!/usr/bin/env Rscript
# =============================================================================
# Fig 6D statistical tests: cross-domain RF performance, n = 8 tissues.
#
# Test 1: WEX-on-WEX (within-domain) vs WGS-on-WEX (cross-applied)
#         Paired Wilcoxon + paired t-test, two-sided.
# Test 2: WGS-on-WGS (within-domain) vs WEX-on-WGS (cross-applied)
#         Paired Wilcoxon + paired t-test, two-sided.
#
# Output:
#   Fig6D_crossDomain_AUCs.csv                    (per-tissue AUCs, 4 columns)
#   Fig6D_crossDomain_paired_tests.csv            (2 rows, one per test)
# to run script: cd /cellfile/cellnet/MutationModel    ##also start the conda R env
#  Rscript scripts/07_statisticalTesting/WholeGenome/Fig6D_compareWEX_to_WGS_viceVersa_and_Self.R
# =============================================================================

suppressPackageStartupMessages({
  .libPaths(new = "/data/public/cschmalo/R-4.1.2/")
})

# ---- paths ------------------------------------------------------------------
out_dir <- "data/Modeling/WholeGenomeData/statisticalTestingOfPredictions"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

f_wgsOnWex <- "data/Modeling/WGSmodelOnWEX.RData"
f_wexOnWgs <- "data/Modeling/WEXmodelOnWGS.RData"
f_wexSelf  <- "data/Modeling/exomeTrainData/RF/ROC_PR_RF_concat.RData"
f_wgsSelf  <- "data/Modeling/WholeGenomeData/RF/ROC_PR_RF_concat.RData"

tissues_6D <- c("brain", "breast", "esophagus", "kidney", "liver",
                "ovary", "prostate", "skin")

# ---- load -------------------------------------------------------------------
cat("Loading cross-domain AUC objects...\n")
load(f_wgsOnWex)   # WGSmodelOnWEX (list of S4 performance objects)
load(f_wexOnWgs)   # WEXmodelOnWGS

cat("Loading self (within-domain) AUC objects...\n")
load(f_wexSelf)    # ROC_PR_RF_concat (exome)
WEXself <- ROC_PR_RF_concat
load(f_wgsSelf)    # ROC_PR_RF_concat (WG, overrides)
WGSself <- ROC_PR_RF_concat
rm(ROC_PR_RF_concat)

# ---- extract AUC scalars ----------------------------------------------------
auc_WEXself     <- vapply(WEXself[tissues_6D],     function(x) x$auc@y.values[[1]], numeric(1))
auc_WGSself     <- vapply(WGSself[tissues_6D],     function(x) x$auc@y.values[[1]], numeric(1))
auc_WGSonWEX    <- vapply(WGSmodelOnWEX[tissues_6D], function(x) x@y.values[[1]],     numeric(1))
auc_WEXonWGS    <- vapply(WEXmodelOnWGS[tissues_6D], function(x) x@y.values[[1]],     numeric(1))

# ---- write per-tissue AUC table ---------------------------------------------
auc_table <- data.frame(
  tissue       = tissues_6D,
  WEX_on_WEX   = auc_WEXself,
  WGS_on_WEX   = auc_WGSonWEX,
  WGS_on_WGS   = auc_WGSself,
  WEX_on_WGS   = auc_WEXonWGS,
  delta_WEX    = auc_WEXself - auc_WGSonWEX,   # within-WEX minus cross
  delta_WGS    = auc_WGSself - auc_WEXonWGS,   # within-WGS minus cross
  stringsAsFactors = FALSE
)
out_aucs <- file.path(out_dir, "Fig6D_crossDomain_AUCs.csv")
write.csv(auc_table, out_aucs, row.names = FALSE)
cat(sprintf("Wrote per-tissue AUC table: %s\n", out_aucs))

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

# ---- paired tests -----------------------------------------------------------
run_paired <- function(within_auc, cross_auc, label) {
  diffs <- within_auc - cross_auc
  wt <- wilcox.test(within_auc, cross_auc, paired = TRUE,
                    alternative = "two.sided", exact = FALSE)
  tt <- t.test(within_auc, cross_auc, paired = TRUE,
               alternative = "two.sided")
  data.frame(
    test_description = label,
    n_tissues        = length(within_auc),
    mean_AUC_within  = mean(within_auc),
    mean_AUC_cross   = mean(cross_auc),
    mean_delta       = mean(diffs),
    median_delta     = median(diffs),
    sd_delta         = sd(diffs),
    n_pos_diffs      = sum(diffs > 0),
    n_neg_diffs      = sum(diffs < 0),
    wilcox_stat      = as.numeric(wt$statistic),
    wilcox_p         = wt$p.value,
    t_stat           = as.numeric(tt$statistic),
    t_p              = tt$p.value,
    stringsAsFactors = FALSE
  )
}

cat("\n=== Test 1: WEX-on-WEX vs WGS-on-WEX (paired, n=8) ===\n")
test1 <- run_paired(auc_WEXself, auc_WGSonWEX,
                    "WEX-on-WEX (within) vs WGS-on-WEX (cross)")

cat("=== Test 2: WGS-on-WGS vs WEX-on-WGS (paired, n=8) ===\n")
test2 <- run_paired(auc_WGSself, auc_WEXonWGS,
                    "WGS-on-WGS (within) vs WEX-on-WGS (cross)")

results <- rbind(test1, test2)
results$wilcox_p_fmt  <- format_p(results$wilcox_p)
results$wilcox_signif <- signif_stars(results$wilcox_p)
results$t_p_fmt       <- format_p(results$t_p)
results$t_signif      <- signif_stars(results$t_p)

out_tests <- file.path(out_dir, "Fig6D_crossDomain_paired_tests.csv")
write.csv(results, out_tests, row.names = FALSE)
cat(sprintf("Wrote test results: %s\n", out_tests))

# ---- console summary --------------------------------------------------------
cat("\n=== Fig 6D summary ===\n")
cat("Significance: *** p<0.001, ** p<0.01, * p<0.05, ns >=0.05\n")
cat("Sign of mean_delta: positive = within-domain > cross-applied\n\n")

cat("--- Per-tissue AUCs ---\n")
print(auc_table, digits = 4, row.names = FALSE)

cat("\n--- Paired tests ---\n")
print(results[, c("test_description","n_tissues","mean_AUC_within","mean_AUC_cross",
                  "mean_delta","n_pos_diffs","n_neg_diffs",
                  "wilcox_p_fmt","wilcox_signif","t_p_fmt","t_signif")],
      digits = 4, row.names = FALSE)

cat("\nDone.\n")

####addition later for better than random text in the manuscript:
##Even though the whole-genome model predicted 
##mutation propensities in exonic regions better
##than random, its performance was worse than the exome-specific 
##model for all tissues that we tested.

auc_WGS_on_WEX <- c(0.582395, 0.551772, 0.580024, 0.538911,
                    0.564300, 0.542229, 0.575928, 0.599756)
t.test(auc_WGS_on_WEX, mu = 0.5)