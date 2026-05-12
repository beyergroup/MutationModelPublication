#!/usr/bin/env Rscript
# =============================================================================
# Fig 6A statistical tests: pairwise model comparison (DeLong) and
# AUC vs random (one-sample Wilcoxon + t-test) for RF and GLM on whole-genome
# CWCV predictions across 8 tissues.
#
# Output:
#   Fig6A_delong_RF_vs_GLM_WG.csv
#   Fig6A_auc_vs_random_WG.csv
# =============================================================================

suppressPackageStartupMessages({
  .libPaths(new = "/data/public/cschmalo/R-4.1.2/")
  library(pROC)
})

# ---- paths ------------------------------------------------------------------
base_dir <- "data/Modeling/WholeGenomeData"
out_dir  <- "data/Modeling/WholeGenomeData/statisticalTestingOfPredictions"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

f_rf_pred  <- file.path(base_dir, "RF",  "predPerTissueRF.RData")
f_glm_pred <- file.path(base_dir, "GLM", "predPerTissueGLM_sig.RData")
f_rf_auc   <- file.path(base_dir, "RF",  "ROC_PR_RF_perChr.RData")
f_glm_auc  <- file.path(base_dir, "GLM", "ROC_PR_glm_perChr_sig.RData")

# ---- load -------------------------------------------------------------------
cat("Loading prediction objects...\n")
load(f_rf_pred)   # predPerTissueRF
load(f_glm_pred)  # predPerTissueGLMsig

cat("Loading per-chromosome AUC objects...\n")
load(f_rf_auc)    # ROC_PR_RF_perChr
load(f_glm_auc)   # ROC_PR_glm_perChr_sig

# ---- alignment sanity checks ------------------------------------------------
cat("\n=== Alignment checks ===\n")
stopifnot(identical(names(predPerTissueRF), names(predPerTissueGLMsig)))
tissues <- names(predPerTissueRF)
cat("Tissues:", paste(tissues, collapse = ", "), "\n")

for (tis in tissues) {
  stopifnot(identical(names(predPerTissueRF[[tis]]),
                      names(predPerTissueGLMsig[[tis]])))
}
cat("Tissue/chromosome name alignment: OK\n")

# ---- assemble per-tissue concatenated frames + label sanity check -----------
cat("\n=== Assembling concatenated predictions and verifying label alignment ===\n")

assemble_tissue <- function(tis) {
  chroms <- names(predPerTissueRF[[tis]])
  rf_list  <- predPerTissueRF[[tis]]
  glm_list <- predPerTissueGLMsig[[tis]]
  
  # row count check per chromosome
  for (cr in chroms) {
    n_rf  <- nrow(rf_list[[cr]])
    n_glm <- nrow(glm_list[[cr]])
    if (n_rf != n_glm) {
      stop(sprintf("Row count mismatch in %s %s: RF=%d GLM=%d",
                   tis, cr, n_rf, n_glm))
    }
  }
  
  rf_df  <- do.call(rbind, rf_list)
  glm_df <- do.call(rbind, glm_list)
  
  if (!identical(rf_df$label, glm_df$label)) {
    stop(sprintf("Label mismatch RF vs GLM in tissue %s", tis))
  }
  
  data.frame(
    label    = rf_df$label,
    pred_RF  = rf_df$pred,
    pred_GLM = glm_df$pred,
    stringsAsFactors = FALSE
  )
}

tissue_data <- lapply(tissues, assemble_tissue)
names(tissue_data) <- tissues
cat("Label alignment across all tissues: OK\n")
for (tis in tissues) {
  cat(sprintf("  %s: %d positions\n", tis, nrow(tissue_data[[tis]])))
}

# ---- spot check: stored chr-level AUC vs recomputed (brain chr1) -----------
cat("\n=== Spot check: stored vs recomputed per-chromosome AUC (brain chr1) ===\n")
spot_check <- function(model_label, pred_obj, auc_obj) {
  df <- pred_obj$brain$chr1
  recomputed <- as.numeric(pROC::auc(df$label, df$pred,
                                     direction = "<", quiet = TRUE))
  stored <- as.numeric(auc_obj$brain$chr1$auc)
  cat(sprintf("  %-4s  stored = %.6f   recomputed = %.6f   diff = %.2e\n",
              model_label, stored, recomputed, abs(stored - recomputed)))
}
spot_check("RF",  predPerTissueRF,    ROC_PR_RF_perChr)
spot_check("GLM", predPerTissueGLMsig, ROC_PR_glm_perChr_sig)

# ---- DeLong: RF vs GLM per tissue ------------------------------------------
cat("\n=== DeLong tests: RF vs GLM per tissue ===\n")

run_delong <- function(tis) {
  df <- tissue_data[[tis]]
  cat(sprintf("  %s (n=%d)... ", tis, nrow(df)))
  t_start <- Sys.time()
  rA <- pROC::roc(df$label, df$pred_RF,  direction = "<", quiet = TRUE)
  rB <- pROC::roc(df$label, df$pred_GLM, direction = "<", quiet = TRUE)
  tt <- pROC::roc.test(rA, rB, method = "delong",
                       paired = TRUE, alternative = "two.sided")
  aucA <- as.numeric(pROC::auc(rA))
  aucB <- as.numeric(pROC::auc(rB))
  cat(sprintf("done (%.1fs)\n", as.numeric(Sys.time() - t_start, units = "secs")))
  data.frame(
    tissue       = tis,
    model_A      = "RF",
    model_B      = "GLM",
    n_positions  = nrow(df),
    AUC_A        = aucA,
    AUC_B        = aucB,
    delta_AUC    = aucA - aucB,
    Z            = as.numeric(tt$statistic),
    p_value      = tt$p.value,
    stringsAsFactors = FALSE
  )
}

delong_results <- do.call(rbind, lapply(tissues, run_delong))
delong_results$p_adj_BH <- p.adjust(delong_results$p_value, method = "BH")

# ---- AUC > 0.5: one-sample Wilcoxon + t-test on 22 per-chromosome AUCs ------
cat("\n=== AUC vs 0.5: one-sample Wilcoxon + t-test on per-chromosome AUCs ===\n")

extract_chr_aucs <- function(auc_obj, tis) {
  chr_list <- auc_obj[[tis]]
  vapply(chr_list, function(x) as.numeric(x$auc), numeric(1))
}

run_auc_vs_random <- function(tis, model_label, auc_obj) {
  aucs <- extract_chr_aucs(auc_obj, tis)
  wt <- wilcox.test(aucs, mu = 0.5, alternative = "two.sided", exact = FALSE)
  tt <- t.test(aucs, mu = 0.5, alternative = "two.sided")
  data.frame(
    tissue       = tis,
    model        = model_label,
    n_chroms     = length(aucs),
    mean_AUC     = mean(aucs),
    median_AUC   = median(aucs),
    sd_AUC       = sd(aucs),
    min_AUC      = min(aucs),
    max_AUC      = max(aucs),
    wilcox_stat  = as.numeric(wt$statistic),
    wilcox_p     = wt$p.value,
    t_stat       = as.numeric(tt$statistic),
    t_p          = tt$p.value,
    stringsAsFactors = FALSE
  )
}

vs_random_results <- do.call(rbind, lapply(tissues, function(tis) {
  rbind(
    run_auc_vs_random(tis, "RF",  ROC_PR_RF_perChr),
    run_auc_vs_random(tis, "GLM", ROC_PR_glm_perChr_sig)
  )
}))
vs_random_results$wilcox_p_adj_BH <- p.adjust(vs_random_results$wilcox_p, method = "BH")
vs_random_results$t_p_adj_BH      <- p.adjust(vs_random_results$t_p,      method = "BH")

# ---- formatting helpers -----------------------------------------------------
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

delong_results$p_value_fmt   <- format_p(delong_results$p_value)
delong_results$p_adj_BH_fmt  <- format_p(delong_results$p_adj_BH)
delong_results$signif        <- signif_stars(delong_results$p_adj_BH)

vs_random_results$wilcox_p_fmt        <- format_p(vs_random_results$wilcox_p)
vs_random_results$wilcox_p_adj_BH_fmt <- format_p(vs_random_results$wilcox_p_adj_BH)
vs_random_results$wilcox_signif       <- signif_stars(vs_random_results$wilcox_p_adj_BH)
vs_random_results$t_p_fmt             <- format_p(vs_random_results$t_p)
vs_random_results$t_p_adj_BH_fmt      <- format_p(vs_random_results$t_p_adj_BH)
vs_random_results$t_signif            <- signif_stars(vs_random_results$t_p_adj_BH)

# ---- write CSVs -------------------------------------------------------------
out_delong  <- file.path(out_dir, "Fig6A_delong_RF_vs_GLM_WG.csv")
out_vs_rand <- file.path(out_dir, "Fig6A_auc_vs_random_WG.csv")
write.csv(delong_results,    out_delong,  row.names = FALSE)
write.csv(vs_random_results, out_vs_rand, row.names = FALSE)

cat("\n=== Wrote outputs ===\n")
cat("  ", out_delong,  "\n")
cat("  ", out_vs_rand, "\n")

# ---- console summary --------------------------------------------------------
cat("\n=== DeLong summary (BH-adjusted p < 0.05 = significant) ===\n")
print(delong_results[, c("tissue","AUC_A","AUC_B","delta_AUC",
                         "p_value_fmt","p_adj_BH_fmt","signif")],
      digits = 4, row.names = FALSE)

cat("\n=== AUC vs 0.5 summary ===\n")
print(vs_random_results[, c("tissue","model","mean_AUC","sd_AUC",
                            "wilcox_p_fmt","wilcox_p_adj_BH_fmt","wilcox_signif",
                            "t_p_fmt","t_p_adj_BH_fmt","t_signif")],
      digits = 4, row.names = FALSE)

cat("\nDone.\n")