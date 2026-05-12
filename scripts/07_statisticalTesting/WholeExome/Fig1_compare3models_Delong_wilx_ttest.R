#!/usr/bin/env Rscript
# =============================================================================
# Fig 1B statistical tests: pairwise model comparisons (DeLong) and
# AUC vs. random (one-sample Wilcoxon + t-test) for RF, GLM, Lasso (SL) on
# whole-exome CWCV predictions across 10 tissues.
#
# Output:
#   Fig1_delong_pairwise_models_exome.csv
#   Fig1_auc_vs_random_exome.csv
# =============================================================================

suppressPackageStartupMessages({
  library(pROC)
})

# ---- paths ------------------------------------------------------------------
base_dir   <- "data/Modeling/exomeTrainData"
out_dir    <- "data/Modeling/WholeExomeData/statisticalTestingOfPredictions"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

f_glm_pred <- file.path(base_dir, "GLM",   "predPerTissueGLM_sig.RData")
f_las_pred <- file.path(base_dir, "Lasso", "predPerTissueLasso.RData")
f_rf_pred  <- file.path(base_dir, "RF",    "predPerTissueRF.RData")

f_glm_auc  <- file.path(base_dir, "GLM",   "ROC_PR_glm_perChr_sig.RData")
f_las_auc  <- file.path(base_dir, "Lasso", "ROC_PR_lasso_perChr.RData")
f_rf_auc   <- file.path(base_dir, "RF",    "ROC_PR_RF_perChr.RData")

# ---- load -------------------------------------------------------------------
cat("Loading prediction objects...\n")
load(f_glm_pred)   # predPerTissueGLMsig
load(f_las_pred)   # predPerTissueLasso
load(f_rf_pred)    # predPerTissueRF

cat("Loading per-chromosome AUC objects...\n")
load(f_glm_auc)    # ROC_PR_glm_perChr_sig
load(f_las_auc)    # ROC_PR_lasso_perChr
load(f_rf_auc)     # ROC_PR_RF_perChr

# ---- alignment sanity checks ------------------------------------------------
cat("\n=== Alignment checks ===\n")
stopifnot(identical(names(predPerTissueGLMsig), names(predPerTissueLasso)))
stopifnot(identical(names(predPerTissueGLMsig), names(predPerTissueRF)))
tissues <- names(predPerTissueRF)
cat("Tissues:", paste(tissues, collapse = ", "), "\n")

for (tis in tissues) {
  stopifnot(identical(names(predPerTissueGLMsig[[tis]]),
                      names(predPerTissueRF[[tis]])))
  stopifnot(identical(names(predPerTissueLasso[[tis]]),
                      names(predPerTissueRF[[tis]])))
}
cat("Tissue/chromosome name alignment: OK\n")

# ---- assemble per-tissue concatenated frames + label sanity check -----------
cat("\n=== Assembling concatenated predictions and verifying label alignment ===\n")

assemble_tissue <- function(tis) {
  chroms <- names(predPerTissueRF[[tis]])
  rf_list  <- predPerTissueRF[[tis]]
  glm_list <- predPerTissueGLMsig[[tis]]
  las_list <- predPerTissueLasso[[tis]]
  
  # row count check per chromosome
  for (cr in chroms) {
    n_rf  <- nrow(rf_list[[cr]])
    n_glm <- nrow(glm_list[[cr]])
    n_las <- nrow(las_list[[cr]])
    if (!(n_rf == n_glm && n_rf == n_las)) {
      stop(sprintf("Row count mismatch in %s %s: RF=%d GLM=%d Lasso=%d",
                   tis, cr, n_rf, n_glm, n_las))
    }
  }
  
  # concatenate
  rf_df  <- do.call(rbind, rf_list)
  glm_df <- do.call(rbind, glm_list)
  las_df <- do.call(rbind, las_list)
  
  # label alignment check
  if (!identical(rf_df$label, glm_df$label)) {
    stop(sprintf("Label mismatch RF vs GLM in tissue %s", tis))
  }
  if (!identical(rf_df$label, las_df$label)) {
    stop(sprintf("Label mismatch RF vs Lasso in tissue %s", tis))
  }
  
  data.frame(
    label    = rf_df$label,
    pred_RF  = rf_df$pred,
    pred_GLM = glm_df$pred,
    pred_SL  = las_df$pred
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
  cat(sprintf("  %-6s  stored = %.6f   recomputed = %.6f   diff = %.2e\n",
              model_label, stored, recomputed, abs(stored - recomputed)))
}
spot_check("RF",  predPerTissueRF,    ROC_PR_RF_perChr)
spot_check("GLM", predPerTissueGLMsig, ROC_PR_glm_perChr_sig)
spot_check("SL",  predPerTissueLasso, ROC_PR_lasso_perChr)

# ---- DeLong: pairwise model comparisons per tissue --------------------------
cat("\n=== DeLong tests: RF vs GLM, RF vs SL, GLM vs SL per tissue ===\n")

run_delong <- function(tis) {
  df <- tissue_data[[tis]]
  pairs <- list(
    c("RF",  "GLM"),
    c("RF",  "SL"),
    c("GLM", "SL")
  )
  out <- lapply(pairs, function(pr) {
    cA <- paste0("pred_", pr[1])
    cB <- paste0("pred_", pr[2])
    rA <- pROC::roc(df$label, df[[cA]], direction = "<", quiet = TRUE)
    rB <- pROC::roc(df$label, df[[cB]], direction = "<", quiet = TRUE)
    tt <- pROC::roc.test(rA, rB, method = "delong",
                         paired = TRUE, alternative = "two.sided")
    aucA <- as.numeric(pROC::auc(rA))
    aucB <- as.numeric(pROC::auc(rB))
    data.frame(
      tissue       = tis,
      model_A      = pr[1],
      model_B      = pr[2],
      n_positions  = nrow(df),
      AUC_A        = aucA,
      AUC_B        = aucB,
      delta_AUC    = aucA - aucB,
      Z            = as.numeric(tt$statistic),
      p_value      = tt$p.value
    )
  })
  do.call(rbind, out)
}

delong_results <- do.call(rbind, lapply(tissues, function(tis) {
  cat(sprintf("  %s...\n", tis))
  run_delong(tis)
}))
delong_results$p_adj_BH <- p.adjust(delong_results$p_value, method = "BH")

# ---- AUC > 0.5: one-sample Wilcoxon + t-test on 22 per-chromosome AUCs ------
cat("\n=== AUC vs 0.5: one-sample Wilcoxon + t-test on 22 per-chromosome AUCs ===\n")

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
    t_p          = tt$p.value
  )
}

vs_random_results <- do.call(rbind, lapply(tissues, function(tis) {
  rbind(
    run_auc_vs_random(tis, "RF",  ROC_PR_RF_perChr),
    run_auc_vs_random(tis, "GLM", ROC_PR_glm_perChr_sig),
    run_auc_vs_random(tis, "SL",  ROC_PR_lasso_perChr)
  )
}))
vs_random_results$wilcox_p_adj_BH <- p.adjust(vs_random_results$wilcox_p, method = "BH")
vs_random_results$t_p_adj_BH      <- p.adjust(vs_random_results$t_p,      method = "BH")

# ---- readability helpers ----------------------------------------------------
# Format p-values: plain decimal for moderate p (>= 1e-4), scientific for tiny.
format_p <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p >= 1e-4,
                formatC(p, format = "f", digits = 4),
                formatC(p, format = "e", digits = 2)))
}

# Significance stars on BH-adjusted p
signif_stars <- function(p_adj) {
  cut(p_adj,
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c("***", "**", "*", "ns"),
      right  = TRUE)
}

# Apply to DeLong table
delong_results$p_value_fmt  <- format_p(delong_results$p_value)
delong_results$p_adj_BH_fmt <- format_p(delong_results$p_adj_BH)
delong_results$signif       <- signif_stars(delong_results$p_adj_BH)

# Apply to AUC-vs-random table (both Wilcoxon and t-test)
vs_random_results$wilcox_p_fmt        <- format_p(vs_random_results$wilcox_p)
vs_random_results$wilcox_p_adj_BH_fmt <- format_p(vs_random_results$wilcox_p_adj_BH)
vs_random_results$wilcox_signif       <- signif_stars(vs_random_results$wilcox_p_adj_BH)
vs_random_results$t_p_fmt             <- format_p(vs_random_results$t_p)
vs_random_results$t_p_adj_BH_fmt      <- format_p(vs_random_results$t_p_adj_BH)
vs_random_results$t_signif            <- signif_stars(vs_random_results$t_p_adj_BH)

# ---- write CSVs -------------------------------------------------------------
out_delong  <- file.path(out_dir, "Fig1_delong_pairwise_models_exome.csv")
out_vs_rand <- file.path(out_dir, "Fig1_auc_vs_random_exome.csv")
write.csv(delong_results,    out_delong,  row.names = FALSE)
write.csv(vs_random_results, out_vs_rand, row.names = FALSE)

cat("\n=== Wrote outputs ===\n")
cat("  ", out_delong,  "\n")
cat("  ", out_vs_rand, "\n")

# ---- console summary --------------------------------------------------------
cat("\n=== DeLong summary (signif: *** p<0.001, ** p<0.01, * p<0.05, ns >=0.05; BH-adjusted) ===\n")
print(delong_results[, c("tissue","model_A","model_B","AUC_A","AUC_B",
                         "delta_AUC","p_value_fmt","p_adj_BH_fmt","signif")],
      digits = 4, row.names = FALSE)

cat("\n=== AUC vs 0.5 summary (Wilcoxon & t-test, BH-adjusted) ===\n")
print(vs_random_results[, c("tissue","model","mean_AUC","sd_AUC",
                            "wilcox_p_fmt","wilcox_p_adj_BH_fmt","wilcox_signif",
                            "t_p_fmt","t_p_adj_BH_fmt","t_signif")],
      digits = 4, row.names = FALSE)

cat("\nDone.\n")