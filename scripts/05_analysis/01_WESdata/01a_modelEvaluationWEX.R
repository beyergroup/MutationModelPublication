library(readxl)
source("scripts/05_analysis/00_NamesAndColors.R")
library(ROCR)
library(ranger) 
library(pROC)

# get mutation rate and nPositions, percTP, correlation etc. (dataInfos) #####
dataInfos = sapply(tissues, function(tissue){
  load(paste0("data/MutTables/exomeTrainData/", 
              tissue, "_Muts_mapped_processed.RData"))
  chrCount = sapply(table(datchroms), as.integer)
  chrPerc = sapply(split(dat$mutated, datchroms),function(x){mean(x==1)})
  cors = cor(dat[sapply(dat, is.numeric)], use = "pair")
  return(list(nMuts = nrow(dat), 
              percTP = mean(dat$mutated == 1), 
              nMutsPerChr = chrCount,
              percTPperChr = chrPerc, 
              cors = cors))
}, simplify = F)
save(dataInfos, file = "data/processedData/dataInfos.RData")
#####

# load predictions for each model and each tissue #####
# random forest
predPerTissueRF = sapply(tissues, function(tissue){
  load(paste0("data/Modeling/exomeTrainData/RF/", tissue,
              "_predictions.RData"))
  return(predictions)
}, simplify=F)
save(predPerTissueRF, 
     file = "data/Modeling/exomeTrainData/RF/predPerTissueRF.RData")
# load sig. GLM predictions for each tissue
predPerTissueGLMsig = sapply(tissues, function(tissue){
  load(file = paste0("data/Modeling/exomeTrainData/GLM/", tissue, 
                     "_predictions_sig.RData"))
  return(predictions)
}, simplify=F)
save(predPerTissueGLMsig, 
     file = "data/Modeling/exomeTrainData/GLM/predPerTissueGLM_sig.RData")
predPerTissueLasso = sapply(tissues, function(tissue){
  load(file = paste0("data/Modeling/exomeTrainData/Lasso/", 
                     tissue, "_predictions_sig.RData"))
  return(predictions)
}, simplify=F)
save(predPerTissueLasso, 
     file = "data/Modeling/exomeTrainData/Lasso/predPerTissueLasso.RData")
#####


# compute ROC, PR, and AUROC for each chromosome  #####
# RF
ROC_PR_RF_perChr = sapply(names(predPerTissueRF), function(tissue){ #iterate through tissues
  pred = predPerTissueRF[[tissue]]
  # get performances
  res = lapply(pred, function(x){ #iterate through chromosomes
    perf = prediction(x$pred, x$label)
    roc = performance(perf, "tpr", "fpr")
    pr = performance(perf,"prec", "rec")
    auc = performance(perf,"auc")@y.values[[1]]
    return(list(roc = roc, pr = pr, auc = auc))
  })
  return(res)
}, simplify=F)
save(ROC_PR_RF_perChr, 
     file = "data/Modeling/exomeTrainData/RF/ROC_PR_RF_perChr.RData")
# glm
ROC_PR_glm_perChr_sig = sapply(names(predPerTissueGLMsig), function(tissue){
  pred = predPerTissueGLMsig[[tissue]]
  # get performances
  res = lapply(pred, function(x){ #iterate through chromosomes
    perf = prediction(x$pred, x$label)
    roc = performance(perf, "tpr", "fpr")
    pr = performance(perf,"prec", "rec")
    auc = performance(perf,"auc")@y.values[[1]]
    return(list(roc = roc, pr = pr, auc = auc))
  })
  return(res)
}, simplify=F)
save(ROC_PR_glm_perChr_sig, 
     file = "data/Modeling/exomeTrainData/GLM/ROC_PR_glm_perChr_sig.RData")
# lasso
ROC_PR_lasso_perChr = sapply(names(predPerTissueLasso), function(tissue){
  pred = predPerTissueLasso[[tissue]]
  # get performances
  res = lapply(pred, function(x){ #iterate through chromosomes
    perf = prediction(x$pred, x$label)
    roc = performance(perf, "tpr", "fpr")
    pr = performance(perf,"prec", "rec")
    auc = performance(perf,"auc")@y.values[[1]]
    return(list(roc = roc, pr = pr, auc = auc))
  })
  return(res)
}, simplify=F)
save(ROC_PR_lasso_perChr, 
     file = "data/Modeling/exomeTrainData/Lasso/ROC_PR_lasso_perChr.RData")
#####


# compute ROC, PR, and AUROC for all chromosomes concatenated #####
ROC_PR_RF_concat = sapply(names(predPerTissueRF), function(tissue){
  predConcat = do.call(rbind,predPerTissueRF[[tissue]])
  ROC_rf = performance(prediction(predConcat$pred, 
                                  predConcat$label), 
                       "tpr", "fpr")
  AUC_rf = performance(prediction(predConcat$pred, 
                                  predConcat$label), 
                       "auc")
  PR_rf = performance(prediction(predConcat$pred, 
                                 predConcat$label), 
                      "prec", "rec")
  return(list(roc = ROC_rf, pr = PR_rf, auc = AUC_rf))
}, simplify=F)
save(ROC_PR_RF_concat, 
     file = "data/Modeling/exomeTrainData/RF/ROC_PR_RF_concat.RData")
#  glm
ROC_PR_glm_concat_sig = sapply(names(predPerTissueGLMsig), function(tissue){
  predConcat = do.call(rbind,predPerTissueGLMsig[[tissue]])
  ROC_rf = performance(prediction(predConcat$pred, 
                                  predConcat$label), 
                       "tpr", "fpr")
  AUC_rf = performance(prediction(predConcat$pred, 
                                  predConcat$label), 
                       "auc")
  PR_rf = performance(prediction(predConcat$pred, 
                                 predConcat$label), 
                      "prec", "rec")
  return(list(roc = ROC_rf, pr = PR_rf, auc = AUC_rf))
}, simplify=F)
save(ROC_PR_glm_concat_sig, 
     file = "data/Modeling/exomeTrainData/GLM/ROC_PR_glm_concat_sig.RData")
# lasso
ROC_PR_lasso_concat = sapply(names(predPerTissueLasso), function(tissue){
  predConcat = do.call(rbind,predPerTissueLasso[[tissue]])
  ROC_rf = performance(prediction(predConcat$pred, 
                                  predConcat$label), 
                       "tpr", "fpr")
  AUC_rf = performance(prediction(predConcat$pred, 
                                  predConcat$label), 
                       "auc")
  PR_rf = performance(prediction(predConcat$pred, 
                                 predConcat$label), 
                      "prec", "rec")
  return(list(roc = ROC_rf, pr = PR_rf, auc = AUC_rf))
}, simplify=F)
save(ROC_PR_lasso_concat, 
     file = "data/Modeling/exomeTrainData/Lasso/ROC_PR_lasso_concat.RData")
#####


# compute concat train performance of RF, GLM and LASSO #####
ROC_PR_RF_train = sapply(tissues, function(tissue){
  load(paste0("data/MutTables/exomeTrainData/",
              tissue, "_Muts_mapped_processed.RData"))
  chroms = unique(datchroms)
  pred = sapply(chroms, function(cr){
    trainData = dat[datchroms != cr,]
    load(paste0("data/Modeling/exomeTrainData/RF/", tissue, "_", 
                cr, "_forPrediction.RData"))
    temp = data.frame(pred = rf$predictions[,2],  label = trainData$mutated)
  }, simplify = F)
  predConcat = do.call(rbind,pred)
  ROC = performance(prediction(predConcat$pred, 
                                  predConcat$label), 
                       "tpr", "fpr")
  AUC = performance(prediction(predConcat$pred, 
                                  predConcat$label), 
                       "auc")
  PR = performance(prediction(predConcat$pred, 
                                 predConcat$label), 
                      "prec", "rec")
  return(list(roc = ROC, pr = PR, auc = AUC))
}, simplify = F)
save(ROC_PR_RF_train, 
     file = "data/Modeling/exomeTrainData/RF/ROC_PR_RF_train.RData")
ROC_PR_glm_sig_train = sapply(tissues, function(tissue){
  load(paste0("data/MutTables/exomeTrainData/",
              tissue, "_Muts_mapped_processed.RData"))
  chroms = unique(datchroms)
  pred = sapply(chroms, function(cr){
    trainData = dat[datchroms != cr,]
    load(paste0("data/Modeling/exomeTrainData/GLM/", 
                tissue, "_", cr, "_sig.RData"))
    temp = data.frame(pred = logR$fitted.values,  label = trainData$mutated)
  }, simplify = F)
  predConcat = do.call(rbind,pred)
  ROC = performance(prediction(predConcat$pred, 
                               predConcat$label), 
                    "tpr", "fpr")
  AUC = performance(prediction(predConcat$pred, 
                               predConcat$label), 
                    "auc")
  PR = performance(prediction(predConcat$pred, 
                              predConcat$label), 
                   "prec", "rec")
  return(list(roc = ROC, pr = PR, auc = AUC))
}, simplify = F)
save(ROC_PR_glm_sig_train, 
     file = "data/Modeling/exomeTrainData/RF/ROC_PR_glm_sig_train.RData")
ROC_PR_lasso_sig_train = sapply(tissues, function(tissue){
  load(paste0("data/MutTables/exomeTrainData/",
              tissue, "_Muts_mapped_processed.RData"))
  chroms = unique(datchroms)
  pred = sapply(chroms, function(cr){
    trainData = dat[datchroms != cr,]
    load(paste0("data/Modeling/exomeTrainData/Lasso/",
                tissue, "_", cr, "_sig.RData"))
    temp = data.frame(pred = logR$fitted.values,  label = trainData$mutated)
  }, simplify = F)
  predConcat = do.call(rbind,pred)
  ROC = performance(prediction(predConcat$pred, 
                               predConcat$label), 
                    "tpr", "fpr")
  AUC = performance(prediction(predConcat$pred, 
                               predConcat$label), 
                    "auc")
  PR = performance(prediction(predConcat$pred, 
                              predConcat$label), 
                   "prec", "rec")
  return(list(roc = ROC, pr = PR, auc = AUC))
}, simplify = F)
save(ROC_PR_lasso_sig_train, 
     file = "data/Modeling/exomeTrainData/RF/ROC_PR_lasso_sig_train.RData")
#####


# get predictor importances #####
# rf
rf_gini = sapply(tissues, function(tissue){
  load(paste0("data/Modeling/exomeTrainData/RF/", tissue,
              "_importances_gini.RData"))
  names(imp) = names(chrCols)
  imp = as.data.frame(imp)
  temp = imp[names(p2P),]
  rownames(temp) = names(p2P)
  return(temp)
}, simplify=F)
save(rf_gini, file = "data/Modeling/exomeTrainData/RF/RF_imps.RData")
# glm predictor importances
glm_imps = sapply(tissues, function(tissue){
  temp = sapply(names(chrCols),function(cr){
    load(paste0("data/Modeling/exomeTrainData/GLM/",
                tissue, "_", cr, ".RData"))
    logR$coefficients[names(p2P)]
  })
  rownames(temp) = names(p2P)
  return(temp)
}, simplify=F)
# glm pvals
glm_pvals = sapply(tissues, function(tissue){
  temp = sapply(names(chrCols),function(cr){
    load(paste0("data/Modeling/exomeTrainData/GLM/",
                tissue, "_", cr, ".RData"))
    pvals = coef(summary(logR))[-1,4][names(p2P)]
  })
  temp[temp == 0] = 2e-16
  rownames(temp) = names(p2P)
  return(temp)
}, simplify=F)
save(glm_imps, glm_pvals, 
     file = "data/Modeling/exomeTrainData/GLM/GLM_impsAndPvals.RData")

# glm significant coefficients
glm_imps_sig = sapply(tissues, function(tissue){
  temp = sapply(names(chrCols),function(cr){
    load(paste0("data/Modeling/exomeTrainData/GLM/", 
                tissue, "_", cr, "_sig.RData"))
    logR$coefficients[names(p2P)]
  })
  rownames(temp) = names(p2P)
  return(temp)
}, simplify=F)
# glm significant pvals
glm_pvals_sig = sapply(tissues, function(tissue){
  temp = sapply(names(chrCols),function(cr){
    load(paste0("data/Modeling/exomeTrainData/GLM/",
                tissue, "_", cr, "_sig.RData"))
    pvals = coef(summary(logR))[-1,4][names(p2P)]
  })
  temp[temp == 0] = 2e-16
  rownames(temp) = names(p2P)
  
  return(temp)
}, simplify=F)
save(glm_imps_sig, glm_pvals_sig, 
     file = "data/Modeling/exomeTrainData/GLM/GLM_impsAndPvals_sig.RData")
# lasso
lasso_stability = sapply(tissues, function(tissue){
  temp = sapply(names(chrCols),function(cr){
    load(paste0("data/Modeling/exomeTrainData/Lasso/",
                tissue, "_", cr, ".RData")) # sp and stab
    sp$x[,stab$lpos][names(p2P)]
  })
  rownames(temp) = names(p2P)
  return(temp)
}, simplify=F)
lasso_imp = sapply(tissues, function(tissue){
  temp = sapply(names(chrCols),function(cr){
    load(paste0("data/Modeling/exomeTrainData/Lasso/",
                tissue, "_", cr, "_sig.RData")) # logR, sigFeatures
    logR$coefficients[names(p2P)]
  })
  rownames(temp) = names(p2P)
  return(temp)
}, simplify=F)
save(lasso_stability, lasso_imp, 
     file = "data/Modeling/exomeTrainData/Lasso/Lasso_impAndStab.RData")
#####
