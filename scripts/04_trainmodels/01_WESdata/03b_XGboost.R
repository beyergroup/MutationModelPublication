args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("brain","breast", "colon","esophagus","kidney","liver", "lung","ovary",
            "prostate", "skin")
tissue = tissues[args]
# .libPaths(new = "/data/public/cschmalo/R-4.1.2/")
library(xgboost)
nThreads = 8
dir.create("data/Modeling/exomeTrainData/",showWarnings=F)
dir.create("data/Modeling/exomeTrainData/xgboost",showWarnings=F)

# define parameter grid
grid = expand.grid(
  eta = c(0.01, 0.1, 0.3),
  max_depth = c(4, 6),
  subsample = c(0.7, 0.8, 0.9),
  colsample_bytree = c(0.7, 0.8, 0.9)
)

print(tissue)

# load data for this tissue######
load(paste0("data/MutTables/exomeTrainData/", 
            tissue, "_Muts_mapped_processed.RData"))
chroms = unique(datchroms)
#####

# CWCV for xgboost #####
print("doing CWCV for xgboost")
cwcv = sapply(chroms, function(cr){
  cat(cr, ' ')
  dtrain = xgb.DMatrix(data = as.matrix(dat[datchroms != cr,-which(colnames(dat) == "mutated")]), 
                       label = as.integer(as.character(dat[datchroms != cr,"mutated"])))
  dtest = xgb.DMatrix(data = as.matrix(dat[datchroms == cr,-which(colnames(dat) == "mutated")]), 
                      label = as.integer(as.character(dat[datchroms == cr,"mutated"])))
  # do grid search for parameter optimization
  results <- apply(grid, 1, function(params) {
    cv = xgb.cv(
      params = as.list(c(params, objective = "binary:logistic", eval_metric = "auc", nthread = nThreads)),
      data = dtrain,
      nrounds = 750,
      nfold = 5,
      early_stopping_rounds = 20,
      verbose = 0
    )
    best_iter <- which.max(cv$evaluation_log$test_auc_mean)
    return(c(best_nrounds = best_iter, best_auc = max(cv$evaluation_log$test_auc_mean)))
  })
  cvResults = cbind(grid, t(results))
  best = cvResults[which.max(cvResults$best_auc), ]
  xgbModel = xgb.train(data = dtrain,nrounds = best$best_nrounds,
                       params = c(best[1:4],objective = "binary:logistic", nthread = nThreads),
                       verbose = 0)
  
  # save model
  xgb.save(xgbModel, paste0("data/Modeling/exomeTrainData/xgboost/", tissue, "_", cr, "_model.ubj"))
  # get importance
  imp = xgb.importance(model = xgbModel)
  # get testData predictions
  pred_prob <- predict(xgbModel, dtest)
  preds = cbind(pred = pred_prob, label = as.integer(as.character(dat[datchroms == cr,"mutated"])))
  return(list(imp = imp$Gain, preds = preds))
}, simplify = F)
save(cwcv, file = paste0("data/Modeling/exomeTrainData/xgboost/", tissue,
                         "_cwcv.RData"))
print("saving importances")
imp = sapply(cwcv, function(x){x$imp})
save(imp, file = paste0("data/Modeling/exomeTrainData/xgboost/", tissue,
                        "_importances_gain.RData"))
print("saving predictions")
predictions = sapply(cwcv, function(x){x$preds}, simplify = F)
save(predictions, file = paste0("data/Modeling/exomeTrainData/xgboost/", tissue,
                                "_predictions.RData"))
cat('\n')
#####


# create final model #####
dtrain = xgb.DMatrix(data = as.matrix(dat[,-which(colnames(dat) == "mutated")]), 
                     label = as.integer(as.character(dat[,"mutated"])))
results <- apply(grid, 1, function(params) {
  cv = xgb.cv(
    params = as.list(c(params, objective = "binary:logistic", eval_metric = "auc", nthread = nThreads)),
    data = dtrain,
    nrounds = 750,
    nfold = 5,
    early_stopping_rounds = 20,
    verbose = 0
  )
  best_iter <- which.max(cv$evaluation_log$test_auc_mean)
  return(c(best_nrounds = best_iter, best_auc = max(cv$evaluation_log$test_auc_mean)))
})
cvResults = cbind(grid, results)
cvResults = cbind(grid, t(results))
best = cvResults[which.max(cvResults$best_auc), ]
xgbModel = xgb.train(data = dtrain,nrounds = best$best_nrounds,
                     params = c(best[1:4],objective = "binary:logistic", nthread = nThreads),
                     verbose = 0)

importance = xgb.importance(model = xgbModel)
# save model
xgb.save(xgbModel, paste0("data/Modeling/exomeTrainData/xgboost/", tissue, "_finalModel", ".ubj"))
# save importance
save(importance, file = paste0("data/Modeling/exomeTrainData/xgboost/", tissue, "_", 
                               "finalModel_importance.RData"))
#####

print("done")


