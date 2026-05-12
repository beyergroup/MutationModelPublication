args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("brain","breast", "colon","esophagus","kidney","liver", "lung","ovary",
            "prostate", "skin")
tissue = tissues[args]
.libPaths(new = "/data/public/cschmalo/R-4.1.2/")
library(ranger)
nThreads = 8
# nPermutations = 10000
# maxData = 50000
dir.create("data/Modeling/exomeTrainData_distancePredictors/",showWarnings=F)
dir.create("data/Modeling/exomeTrainData_distancePredictors/RF",showWarnings=F)


print(tissue)

# load data for this tissue######
load(paste0("data/MutTables/exomeTrainData_distancePredictors/", 
            tissue, "_Muts_mapped_processed.RData"))
chroms = unique(datchroms)
#####

# grow forest with impurity_corrected #####
print("growing forests with impurity_corrected")
imp = sapply(chroms, function(cr){
  cat(cr, ' ')
  trainData = dat[datchroms != cr,]
  rf = ranger(mutated ~ ., data = trainData, importance = "impurity_corrected",
              write.forest = F, seed = 1234, num.threads =  nThreads,
              respect.unordered.factors = 'partition',
              probability = T, verbose=F)
  return(rf$variable.importance)
})
save(imp, file = paste0("data/Modeling/exomeTrainData_distancePredictors/RF/", tissue,
                        "_importances_gini.RData"))
######


# create rf for prediction and test on held-out chromosomes #####
print("predictions")
predictions = lapply(chroms, function(cr){
  cat(cr, ' ')
  trainData = dat[datchroms != cr,]
  rf = ranger(mutated ~ ., data = trainData,
              write.forest = T, seed = 1234, num.threads =  nThreads,
              respect.unordered.factors = 'partition',
              probability = T, verbose=F, importance = "permutation")
  save(rf, file = paste0("data/Modeling/exomeTrainData_distancePredictors/RF/", tissue, "_", 
                         cr, "_forPrediction.RData"))
  testData = dat[datchroms == cr,]
  p = predict(rf, data = testData, num.threads=nThreads, verbose=F)
  temp = data.frame(pred = p$predictions[,2],  label = testData$mutated)
  return(temp)
})
names(predictions) = chroms
save(predictions, file = paste0("data/Modeling/exomeTrainData_distancePredictors/RF/", tissue,
                                "_predictions.RData"))
cat('\n')
#####


# create final model #####
rf = ranger(mutated ~ ., data = dat, importance = "impurity_corrected",
            write.forest = F, seed = 1234, num.threads =  nThreads,
            respect.unordered.factors = 'partition',
            probability = T, verbose=F)
importance = rf$variable.importance
rf = ranger(mutated ~ ., data = dat,
            write.forest = T, seed = 1234, num.threads =  nThreads,
            respect.unordered.factors = 'partition',
            probability = T, verbose=F)
save(rf, importance, file = paste0("data/Modeling/exomeTrainData_distancePredictors/RF/", tissue, "_", 
                                   "finalModel.RData"))

#####

print("done")




randomVar = "finished"
save(randomVar, file = "temp/finishedWithRF.R")