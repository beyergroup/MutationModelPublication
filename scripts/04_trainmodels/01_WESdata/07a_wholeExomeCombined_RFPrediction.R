##steps:
##Load WE RF model
##Load liver predictors, we tried colon first, then switched to liver  to stay consistent with genome results
##make predictions for WE positions
library(ranger)

args = as.numeric(commandArgs(trailingOnly=T))
tissue = c("liver")
#tissue = tissues[args]
print(tissue)

nThreads = 28
#.libPaths(new = "/data/user/ypaul1/.cache/R/renv/source/repository")
#library(ranger)


#dir.create("/cellfile/datapublic/ypaul1/Mutations/results/WholeExomePredictions", showWarnings = F)

#load(paste0("/cellfile/datapublic/ypaul1/Mutations/data/fromCorinna/generalModel_","finalModel.RData"))

load(paste0("data/Modeling/exomeTrainData/RF/generalModel_","finalModel.RData"))

dumpVar = sapply(1:50, function(i){
  print(i)
  # load data
  load(paste0("data/MutTables/WholeExomeData/exomeMuts_part",i,"_", tissue, "_mapped.RData"))
  testDat = data$pred
  yhat = predict(rf, testDat, num.threads = nThreads)
  
  predictions = data.frame(data$muts,
                           prediction = yhat$predictions[,2])
  save(predictions, 
       file = paste0("data/Modeling/exomeTrainData/RF/exomeMuts_part",
                     i,"_",tissue,"_RFpredictionsExomeGeneric.RData"))
  return(NA)
})
print("done")
