

library(GenomicRanges)
tissues = c("brain","breast", "colon","esophagus","kidney","liver", "lung", "ovary", 
            "prostate", "skin") # "allTissues"
methods = c("RFpreds", "contextPreds", "mult", "mean", "combOdds", "LMcombination")
source("scripts/05_analysis/00_NamesAndColors.R")
# load the bin definitions
load("data/processedData/exomeBins.RData") # bins100, bins1kb
# includes bins that are smaller than defined range
# this was done because we otherwise would have very few complete bins
bins1kb = bins1kb[width(bins1kb) == 1000]
bins100 =  bins100[width(bins100) == 100]

# 1kb bins #####
perf1kb = sapply(tissues, function(tissue){
  print(tissue)
  # load exome-wide prediction table
  load(paste0("data/Modeling/WholeExomeData/combinedPredictions/combinedPredictions_",
              tissue, ".RData")) #data
  data = GRanges(seqnames=data$chr, 
                 ranges=IRanges(start=data$pos, width = 1),
                 mcols = data[,-(1:2)])
  # assign positions to bins
  hits  =  findOverlaps(query = bins1kb, subject = data)
  bin_idx <- queryHits(hits)
  
  # compute mean score for each bin
  avgScores = sapply(paste0("mcols.", methods), function(varname){
    scores <- mcols(data)[,varname][subjectHits(hits)]
    # Compute mean score per bin using tapply
    meanPerBin <- tapply(scores, bin_idx, mean)
    return(meanPerBin)
  })
  colnames(avgScores) = methods
  # compute sum of score for each bin
  sumScores = sapply(paste0("mcols.", methods), function(varname){
    scores <- mcols(data)[,varname][subjectHits(hits)]
    # Compute mean score per bin using tapply
    meanPerBin <- tapply(scores, bin_idx, sum)
    return(meanPerBin)
  })
  colnames(sumScores) = methods
  # compute number of mutations for each bin
  binnedMutations = sapply(c("mcols.cancerMuts","mcols.healthyMuts"), function(varname){
    scores <- mcols(data)[,varname][subjectHits(hits)]
    # Compute mean score per bin using tapply
    meanPerBin <- tapply(scores, bin_idx, sum)
    return(meanPerBin)
  })
  colnames(binnedMutations) = c("n_cancerMuts", "n_healthyMuts")
  if(tissue == "ovary"){binnedMutations = binnedMutations[,"n_cancerMuts", drop = F]}
  # binned = cbind(as.data.frame(bins1kb), avgScores, binnedMutations)
  perf = list(avgScores = cor(avgScores, binnedMutations, method = "spearman"),
              avgScores_pval = apply(binnedMutations,2,function(x){
                apply(avgScores,2,function(y){
                  cor.test(x,y, method = "spearman", exact = F)$p.value
                })
              }),
              sumScores = cor(sumScores, binnedMutations, method = "spearman"),
              sumScores_pval = apply(binnedMutations,2,function(x){
                apply(sumScores,2,function(y){
                  cor.test(x,y, method = "spearman", exact = F)$p.value
                })
              }))
  return(perf)
}, simplify = F)
save(perf1kb, file = "data/Modeling/WholeExomeData/perf_Allmethods_binned1kb.RData")
#####

# 100bp bins #####
perf100 = sapply(tissues, function(tissue){
  print(tissue)
  # load exome-wide prediction table
  load(paste0("data/Modeling/WholeExomeData/combinedPredictions/combinedPredictions_",
              tissue, ".RData")) #data
  data = GRanges(seqnames=data$chr,
                 ranges=IRanges(start=data$pos, width = 1),
                 mcols = data[,-(1:2)])
  # assign positions to bins
  hits  =  findOverlaps(query = bins100, subject = data)
  bin_idx <- queryHits(hits)
  
  # compute mean score for each bin
  avgScores = sapply(paste0("mcols.", methods), function(varname){
    scores <- mcols(data)[,varname][subjectHits(hits)]
    # Compute mean score per bin using tapply
    meanPerBin <- tapply(scores, bin_idx, mean)
    return(meanPerBin)
  })
  colnames(avgScores) = methods
  # compute sum of score for each bin
  sumScores = sapply(paste0("mcols.", methods), function(varname){
    scores <- mcols(data)[,varname][subjectHits(hits)]
    # Compute mean score per bin using tapply
    meanPerBin <- tapply(scores, bin_idx, sum)
    return(meanPerBin)
  })
  colnames(sumScores) = methods
  # compute number of mutations for each bin
  binnedMutations = sapply(c("mcols.cancerMuts","mcols.healthyMuts"), function(varname){
    scores <- mcols(data)[,varname][subjectHits(hits)]
    # Compute mean score per bin using tapply
    meanPerBin <- tapply(scores, bin_idx, sum)
    return(meanPerBin)
  })
  colnames(binnedMutations) = c("n_cancerMuts", "n_healthyMuts")
  if(tissue == "ovary"){binnedMutations = binnedMutations[,"n_cancerMuts", drop = F]}
  
  # binned = cbind(as.data.frame(bins1kb), avgScores, binnedMutations)
  perf = list(avgScores = cor(avgScores, binnedMutations, method = "spearman"),
              avgScores_pval = apply(binnedMutations,2,function(x){
                apply(avgScores,2,function(y){
                  cor.test(x,y, method = "spearman", exact = F)$p.value
                })
              }),
              sumScores = cor(sumScores, binnedMutations, method = "spearman"),
              sumScores_pval = apply(binnedMutations,2,function(x){
                apply(sumScores,2,function(y){
                  cor.test(x,y, method = "spearman", exact = F)$p.value
                })
              }))
  return(perf)
}, simplify = F)
save(perf100, file = "data/Modeling/WholeExomeData/perf_Allmethods_binned100bp.RData")

#####
# 1bp #####
perf1bp = sapply(tissues, function(tissue){
  print(tissue)
  # load exome-wide prediction table
  load(paste0("data/Modeling/WholeExomeData/combinedPredictions/combinedPredictions_",
              tissue, ".RData")) #data
  if(tissue == "ovary"){
    data$healthyMuts = NULL
    data = na.omit(data)
    samp = sample(1:nrow(data), size = 5000000)
    mutations = data[samp,"cancerMuts", drop = F]
  }else{
    data = na.omit(data)
    samp = sample(1:nrow(data), size = 5000000)
    mutations = data[samp,c("cancerMuts", "healthyMuts"), drop = F]
  }
  
  scores = data[samp,methods]
  
  rm(data, samp);gc()
  perf = list(scores = cor(scores, mutations, method = "spearman"),
              avgScores_pval = apply(mutations,2,function(x){
                apply(scores,2,function(y){
                  cor.test(x,y, method = "spearman", exact = F)$p.value
                })
              }))
  return(perf)
}, simplify = F)
save(perf1bp, file = "data/Modeling/WholeExomeData/perf_Allmethods_1bp.RData")

######

# visualize everything #####
load("data/Modeling/WholeExomeData/perf_Allmethods_binned100bp.RData")
load("data/Modeling/WholeExomeData/perf_Allmethods_binned1kb.RData")
load("data/Modeling/WholeExomeData/perf_Allmethods_1bp.RData")
plotPerf = function(x, tissue, names.arg = NA,p){
  dens = ifelse(p<=0.05, yes = -1, no = 20)
  barplot(x, ylim = c(0,max(x)),density = dens,
          las = 1, names.arg = names.arg, col = tissueCols[tissue], xpd = NA) #, horiz = T
}
png("fig/wholeExomePredictions/perf_cancer_binned_overview.png", 
    width = 1200, height = 1200, pointsize = 30)
par(mfrow = c(10,4), mar = c(1,2,1,1), oma = c(6,8,2,0))
dumpVar = sapply(tissues, function(tissue){
  # cancer, 1kb, avg
  temp = plotPerf(perf1kb[[tissue]]$avgScores[,"n_cancerMuts"], 
                  tissue = tissue,
                  p = perf1kb[[tissue]]$avgScores_pval[,"n_cancerMuts"])
  text(x = -3, y = max(perf1kb[[tissue]]$avgScores[,"n_cancerMuts"])/2,
       labels = t2T[tissue], adj = 1, xpd = NA)
  if(tissue == tissues[1])
    title(main="1kb, avg", line = 1, xpd = NA)
  if(tissue == tail(tissues,1))
    text(x = temp[,1], y = -0.1, labels  = methods, srt = 45, xpd = NA, adj = 1, offset = 5)
  # cancer, 1kb, sum
  plotPerf(perf1kb[[tissue]]$sumScores[,"n_cancerMuts"], tissue = tissue,
           p = perf1kb[[tissue]]$sumScores_pval[,"n_cancerMuts"])
  if(tissue == tissues[1])
    title(main="1kb, sum", line = 1, xpd = NA)
  if(tissue == tail(tissues,1))
    text(x = temp[,1], y = -0.1, labels  = methods, srt = 45, xpd = NA, adj = 1)
  # cancer, 100bp, avg
  plotPerf(perf100[[tissue]]$avgScores[,"n_cancerMuts"], tissue = tissue,
           p = perf100[[tissue]]$avgScores_pval[,"n_cancerMuts"])
  if(tissue == tissues[1])
    title(main="100bp, avg", line = 1, xpd = NA)
  if(tissue == tail(tissues,1))
    text(x = temp[,1], y = -0.1, labels  = methods, srt = 45, xpd = NA, adj = 1)
  # cancer, 100bp, sum
  plotPerf(perf100[[tissue]]$sumScores[,"n_cancerMuts"], tissue = tissue,
           p = perf100[[tissue]]$sumScores_pval[,"n_cancerMuts"])
  if(tissue == tissues[1])
    title(main="100bp, sum", line = 1, xpd = NA)
  if(tissue == tail(tissues,1))
    text(x = temp[,1], y = -0.1, labels  = methods, srt = 45, xpd = NA, adj = 1)
})
mtext("Correlation with binned mutation count", side = 2, outer = T, line = 5.5)
dev.off()

png("fig/wholeExomePredictions/perf_healthy_binned_overview.png", 
    width = 1200, height = 1200, pointsize = 30)
par(mfrow = c(10,4), mar = c(1,2,1,1), oma = c(6,8,2,0))
dumpVar = sapply(tissues, function(tissue){
  if(tissue == "ovary"){
    plot(NULL,xlim = c(0,7), ylim = c(0,1), bty = "n", xaxt = "n", yaxt = "n")
    text(x = -3, y = 0.5,
         labels = t2T[tissue], adj = 1, xpd = NA)
    plot.new();plot.new();plot.new()
    return(NA)
  }
  # cancer, 1kb, avg
  temp = plotPerf(perf1kb[[tissue]]$avgScores[,"n_healthyMuts"], 
                  tissue = tissue,
                  p = perf1kb[[tissue]]$avgScores_pval[,"n_healthyMuts"])
  text(x = -3, y = max(perf1kb[[tissue]]$avgScores[,"n_healthyMuts"])/2,
       labels = t2T[tissue], adj = 1, xpd = NA)
  if(tissue == tissues[1])
    title(main="1kb, avg", line = 1, xpd = NA)
  if(tissue == tail(tissues,1))
    text(x = temp[,1], y = -0.1, labels  = methods, srt = 45, xpd = NA, adj = 1, offset = 5)
  # cancer, 1kb, sum
  plotPerf(perf1kb[[tissue]]$sumScores[,"n_healthyMuts"], tissue = tissue,
           p = perf1kb[[tissue]]$sumScores_pval[,"n_healthyMuts"])
  if(tissue == tissues[1])
    title(main="1kb, sum", line = 1, xpd = NA)
  if(tissue == tail(tissues,1))
    text(x = temp[,1], y = -0.1, labels  = methods, srt = 45, xpd = NA, adj = 1)
  # cancer, 100bp, avg
  plotPerf(perf100[[tissue]]$avgScores[,"n_healthyMuts"], tissue = tissue,
           p = perf100[[tissue]]$avgScores_pval[,"n_healthyMuts"])
  if(tissue == tissues[1])
    title(main="100bp, avg", line = 1, xpd = NA)
  if(tissue == tail(tissues,1))
    text(x = temp[,1], y = -0.1, labels  = methods, srt = 45, xpd = NA, adj = 1)
  # cancer, 100bp, sum
  plotPerf(perf100[[tissue]]$sumScores[,"n_healthyMuts"], tissue = tissue,
           p = perf100[[tissue]]$sumScores_pval[,"n_healthyMuts"])
  if(tissue == tissues[1])
    title(main="100bp, sum", line = 1, xpd = NA)
  if(tissue == tail(tissues,1))
    text(x = temp[,1], y = -0.1, labels  = methods, srt = 45, xpd = NA, adj = 1)
})
mtext("Correlation with binned mutation count", side = 2, outer = T, line = 5.5)
dev.off()


png("fig/wholeExomePredictions/perf_cancer_binned_with1bp.png", 
    width = 1200, height = 1200, pointsize = 30)
par(mfrow = c(10,3), mar = c(1,2,1,1), oma = c(6,8,2,0))
dumpVar = sapply(tissues, function(tissue){
  print(tissue)
  # cancer, 1bp, avg
  temp = plotPerf(perf1bp[[tissue]]$scores[,"cancerMuts"], tissue = tissue,
                  p = perf100[[tissue]]$avgScores_pval[,"n_cancerMuts"])
  if(tissue == tissues[1])
    title(main="1bp", line = 1, xpd = NA)
  if(tissue == tail(tissues,1))
    text(x = temp[,1], y = par("usr")[3]-(par("usr")[4]*0.2), labels  = methods, 
         srt = 45, xpd = NA, adj = 1)
  text(x = -2.5, y = max(perf1bp[[tissue]]$scores[,"cancerMuts"])/2,
       labels = t2T[tissue], adj = 1, xpd = NA)
  # cancer, 100bp, avg
  plotPerf(perf100[[tissue]]$avgScores[,"n_cancerMuts"], tissue = tissue,
           p = perf100[[tissue]]$avgScores_pval[,"n_cancerMuts"])
  if(tissue == tissues[1])
    title(main="100bp", line = 1, xpd = NA)
  if(tissue == tail(tissues,1))
    text(x = temp[,1], y = par("usr")[3]-(par("usr")[4]*0.2), labels  = methods, 
         srt = 45, xpd = NA, adj = 1)
  # cancer, 1kb, avg
  plotPerf(perf1kb[[tissue]]$avgScores[,"n_cancerMuts"], 
           tissue = tissue,
           p = perf1kb[[tissue]]$avgScores_pval[,"n_cancerMuts"])
  if(tissue == tissues[1])
    title(main="1kb", line = 1, xpd = NA)
  if(tissue == tail(tissues,1))
    text(x = temp[,1], y = par("usr")[3]-(par("usr")[4]*0.2), labels  = methods, 
         srt = 45, xpd = NA, adj = 1, offset = 5)
})
mtext("Correlation with binned mutation count", side = 2, outer = T, line = 6)
dev.off()

png("fig/wholeExomePredictions/perf_cancer_binned_with1bp_onlyfinalMethods.png", 
    width = 1200, height = 1200, pointsize = 30)
par(mfrow = c(10,3), mar = c(1,2,1,1), oma = c(6,8,2,0))
methods2Use = c("RFpreds", "contextPreds", "mult", "combOdds")
dumpVar = sapply(tissues, function(tissue){
  print(tissue)
  # cancer, 1bp, avg
  temp = plotPerf(perf1bp[[tissue]]$scores[methods2Use,"cancerMuts"], tissue = tissue,
                  p = perf100[[tissue]]$avgScores_pval[methods2Use,"n_cancerMuts"])
  if(tissue == tissues[1])
    title(main="1bp", line = 1, xpd = NA)
  if(tissue == tail(tissues,1))
    text(x = temp[,1], y = par("usr")[3]-(par("usr")[4]*0.2), labels  = methods2Use, 
         srt = 45, xpd = NA, adj = 1)
  text(x = -1.5, y = max(perf1bp[[tissue]]$scores[,"cancerMuts"])/2,
       labels = t2T[tissue], adj = 1, xpd = NA)
  # cancer, 100bp, avg
  plotPerf(perf100[[tissue]]$avgScores[methods2Use,"n_cancerMuts"], tissue = tissue,
           p = perf100[[tissue]]$avgScores_pval[methods2Use,"n_cancerMuts"])
  if(tissue == tissues[1])
    title(main="100bp", line = 1, xpd = NA)
  if(tissue == tail(tissues,1))
    text(x = temp[,1], y = par("usr")[3]-(par("usr")[4]*0.2), labels  = methods2Use, 
         srt = 45, xpd = NA, adj = 1)
  # cancer, 1kb, avg
  plotPerf(perf1kb[[tissue]]$avgScores[methods2Use,"n_cancerMuts"], 
           tissue = tissue,
           p = perf1kb[[tissue]]$avgScores_pval[methods2Use,"n_cancerMuts"])
  if(tissue == tissues[1])
    title(main="1kb", line = 1, xpd = NA)
  if(tissue == tail(tissues,1))
    text(x = temp[,1], y = par("usr")[3]-(par("usr")[4]*0.2), labels  = methods2Use, 
         srt = 45, xpd = NA, adj = 1, offset = 5)
})
mtext("Correlation with binned mutation count", side = 2, outer = T, line = 6)
dev.off()
#####



# 