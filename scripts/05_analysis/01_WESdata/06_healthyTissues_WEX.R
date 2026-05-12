# preparation #####
.libPaths(new = "/data/public/cschmalo/R-4.1.2/")
library(ROCR)
library(ggplot2)
library(sinaplot)
library(ranger)
library(plotrix)
source("./scripts/05_analysis/00_NamesAndColors.R")
modelTissues = c("brain","breast", "colon","esophagus",
                 "kidney", "liver", "lung",# "ovary",
                 "prostate", "skin")
load("data/MutTables/SomamutDB/WEStissues.RData") #WEStissues

nThread = 14
capitalize <- function(x) {
  s <- strsplit(x, "_")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " ")
}
lightCol = function(x,alpha){
  x = col2rgb(x)
  rgb(red=x[1,], green=x[2,], blue=x[3,], alpha=alpha, maxColorValue=255)
}
dir.create("fig/healthyTissuesWEX", showWarnings = F)
plotEnding = "_20251124"
source("lib/general_function.R")
#####

# number of mutations #####
nMuts = sapply(WEStissues, function(tissue){
  print(tissue)
  # load healthy tissue trainingData
  load(paste0("data/MutTables/SomamutDB/",
              tissue, "_WES_mapped_processed.RData")) #dat,datchroms
  return(sum(dat$mutated == 1))
})  
sum(nMuts)
tissues2Remove = names(which(nMuts<100))
png(paste0("fig/healthyTissuesWEX/nMuts", plotEnding, ".png"),
    width = 1000, height = 1000, pointsize = 20)
par(mar= c(4,10,1,1))
barplot(rev(nMuts), horiz = T, las = 1, xlab = "n Mutations",
        border = (names(rev(nMuts)) %in% tissues2Remove)+1)
abline(v = 100)
dev.off()
sum(nMuts[!names(nMuts) %in% tissues2Remove])
#####

# number of mutations #####for rebuttal - removing tissues2Remove = names(which(nMuts<100)) from the plot

nMuts2Plot <- nMuts[!(names(nMuts) %in% tissues2Remove)]
png(paste0("fig/healthyTissuesWEX/nMuts_rebuttal_suppFig13A_20260422_YP.png"),
    width = 2000, height = 1000, pointsize = 20)
par(mar = c(10, 4, 1, 1), mgp = c(4, 1, 0))
bp <- barplot(nMuts2Plot, las = 2, ylab = "Number of mutations", xaxt = "n", cex.axis = 1.2, cex.lab = 1.6)
abline(h = 100)
text(x = bp, y = -max(nMuts2Plot)*0.02, labels = names(nMuts2Plot),
     srt = 45, adj = 1, xpd = TRUE, cex = 1.4)
dev.off()






# tissue-specific models #####
perfsTissueSpecific = sapply(modelTissues, function(tissue){
  print(tissue)
  # load healthy tissue trainingData
  load(paste0("data/MutTables/SomamutDB/",
              capitalize(tissue), "_WES_mapped_processed.RData")) #dat,datchroms
  # load model
  load(paste0("data/Modeling/exomeTrainData/RF/", tissue, "_finalModel.RData"))
  # predict
  yhat = predict(rf, data = dat, num.threads = nThread, type = "response")
  # compute performance
  temp = prediction(pred = yhat$predictions[,2],labels = dat$mutated)
  roc = performance(temp,  "tpr", "fpr")
  auc = performance(temp, "auc")
  pr = performance(temp, "prec", "rec")
  return(list(roc = roc, pr = pr, auc = auc))
}, simplify = F)
save(perfsTissueSpecific, file="data/Modeling/healthyTissues/perfsTissueSpecific.RData")
#####

# allTissue model #####
# load model
load(paste0("data/Modeling/exomeTrainData/RF/generalModel_", 
            "finalModel.RData"))
perfsAllTissueModel = sapply(WEStissues, function(tissue){
  print(tissue)
  # load healthy tissue trainingData
  load(paste0("data/MutTables/SomamutDB/",
              tissue, "_WES_mapped_processed.RData")) #dat,datchroms
  
  # predict
  yhat = predict(rf, data = dat, num.threads = nThread, type = "response")
  # compute performance
  temp = prediction(pred = yhat$predictions[,2],labels = dat$mutated)
  roc = performance(temp,  "tpr", "fpr")
  auc = performance(temp, "auc")
  pr = performance(temp, "prec", "rec")
  return(list(roc = roc, pr = pr, auc = auc))
}, simplify = F)
save(perfsAllTissueModel, file="data/Modeling/healthyTissues/perfsAllTissueModel.RData")
#####


# load/compute training performance of general and tissue-specific models #####
load("data/Modeling/exomeTrainData/RF/generalModel__predictions.RData") # predictions
GeneralPredConcat = do.call(rbind,predictions)
temp = prediction(GeneralPredConcat$pred, GeneralPredConcat$label)
GeneralSelfPerf = list(roc = performance(temp, "tpr", "fpr"),
                  pr = performance(temp,"prec", "rec"),
                  auc = performance(temp,"auc")@y.values[[1]])

load("data/Modeling/exomeTrainData/RF/ROC_PR_RF_concat.RData")

#####

# plotting #####
load("data/Modeling/healthyTissues/perfsTissueSpecific.RData")
load("data/Modeling/healthyTissues/perfsAllTissueModel.RData")
# ROC, and PR comparing for each tissue: 
# general model, training perf of general model, 
# tissue-specific model (when available), perf on tissue-specific training
# ROC
png(paste0("fig/healthyTissuesWEX/ROC", plotEnding, ".png"),
    width = 2000, height = 780, pointsize = 20)
par(mfrow = c(3,10), mar = c(1,1,2,1), oma = c(3,3,0,0))
dumpVar = sapply(1:length(WEStissues), function(i){
  tissue = WEStissues[i]
  plot(NA, xlim = c(0,1), ylim = c(0,1),
       xlab = "", ylab = "", las = 1, xaxt = "n", yaxt = "n", main = tissue)
  abline(0,1, col = "dark grey", lty = 2, lwd = 1.5)
  if((i-1)%%(par()$mfrow[2]) == 0){
    axis(2, las = 1)
  }
  if(i>((par()$mfrow[1]-1)*(par()$mfrow[2]))){
    axis(1, las = 1)
  }
  # general model
  lines(perfsAllTissueModel[[tissue]]$roc@x.values[[1]],
        perfsAllTissueModel[[tissue]]$roc@y.values[[1]], lwd = 1.5)
  # general model training performance
  lines(GeneralSelfPerf$roc@x.values[[1]],
        GeneralSelfPerf$roc@y.values[[1]], lty = 2, lwd = 1.5)
  # tissue-specific model
  if(tolower(tissue) %in% modelTissues){
    lines(perfsTissueSpecific[[tolower(tissue)]]$roc@x.values[[1]],
          perfsTissueSpecific[[tolower(tissue)]]$roc@y.values[[1]],
          col = tissueCols[tolower(tissue)], lwd = 1.5)
    lines(ROC_PR_RF_concat[[tolower(tissue)]]$roc@x.values[[1]],
          ROC_PR_RF_concat[[tolower(tissue)]]$roc@y.values[[1]], 
          col = tissueCols[tolower(tissue)], lty = 2, lwd = 1.5)
  }
})
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
legend("topleft", lty = c(1,2,1,2), col = c("black", "black", "red", "red"),xpd = NA, lwd = 2,
      legend = c("General model", "Training performance\ngeneral model", 
                 "Tissue-specific model", "Training performance\ntissue-specific model"))
mtext(text = "FPR", side = 1, outer = T, line = 1.9)
mtext(text = "TPR", side = 2, outer = T, line = 1.9)
dev.off()

# PR
png(paste0("fig/healthyTissuesWEX/PR", plotEnding, ".png"),
    width = 2000, height = 800, pointsize = 20)
par(mfrow = c(3,10), mar = c(1,1,2,1), oma = c(3,3,0,0))
dumpVar = sapply(1:length(WEStissues), function(i){
  tissue = WEStissues[i]
  plot(NA, xlim = c(0,1), ylim = c(0,1),
       xlab = "", ylab = "", las = 1, xaxt = "n", yaxt = "n", main = tissue)
  if((i-1)%%(par()$mfrow[2]) == 0){
    axis(2, las = 1)
  }
  if(i>((par()$mfrow[1]-1)*(par()$mfrow[2]))){
    axis(1, las = 1)
  }
  # general model
  lines(perfsAllTissueModel[[tissue]]$pr@x.values[[1]],
        perfsAllTissueModel[[tissue]]$pr@y.values[[1]], lwd = 1.5)
  # general model training performance
  lines(GeneralSelfPerf$pr@x.values[[1]],
        GeneralSelfPerf$pr@y.values[[1]], lty = 2, lwd = 1.5)
  # tissue-specific model
  if(tolower(tissue) %in% modelTissues){
    lines(perfsTissueSpecific[[tolower(tissue)]]$pr@x.values[[1]], lwd = 1.5,
          perfsTissueSpecific[[tolower(tissue)]]$pr@y.values[[1]], 
          col = tissueCols[tolower(tissue)])
    lines(ROC_PR_RF_concat[[tolower(tissue)]]$pr@x.values[[1]], lwd = 1.5,
          ROC_PR_RF_concat[[tolower(tissue)]]$pr@y.values[[1]], 
          col = tissueCols[tolower(tissue)], lty = 2)
  }
})
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
legend("topleft", lty = c(1,2,1,2), col = c("black", "black", "red", "red"),xpd = NA, lwd = 2,
       legend = c("General model", "Training performance general model", 
                  "Tissue-specific model", "Training performance tissue-specific model"))
mtext(text = "Recall", side = 1, outer = T, line = 1.9)
mtext(text = "Precision", side = 2, outer = T, line = 1.9)
dev.off()


# AUC general model
pdfAndPng(file = paste0("fig/healthyTissuesWEX/AUC", plotEnding), width = 4, 
          height = 10, pngArgs = list(pointsize=15), pdfArgs = list(pointsize=12),
          expr = expression({
            AUCs = rbind("General" = sapply(perfsAllTissueModel, function(x)x$auc@y.values[[1]]),
                         "Tissue-specific" = sapply(perfsTissueSpecific, function(x)x$auc@y.values[[1]])[tolower(WEStissues)])
            par(mar = c(3,7,3,0.5))
            AUCs = AUCs[,ncol(AUCs):1]
            offset = 0.4
            AUCs[is.na(AUCs)] = offset
            AUCs = AUCs[,!colnames(AUCs) %in% tissues2Remove]
            temp = barplot(AUCs-offset, 
                           beside = T, horiz = T, las = 1, legend.text = F, 
                           col = rbind("grey20",  tissueCols[tolower(colnames(AUCs))]), 
                           density = c(NA, 40),border = T,
                           names.arg = sapply((colnames(AUCs)), capitalize), xaxt = "n", mgp = c(1.9,0.4,0))
            
            abline(v=0.5-offset, lty = 2)
            tickPos =axTicks(1)
            axis(1,at = tickPos, labels = c(0,tickPos[-1]+offset))
            axis.break(axis = 1, breakpos = 0.01, style = "gap") #, bgcol = "grey"
            abline(v=0)
            axis(2, at = colMeans(temp), labels = rep("", ncol(temp)),lwd = 0, lwd.ticks = 1, pos = 0, tck=-0.02)
            title(xlab = "AUC", mgp = c(1.8,0,0))
            legend(x = par()$usr[1], y = max(temp)*1.05, xjust = 0.3, yjust = 0,xpd = NA, ncol = 2, 
                   legend = c("Tissue-specific model", "Universal all-tissue model"), 
                   fill = c("grey20", "grey20"), density = c(40, NA))
          }))
pdfAndPng(file = paste0("fig/healthyTissuesWEX/AUC_vertical", plotEnding), width = 10, 
          height = 4, pngArgs = list(pointsize=15), pdfArgs = list(pointsize=12),
          expr = expression({
            AUCs = rbind("General" = sapply(perfsAllTissueModel, function(x)x$auc@y.values[[1]]),
                         "Tissue-specific" = sapply(perfsTissueSpecific, function(x)x$auc@y.values[[1]])[tolower(WEStissues)])
            par(mar = c(7,4,4,0.5))
            # AUCs = AUCs[,ncol(AUCs):1]
            offset = 0.4
            AUCs[is.na(AUCs)] = offset
            AUCs = AUCs[,!colnames(AUCs) %in% tissues2Remove]
            
            temp = barplot(AUCs-offset, 
                           beside = T,  las = 2, legend.text = F, 
                           col = rbind("grey20",  tissueCols[tolower(colnames(AUCs))]), 
                           density = c(40,NA),border = T,
                           names.arg = sapply((colnames(AUCs)), capitalize), yaxt = "n", mgp = c(1.9,0.4,0))
            abline(h=0.5-offset, lty = 2)
            tickPos =axTicks(2)
            axis(2,at = tickPos, labels = c(0,tickPos[-1]+offset), las = 1)
            axis.break(axis =2, breakpos = mean(tickPos[1:2]), style = "gap") #, bgcol = "grey"
            abline(h=0)
            axis(1, at = colMeans(temp), labels = rep("", ncol(temp)),lwd = 0, lwd.ticks = 1, pos = 0, tck=-0.02)
            title(ylab = "AUC", mgp = c(2,0,0))
            legend(x = 0, y = 0.675-offset, xpd = NA,
                   legend = c("Universal all-tissue model", "Tissue-specific model"), 
                   fill = c("grey20", "grey20"), density = c( 40, NA))
          }))

# AUCs vs nMuts
nMutsHealthy = sapply(WEStissues, function(tissue){
  load(paste0("data/MutTables/SomamutDB/",
              tissue, "_WES_mapped_processed.RData")) #dat,datchroms
  return(nrow(dat))
})
pdfAndPng(paste0("fig/healthyTissuesWEX/AUCvsNmuts", plotEnding),
          width = 6, height = 6, 
          pngArgs = list(pointsize = 15), pdfArgs = list(pointsize=10),
          expr = expression({
            plot(nMutsHealthy, AUCs["General", names(nMutsHealthy)], 
                 ylab = "AUC", xlab = "Test data size", las = 1)
            points(nMutsHealthy, AUCs[2, names(nMutsHealthy)], 
                   col = tissueCols[tolower(names(nMutsHealthy))], pch = 19)
            legend("bottomright", c("General model", "Tissue-specific model"), col = c("black", "grey30"), pch = c(1,19))
          }))
# png(paste0("fig/healthyTissuesWEX/AUCvsNmuts", plotEnding, ".png"),
#     width = 600, height = 600, pointsize = 15)
# 
# dev.off()
####
