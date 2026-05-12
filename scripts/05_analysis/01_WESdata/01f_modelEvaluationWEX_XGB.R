library(readxl)
source("scripts/05_analysis/00_NamesAndColors.R")
library(ROCR)
library(ranger) 
library(sinaplot)
library(berryFunctions) # for smallPlot

plotEnding = "20260412"
methodCols = setNames(RColorBrewer::brewer.pal(4, "Dark2"), c("RF", "GLM", "SL", "XGB"))

load("data/Modeling/exomeTrainData/RF/predPerTissueRF.RData")
load("data/Modeling/exomeTrainData/GLM/predPerTissueGLM_sig.RData")
load("data/Modeling/exomeTrainData/Lasso/predPerTissueLasso.RData")
load("data/Modeling/exomeTrainData/RF/ROC_PR_RF_perChr.RData")
load("data/Modeling/exomeTrainData/GLM/ROC_PR_glm_perChr_sig.RData")
load("data/Modeling/exomeTrainData/Lasso/ROC_PR_lasso_perChr.RData")
load("data/Modeling/exomeTrainData/RF/ROC_PR_RF_concat.RData")
load("data/Modeling/exomeTrainData/GLM/ROC_PR_glm_concat_sig.RData")
load("data/Modeling/exomeTrainData/Lasso/ROC_PR_lasso_concat.RData")


# 
predPerTissueXGB = sapply(tissues, function(tissue){
  load( paste0("data/Modeling/exomeTrainData/xgboost/", tissue,
               "_predictions.RData"))
  return(predictions)
}, simplify=F)
ROC_PR_XGB_perChr = sapply(names(predPerTissueXGB), function(tissue){ #iterate through tissues
  pred = predPerTissueXGB[[tissue]]
  # get performances
  res = lapply(pred, function(x){ #iterate through chromosomes
    perf = prediction(x[,"pred"], x[,"label"])
    roc = performance(perf, "tpr", "fpr")
    pr = performance(perf,"prec", "rec")
    auc = performance(perf,"auc")@y.values[[1]]
    return(list(roc = roc, pr = pr, auc = auc))
  })
  return(res)
}, simplify=F)
ROC_PR_XGB_concat = sapply(names(predPerTissueXGB), function(tissue){
  predConcat = do.call(rbind,predPerTissueXGB[[tissue]])
  ROC_rf = performance(prediction(predConcat[,"pred" ], 
                                  predConcat[,"label"]), 
                       "tpr", "fpr")
  AUC_rf = performance(prediction(predConcat[,"pred" ], 
                                  predConcat[,"label"]), 
                       "auc")
  PR_rf = performance(prediction(predConcat[,"pred" ], 
                                 predConcat[,"label"]), 
                      "prec", "rec")
  return(list(roc = ROC_rf, pr = PR_rf, auc = AUC_rf))
}, simplify=F)

png(paste0("fig/modelEvaluation/compareMethodperformance_withXGB_", plotEnding,".png"), 
    height=1600, width=1300, pointsize=30)#, res = 200)
par(mfrow = c(length(tissues),3), mar = c(0.5,6,0.1,0.1), 
    mgp = c(2.5,1,0), oma = c(3,6,0.05,0.05))
lwd = 2
plotDump = sapply(tissues, function(tissue){
  print(tissue)
  # auroc
  AUCs = data.frame(RF = sapply(ROC_PR_RF_perChr[[tissue]],function(x){x$auc}),
                    GLM = sapply(ROC_PR_glm_perChr_sig[[tissue]],function(x){x$auc}),
                    SL = sapply(ROC_PR_lasso_perChr[[tissue]],function(x){x$auc}),
                    XBG = sapply(ROC_PR_XGB_perChr[[tissue]],function(x){x$auc}))
  # boxplot(AUCs, las = 1, outwex = 1.5,
  #         border = methodCols, xaxt = "n", yaxt = "n")
  sinaplot(AUCs, adjust =0.1, maxwidth = 0.5, col = methodCols, xaxt = "n", yaxt = "n")
  axis(2, mgp = c(2,0.7,0), las = 1,cex.axis = 1.4,
       at = par("yaxp")[1:2], lwd = 0, lwd.ticks=1)
  AUCs = c(RF = ROC_PR_RF_concat[[tissue]]$auc@y.values[[1]],
           GLM = ROC_PR_glm_concat_sig[[tissue]]$auc@y.values[[1]],
           SL = ROC_PR_lasso_concat[[tissue]]$auc@y.values[[1]],
           XGB = ROC_PR_XGB_concat[[tissue]]$auc@y.values[[1]])
  points(AUCs, col = methodCols, pch = 19, cex = 2)
  mtext(side=2, text=t2T[tissue], line = 5, las = 1)
  # add auroc axis labels
  if(tissue == tail(tissues,1)){   
    # mtext(text=names(methodCols), side=1,  at=1:3, cex.lab = 1.4)
    axis(1,at = 1:4, labels=names(methodCols), mgp = c(2,0.7,0),
         cex.axis=1.2)
  }
  if(tissue == tissues[ceiling(length(tissues)/2)]){
    # title(ylab = "AUC", xpd = NA, mgp = c(3,1,0), cex.lab = 1.4)
    mtext(text = "AUC", side = 2, adj = -0.5,cex.lab = 1.4, line = 3)
  }
  #prc
  plot(ROC_PR_RF_concat[[tissue]]$pr@x.values[[1]],
       ROC_PR_RF_concat[[tissue]]$pr@y.values[[1]], lwd = 3, 
       las = 1,xlab = "", xaxt = "n", type = "l", col= methodCols["RF"],
       ylab = "",  ylim = c(0,1), mgp =  c(2,0.7,0),  xaxt = "n",yaxt = "n")
  axis(2, mgp = c(2,0.7,0), las = 1,
       at = c(0,0.5), lwd = 0, lwd.ticks=1, cex.axis = 1.4)
  axis(2, mgp = c(2,0.7,0), las = 1,
       at = 1, padj = 0.7, cex.axis = 1.4)
  plot(ROC_PR_glm_concat_sig[[tissue]]$pr, lwd = 3, 
       col = methodCols["GLM"], add = T, lty = 2)
  plot(ROC_PR_lasso_concat[[tissue]]$pr, lwd = 3, 
       col = methodCols["SL"], add = T, lty = "4414")
  plot(ROC_PR_XGB_concat[[tissue]]$pr, lwd = 3, 
       col = methodCols["XGB"], add = T, lty = "4414")
  rect(xleft = -0.005, xright=0.02, ybottom=0.5, ytop=1, lwd = 1.5)
  # add pr axis labels
  if(tissue == tail(tissues,1)){   
    axis(1, mgp = c(2,0.7,0), las = 1,
         at = c(0,0.5), lwd = 0, lwd.ticks=1, cex.axis = 1.4)
    axis(1, mgp = c(2,0.7,0), las = 1,at = 1, cex.axis = 1.4)
    title(xlab = "Recall", line = 2, xpd = NA, cex.lab = 1.4)
  }  
  if(tissue == tissues[ceiling(length(tissues)/2)]){
    mtext(text = "Precision", side = 2, adj = -5,cex.lab = 1.4, line = 3)
  }
  
  # create inset with zoom
  lim = par("plt")
  xlims = lim[2]-lim[1]
  ylims = lim[4]-lim[3]
  inset = 0.05
  smallPlot(expr={
    plot(x = ROC_PR_RF_concat[[tissue]]$pr@x.values[[1]],
         y = ROC_PR_RF_concat[[tissue]]$pr@y.values[[1]],
         xlab = "", ylab = "", xaxt = "n", yaxt = "n",
         xlim = c(0,0.02),ylim = c(0.5,1), 
         lwd = 3, type = "l", col = methodCols["RF"])
    lines(x = ROC_PR_glm_concat_sig[[tissue]]$pr@x.values[[1]],
          y = ROC_PR_glm_concat_sig[[tissue]]$pr@y.values[[1]],
          col = methodCols["GLM"], lty = 2, lwd = 3)
    lines(x = ROC_PR_lasso_concat[[tissue]]$pr@x.values[[1]],
          y = ROC_PR_lasso_concat[[tissue]]$pr@y.values[[1]],
          col = methodCols["SL"], lty = 2, lwd = 3)
    lines(x = ROC_PR_XGB_concat[[tissue]]$pr@x.values[[1]],
          y = ROC_PR_XGB_concat[[tissue]]$pr@y.values[[1]],
          col = methodCols["XGB"], lty = 2, lwd = 3)
    box(lwd = 2)},  
    x1 = lim[1]+0.3*xlims, x2 = lim[2]-0.05*xlims,
    y1 = lim[3]+0.05*ylims, y2 = lim[4]-0.3*ylims, xpd = F,
    mar = c(0,0,0,0), border = "transparent")
  
  
  # violin of predictions
  rfpreds = do.call(rbind,predPerTissueRF[[tissue]])
  glmpreds = do.call(rbind,predPerTissueGLMsig[[tissue]])
  lassopreds =  do.call(rbind,predPerTissueLasso[[tissue]])
  xgbpreds = do.call(rbind,predPerTissueXGB[[tissue]])
  rfpredssplit = split(rfpreds$pred, rfpreds$label)
  glmpredssplit = split(glmpreds$pred, glmpreds$label)
  lassopredssplit = split(lassopreds$pred, lassopreds$label)
  xgbpredssplit = split(xgbpreds[,"pred"], xgbpreds[,"label"])
  
  plotDat = list("TN.RF" = rfpredssplit$`0`,
                 "TP.RF" = rfpredssplit$`1`,
                 "TN.GLM" = glmpredssplit$`0`,
                 "TP.GLM" = glmpredssplit$`1`,
                 "TN.SL" = lassopredssplit$`0`,
                 "TP.SL" = lassopredssplit$`1`,
                 "TN.XGB" = xgbpredssplit$`0`,
                 "TP.XGB" = xgbpredssplit$`1`)
  tempCols = RColorBrewer::brewer.pal(4,"Set2")
  sinaplot(plotDat, col = c(tempCols,methodCols)[c(1,5,2,6,3,7,4,8)], 
           xaxt = "n",yaxt = "n",las = 1,  ylab = "", cex = 0.8, ylim= c(0,1))
  axis(2, mgp = c(2,0.7,0), las = 1,
       at = c(0,0.5), lwd = 0, lwd.ticks=1, cex.axis = 1.4)
  axis(2, mgp = c(2,0.7,0), las = 1,
       at = 1, padj = 0.7, cex.axis = 1.4)
  abline(h=0.5, col = "grey", lty = 2, lwd = 2)
  abline(v=c(2.5,4.5))
  boxplot(plotDat, col = c("grey80", "grey40"),
          add = T, boxwex = 0.2, outline = F, #col = rgb(0,0,0,alpha = 0),
          ann = F, yaxt = "n", xaxt = "n", mgp = c(2,0.7,0))
  if(tissue == tail(tissues,1)){   
    axis(1,at = 1:8, labels=rep(c("0","1"),4), mgp = c(2,0.7,0),
         cex.axis=1.2)
    title(xlab = "True labels", line = 2, xpd = NA, 
          mgp = c(2,0.7,0), cex.lab = 1.4)
    # axis(1,at = c(1.5,3.5,5.5), labels=names(methodCols), mgp = c(2,0.7,0),
    #      cex.axis=1.2, tick = F)
    text(x= c(1.5,3.5,5.5,7.5), y=0.06,labels=names(methodCols), cex = 1.2, font = 2)
  }  
  if(tissue == tissues[ceiling(length(tissues)/2)]){
    mtext(text = "Prediction", side = 2, adj = -15, cex.lab = 1.4, line = 3)
  }
})
dev.off() 


# xgb_imp_cwcv = sapply(tissues, function(tissue){
#   load(paste0("data/Modeling/exomeTrainData/xgboost/", tissue,
#               "_importances_gain.RData"))
#   names(imp) = names(chrCols)
#   imp = as.data.frame(imp)
#   temp = imp[names(p2P),]
#   rownames(temp) = names(p2P)
#   return(temp)
# }, simplify=F)
xgb_imps = do.call(rbind,sapply(tissues, function(tissue){
  xgbModel = xgb.load(paste0("data/Modeling/exomeTrainData/xgboost/", tissue, "_finalModel", ".ubj"))
  imp = as.data.frame(xgb.importance(model = xgbModel))
  rownames(imp) = imp$Feature
  # imp = imp[names(p2P),]
  # imp$Feature = factor(imp$Feature, levels = imp$Feature)
  
  imp$tissue = tissue 
  imp$tissue = t2T[tissue]
  imp$group = p2G[imp$Feature]
  # importance$Feature = factor(imp$Feature, levels = imp$Feature)
  return(imp)
}, simplify=F))
xgb_imps$group = factor(p2G[xgb_imps$Feature], levels = rev(unique(p2G)))
xgb_imps$Feature = factor(xgb_imps$Feature, levels = names(p2P))

ggplot(xgb_imps, aes(x = tissue, y = Feature, fill = Gain)) + 
  geom_raster() +
  scale_fill_gradient(low="grey90", high="red",na.value="grey") +
  facet_grid(rows = vars(group), space = "free_y", scales = "free_y", switch = "y") +
  # labs(y = "Predictor", x = "Tissue") + 
  labs(fill = "Gain")+
  scale_x_discrete(labels=t2T) +
  scale_y_discrete(labels = p2P)+
  theme(axis.text.y=element_text(size=4.5),
        axis.text.x=element_text(size=8 , angle = 45,vjust = 1, hjust=1),
        axis.title=element_text(size=4,face="bold"),
        strip.text.y.left = element_text(angle = 0, size=5.5),
        strip.placement = "outside",
        panel.spacing = unit(0.1, "lines")) 
ggsave(paste0("fig/modelEvaluation/XGBgain_Tissues", plotEnding, ".png"), 
       height=8, width=6)
ggsave(paste0("fig/modelEvaluation/XGBgain_Tissues", plotEnding, ".pdf"), 
       height=8, width=6)
