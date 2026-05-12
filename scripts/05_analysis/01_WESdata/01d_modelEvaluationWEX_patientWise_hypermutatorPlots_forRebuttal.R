# Modified mutation outlier plots: 96-type spectrum with Plot 1 styling
# Complete self-contained script including data preparation
setwd('/cellfile/cellnet/MutationModel/')
graphics.off()
# preparation #####
library(ROCR)
source("scripts/05_analysis/00_NamesAndColors.R")
#plotEnding = "_20230916"
source("lib/general_function.R")
basesCol = c("#14baeb","#030303","#df2b27","#999999","#a2ca61","#ebc7c3") # Same as used by Alexandrov et al.
tissues <- c("brain", "breast", "colon", "esophagus", "kidney", 
             "liver", "lung", "ovary", "prostate", "skin")

bases = c("A", "C", "G", "T")
mutTypes = cbind(paste0(paste0(rep(rep(bases, each = 4), 6),"[",
                               c(rep("C", 48), rep("T", 48)), ">",
                               c(rep(c("A", "G", "T"), each =16),
                                 rep(c("A", "C", "G"), each =16)), "]",
                               rep(bases, 24))),
                 paste0(paste0(rep(rep(rev(bases), each = 4),  6),"[",
                               c(rep("G", 48), rep("A", 48)), ">",
                               c(rep(c("T", "C", "A"), each =16),
                                 rep(c("T", "G", "C"), each =16)), "]",
                               rep(rev(bases), 24))))
rownames(mutTypes) = paste(mutTypes[,1])
mutTypeTranslator= setNames(nm=c(mutTypes[,1], mutTypes[,2]),
                            object=c(rownames(mutTypes), rownames(mutTypes)))

tissue2Cancer = list("lung" = "Lung adenocarcinoma",
                     "breast" = "Breast invasive carcinoma",
                     "skin" = "Skin Cutaneous Melanoma",
                     "colon" = c("Colon adenocarcinoma", "Rectum adenocarcinoma"),
                     "ovary" = "Ovarian serous cystadenocarcinoma",
                     "kidney" = "Kidney renal clear cell carcinoma",
                     "prostate" = "Prostate adenocarcinoma",
                     "esophagus" = "Esophageal carcinoma",
                     "liver" = "Liver hepatocellular carcinoma",
                     "brain" = "Brain Lower Grade Glioma")
cancer2Tissue = do.call(rbind,sapply(names(tissue2Cancer), function(tissue){
  cbind(tissue, tissue2Cancer[[tissue]])
}))
cancer2Tissue = setNames(object = cancer2Tissue[,1], nm = cancer2Tissue[,2])
chrs = paste0("chr", c(1:22))
signatures = read.table("data/rawdata/COSMIC_catalogue-signatures_SBS96_v3.4/COSMIC_v3.4_SBS_GRCh37.txt", 
                        header = T, row.names = 1)
signatures = signatures[mutTypes[,1],]
#####


# get table with relevant mutation data #####
load("data/MutTables/exomeTrainData/muts.RData")
mutData = muts[,c("Chromosome", "Start_Position","Reference_Allele", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode", "CONTEXT","cancerType")]
mutData$patientBarcodes = substr(mutData$Tumor_Sample_Barcode, 1,12)
mutData$mutType = paste0(substr(mutData$CONTEXT,5,5),
                         "[",mutData$Reference_Allele, ">",mutData$Tumor_Seq_Allele2, "]",
                         substr(mutData$CONTEXT,7,7))
mutData$tissue = cancer2Tissue[mutData$cancerType]
mutData$position = paste0(mutData$Chromosome, "_",mutData$Start_Position)
#####


# get patient information #####
# meta = read.table("~/Documents/clinical_PANCAN_patient_with_followup.tsv",
#                   header = T, sep = "\t", quote = "")
# toRemove = sapply(meta, function(x){
#   all(x =="" | x == "[Not Applicable]" | x == "[Not Available]", na.rm = T) })
# meta = meta[,!toRemove]
# subMeta = meta[meta$bcr_patient_barcode %in% mutData$patientBarcodes,c("bcr_patient_barcode", "days_to_birth", "radiation_therapy")]
# rownames(subMeta) = subMeta[,1]
#####


# identify outlier samples #####
topPercent = sapply(tissues, function(tissue){
  temp = table(mutData$patientBarcodes[mutData$tissue == tissue])
  temp = sort(temp, decreasing = T)
  topPerc = temp[1:5]/sum(temp)*100
  topPerc
}, simplify = F)

outliers = do.call(rbind,sapply(tissues, function(tissue){
  temp = topPercent[[tissue]]
  res = temp[temp>=5]
  if(length(res)>0){
    return(cbind(tissue, names(res), res))
  } else{
    NULL
  }
}))
#outlierInfo  = cbind(outliers, subMeta[outliers[,2],])
#outlierInfo
#####


# Plot: mutation outliers with Plot 1-style 96-type spectrum #####
barCols <- rep(rev(basesCol), each = 16)

dumpVar = sapply(unique(outliers[,"tissue"]), function(tissue){
  subDat = mutData[mutData$tissue == tissue,]
  subSig = factor(mutTypeTranslator[subDat$mutType], levels = mutTypes[,1])
  
  patients = outliers[outliers[,1] == tissue, 2]
  nPanels = length(patients) * 3 + 1
  
  png(paste0("fig/For_rebuttal_mutationOutliers_", tissue, ".png"), 
      width = 200 * nPanels, height = 1000, pointsize = 15)
  
  #png(paste0("fig/For_rebuttal_mutationOutliers_", tissue, ".png"), 
  #   width = 1400, height = 1000, pointsize = 15)
  
  #par(mfrow = c(1, nPanels), mar = c(2, 0.5, 2, 0), oma = c(2, 6, 0, 0))
  
  
  # Tissue-specific cex to compensate for different nPanels
  tissueCex = list(
    "brain"     = 0.9,
    "breast"    = 0.9,
    "colon"     = 0.6,
    "esophagus" = 0.6,
    "ovary"     = 0.6,
    "prostate"  = 0.6,
    "skin"      = 0.6
  )
  cexVal = if(tissue %in% names(tissueCex)) tissueCex[[tissue]] else 1.0
  
  par(mfrow = c(1, nPanels), mar = c(2, 0.5, 2, 0), oma = c(2, 6, 0, 0), cex = cexVal)
  
  
  
  # ---- Panel 1: All samples ----
  par(mar = c(2, 0.5, 2, 0))
  barplotTemp = barplot(rev(table(subSig)), las = 1, horiz = T,
                        main = "All samples", names.arg = NA, 
                        col = barCols, space = 0, yaxs = "i",
                        border = NA)
  abline(h = seq(16, 80, by = 16), lty = 2, lwd = 1.5)
  abline(h = 0, lty = 1, lwd = 2)
  abline(h = 96, lty = 1, lwd = 1)
  
  # Formatted trinucleotide labels (matching Plot 1 style)
  classes = mutTypes[, 1]  # e.g., "A[C>A]C"
  names_3letter = rev(paste0(substr(classes, 1, 1), 
                             substr(classes, 3, 3), 
                             substr(classes, 7, 7)))
  names_center = rev(paste0(" ", substr(classes, 3, 3), " "))
  
  ypos = seq(0.5, 95.5, by = 1)
  
  text(x = 0, y = ypos, labels = names_3letter, family = "mono",
       xpd = NA, adj = 1.1, cex = 0.95)
  text(x = 0, y = ypos, labels = names_center, family = "mono",
       xpd = NA, adj = 1.1, col = barCols, font = 2, cex = 0.95)
  
  # Mutation class labels on far left
  classLabels = c("T>G", "T>C", "T>A", "C>T", "C>G", "C>A")
  classYpos = seq(8, 88, by = 16)
  
  
  # Tissue-specific NDC positions (adjust values as needed per tissue)
  ndcOffsets = list(
    "colon"    = c(0.05, 0.06),
    "skin"     = c(0.05, 0.06),
    "breast"   = c(0.04, 0.05),
    "prostate" = c(0.05, 0.06),
    "brain"    = c(0.04, 0.05),
    "ovary"    = c(0.05, 0.06),
    "esophagus"= c(0.05, 0.06)
    # add more tissues as needed
  )
  # Fallback default
  offsets = if(tissue %in% names(ndcOffsets)) ndcOffsets[[tissue]] else c(0.04, 0.06)
  
  xLeft = grconvertX(offsets[1], from = "ndc", to = "user")
  xSeg  = grconvertX(offsets[2], from = "ndc", to = "user")
  
  text(x = xLeft, y = classYpos, labels = classLabels,
       xpd = NA, adj = 1, col = rev(basesCol), font = 2, cex = 1.2/cexVal)
  segments(y0 = seq(0, 80, by = 16), y1 = seq(16, 96, by = 16), 
           x0 = xSeg, lwd = 4, col = rev(basesCol), xpd = NA, lend = 1)
  
  
  
  
  
  title(ylab = "Mutation type", mgp = c(5, 1, 0), xpd = NA, cex.lab = 2)
  
  # ---- Patient panels ----
  patientMutCounts = sapply(patients, function(patient){
    patientMuts = factor(
      mutTypeTranslator[subDat$mutType[subDat$patientBarcodes == patient]], 
      levels = mutTypes[,1])
    
    barplot(rev(table(patientMuts)), las = 1, horiz = T, yaxt = "n",
            main = patient, col = barCols, space = 0, yaxs = "i",
            border = NA, names.arg = NA)
    abline(h = seq(16, 80, by = 16), lty = 2, lwd = 1.5)
    abline(h = 0, lty = 1, lwd = 2)
    
    # Most correlated COSMIC signature
    mostCorrSign = cor(table(patientMuts), signatures)
    maxCorr = which.max(mostCorrSign)
    
    barplot(rev(signatures[, maxCorr]), las = 1, horiz = T, yaxt = "n",
            main = paste0("\n", colnames(signatures)[maxCorr], 
                          "\nr = ", format(mostCorrSign[maxCorr], digits = 2)),
            col = barCols, space = 0, yaxs = "i",
            border = NA, names.arg = NA)
    abline(h = seq(16, 80, by = 16), lty = 2, lwd = 1.5)
    abline(h = 0, lty = 1, lwd = 2)
    
    # Without this patient
    nonPatientMuts = factor(
      mutTypeTranslator[subDat$mutType[subDat$patientBarcodes != patient]], 
      levels = mutTypes[,1])
    
    fadedCols = sapply(barCols, function(col) adjustcolor(col, alpha.f = 0.5))
    
    barplot(rev(table(nonPatientMuts)), las = 1, horiz = T, yaxt = "n",
            main = paste0("\nWithout\n", patient), 
            col = fadedCols, space = 0, yaxs = "i",
            border = NA, names.arg = NA)
    abline(h = seq(16, 80, by = 16), lty = 2, lwd = 1.5)
    abline(h = 0, lty = 1, lwd = 2)
    
    return(table(patientMuts))
  })
  
  mtext("Count", side = 1, outer = T, line = 1, cex = 1.5)
  dev.off()
})
#####