source("scripts/05_analysis/01_WESdata/prep_data_for_python_predictorPlots_forRebuttal.R")
source("scripts/05_analysis/00_NamesAndColors.R") #  p2P mapping
load("data/processedData/dataInfos.RData")

# Export all tissues to CSV
export_cors_for_python(
  dataInfos = dataInfos,
  tissues = tissues,
  p2P = p2P,
  output_dir = "fig/modelEvaluation/WES_predictor_cors_csv"
)


##from command line then run for each tissue:
#python3 /cellfile/cellnet/MutationModel/scripts/05_analysis/01_WESdata/corrplot_complexheatmap_forRebuttal.py fig/modelEvaluation/WES_predictor_cors_csv/cors_skin.csv --tissue 'Skin'
#--output fig/modelEvaluation/