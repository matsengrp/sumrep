devtools::load_all()

# Run partis for each data set and save annotation R objects
source("Examples/BCRAnalyses/RunPartisAnnotations.R")

# Run IgBlast for each data set and save annotation R objects
source("Examples/BCRAnalyses/RunIgBlastAnnotations.R")

# Run partis using IgBlast's germline database for each data set and save 
#   annotation R objects
source("Examples/BCRAnalyses/RunPartisIgBlastAnnotations.R")

# Load each annotation object and run compareRepertoires for each pair
source("Examples/BCRAnalyses/CompareBCRDatasets.R")

# Make some plots based on the above comparisons
source("Examples/BCRAnalyses/BCRSummaryPlots.R")
