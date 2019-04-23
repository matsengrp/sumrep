devtools::load_all()

# Run partis for each data set and save annotation R objects
source("Examples/BCRAnalyses/RunPartisAnnotations.R")

# Load each annotation object and run compareRepertoires for each pair
source("Examples/BCRAnalyses/CompareBCRDatasets.R")

# Make some plots based on the above comparisons
source("Examples/BCRAnalyses/BCRSummaryPlots.R")
