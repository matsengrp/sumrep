source("R/SummaryStats.R")

#' Run through a full repertoire comparison
#'
#' \code{compareRepertoires} iterates through the various comparison functions
#'   in \code{sumrep}, printing the distance or divergence value of each one.
#'   Both repertoires are assumed to be both annotated and partitioned.
#' @param repertoire_1 First repertoire
#' @param repertoire_2 Second repertoire
compareRepertoires <- function(repertoire_1, repertoire_2) {
    function_strings <- list("compareDistancesFromNaiveToMature",
                             "compareCDR3Lengths",
                             "compareCDR3Distributions",
                             "compareDistanceBetweenMutationsDistributions",
                             "compareSubstitutionModels",
                             "compareVGene3PrimeDeletionLengths",
                             "compareDGene3PrimeDeletionLengths",
                             "compareDGene5PrimeDeletionLengths",
                             "compareJGene5PrimeDeletionLengths", 
                             "compareVDInsertionLengths",
                             "compareDJInsertionLengths",
                             "compareClusterSizes",
                             "compareHillNumbers")
    for(f_string in function_strings) {
        f <- eval(parse(text=f_string)) 
        cat("Result of", f_string, ": ", f(repertoire_1, repertoire_2), '\n')
    }
}
