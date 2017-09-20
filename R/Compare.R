source("R/SummaryStats.R")

#' Resample a dataset with replacement
#'
#' \code{resample_data} returns a data.table (or data.frame) object with
#'   resampled rows from the input dataset. This yields a new dataset with the
#'   same number of rows and similar characteristics.
#' @param dat data.table or data.frame object to be resampled
#' @return A new, resampled data.table or data.frame object, depending on 
#'   the type of \code{dat}
resample_data <- function(dat) {
    return( dat[sample(nrow(dat), replace=TRUE), ] )
}

#' Run through a full repertoire comparison
#'
#' \code{compareRepertoires} iterates through the various comparison functions
#'   in \code{sumrep}, printing the distance or divergence value of each one.
#'   Both repertoires are assumed to be both annotated and partitioned.
#' @param repertoire_1 First repertoire
#' @param repertoire_2 Second repertoire
compareRepertoires <- function(repertoire_1, repertoire_2) {
    sig_digs <- 4
    function_strings <- list("compareGCContents",
                             "compareHotspotCounts",
                             "compareColdspotCounts",
                             "compareDistancesFromNaiveToMature",
                             "compareCDR3Lengths",
                             "compareVGeneDistributions",
                             "compareDGeneDistributions",
                             "compareJGeneDistributions",
                             "compareVDJDistributions",
                             "compareHydrophobicityDistributions",
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
        pt <- proc.time()
        comparison <- f(repertoire_1, repertoire_2)
        elapsed <- (proc.time() - pt)[3]
        cat("Result of ", f_string, ": ", 
            crayon::green(comparison %>% signif(4) %>% toString),
            ' (', elapsed, 's)', '\n', 
            sep='')
    }
}
