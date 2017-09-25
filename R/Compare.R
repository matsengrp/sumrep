#' Resample a dataset with replacement
#'
#' \code{resample_data} returns a data.table (or data.frame) object with
#'   resampled rows from the input dataset. This yields a new dataset with the
#'   same number of rows and similar characteristics.
#' @param dat data.table or data.frame object to be resampled
#' @return A new, resampled data.table or data.frame object, depending on 
#'   the type of \code{dat}
resampleData <- function(dat) {
    return( dat[sample(nrow(dat), replace=TRUE), ] )
}

#' Apply \code{function_string} to inputs \code{input_1} and \code{input_2},
#'   and print the results
#'
#' @param function_string The name of the function to apply. Must be a 
#'   comparison that takes two inputs
#' @param input_1 First input to \code{function_string}
#' @param input_2 Second input to \code{function_string}
doComparison <- function(function_string, input_1, input_2) {
    f <- eval(parse(text=function_string)) 
    pt <- proc.time()
    comparison <- f(input_1, input_2)
    elapsed <- (proc.time() - pt)[3]
    cat("Result of ", function_string, ": ", 
        crayon::green(comparison %>% signif(4) %>% toString),
        ' (', elapsed, 's)', '\n', 
        sep='')
}

#' Run through a full repertoire comparison
#'
#' \code{compareRepertoires} iterates through the various comparison functions
#'   in \code{sumrep}, printing the distance or divergence value of each one.
#'   Both repertoires are assumed to be both annotated.
#'   Partition-based comparisons are available if the repertoires are 
#'   partitioned.
#'   Comparisons based on per-gene and per-gene-per-position mutation rates are 
#'   also available given mutation rate information.
#' @param repertoire_1 First repertoire
#' @param repertoire_2 Second repertoire
compareRepertoires <- function(repertoire_1, repertoire_2) {
    annotations_1 <- repertoire_1$annotations
    annotations_2 <- repertoire_2$annotations

    sig_digs <- 4
    function_strings <- list("comparePairwiseDistanceDistributions",
                             "compareNNDistanceDistributions",
                             "compareGCContents",
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
                             "compareDJInsertionLengths"
                             )

    partition_function_strings <- list("compareClusterSizes",
                                       "compareHillNumbers")
    if("clone" %in% intersect(names(annotations_1),
                              names(annotations_2))) {
        function_strings <- c(function_strings, partition_function_strings)
    }

    for(f_string in function_strings) {
        doComparison(f_string, annotations_1, annotations_2)
    }

    mutation_rates_1 <- repertoire_1$mutation_rates
    mutation_rates_2 <- repertoire_2$mutation_rates
    if(!is.null(mutation_rates_1) && !is.null(mutation_rates_2)) {
        mutation_function_strings <- 
            list("comparePerGeneMutationRates",
                 "comparePerGenePerPositionMutationRates")
        for(f_string in mutation_function_strings) {
            doComparison(f_string, mutation_rates_1, mutation_rates_2)
        }
    }
}
