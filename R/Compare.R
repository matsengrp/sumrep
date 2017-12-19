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

#' Resample DNA sequences from a fasta file with replacement, and save the
#'   resultant sequences to a new fasta file
#'
#' This acts as a bootstrap on a DNA sequence dataset as an attempt to get 
#'   a benchmark 'control' for assessing similarity in comparison functions.
#'   See \code{compareRepertoires} for more details.
#' @param fasta_file Fasta file from which to resample sequences
#' @param output_filename Desired output_filename
bootstrapFasta <- function(fasta_file, output_filename) {
    if(output_filename %>% file.exists) {
        stop(paste("File", output_filename, "already exists.",
                   "Please choose a different output filename."))
    }

    sequences <- fasta_file %>% 
        seqinr::read.fasta(as.string=TRUE) %>%
        sample(replace=TRUE) %>%
        seqinr::write.fasta(names=seq(1, length(.)), 
                    file.out=output_filename)
}

#' Print the result and elapsed time of a comparison
#'
#' @param f Comparison function 
#' @param input_1 First input to \code{f}
#' @param input_2 Second input to \code{f}
#' @param string_header Text to precede the comparison value and elapsed time
#' @param color Color for text. Green for usual comparison, yellow for 
#'   bootstrap
getAndPrintComparison <- function(f, input_1, input_2, string_header, color) {
    pt <- proc.time()
    comparison <- f(input_1, input_2)
    elapsed_time <- (proc.time() - pt)[3]
    cat(string_header,
        color(comparison %>% signif(4) %>% toString),
        ' (',
        elapsed_time,
        's)',
        '\n', sep='')
    return(comparison)
}

#' Apply \code{function_string} to inputs \code{input_1} and \code{input_2},
#'   and print the results
#'
#' @param function_string The name of the function to apply. Must be a 
#'   comparison that takes two inputs
#' @param input_list List of inputs on which to do comparisons. If 
#'   \code{length(input_list) == 2}, do a usual comparison on the two 
#'   repertoire inputs. If \code{length(input_list) == 3}, a comparison of
#'   the first input with a bootstrapped version will be included.
doComparison <- function(function_string, input_list) {
    f <- eval(parse(text=function_string)) 
    input_1 <- input_list[[1]]
    input_2 <- input_list[[2]]
    comparisons <- getAndPrintComparison(f, input_1, input_2, 
                    string_header=paste0("Result of ", function_string, ": "),
                    color=crayon::green)
    if(length(input_list) == 3) {
        input_1_boot <- input_list[[3]]
        comparisons[2] <- getAndPrintComparison(f, input_1, input_1_boot, 
                        string_header="    Bootstrapped result: ",
                        color=crayon::yellow)
    }
    return(comparisons)
}

#' Run through a full repertoire comparison
#'
#' \code{compareRepertoires} iterates through the various comparison functions
#'   in \code{sumrep}, printing the distance or divergence value of each one.
#'   Both repertoires are assumed to be annotated.
#'   Partition-based comparisons are available if the repertoires are 
#'   partitioned.
#'   Comparisons based on per-gene and per-gene-per-position mutation rates are 
#'   also available given mutation rate information.
#' @param repertoire_1 First repertoire
#' @param repertoire_2 Second repertoire
#' @param rep_1_bootstrap Annotated repertoire based on bootstrapping the DNA
#'   sequences from the first repertoire
compareRepertoires <- function(repertoire_1, repertoire_2, rep_1_bootstrap) {
    annotations_1 <- repertoire_1$annotations
    annotations_2 <- repertoire_2$annotations
    annotations_list <- list(annotations_1, annotations_2)

    if(hasArg(rep_1_bootstrap)) {
        annotations_1_boot <- rep_1_bootstrap$annotations

        if(!identical(names(annotations_1), names(annotations_1_boot))) {
            stop("Bootstrapped repertoire annotations do not match the original's.")
        }

        annotations_list[[3]] <- annotations_1_boot
    }

    sig_digs <- 4
    function_strings <- list(
                             # Distance-based metrics
                             "comparePairwiseDistanceDistributions",
                             "compareNNDistanceDistributions",
                             "compareCDR3Distributions",
                             "compareDistancesFromNaiveToMature",
                             # SHM metrics
                             "compareDistanceBetweenMutationsDistributions",
                             "compareSubstitutionAndMutabilityModels",
                             "compareSelectionStrengths",
                             # Sequence-based metrics
                             "compareGCContents",
                             "compareHotspotCounts",
                             "compareColdspotCounts",
                             "compareHydrophobicityDistributions",
                             "compareAtchleyFactorDistributions",
                             "compareAliphaticIndexDistributions",
                             "compareGRAVYDistributions",
                             "compareInFramePercentages",
                             # Recombination metrics
                             "compareCDR3Lengths",
                             "compareVGeneDistributions",
                             "compareDGeneDistributions",
                             "compareJGeneDistributions",
                             "compareVDJDistributions",
                             "compareVGene3PrimeDeletionLengths",
                             "compareDGene3PrimeDeletionLengths",
                             "compareDGene5PrimeDeletionLengths",
                             "compareJGene5PrimeDeletionLengths", 
                             "compareVDInsertionLengths",
                             "compareDJInsertionLengths",
                             "compareVDInsertionMatrices",
                             "compareDJInsertionMatrices"
                             )
    
    comparison_dat <- matrix(NA, nrow=0, ncol=2) %>% data.table %>%
        setNames(c("Comparison", "Value"))
    boot_comparisons <- {}

    for(f_string in function_strings) {
        comparisons <- doComparison(f_string, annotations_list)
        comparison_dat <- rbind(comparison_dat,
                                data.table(Comparison=f_string,
                                Value=comparisons[1]))
    }

    mutation_rates_1 <- repertoire_1$mutation_rates
    mutation_rates_2 <- repertoire_2$mutation_rates
    if(!is.null(mutation_rates_1) && !is.null(mutation_rates_2)) {
        mutation_function_strings <- 
            list(
                 # More SHM metrics
                 "comparePerGeneMutationRates",
                 "comparePerGenePerPositionMutationRates"
                 )
        for(f_string in mutation_function_strings) {
            comparisons <- doComparison(f_string, 
                                        list(mutation_rates_1, 
                                             mutation_rates_2))
            comparison_dat <- rbind(comparison_dat,
                                data.table(Comparison=f_string,
                                Value=comparisons[1]))
        }
    }


    partition_function_strings <- list(
                                       # Clonal family metrics
                                       "compareClusterSizes",
                                       "compareHillNumbers"
                                       )
    if("clone" %in% intersect(names(annotations_1),
                              names(annotations_2))) {
        function_strings <- c(function_strings, partition_function_strings)
    }



    return(comparison_dat)
}
