library(alakazam)
library(ape)
library(Biostrings)
library(data.table)
library(dplyr)
library(HDMD)
library(jsonlite)
library(magrittr)
library(pegas)
library(Peptides)
library(RecordLinkage)
library(shazam)
library(seqinr)
library(stringdist)
library(textmineR)
library(yaml)

source("R/PartisFunctions.R")

#' Similarity metric for joint categorical distributions
#'
#' Currently unused. This metric attempts to compare contingency tables of two
#'   distributions that accounts for joint and all marginal usage rates.
#' @param table_a First contingency table
#' @param table_b Second contingency table
#' @return A positive-valued distance metric reflecting the similarity of the
#'   distributions manifest in table_a and table_b
compareCategoricalDistributions <- function(table_a, table_b) {
    table_a_dims <- table_a %>% dim
    table_b_dims <- table_b %>% dim
    dim_count_a <- table_a_dims %>% length
    dim_count_b <- table_b_dims %>% length
    if(dim_count_a != dim_count_b || !all(table_a_dims == table_b_dims)) {
        stop("Table dimensions must agree")
    }
    dims <- table_a_dims
    entrywise_sum_of_absolute_differences <- (table_a - table_b) %>% abs %>% sum
    entrywise_comparison <- entrywise_sum_of_absolute_differences/prod(dims)
    dimension_count <- dims %>% length
    dimension_comparisons <- rep(NA, dimension_count)
    for(i in 1:dimension_count) {
        dimension_sums_a <- apply(table_a, i, sum)
        dimension_sums_b <- apply(table_b, i, sum)
        dimension_sum_of_absolute_differences <- 
            (dimension_sums_a - dimension_sums_b) %>% abs %>% sum
        other_dims <- dims[-i]
        dimension_comparisons[[i]] <- 
            dimension_sum_of_absolute_differences/prod(other_dims)
    }
    total_comparison <- entrywise_comparison + sum(dimension_comparisons)
    return(total_comparison)
}

#' Discretize two lists of continuous data into mutual, well-defined bins.
#' 
#' @param list_a First list to bin
#' @param list_b Second list to bin
#' @return a list of vectors containing the corresponding counts for each bin 
#' @examples
#' list_a <- runif(100)
#' list_b <- runif(100)
#' binned <- binContinuousListsAsDiscrete(list_a, list_b) %>% 
#'     melt(value.factor=TRUE) %>% 
#'     mutate(bin_number=rep(seq(1:(nrow(.)/2)), 2))
#' names(binned)[1:2] <- c("count", "list_id")
#' ggplot(binned, aes(x=bin_number, y=count, fill=as.factor(list_id))) + 
#'     geom_bar(stat="identity", position="dodge")

binContinuousListsAsDiscrete <- function(list_a, list_b) {
    a_length <- list_a %>% length
    b_length <- list_b %>% length
    bin_count <- min(a_length, b_length) %>% sqrt %>% ceiling
    bin_count <- ifelse(bin_count < 2, 2, bin_count)
    bins <- c(list_a, list_b) %>% cut(breaks=bin_count, labels=1:bin_count)
    table_a <- bins[1:a_length] %>% table %>% unname %>% as.vector
    table_b <- bins[-(1:a_length)] %>% table %>% unname %>% as.vector
    return(list(table_a, table_b))
}

#' Compute JS divergence for continuous samples
#'
#' This is currently not the default method for continuous data, as it requires
#'   constructing approximating functions of the two densities, and then using
#'   numerical integrals for the expectation steps. 
#' @param sample_1 First data sample
#' @param sample_2 Second data sample
#' @return Approximate JS divergence of the distributions induced from sample_1
#'   and sample_2
getContinuousJSDivergence <- function(sample_1, sample_2) {
    m <- function(x) {
        result <- 0.5*(p(x) + q(x))
        return(result)
    }

    integrand <- function(f, g) {
        func <- function(x) {
            result <- ifelse(f(x) == 0, 0, ifelse(g(x) == 0, Inf, 
                                                  f(x)*(log(f(x)) - log(g(x)))))
            return(result)
        }
        return(func)
    }

    p <- sample_1 %>% density %>% approxfun
    q <- sample_2 %>% density %>% approxfun
    lower <- max(min(sample_1), min(sample_2))
    upper <- min(max(sample_1), max(sample_2))
    KL_div_1 <- integrate(integrand(p, m), lower, upper)$value
    KL_div_2 <- integrate(integrand(q, m), lower, upper)$value
    JS_divergence <- 0.5*(KL_div_1 + KL_div_2)
    return(JS_divergence)
}

#' Compute the JS Divergence of two samples
#'
#' The samples are assumed to contain discrete data by
#' default. If continuous, the lists are passed to 
#' \link{binContinuousListsAsDiscrete}
#' to be discretized into bins commensurate to the list sizes. 
#' Note: This function is symmetric in \code{list_a} and \code{list_b}, 
#' since JS-divergence is symmetric.
#' @param list_a First sample
#' @param list_b Second sample
#' @return The positive-valued JS-divergence of the distributions induced from 
#'   list_a and list_b
#' @examples
#' l1 <- sample.int(100, replace=TRUE)
#' l2 <- sample.int(100, replace=TRUE)
#' getJSDivergence(l1, l2)
#' getJSDivergence(l2, l1)
#' getJSDivergence(l1, l1)
getJSDivergence <- function(list_a, list_b, continuous=FALSE) {
    if(continuous) {
        binned <- binContinuousListsAsDiscrete(list_a, list_b)
        divergence <- textmineR::CalcJSDivergence(binned[[1]], binned[[2]])
    } else {
        max_val <- max(list_a, list_b)
        table_a <- list_a %>% factor(levels=0:max_val) %>% table %>% as.vector
        table_b <- list_b %>% factor(levels=0:max_val) %>% table %>% as.vector
        divergence <- textmineR::CalcJSDivergence(table_a, table_b)
    }
    return(divergence)
}

#' Remove empty strings from a list or vector
#'
#' @param l List or vector of strings
#' @return The given list or vector with all empty strings removed
removeEmptyStrings <- function(l) {
    return(l[l != ""])
}

#' Convert lists of factors or vectors of characters into a vector of strings
#' 
#' @param List or vector of either strings or char vectors of DNA sequences
#' @return Vector of strings of DNA sequences
standardizeList <- function(l) {
    new_list <- l %>% 
                sapply(toString) %>%
                gsub(pattern=", *", replace="") %>%
                sapply(paste, collapse='') %>% unname
    return(new_list)
}

#' Determine the comparison method for stringdistmatrix 
#'
#' The comparison method is determined based on whether 
#' all sequences are of the same length. If so, use hamming distance; 
#' else use Levenshtein.
#' @param List/vector of DNA sequences
#' @return String-valued comparison method for use in stringdistmatrix within 
#'   the getDistanceMatrix function
determineComparisonMethod <- function(sequence_list) {
    length_count <- sequence_list %>% standardizeList %>% sapply(nchar) %>% 
        table %>% length 
    comparison_method <- ifelse(length_count > 1, "lv", "hamming")
    return(comparison_method)
}

#' Get distance matrix from stringdist
#' 
#' \code{getDistanceMatrix} determines 
#' @param raw_sequences List or vector of DNA sequences
#' @return Distance matrix of the sequences, using hamming distances if all
#'   sequences are the same length, and levenshtein otherwise
getDistanceMatrix <- function(raw_sequences) {
    sequence_list <- raw_sequences %>% standardizeList
    comparison_method <- sequence_list %>% determineComparisonMethod
    mat <- sequence_list %>% 
        stringdist::stringdistmatrix(method=comparison_method) %>% as.matrix
    return(mat)
}

#' Get sorted vector of pairwise distances
#'
#' Convert the matrix given by \link{getDistanceMatrix} into a sorted
#'   vector of distance values.
#' @param sequence_list List of DNA sequence strings
#' @return Vector of pairwise distances
getDistanceVector <- function(sequence_list) {
    mat <- sequence_list %>% getDistanceMatrix
    vec <- mat[mat %>% lower.tri] %>% as.vector %>% sort
    return(vec)
}

#' Subsample a vector
#' 
#' @inheritParams subsample
#' @param v Vector to subsample
#' @return A subsampled vector 
subsampleVector <- function(v, sample_count) {
    new_vector <- {}
    if(length(v) > sample_count) {
        new_vector <- sample(v, sample_count)
    } else {
        new_vector <- v
    }
    return(new_vector)
}

#' Subsample a dataset
#'
#' @param dataset A data.table, data.frame, or vector object
#' @param sample_count Number of samples to retain in the subsampled data.
#'   Samples refers to elements in a vector or rows in a data.table/data.frame
#' @return A subsampled dataset of the same type given by \code{dataset}
subsample <- function(dataset, sample_count) {
    new_dataset <- {}
    if(is.null(dim(dataset))) {
        new_dataset <- subsampleVector(dataset, sample_count)
    } else if(nrow(dataset) > sample_count) {
        new_dataset <- dataset[sample(1:nrow(dataset), sample_count), ]
    } else {
        new_dataset <- dataset
    }
    return(new_dataset)
}

#' Get the mean absolute difference of two sequences 
#' 
#' @param sequence_a First sequence, either a vector or matrix
#' @param sequence_b Second sequence
#' @param ignore_na Should we ignore na values when computing the mean?
getMeanAbsoluteDifference <- function(sequence_a, sequence_b, na_rm=TRUE) {
    difference <- (sequence_a - sequence_b) %>% abs %>% mean(na.rm=ignore_na)
    return(difference)
}

#' Estimate a divergence by subsampling and averaging over the comparison
#'   function output
#'
#' Arguemnts to \code{func} may be supplied after the other required 
#'   parameters.
#' @param dataset_a First dataset (vector, matrix, or data.table) input to 
#'   \code{func}
#' @param dataset_b Second dataset input to \code{func}
#' @param func Comparison function for which to calculate the JS divergence
#' @param subsample_count Number of samples to extract on each trial
#' @param trial_count Number of trials to average over
#' @param divergenceFunction The divergence function to apply to the two 
#'   datasets. JS divergence by default, although this only makes sense for
#'   vectors corresponding to empirical distributions
#' @return The average computed JS divergence across trials, used as an
#'   estimate of the true JS divergence
#' @examples
#' seq_1 <- c("AAAAA", "AATTT", "TTATT", "ACATG", "GGGAA")
#' seq_2 <- c("AAAAC", "ACATG", "CGGGA", "ACATG", "GGACA")
#' getAverageDivergence(seq_1, seq_2, getNearestNeighborDistances, 
#'   subsample_count=3, trial_count=20, k=2)
getAverageDivergence <- function(dataset_a, dataset_b, func, 
                                   subsample_count, 
                                   trial_count, 
                                   divergenceFunction=getJSDivergence,
                                   ...) {
    divergences <- {}
    for(trial in 1:trial_count) {
        distances_a <- dataset_a %>% 
            subsample(sample_count=subsample_count) %>%
            func
        distances_b <- dataset_b %>% 
            subsample(sample_count=subsample_count) %>%
            func
        divergences[trial] <- divergenceFunction(distances_a, distances_b)
    }
    return(divergences %>% mean)
}

#' Compare pairwise distance distributions of two lists of sequences
#'
#' \code{comparePairwiseDistanceDistributions} computes the JS
#'   divergence of the pairwise distance distribution of two lists of
#'   DNA sequences. This function iterates through a number of trials (given
#'   by \code{trial_count}, subsampling the full datasets by the amount given
#'   by \code{subsample_count}, and returns the mean divergence.
#' @param dat_a First dataset, containing mature sequences
#' @param dat_b Second dataset, containing mature sequences
#' @inheritParams getAverageDivergence 
#' @return Estimated JS divergence of the distributions inferred from list_a
#'   and list_b
comparePairwiseDistanceDistributions <- function(dat_a, dat_b,
                                                 subsample_count=100,
                                                 trial_count=10) {
    average_divergence <- getAverageDivergence(dat_a %$% mature_seq,
                                                dat_b %$% mature_seq,
                                                getDistanceVector,
                                                subsample_count,
                                                trial_count)
    return(average_divergence)
}

#' Get sorted list of nearest neighbor distances
#'
#' \code{getNearestNeighborDistances} returns a list of distances of
#'   the kth nearest neighbor.
#' @param sequence_list List of DNA sequence strings
#' @param k The separation depth for the nearest neighbor distances.
#'   k = 1 corresponds to the nearest neighbor, k = 2 corresponds to
#'   the second-nearest neighbor, etc.
#' @return Vector of kth nearest neighbor distances
getNearestNeighborDistances <- function(sequence_list, k=1) {
    mat <- sequence_list %>% getDistanceMatrix
    n <- sequence_list %>% length
    distances <- rep(NA, n)
    for(i in 1:n) {
        distances[i] <- sort(mat[i, -i], partial=k)[k]
    }

    return(distances)
}

#' Compare kth nearest neighbor distance distributions of two lists of
#'   sequences
#' \code{compareNNDistanceDistribution} computes the JS divergence of
#'   the kth nearest neighbor distance distribution of two lists of
#'   DNA sequences
#' @inheritParams getAverageDivergence
#' @param dat_a First dataset, containing mature sequences
#' @param dat_b Second dataset, containing mature sequences
#' @param k The separation depth for the nearest neighbor distances.
#'   k = 1 corresponds to the nearest neighbor, k = 2 corresponds to
#'   the second-nearest neighbor, etc.
#' @return Estimated JS divergence of the distributions inferred from list_a
#'   and list_b
compareNNDistanceDistributions <- function(dat_a, dat_b, k=1,
                                           subsample_count=100,
                                           trial_count=10) {
    average_divergence <- getAverageDivergence(dat_a %$% mature_seq,
                                                 dat_b %$% mature_seq,
                                                 getNearestNeighborDistances,
                                                 subsample_count,
                                                 trial_count,
                                                 k=k)
    return(average_divergence)
}

#' Get the GC content distribution of a list of DNA sequences
#'
#' \code{getGCContentDistribution} returns a list of the GC content of each DNA
#'   sequence in \code{raw_sequences}, given by the \code{ape} library
#' @param raw_sequences List or vector of strings or character sequences
#'   corresponding to DNA sequences
#' @return A vector of GC content values
getGCContentDistribution <- function(raw_sequences) {
    sequence_list <- raw_sequences %>% sapply(paste, collapse='') %>% unname
    dna_list <- sequence_list %>% strsplit(split='') %>% lapply(ape::as.DNAbin)
    gc_dist <- dna_list %>% sapply(ape::GC.content)
    return(gc_dist)
}

#' Compare the GC distributions of two lists of DNA sequences
#'
#' @param list_a First list or vector of DNA sequences
#' @param list_b Second list or vector of DNA sequences
#' @return JS divergence of the GC content distributions inferred from list_a
#'   and list_b
compareGCContents <- function(dat_a, dat_b) {
    density_a <- dat_a %$% mature_seq %>% getGCContentDistribution
    density_b <- dat_b %$% mature_seq %>% getGCContentDistribution
    divergence <- getJSDivergence(density_a, density_b, continuous=TRUE)
    return(divergence)
}

#' Get number of times a motif occurs in a set of one or more reference 
#'   sequences
#'
#' @param motif String representing the motif pattern 
#' @param dna_sequences List, vector of reference sequences
#' @return The number of occurrences of \code{motif} in \code{dna_sequences}
getMotifCount <- function(motif, dna_sequences) {
    dna_strings <- dna_sequences %>% unlist %>% Biostrings::DNAStringSet()
    count <- motif %>% Biostrings::vcountPattern(dna_strings, fixed=FALSE) %>% 
        sum
    return(count)
}

#' Get number of times a set of motif (hot/cold)spots occur in a set of 
#'   reference sequences
#'
#' @inheritParams getMotifCount
#' @param spots Vector of hot or cold spots of interest
#' @return The total number of occurrences of each motif in \code{spots}, in
#'   \code{dna_sequences}
getSpotCount <- function(dna_sequences, spots) {
    count <- spots %>% sapply(getMotifCount, dna_sequences=dna_sequences) %>% sum
    return(count)
}

#' Get the number of occurrences of AID hotspots in a set of reference 
#'   sequences
#' 
#' @inheritParams getMotifCount
#' @return The number of AID hotspot occurrences in \code{dna_sequences}
getHotspotCount <- function(dna_sequences, hotspots=c("WRC", "WA")) {
    return(getSpotCount(dna_sequences, hotspots))
}

#' Get the number of occurrences of AID coldspots in a set of reference 
#'   sequences
#' 
#' @inheritParams getMotifCount
#' @return The number of AID coldhotspot occurrences in \code{dna_sequences}
getColdspotCount <- function(dna_sequences, coldspots=c("SYC")) {
    return(getSpotCount(dna_sequences, coldspots))
}

#' Compare hot or coldspot count distributions of two sets of mature BCR sequences
#'
#' @param dat_a First dataset
#' @param dat_b Second dataset
#' @param count_function The comparison function corresponding to hotspot or 
#'   coldspot counts. Must be either getHotspotCount or getColdspotCount
#' @param subsample_count Number of subsamples for averaging
#' @param trial_count Number of subsample trials
#' @return The average JS divergence of the count distributions inferred from
#'   \code{dat_a$mature_seq} and \code{dat_b$mature_seq}, respectively
compareCounts <- function(dat_a,
                          dat_b,
                          count_function,
                          subsample_count,
                          trial_count) {

    divergence <- getAverageDivergence(dat_a %$% mature_seq,
                                         dat_b %$% mature_seq,
                                         count_function,
                                         subsample_count,
                                         trial_count)
    return(divergence)
}

#' Compare hotspot count distributions of two sets of mature BCR sequences
#'
#' @inheritParams compareCounts
#' @return The JS divergence of the hotspot count distributions inferred from
#'   \code{dat_a$mature_seq} and \code{dat_b$mature_seq}, respectively
compareHotspotCounts <- function(dat_a, 
                                 dat_b,
                                 subsample_count=1000,
                                 trial_count=10) {
    divergence <- compareCounts(dat_a,
                                dat_b,
                                getHotspotCount,
                                subsample_count,
                                trial_count)
    return(divergence)
}

#' Get the number of occurrences of AID coldspots in a set of reference 
#'   sequences
#' 
#' @inheritParams getMotifCount
#' @return The number of AID coldhotspot occurrences in \code{dna_sequences}
getColdspotCount <- function(dna_sequences, coldspots=c("SYC")) {
    return(getSpotCount(dna_sequences, coldspots))
}

#' Compare coldspot count distributions of two sets of mature BCR sequences
#'
#' @inheritParams compareCounts
#' @return The JS divergence of the coldspot count distributions inferred from
#'   \code{dat_a$mature_seq} and \code{dat_b$mature_seq}, respectively
compareColdspotCounts <- function(dat_a, 
                                  dat_b,
                                  subsample_count=1000,
                                  trial_count=10) {
    divergence <- compareCounts(dat_a,
                                dat_b,
                                getColdspotCount,
                                subsample_count,
                                trial_count)
    return(divergence)
}

#' Get distribution of Levenshtein distances from the inferred naive sequences
#' to the observed mature ones, sequence by sequence
#'
#' @param dat Annotated dataset
#' @return Vector of Levenshtein distances from naive to mature
getDistancesFromNaiveToMature <- function(dat) {
    distances <- dat$mature_seq %>% 
        mapply(FUN=stringdist::stringdist, b=dat$naive_seq, method="lv") %>% 
        sort %>% unname
    return(distances)
}

#' Compare Levenshtein distance distributions from naive sequences to 
#  their corresponding mature ones for two repertoires
#' @param dat_a First dataset
#' @param dat_b Second dataset
#' @return JS divergence of the two distance distributions
compareDistancesFromNaiveToMature <- function(dat_a, dat_b) {
    distances_a <- getDistancesFromNaiveToMature(dat_a)
    distances_b <- getDistancesFromNaiveToMature(dat_b)
    divergence <- getJSDivergence(distances_a, distances_b)
    return(divergence)
}

#' Import a vector of DNA sequence strings from a fasta file
#' 
#' @param filename Name of fasta file including the sequences
#' @return A vector of DNA sequence strings
getSequenceListFromFasta <- function(filename) {
    sequences <- filename %>% seqinr::read.fasta() %>% 
        lapply(paste, collapse="") %>% unlist %>% unname
    return(sequences)
}

#' Get the distribution of inferred length of CDR3 regions of each sequence
#' 
#' @param Annotated dataset
#' @return Vector of CDR3 lengths (in nt units)
getCDR3Lengths <- function(dat) {
    CDR3_lengths <- dat$cdr3_length %>% na.omit
    return(CDR3_lengths)
}

#' Compare the distribution of CDR3 lengths for two datasets
#'
#' @param dat_a First dataset
#' @param dat_b Second dataset
#' @return The JS divergence of the two CDR3 length distributions
compareCDR3Lengths <- function(dat_a, dat_b) {
    a_lengths <- getCDR3Lengths(dat_a)
    b_lengths <- getCDR3Lengths(dat_b)
    divergence <- getJSDivergence(a_lengths, b_lengths)
    return(divergence)
}

#' Get gene usage table from full list of genes, setting zero to unused ones
#' 
#' @param gene_list List of genes being used in a given repertoire
#' @param gene_list List of all available genes present in at least one
#'   repertoire under consideration
getGeneUsageTableFromFullList <- function(gene_list, full_gene_list) {
    usage_table <- gene_list %>% table
    usage_table[full_gene_list %>% setdiff(gene_list) %>% unlist] <- 0
    return(usage_table)
}

#' Ignore allelic differences of common genes
#' 
#' This function is used when reporting gene usage, so that two alleles of a
#' common gene are both treated as the one gene, and not separately.
#' For example, IGHV3-13*01 and IGHV3-13*02 both become IGHV3-13.
#' @param gene_list Vector of genes possibly broken down by allele
#' @return Vector of genes with allelic variant information discarded
collapseAlleles <- function(gene_list) {
    collapsed_gene_list <- gene_list %>% 
        sapply(gsub, pattern="\\*\\d+", replace="") %>%
        unname %>%
        sapply(as.factor)
    return(collapsed_gene_list)
}

#' Compare gene usage counts for two lists of genes
#'
#' @param gene_list_a First vector of genes 
#' @param gene_list_b Second vector of genes
#' @param collapse_alleles Should allelic variants be ignored?
#' @return Mean absolute difference of gene counts between \code{gene_list_a}
#'   and \code{gene_list_b}
compareGeneUsage <- function(gene_list_a, gene_list_b, collapse_alleles) {
    if(collapse_alleles) {
        gene_list_a <- gene_list_a %>% collapseAlleles
        gene_list_b <- gene_list_b %>% collapseAlleles
    }
    full_gene_list <- union(gene_list_a, gene_list_b)
    table_a <- gene_list_a %>% getGeneUsageTableFromFullList(full_gene_list)
    table_b <- gene_list_b %>% getGeneUsageTableFromFullList(full_gene_list)
    divergence <- full_gene_list %>% 
        sapply(function(x) { abs(table_a[x] - table_b[x]) }) %>% mean
    return(divergence)
}

#' Compare germline V, D, or J gene usage for two repertoires
#' @param dat_a First annotated dataset
#' @param dat_b Second annotated dataset
#' @param gene_type String for gene type, taken as a column name of the 
#'   annotated datasets. Must be "v_gene", "d_gene", or "j_gene"
#' @inheritParams compareGeneUsage
#' @return Mean absolute difference of gene counts between the two
#'   repertoires
compareGermlineGeneDistributions <- function(dat_a, dat_b, gene_type,
                                             collapse_alleles) {
    if(gene_type %>% missing) {
        stop("gene_type needs to be supplied.")
    }

    gene_list_a <- dat_a[, gene_type, with=FALSE] 
    gene_list_b <- dat_b[, gene_type, with=FALSE]
    divergence <- compareGeneUsage(gene_list_a, gene_list_b, 
                                   collapse_alleles) 
    return(divergence)
}

#' Compare V gene distribution between two repertoires
#'
#' @inheritParams compareGermlineGeneDistributions
#' @return Mean absolute difference of V gene counts between the two
#'   repertoires
compareVGeneDistributions <- function(dat_a, dat_b) {
    return(compareGermlineGeneDistributions(dat_a, dat_b, gene_type="v_gene",
                                            collapse_alleles=TRUE))
}

#' Compare D gene distribution between two repertoires
#'
#' @inheritParams compareGermlineGeneDistributions
#' @return Mean absolute difference of D gene counts between the two
#'   repertoires
compareDGeneDistributions <- function(dat_a, dat_b) {
    return(compareGermlineGeneDistributions(dat_a, dat_b, gene_type="d_gene",
                                            collapse_alleles=TRUE))
}

#' Compare J gene distribution between two repertoires
#'
#' @inheritParams compareGermlineGeneDistributions
#' @return Mean absolute difference of J gene counts between the two
#'   repertoires
compareJGeneDistributions <- function(dat_a, dat_b) {
    return(compareGermlineGeneDistributions(dat_a, dat_b, gene_type="j_gene",
                                            collapse_alleles=TRUE))
}

#' Get table of combined V, D, and J usage frequencies
#'
#' @param dat Annotated dataset
#' @inheritParams compareGermlineGeneDistributions
#' @return A table of joint gene IDs and usage counts
getJointGeneTable <- function(dat, collapseAlleles) {
    if(collapseAlleles) {
        v_genes <- dat %$% v_gene %>% collapseAlleles
        d_genes <- dat %$% d_gene %>% collapseAlleles
        j_genes <- dat %$% j_gene %>% collapseAlleles
        gene_dat <- data.table(v_gene=v_genes,
                                d_gene=d_genes,
                                j_gene=j_genes)
    } else {
        gene_dat <- dat
    }

    gene_type_list <- c("v_gene", "d_gene", "j_gene")
    gene_table <- gene_dat %>% plyr::count(gene_type_list)
    gene_table$concat <- do.call(paste0, gene_table[gene_type_list])

    return(gene_table[, c("concat", "freq")])
}

#' Compare joint V, D, and J gene usage between two annotated repertoires
#' @inheritParams compareGermlineGeneDistributions
#' @return Mean absolute difference between joint gene usage, over all
#'   observed genes between the two repertoires
compareVDJDistributions <- function(dat_a, dat_b, collapseAlleles=TRUE) {
    table_a <- getJointGeneTable(dat_a, collapseAlleles)
    table_b <- getJointGeneTable(dat_b, collapseAlleles)

    full_triple_list <- union(table_a$concat, table_b$concat)

    summands <- rep(NA, length(full_triple_list))
    i <- 1
    for(triple in full_triple_list) {
        a_count <- table_a[table_a$concat == triple, ]$freq
        b_count <- table_b[table_b$concat == triple, ]$freq
        a_count <- ifelse(length(a_count) == 0, 0, a_count)
        b_count <- ifelse(length(b_count) == 0, 0, b_count)
        summands[i] <- abs(a_count - b_count)
        i <- i + 1
    }
    
    divergence <- summands %>% mean
    return(divergence)
}

#' Convert a DNA string to a string of single-letter amino acid codes
#'
#' @param sequence String of DNA bases
#' @return String of single-letter amino acid codes
convertDNAToAminoAcids <- function(sequence) {
    aa_sequence <- sequence %>% sapply(strsplit, '') %>% 
        sapply(seqinr::translate) %>% paste0(collapse='')
    return(aa_sequence)
}

#' Check if an amino acid sequence has unrecognized codons
#' 
#' @param aa_sequence String of single-letter amino acid codons
#' @return A boolean stating whether any unrecognized AA codons are present 
#'   in \code{aa_sequence}
hasUnrecognizedAminoAcids <- function(aa_sequence) {
    return(grepl("X", aa_sequence))  
}

#' Check if an amino acid sequence has a stop codon
#' 
#' @param aa_sequence String of single-letter amino acid codes
#' @return A boolean stating whether any stop codons are present 
#'   in \code{aa_sequence}
hasStopCodon <- function(aa_sequence) {
    return(grepl("\\*", aa_sequence))
}

#' Remove strings in a list of amino acid sequences which contain at least one
#'   unrecognized amino acid or stop codon
#'
#' @param aa_sequences Vector or list of amino acid sequences
#' @return Filtered list with "bad" sequences removed
removeBadAminoAcidSequences <- function(aa_sequences) {
    return(aa_sequences[!hasUnrecognizedAminoAcids(aa_sequences) &
                        !hasStopCodon(aa_sequences)])
}

getKideraFactorsBySequence <- function(sequence) {
    kidera_factors <- sequence %>% 
        Peptides::kideraFactors() %>% 
        first
    return(kidera_factors)
}

getKideraFactors <- function(sequence_list) {
    kidera_factors <- sequence_list %>%
        sapply(convertDNAToAminoAcids) %>%
        removeBadAminoAcidSequences %>%
        lapply(getKideraFactorsBySequence) %>% 
        do.call(what=rbind)
    return(kidera_factors)
}

getHydrophobicityDistribution <- function(sequence_list) {
    hydrophobicity_list <- sequence_list %>% 
        removeEmptyStrings %>%
        getKideraFactors %>%
        data.table %>% select_("KF4") %>% unlist
    return(hydrophobicity_list)
}

compareHydrophobicityDistributions <- function(dat_a, dat_b) {
    dist_a <- dat_a %$% cdr3s %>% getHydrophobicityDistribution
    dist_b <- dat_b %$% cdr3s %>% getHydrophobicityDistribution
    divergence <- getJSDivergence(dist_a, dist_b, continuous=TRUE)
    return(divergence)
}

#' Get the mean Atchley factor of a DNA sequence
#'
#' TODO: Do something more clever than taking the mean. Is there
#'   a set of factors of interest that we can isolate? Or should
#'   we do an average pairwise comparison?
#' @param dna_seq String of a DNA sequence
#' @return The mean of the set of Atchley factors for 
#'   \code{dna_seq}
getMeanAtchleyFactorBySequence <- function(dna_seq) {
    atchley_factors <- dna_seq %>%
        HDMD::FactorTransform() %>%
        first %>%
        mean
    return(atchley_factors)
}

#' Get the distribution of mean Atchley factors for a list of DNA seuqnces
#'
#' @param sequence_list List or vector of DNA sequences
#' @return Vector of mean Atchley factors of \code{sequence_list}
getMeanAtchleyFactorDistribution <- function(sequence_list) {
    atchley_factors <- sequence_list %>% 
        sapply(convertDNAToAminoAcids) %>%
        removeBadAminoAcidSequences %>%
        sapply(getMeanAtchleyFactorBySequence) %>%
        unname %>%
        unlist
    return(atchley_factors)
}

#' Compare the distributions of mean Atchley factors for two datasets
#'
#' @param dat_a First dataset, a data.table or data.frame
#' @param dat_b Second dataset, a data.table or data.frame
#' @return JS divergence of Atchley factor distributions
compareAtchleyFactorDistributions <- function(dat_a, dat_b) {
    dist_a <- dat_a$cdr3s %>% getMeanAtchleyFactorDistribution
    dist_b <- dat_b$cdr3s %>% getMeanAtchleyFactorDistribution
    divergence <- getJSDivergence(dist_a, dist_b, continuous=TRUE)
    return(divergence)
}

#' Get the aliphatic index of a DNA sequence
#'
#' @param dna_sequence String of DNA characters
#' @return The aliphatic index of \code{dna_sequence}, if applicable
getAliphaticIndex <- function(dna_sequence) {
    aliphatic_index <- dna_sequence %>% 
        Peptides::aIndex()
    return(aliphatic_index)
}

#' Get the distribution of aliphatic indices of a list of sequences
#'
#' @param sequence_list List or vector of DNA sequences
#' @return Vector of aliphatic indices
getAliphaticIndexDistribution <- function(sequence_list) {
    a_indices <- sequence_list %>% 
        sapply(convertDNAToAminoAcids) %>%
        removeBadAminoAcidSequences %>%
        sapply(getAliphaticIndex) %>%
        unname %>%
        unlist
    return(a_indices)
}

#' Compare the distributions of aliphatic indices of two datasets
#'
#' @param dat_a First dataset, a data.table or data.frame
#' @param dat_b Second datset, a data.table or data.frame
#' @return The JS divergence of the two distributions
compareAliphaticIndexDistributions <- function(dat_a, dat_b) {
    dist_a <- dat_a$cdr3s %>% getAliphaticIndexDistribution
    dist_b <- dat_b$cdr3s %>% getAliphaticIndexDistribution
    divergence <- getJSDivergence(dist_a, dist_b, continuous=TRUE)
    return(divergence)
}

#' Get the distribution of GRAVY values from a list or vector of sequences
#'
#' @param List or vector of DNA sequence strings
#' @return Vector of GRAVY values for \code{sequence_list}
getGRAVYDistribution <- function(sequence_list) {
    dist <- sequence_list %>% 
            removeEmptyStrings %>% 
            standardizeList %>%
            toupper %>%
            sapply(alakazam::gravy) %>% 
            unname
    return(dist)
}

#' Compare the GRAVY distributions of two datasets
#'
#' @param dat_a First dataset, a data.table or data.frame object
#' @param dat_b Second dataset, a data.table or data.frame object
#' @return The JS divergence of GRAVY distributions
compareGRAVYDistributions <- function(dat_a, dat_b) {
    dist_a <- dat_a %$% cdr3s %>% getGRAVYDistribution
    dist_b <- dat_b %$% cdr3s %>% getGRAVYDistribution
    divergence <- getJSDivergence(dist_a, dist_b, continuous=TRUE)
    return(divergence)
}

#' Parse a python dictionary string into an R list
#'
#' @param dictionary String that represents a dictionary in Python's syntax
#' @return An R list containing the same information as \code{dictionary}
parsePythonDictionary <- function(dictionary) {
    parsed <- dictionary %>% gsub(pattern="'", replacement='"') %>% 
        gsub(pattern="\\(", replacement="\\[") %>% 
        gsub(pattern="\\)", replacement="\\]") %>% jsonlite::fromJSON()
    return(parsed)
}

#' Extract CDR3 codon start positions from partis-returned dictionary strings
#'
#' This is a helper function to \link{getCDR3s}.
#' Note: We add one to \code{positions} because partis returns position values
#'   that are zero-based (an artifact) of Python, whereas R is one-based by
#'   default.
#' @param dictionary_list A vector of strings corresponding to a Python 
#'   dictionary, present in the column \code{codon_positions} of a data.table 
#'   returned by \code{annotateSequences}.
#' @return A vector of positions of the CDR3 start codons of the BCR sequences
extractCDR3CodonStartPositions <- function(dictionary_list) {
    positions <- dictionary_list %>% sapply(toString) %>% 
        lapply(parsePythonDictionary) %>% sapply(extract, "v") %>% unlist %>% 
        as.numeric

    # partis returns zero-based positions, so add one. 
    positions <- positions + 1
    return(positions)
}

#' Get a vector of CDR3 DNA strings of an annotated dataset
#' 
#' Warning: This function is called on the raw dataset from 
#'   \code{doFullAnnotation} before \code{input_seqs} is changed to
#'   \code{mature_seq}.
#' @param dat Dataset, either a data.table or data.frame object
#' @return Vector of CDR3 strings
getCDR3s <- function(dat) {
    codon_starts <- dat$codon_positions %>% extractCDR3CodonStartPositions
    codon_ends <- codon_starts + dat$cdr3_length
    cdr3s <- dat %$% input_seqs %>% 
        substr(codon_starts, codon_ends - 1) %>% 
        unname
    return(cdr3s)
}

#' Compare levenshtein distance distributions of two CDR3 repertoires
#'
#' @param dat_a First dataset, a data.table or data.frame object
#' @param dat_b Second dataset, a data.table or data.frame object
#' @return The JS divergence of the levenshtein distance distributions of the
#'   CDR3s of the two repertoires
compareCDR3Distributions <- function(dat_a, dat_b, subsample=TRUE, 
                                     subsample_count=10000) {
    subsample_count <- min(nrow(dat_a), nrow(dat_b), subsample_count)
    dist_a <- dat_a[sample(nrow(dat_a), subsample_count), ] %$% cdr3s %>% 
        getDistanceVector
    dist_b <- dat_b[sample(nrow(dat_b), subsample_count), ] %$% cdr3s %>%
        getDistanceVector
    divergence <- getJSDivergence(dist_a, dist_b)
    return(divergence)
}

getDistancesBetweenMutationsBySequence <- function(naive, mature) {
    if(nchar(naive) != nchar(mature)) {
        stop(paste0("nchar(naive) [", nchar(naive), "] != ", "nchar(mature) [", 
                    nchar(mature), "]"))
    }
    naive_char_list <- naive %>% strsplit(split='') %>% unlist
    mature_char_list <- mature %>% strsplit(split='') %>% unlist
    mut_indices_list <- which(naive_char_list != mature_char_list)
    if(length(mut_indices_list) > 1) {
        distances <- (mut_indices_list %>% diff - 1) %>% list
    } else {
        # If there are less than two mutations, this statistic is undefined,
        # so return NA and eventually ignore
        distances <- NA
    }
    return(distances)
}

getDistancesBetweenMutations <- function(naive_list, mature_list) {
    dists <- mapply(function(x, y) {
                        if(nchar(x) == nchar(y)) {
                            res <- getDistancesBetweenMutationsBySequence(x, y)
                        } else {
                            res <- NA
                        }
                    },
                    naive_list,
                    mature_list) %>% unname %>% unlist %>% subset(!is.na(.))
    return(dists)
}

compareDistanceBetweenMutationsDistributions <- function(dat_a, dat_b) {
    dists_a <- getDistancesBetweenMutations(dat_a$naive_seq, dat_a$mature_seq)
    dists_b <- getDistancesBetweenMutations(dat_b$naive_seq, dat_b$mature_seq)
    divergence <- getJSDivergence(dists_a, dists_b)
    return(divergence)
}

getPerGeneMutationRates <- function(rate_dat) {
    rates <- rate_dat %>% sapply( function(gene) { 
                                      gene$overall_mut_rate } )
    return(rates)
}

comparePerGeneMutationRates <- function(dat_a, dat_b) {
    rates_a <- dat_a %>% getPerGeneMutationRates
    rates_b <- dat_b %>% getPerGeneMutationRates
    common_genes <- intersect(rates_a %>% names, rates_b %>% names)
    rates_a_common <- rates_a[names(rates_a) %in% 
                              common_genes][common_genes] %>% unname
    rates_b_common <- rates_b[names(rates_b) %in% 
                              common_genes][common_genes] %>% unname
    divergence <- (rates_a_common - rates_b_common) %>% abs %>% mean
    return(divergence/length(common_genes)) 
}

getPerGenePerPositionMutationRates <- function(rate_dat) {
    rates <- rate_dat %>% sapply( function(gene) {
                                      gene$mut_rate_by_position } )
    return(rates)
}

comparePerGenePerPositionMutationRates <- function(dat_a, dat_b) {
    rates_a <- dat_a %>% getPerGenePerPositionMutationRates
    rates_b <- dat_b %>% getPerGenePerPositionMutationRates
    common_genes <- intersect(rates_a %>% names, rates_b %>% names)
    rates_a_common <- rates_a[names(rates_a) %in% 
                              common_genes][common_genes] %>% unname
    rates_b_common <- rates_b[names(rates_b) %in% 
                              common_genes][common_genes] %>% unname
    divergence <- mapply(function(positions_a, positions_b) { 
                            common_positions <- intersect(positions_a %>% names, 
                                                          positions_b %>% names)
                            a_common <- positions_a[names(positions_a) %in% 
                                            common_positions][common_positions]
                            b_common <- positions_b[names(positions_b) %in% 
                                            common_positions][common_positions]
                            abs(a_common - b_common)/length(common_positions)
                         }, 
                         rates_a_common, rates_b_common) %>% unlist %>% mean
    return(divergence/length(common_genes)) 
}

removeSequencesWithDifferentNaiveAndMatureLengths <- function(dat) {
    return(dat %>% subset(nchar(dat$mature_seq) == nchar(dat$naive_seq)))
}

getSubstitutionModel <- function(dat) {
    sub_mat <- dat %>% removeSequencesWithDifferentNaiveAndMatureLengths %>%
        shazam::createSubstitutionMatrix(sequenceColumn="mature_seq",
                                         germlineColumn="naive_seq",
                                         vCallColumn="v_gene")
    return(sub_mat)
}

compareSubstitutionModels <- function(dat_a, dat_b) {
    model_a <- dat_a %>% getSubstitutionModel %>% c %>% subset(!is.na(.))
    model_b <- dat_b %>% getSubstitutionModel %>% c %>% subset(!is.na(.))
    divergence <- (model_a - model_b) %>% abs %>% mean
    return(divergence)
}

getMutabilityModel <- function(dat, 
                               substitution_model=getSubstitutionModel(dat)) {
    mut_mat <- dat %>% removeSequencesWithDifferentNaiveAndMatureLengths %>%
        shazam::createMutabilityMatrix(substitutionModel=substitution_model,
                                       sequenceColumn="mature_seq",
                                       germlineColumn="naive_seq",
                                       vCallColumn="v_gene")
    return(mut_mat)
}

compareMutabilityModels <- function(dat_a, dat_b, 
                                    sub_mod_a=getSubstitutionModel(dat_a),
                                    sub_mod_b=getSubstitutionModel(dat_b)) {
    model_a <- dat_a %>% getMutabilityModel(substitution_model=sub_mod_a)
    model_b <- dat_b %>% getMutabilityModel(substitution_model=sub_mod_b)
    divergence <- (model_a - model_b) %>% abs %>% mean
    return(divergence)
}

getDeletionLengths <- function(dat, column) {
    lengths <- dat %>% dplyr::select_(column) %>% unlist(use.names=FALSE)
    return(lengths)
}

getVGene3PrimeDeletionLengths <- function(dat) {
    return(getDeletionLengths(dat, "v_3p_del"))
}

getVGene5PrimeDeletionLengths <- function(dat) {
    return(getDeletionLengths(dat, "v_5p_del"))
}

getDGene3PrimeDeletionLengths <- function(dat) {
    return(getDeletionLengths(dat, "d_3p_del"))
}

getDGene5PrimeDeletionLengths <- function(dat) {
    return(getDeletionLengths(dat, "d_5p_del"))
}

getJGene3PrimeDeletionLengths <- function(dat) {
    return(getDeletionLengths(dat, "j_3p_del"))
}

getJGene5PrimeDeletionLengths <- function(dat) {
    return(getDeletionLengths(dat, "j_5p_del"))
}

compareDeletionLengths <- function(dat_a, dat_b, gene, end) {
    deletion_length_function <- paste0("get", gene, end, "DeletionLengths") %>%
        get
    dist_a <- dat_a %>% deletion_length_function
    dist_b <- dat_b %>% deletion_length_function
    divergence <- getJSDivergence(dist_a, dist_b)
    return(divergence)
}

compareVGene3PrimeDeletionLengths <- function(dat_a, dat_b) {
    return(compareDeletionLengths(dat_a, dat_b, "VGene", "3Prime"))
}

compareDGene3PrimeDeletionLengths <- function(dat_a, dat_b) {
    return(compareDeletionLengths(dat_a, dat_b, "DGene", "3Prime"))
}

compareDGene5PrimeDeletionLengths <- function(dat_a, dat_b) {
    return(compareDeletionLengths(dat_a, dat_b, "DGene", "5Prime"))
}

compareJGene5PrimeDeletionLengths <- function(dat_a, dat_b) {
    return(compareDeletionLengths(dat_a, dat_b, "JGene", "5Prime"))
}

getInsertionLengths <- function(dat, column) {
    lengths <- dat %>% dplyr::select_(column) %>% unlist %>% 
        sapply(toString) %>% sapply(nchar) %>% unname
    return(lengths)
}

getVDInsertionLengths <- function(dat) {
    return(getInsertionLengths(dat, "vd_insertion"))
}

getDJInsertionLengths <- function(dat) {
    return(getInsertionLengths(dat, "dj_insertion"))
}

compareInsertionLengths <- function(dat_a, dat_b, genes) {
    insertion_length_function <- paste0("get", genes, "InsertionLengths") %>% 
        get
    dist_a <- dat_a %>% insertion_length_function
    dist_b <- dat_b %>% insertion_length_function
    divergence <- getJSDivergence(dist_a, dist_b)
    return(divergence)
}

compareVDInsertionLengths <- function(dat_a, dat_b) {
    return(compareInsertionLengths(dat_a, dat_b, "VD"))
}

compareDJInsertionLengths <- function(dat_a, dat_b) {
    return(compareInsertionLengths(dat_a, dat_b, "DJ"))
}

getCloneList <- function(dat) {
    clone_list <- dat$partition %>% first %>% toString %>% strsplit(";") %>% 
        lapply(strsplit, ":") %>% first %>% lapply(as.numeric)
    return(clone_list)
}

includeClonalMemberships <- function(annotations, partitions) {
    clone_list <- partitions %>% getCloneList
    clone_df <- clone_list %>% data.table::melt()
    names(clone_df) <- c("unique_ids", "clone")
    new_df <- merge(annotations, clone_df, by="unique_ids")
    return(new_df)
}

getClusterSizes <- function(dat) {
    sizes <- dat$clone %>% table %>% unname %>% c %>% sort
    return(sizes)
}

compareClusterSizes <- function(dat_a, dat_b) {
    dist_a <- dat_a %>% getClusterSizes
    dist_b <- dat_b %>% getClusterSizes
    divergence <- getJSDivergence(dist_a, dist_b)
    return(divergence)
}

getHillNumbers <- function(dat, diversity_orders=c(0, 1, 2)) {
    counts <- dat %>% getClusterSizes
    diversity <- alakazam::calcDiversity(counts, diversity_orders)
    return(diversity)
}

#' Compare one or multiple Hill numbers of two datasets. The measure of distance
#' is a sum of absolute differences (or just the absolute difference if there is
#' only one diversity_order value
#' @param dat_a First dataset to compare
#' @param dat_b Second dataset to compare
#' @param diversity_orders Scalar- or vector-valued list of parameters to the
#' Hill diversity index. Can be any real value although nonnegative values are
#' recommended as biologically meaningful.
compareHillNumbers <- function(dat_a, dat_b, diversity_orders=c(0, 1, 2)) {
    hill_numbers_a <- dat_a %>% getHillNumbers(diversity_orders)
    hill_numbers_b <- dat_b %>% getHillNumbers(diversity_orders)
    distance <- (hill_numbers_a - hill_numbers_b) %>% abs %>% mean
    return(distance)
}
