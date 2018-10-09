library(alakazam)
library(ape)
library(Biostrings)
library(CollessLike)
library(data.table)
library(dplyr)
library(HDMD)
library(jsonlite)
library(magrittr)
library(pegas)
library(Peptides)
library(phytools)
library(RecordLinkage)
library(shazam)
library(seqinr)
library(stringdist)
library(textmineR)
library(yaml)

source("R/PartisFunctions.R")

#' Determine the comparison method for stringdistmatrix 
#'
#' The comparison method is determined based on whether 
#' all sequences are of the same length. If so, use hamming distance; 
#' else use Levenshtein.
#' @param sequence_list List/vector of DNA sequences
#' @return String-valued comparison method for use in stringdistmatrix within 
#'   the getDistanceMatrix function
determineComparisonMethod <- function(sequence_list) {
    length_count <- sequence_list %>% 
        standardizeList %>% 
        sapply(nchar) %>% 
        table %>% 
        length 
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
    sequence_list <- raw_sequences %>% 
        standardizeList
    comparison_method <- sequence_list %>% 
        determineComparisonMethod
    mat <- sequence_list %>% 
        stringdist::stringdistmatrix(method=comparison_method) %>% 
        as.matrix
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
    vec <- mat[mat %>% lower.tri] %>% 
        as.vector %>% 
        sort
    return(vec)
}

#' Get an exact or approximate distribution of pairwise distances
#'
#' @param sequence_list vector of DNA sequence strings
#' @param approximate if TRUE, approximate the distribution by subsampling
#'   and averaging
#' @return vector of integer-valued distances
getPairwiseDistanceDistribution <- function(dat,
                                            column="mature_seq",
                                            approximate=TRUE,
                                            ...
                                            ) {
    sequence_list <- dat[[column]]
    if(approximate) {
        distribution <- sequence_list %>%
            getApproximateDistribution(summary_function=getDistanceVector,
                                       divergence_function=getJSDivergence,
                                       ...
                                       )

    } else {
        distribution <- sequence_list %>%
            getDistanceVector
    }

    return(distribution)
}

plotDistribution <- function(dat_list,
                             summary_function,
                             do_exact=FALSE,
                             x_label,
                             names,
                             ...
                            ) {
    distribution_list <- dat_list %>%
        lapply(summary_function,
               ...)
    if(is.null(names)) {
        names <- paste("Dataset",
                       1:length(distribution_list)
                      )
    }
    if(do_exact) {
        # Actually plot the frequency of each value (which can be >= 200 bars)
        distance_table <- table(distribution_list)
        distribution <- distance_table/sum(distance_table) 
        d <- distribution %>%
            data.frame %>%
            setNames(c("Value", "Frequency"))
        p <- ggplot(d) +
            geom_bar(aes(x=as.numeric(Value), y=Frequency), stat="identity") 
    } else {
        # Combine distributions (of various lengths) into a common data.frame
        # with corresponding IDs
        dat <- distribution_list %>%
            Map(function(x, y) {
                    data.frame(Value=x, Dataset=y)
                },
                .,
                names
               ) %>%
            do.call("rbind", .)
        # Just plot a histogram
        p <- ggplot(dat,
                    aes(x=Value, 
                        y=c(..count..[..group..==1]/sum(..count..[..group..==1]),
                            ..count..[..group..==2]/sum(..count..[..group..==2]),
                            ..count..[..group..==3]/sum(..count..[..group..==3])
                           ),
                        group=Dataset,
                        fill=Dataset
                       )
                    ) +
            geom_histogram(alpha=0.6, position="identity")
    }
    p <- p + 
        xlab(x_label) +
        ylab("Frequency") +
        theme(legend.position="none")
    return(p)
}
                                     

plotPairwiseDistanceDistribution <- function(dat_list,
                                             do_exact=FALSE,
                                             names=NULL,
                                             ...
                                             ) {
    p <- plotDistribution(dat_list=dat_list,
                          summary_function=getPairwiseDistanceDistribution,
                          do_exact=do_exact,
                          x_label="Pairwise distance",
                          names=names,
                          ...
                         )
    return(p)
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
#' @param approximate If TRUE, uses approximate pairwise distance distributions
#' @param do_automatic If TRUE, approximate divergence using subsampling
#' @inheritParams getAutomaticAverageDivergence 
#' @return Estimated JS divergence of the distributions inferred from list_a
#'   and list_b
comparePairwiseDistanceDistributions <- function(dat_a, 
                                                 dat_b,
                                                 column="mature_seq",
                                                 do_automatic=TRUE,
                                                 approximate=TRUE,
                                                 subsample_count=100,
                                                 trial_count=10
                                                 ) {
    if(do_automatic) {
        divergence <- getAutomaticAverageDivergence(sequences_a,
                                                    sequences_b,
                                                    getDistanceVector,
                                                    subsample_count)
    } else {
        dist_a <- dat_a %>% 
            getPairwiseDistanceDistribution(approximate=approximate)
        dist_b <- dat_b %>% 
            getPairwiseDistanceDistribution(approximate=approximate)
        divergence <- getJSDivergence(dist_a, dist_b)
    }
    return(divergence)
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
getNearestNeighborDistances <- function(sequence_list, 
                                        k=1
                                        ) {
    mat <- sequence_list %>% 
        getDistanceMatrix
    n <- sequence_list %>% 
        length
    distances <- rep(NA, n)
    for(i in 1:n) {
        distances[i] <- sort(mat[i, -i], partial=k)[k]
    }

    return(distances)
}

#' Get exact or approximate nearest neighbor distribution
#'
#' @param sequence_list vector of DNA sequence strings
#' @param approximate If TRUE, approximate distribution by subsampling.
#'   Note: this inhibits interpretability since any subsample will have a 
#'   different definition of a "nearest neighbor"
#' @inheritParams getNearestNeighborDistances
#' @return vector of integer-value distances
getNearestNeighborDistribution <- function(dat,
                                           column="mature_seq", 
                                           k=1,
                                           approximate=TRUE,
                                           ...
                                           ) {
    sequence_list <- dat[[column]]
    if(approximate) {
        distribution <- sequence_list %>%
            getApproximateDistribution(summary_function=getNearestNeighborDistances,
                                       divergence_function=getJSDivergence,
                                       k=k,
                                       ...
                                       )
    } else {
        distribution <- sequence_list %>%
            getNearestNeighborDistances(k=k)
    }
    return(distribution)
}

plotNearestNeighborDistribution <- function(dat_list,
                                            do_exact=FALSE,
                                            names=NULL,
                                            ...
                                           ) {
    p <- plotDistribution(dat_list=dat_list,
                          summary_function=getNearestNeighborDistribution,
                          do_exact=do_exact,
                          x_label="Nearest neighbor distance",
                          names=names,
                          ...
                         )
    return(p)
}

#' Compare kth nearest neighbor distance distributions of two lists of
#'   sequences
#' \code{compareNNDistanceDistribution} computes the JS divergence of
#'   the kth nearest neighbor distance distribution of two lists of
#'   DNA sequences
#' @inheritParams getAutomaticAverageDivergence
#' @param dat_a First dataset, containing mature sequences
#' @param dat_b Second dataset, containing mature sequences
#' @param k The separation depth for the nearest neighbor distances.
#'   k = 1 corresponds to the nearest neighbor, k = 2 corresponds to
#'   the second-nearest neighbor, etc.
#' @inheritParams getAutomaticAverageDivergence
#' @return Estimated JS divergence of the distributions inferred from list_a
#'   and list_b
compareNNDistanceDistributions <- function(dat_a, 
                                           dat_b, 
                                           column="mature_seq",
                                           k=1,
                                           subsample_count=100,
                                           trial_count=10
                                           ) {
    average_divergence <- 
        getAutomaticAverageDivergence(dat_a[[mature_seq]],
                                      dat_b[[mature_seq]],
                                      getNearestNeighborDistances,
                                      subsample_count,
                                      k=k)
    return(average_divergence)
}

getGCContent <- function(raw_sequences) {
    sequence_list <- raw_sequences %>% 
        sapply(paste, collapse='') %>% 
        unname
    dna_list <- sequence_list %>% 
        strsplit(split='') %>% 
        lapply(ape::as.DNAbin)
    gc_dist <- dna_list %>% 
        sapply(ape::GC.content)
    return(gc_dist)
}

#' Get the GC content distribution of a list of DNA sequences
#'
#' \code{getGCContentDistribution} returns a list of the GC content of each DNA
#'   sequence in \code{raw_sequences}, given by the \code{ape} library
#' @param sequence_list List or vector of strings or character sequences
#'   corresponding to DNA sequences
#' @return A vector of GC content values
getGCContentDistribution <- function(dat,
                                     column="mature_seq",
                                     approximate=FALSE,
                                     ...
                                     ) {
    sequence_list <- dat[[column]]
    if(approximate) {
        distribution <- sequence_list %>%
            getApproximateDistribution(summary_function=getGCContent,
                                       divergence_function=getContinuousJSDivergence,
                                       ...
                                       )

    } else {
        distribution <- sequence_list %>%
            getGCContent
    }

    return(distribution)
}

plotGCContentDistribution <- function(dat_list,
                                      do_exact=FALSE,
                                      names=NULL
                                     ) {
    p <- plotDistribution(dat_list,
                          getGCContentDistribution,
                          do_exact=do_exact,
                          x_label="GC content",
                          names=names
                         )
    return(p)
}

#' Compare the GC distributions of two lists of DNA sequences
#'
#' @param list_a First list or vector of DNA sequences
#' @param list_b Second list or vector of DNA sequences
#' @return JS divergence of the GC content distributions inferred from list_a
#'   and list_b
compareGCContents <- function(dat_a, 
                              dat_b, 
                              column="mature_seq"
                             ) {
    density_a <- dat_a %>%
        getGCContentDistribution(column=column)
    density_b <- dat_b %>%
        getGCContentDistribution(column=column)
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
    dna_strings <- dna_sequences %>% 
        unlist %>% 
        Biostrings::DNAStringSet()
    count <- motif %>% 
        Biostrings::vcountPattern(dna_strings, fixed=FALSE)
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
    count <- spots %>% 
        # Get a length(dna_sequences) x length(spots) matrix of counts
        sapply(getMotifCount, dna_sequences=dna_sequences) %>%
        # Sum over each spot count for each sequence
        apply(1, sum)     
    return(count)
}

#' Get the number of occurrences of AID hotspots in a set of reference 
#'   sequences
#' 
#' @inheritParams getMotifCount
#' @return The number of AID hotspot occurrences in \code{dna_sequences}
getHotspotCount <- function(dat,
                            column,
                            hotspots
                           ) {
    return(getSpotCount(dna_sequences=dat[[column]], 
                        spots=hotspots))
}

#' Get the distribution of hotspot counts of sequences in column \code{column}
#'   of \code{dat}
getHotspotCountDistribution <- function(dat,
                                        column="mature_seq",
                                        hotspots=c("WRC", "WA"),
                                        approximate=TRUE,
                                        ...
                                       ) {
    if(approximate) {
        counts <- dat %>% 
            getApproximateDistribution(summary_function=getHotspotCount,
                                       divergence_function=getJSDivergence,
                                       column=column,
                                       hotspots=hotspots,
                                       ...
                                      )
    } else {
        counts <- dat %>% 
            getHotspotCount(column=column,
                            hotspots=hotspots
                           )
    }

    return(counts)
}

plotHotspotCountDistribution <- function(dat_list,
                                         do_exact=FALSE,
                                         names=NULL
                                        ) {
    p <- plotDistribution(dat_list,
                          getHotspotCountDistribution,
                          do_exact=do_exact,
                          x_label="Hotspot count",
                          names=names
                         )
    return(p)
}

#' Get the number of occurrences of AID coldspots in a set of reference 
#'   sequences
#' 
#' @inheritParams getMotifCount
#' @return The number of AID coldhotspot occurrences in \code{dna_sequences}
getColdspotCount <- function(dat,
                             column, 
                             coldspots
                            ) {
    return(getSpotCount(dna_sequences=dat[[column]], 
                        spots=coldspots))
}

getColdspotCountDistribution <- function(dat,
                                         column="mature_seq",
                                         coldspots="SYC",
                                         approximate=TRUE,
                                         ...
                                        ) {
    if(approximate) {
        counts <- dat %>% 
            getApproximateDistribution(summary_function=getColdspotCount,
                                       divergence_function=getJSDivergence,
                                       column=column,
                                       coldspots=coldspots,
                                       ...
                                      )
    } else {
        counts <- dat %>% 
            getColdspotCount(column=column,
                             coldspots=coldspots
                            )
    }

    return(counts)
}

plotColdspotCountDistribution <- function(dat_list,
                                          do_exact=FALSE,
                                          names=NULL
                                         ) {
    p <- plotDistribution(dat_list,
                          getColdspotCountDistribution,
                          do_exact=do_exact,
                          x_label="Coldspot count",
                          names=names
                         )
    return(p)
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

    divergence <- getAutomaticAverageDivergence(dat_a %$% mature_seq,
                                                dat_b %$% mature_seq,
                                                count_function,
                                                subsample_count,
                                                getSumOfAbsoluteDifferences,
                                                tolerance=1e-4)
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
getDistancesFromNaiveToMature <- function(dat,
                                          v_gene_only=TRUE
                                          ) {
    mature_column <- ifelse(v_gene_only, "v_qr_seqs", "mature_seq")
    naive_column <- ifelse(v_gene_only, "v_gl_seq", "naive_seq")
    distances <- dat[[mature_column]] %>% 
        mapply(FUN=stringdist::stringdist, 
               b=dat[[naive_column]], 
               method="lv") %>% 
        sort %>% 
        unname
    return(distances)
}

getDistanceFromNaiveToMatureDistribution <- function(dat,
                                                     approximate=TRUE,
                                                     v_gene_only=FALSE,
                                                     ...
                                                     ) {
    if(approximate) {
        distribution <- dat %>%
            getApproximateDistribution(summary_function=getDistancesFromNaiveToMature,
                                       divergence_function=getJSDivergence,
                                       ...
                                       )
    } else {
        distribution <- dat %>%
            getDistancesFromNaiveToMature(v_gene_only=v_gene_only)
    }
}

plotDistanceFromNaiveToMatureDistribution <- function(dat_list,
                                                      do_exact=FALSE,
                                                      names=NULL
                                                     ) {
    p <- plotDistribution(dat_list,
                          getDistanceFromNaiveToMatureDistribution,
                          do_exact=do_exact,
                          x_label="Distance from naive to mature",
                          names=names
                         )
    return(p)
}

#' Compare Levenshtein distance distributions from naive sequences to 
#  their corresponding mature ones for two repertoires
#' @param dat_a First dataset
#' @param dat_b Second dataset
#' @return JS divergence of the two distance distributions
compareDistancesFromNaiveToMature <- function(dat_a, 
                                              dat_b, 
                                              do_automatic=TRUE,
                                              approximate=TRUE,
                                              ...
                                              ) {
    if(do_automatic) {
        divergence <- 
            getAutomaticAverageDivergence(
                dat_a,
                dat_b,
                getDistanceFromNaiveToMatureDistribution,
                approximate=approximate,
                subsample_count=100,
                tolerance=1e-3,
                ...
            )
    } else {
        dist_a <- dat_a %>% 
            getDistanceFromNaiveToMatureDistribution(approximate=approximate,
                                                     ...
                                                     )
        dist_b <- dat_b %>% 
            getDistanceFromNaiveToMatureDistribution(approximate=approximate,
                                                     ...
                                                     )
        divergence <- getJSDivergence(dist_a, dist_b)
    }

    return(divergence)
}

#' Get the distribution of inferred length of CDR3 regions of each sequence
#' 
#' @param Annotated dataset
#' @return Vector of CDR3 lengths (in nt units)
getCDR3Lengths <- function(dat) {
    CDR3_lengths <- dat$cdr3_length %>% 
        na.omit
    return(CDR3_lengths)
}

plotCDR3Lengths <- function(dat_list,
                            do_exact=FALSE,
                            names=NULL
                           ) {
    p <- plotDistribution(dat_list,
                          getCDR3Lengths,
                          do_exact=do_exact,
                          x_label="CDR3 length",
                          names=names
                         )
    return(p)
}

#' Compare the distribution of CDR3 lengths for two datasets
#'
#' @param dat_a First dataset
#' @param dat_b Second dataset
#' @return The JS divergence of the two CDR3 length distributions
compareCDR3Lengths <- function(dat_a, dat_b) {
    a_lengths <- dat_a %>% 
        getCDR3Lengths
    b_lengths <- dat_b %>% 
        getCDR3Lengths
    divergence <- getJSDivergence(a_lengths, b_lengths)
    return(divergence)
}

#' Get usage table from full list of categorical variables, such as
#'   genes or amino acids, setting zero to unused ones
#' 
#' @param factor_list List of factors
#' @param full_factor_list Full list of reference factors under consideration
#' @return A table of usage counts of \code{factor_list}
getUsageTableFromFullList <- function(factor_list, full_factor_list) {
    usage_table <- factor_list %>% 
        table
    usage_table[full_factor_list %>% 
                setdiff(factor_list) %>% 
                unlist] <- 0
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
        gene_list_a <- gene_list_a %>% 
            collapseAlleles
        gene_list_b <- gene_list_b %>% 
            collapseAlleles
    }
    full_gene_list <- union(gene_list_a, gene_list_b)
    table_a <- gene_list_a %>% 
        getUsageTableFromFullList(full_gene_list)
    table_b <- gene_list_b %>% 
        getUsageTableFromFullList(full_gene_list)
    divergence <- full_gene_list %>% 
        sapply(function(x) { abs(table_a[x] - table_b[x]) }) %>% 
        sum
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
        v_genes <- dat %$% 
            v_gene %>% 
            collapseAlleles
        d_genes <- dat %$% 
            d_gene %>% 
            collapseAlleles
        j_genes <- dat %$% 
            j_gene %>% 
            collapseAlleles
        gene_dat <- data.table(v_gene=v_genes,
                               d_gene=d_genes,
                               j_gene=j_genes)
    } else {
        gene_dat <- dat
    }

    gene_type_list <- c("v_gene", "d_gene", "j_gene")
    gene_table <- gene_dat %>% 
        plyr::count(gene_type_list)
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
    
    divergence <- summands %>% 
        sum
    return(divergence)
}

getKideraFactorsBySequence <- function(sequence) {
    kidera_factors <- sequence %>% 
        Peptides::kideraFactors() %>% 
        first
    return(kidera_factors)
}

#' Get hydrophobicity from amino acid sequence, corresponding to
#' the fourth kidera factor
#'
#' @param sequence Amino acid sequence string
#' @return real-valued hydrophobicity estimate, computed as the average
#'   hydrophobicity of each amino acid in the sequence
getHydrophobicityFromAASequence <- function(sequence) {
    if(!is.na(sequence)) {
        hydrophobicity <- sequence %>%
            getKideraFactorsBySequence %>%
            t %>%
            as.data.table %$%
            KF4
    } else {
        hydrophobicity <- NA
    }

    return(hydrophobicity)
}

getHydrophobicityDistribution <- function(dat,
                                          include_NA=FALSE
                                          ) {
    hydrophobicity_list <- dat %$%
        cdr3_aa %>% 
        filterAminoAcidSequences %>%
        sapply(getHydrophobicityFromAASequence)

    return(hydrophobicity_list)
}

plotHydrophobicityDistribution <- function(dat_list,
                                           do_exact=FALSE,
                                           names=NULL
                                          ) {
    p <- plotDistribution(dat_list,
                          getHydrophobicityDistribution,
                          do_exact=do_exact,
                          x_label="Hydrophobicity",
                          names=names
                         )
    return(p)
}

compareHydrophobicityDistributions <- function(dat_a, dat_b) {
    divergence <- 
        getAutomaticAverageDivergence(
            dat_a,
            dat_b,
            getHydrophobicityDistribution,
            subsample_count=100,
            divergenceFunction=getContinuousJSDivergence,
            tolerance=1e-3
        )
    return(divergence)
}

#' Get the mean Atchley factor of an amino acid sequence
#'
#' @param aa_seq String of an amino acid sequence
#' @param factor_number The Atchley factor to be applied to \code{aa_seq}. 
#'  Must be 1, 2, 3, 4, or 5.
#' @return The mean of the set of Atchley factors for 
#'   \code{aa_seq}
getAtchleyFactorList <- function(aa_seq, factor_number) {
    factor_list <- aa_seq %>% 
        HDMD::FactorTransform(Factor=factor_number)
    return(factor_list)
}

#' Get a list of the means of each of the five Atchley factors for 
#'   a list of DNA seuqnces
#'
#' @param sequence_list List or vector of DNA sequences
#' @return Vector of means of each of the five Atchley factors of \code{sequence_list}
getMeanAtchleyFactorDistribution <- function(sequence_list) {
    collapsed_sequence <- sequence_list %>% 
        filterAminoAcidSequences %>%
        paste(collapse='')
    factor_numbers <- 1:5
    mean_atchley_factors <- rep(NA, length(factor_numbers))
    for(factor_number in factor_numbers) {
        mean_atchley_factors[factor_number] <- collapsed_sequence %>% 
            getAtchleyFactorList(factor_number=factor_number) %>%
            first %>%
            mean
    }

    return(mean_atchley_factors)
}

#' Compare the distributions of mean Atchley factors for two datasets
#'
#' @param dat_a First dataset, a data.table or data.frame
#' @param dat_b Second dataset, a data.table or data.frame
#' @return Estimated mean absolute difference of the mean of each Atchley factor
#'   of \code{dat_a} and \code{dat_b}
compareAtchleyFactorDistributions <- function(dat_a, dat_b) {
    divergence <- getAutomaticAverageDivergence(
        dat_a %$% cdr3_aa,
        dat_b %$% cdr3_aa,
        getMeanAtchleyFactorDistribution,
        subsample_count=100,
        divergenceFunction=getSumOfAbsoluteDifferences)
    return(divergence)
}

#' Get the aliphatic index of a DNA sequence
#'
#' @param dna_sequence String of DNA characters
#' @return The aliphatic index of \code{dna_sequence}, if applicable
getAliphaticIndex <- function(aa_sequence) {
    aliphatic_index <- aa_sequence %>% 
               Peptides::aIndex() %>%
               ifelse(. < 125, ., NA) # Don't bother with outliers

    return(aliphatic_index)
}

#' Get the distribution of aliphatic indices of a list of sequences
#'
#' @param sequence_list List or vector of DNA sequences
#' @return Vector of aliphatic indices
getAliphaticIndexDistribution <- function(dat) {
    a_indices <- dat %$%
        cdr3_aa %>% 
        filterAminoAcidSequences %>%
        sapply(function(x) {
                   ifelse(!is.na(x),
                          getAliphaticIndex(x),
                          NA)
        }
        )

    return(a_indices)
}

plotAliphaticIndexDistribution <- function(dat_list,
                                           do_exact=FALSE,
                                           names=NULL
                                          ) {
    p <- plotDistribution(dat_list,
                          getAliphaticIndexDistribution,
                          do_exact=do_exact,
                          x_label="Aliphatic index",
                          names=names
                         )
    return(p)
}

#' Compare the distributions of aliphatic indices of two datasets
#'
#' @param dat_a First dataset, a data.table or data.frame
#' @param dat_b Second datset, a data.table or data.frame
#' @return The JS divergence of the two distributions
compareAliphaticIndexDistributions <- function(dat_a, dat_b) {
    divergence <- getAutomaticAverageDivergence(
        dat_a,
        dat_b,
        getAliphaticIndexDistribution,
        subsample_count=100,
        divergenceFunction=getContinuousJSDivergence)
    return(divergence)
}

#' Get the distribution of GRAVY values from a list or vector of sequences
#'
#' @param List or vector of DNA sequence strings
#' @return Vector of GRAVY values for \code{sequence_list}
getGRAVYDistribution <- function(dat) {
    dist <- dat %$%
            cdr3_aa %>% 
            filterAminoAcidSequences %>%
            sapply(function(x) {
                   ifelse(!is.na(x),
                          alakazam::gravy(x),
                          NA)
                   }
            ) %>% 
            unname
    return(dist)
}

plotGRAVYDistribution <- function(dat_list,
                                  do_exact=FALSE,
                                  names=NULL
                                 ) { 
    p <- plotDistribution(dat_list,
                          getGRAVYDistribution,
                          do_exact=do_exact,
                          x_label="GRAVY index",
                          names=names
                         )
    return(p)
}


#' Compare the GRAVY distributions of two datasets
#'
#' @param dat_a First dataset, a data.table or data.frame object
#' @param dat_b Second dataset, a data.table or data.frame object
#' @return The JS divergence of GRAVY distributions
compareGRAVYDistributions <- function(dat_a, dat_b) {
    dist_a <- dat_a %>% 
        getGRAVYDistribution
    dist_b <- dat_b %>%
        getGRAVYDistribution
    divergence <- getJSDivergence(dist_a, dist_b, continuous=TRUE)
    return(divergence)
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
    positions <- dictionary_list %>% 
        sapply(toString) %>% 
        lapply(parsePythonDictionary) %>% 
        sapply(extract, "v") %>% 
        unlist %>% 
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
    codon_starts <- dat %$%
        codon_positions %>% 
        extractCDR3CodonStartPositions
    codon_ends <- codon_starts + dat$cdr3_length
    cdr3s <- dat %$% 
        input_seqs %>% 
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
compareCDR3Distributions <- function(dat_a, 
                                     dat_b, 
                                     subsample=TRUE, 
                                     subsample_count=100,
                                     trial_count=10) {
    divergence <- getAutomaticAverageDivergence(dat_a %$% cdr3s,
                                                dat_b %$% cdr3s,
                                                getDistanceVector,
                                                subsample_count)
    return(divergence)
}

getDistancesBetweenMutationsBySequence <- function(naive, mature) {
    if(nchar(naive) != nchar(mature)) {
        stop(paste0("nchar(naive) [", nchar(naive), "] != ", "nchar(mature) [", 
                    nchar(mature), "]"))
    }
    naive_char_list <- naive %>% 
        strsplit(split='') %>% 
        unlist
    mature_char_list <- mature %>% 
        strsplit(split='') %>% 
        unlist
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

getDistancesBetweenMutations <- function(dat) {
    dists <- mapply(function(x, y) {
                        if(nchar(x) == nchar(y)) {
                            res <- getDistancesBetweenMutationsBySequence(x, y)
                        } else {
                            res <- NA
                        }
                    },
                    dat$naive_seq,
                    dat$mature_seq
                   ) %>% 
        unname %>% 
        unlist %>% 
        subset(!is.na(.))
    return(dists)
}

plotDistanceBetweenMutationsDistribution <- function(dat_list,
                                                     do_exact=FALSE,
                                                     names=NULL
                                                    ) { 
    p <- plotDistribution(dat_list,
                          getDistancesBetweenMutations,
                          do_exact=do_exact,
                          x_label="Positional distance between mutations",
                          names=names
                         )
    return(p)
}


compareDistanceBetweenMutationsDistributions <- function(dat_a, dat_b) {
    dists_a <- getDistancesBetweenMutations(dat_a)
    dists_b <- getDistancesBetweenMutations(dat_b)
    divergence <- getJSDivergence(dists_a, dists_b)
    return(divergence)
}

getPerGeneMutationRates <- function(rate_dat) {
    rates <- rate_dat %>% 
        sapply( function(gene) { 
                    gene$overall_mut_rate 
                }
        )
    return(rates)
}

comparePerGeneMutationRates <- function(dat_a, dat_b) {
    rates_a <- dat_a %>% 
        getPerGeneMutationRates
    rates_b <- dat_b %>% 
        getPerGeneMutationRates
    common_genes <- intersect(rates_a %>% names, 
                              rates_b %>% names)
    rates_a_common <- rates_a[names(rates_a) %in% 
                              common_genes][common_genes] %>% 
        unname
    rates_b_common <- rates_b[names(rates_b) %in% 
                              common_genes][common_genes] %>% 
        unname
    divergence <- (rates_a_common - rates_b_common) %>% 
        abs %>% 
        mean
    return(divergence/length(common_genes)) 
}

getPerGenePerPositionMutationRates <- function(rate_dat) {
    rates <- rate_dat %>% 
        sapply( function(gene) {
                    gene$mut_rate_by_position 
                } 
        )
    return(rates)
}

comparePerGenePerPositionMutationRates <- function(dat_a, dat_b) {
    rates_a <- dat_a %>% 
        getPerGenePerPositionMutationRates
    rates_b <- dat_b %>% 
        getPerGenePerPositionMutationRates
    common_genes <- intersect(rates_a %>% names, 
                              rates_b %>% names)
    rates_a_common <- rates_a[names(rates_a) %in% 
                              common_genes][common_genes] %>% 
        unname
    rates_b_common <- rates_b[names(rates_b) %in% 
                              common_genes][common_genes] %>% 
        unname
    divergence <- mapply(function(positions_a, positions_b) { 
                            common_positions <- intersect(positions_a %>% names, 
                                                          positions_b %>% names)
                            a_common <- positions_a[names(positions_a) %in% 
                                            common_positions][common_positions]
                            b_common <- positions_b[names(positions_b) %in% 
                                            common_positions][common_positions]
                            abs(a_common - b_common)/length(common_positions)
                         }, 
                         rates_a_common, rates_b_common) %>% 
        unlist %>% 
        mean
    return(divergence/length(common_genes)) 
}

#' Get the inferred substitution model of a repertoire
#'
#' @param dat Annotated dataset
#' @return The default shazam substitution model 
getSubstitutionModel <- function(dat) {
    sub_mat <- dat %>% 
        removeSequencesWithDifferentNaiveAndMatureLengths %>%
        shazam::createSubstitutionMatrix(sequenceColumn="mature_seq",
                                         germlineColumn="naive_seq",
                                         vCallColumn="v_gene") 
    return(sub_mat)
}

#' Get the inferred mutability model of a repetoire
#'
#' This requires the inferred substitution model of \code{dat} as input
#' @param dat Annotated dataset
#' @param substitution_model The inferred substitution model, which can be
#'   obtained via \code{getSubstitutionModel}
#' @return The default shazam mutability model
getMutabilityModel <- function(dat, 
                               substitution_model) {
    mut_mat <- dat %>% 
        removeSequencesWithDifferentNaiveAndMatureLengths %>%
        shazam::createMutabilityMatrix(substitutionModel=substitution_model,
                                       sequenceColumn="mature_seq",
                                       germlineColumn="naive_seq",
                                       vCallColumn="v_gene")
    return(mut_mat)
}

compareMutabilityModels <- function(dat_a, 
                                    dat_b, 
                                    sub_mod_a=getSubstitutionModel(dat_a),
                                    sub_mod_b=getSubstitutionModel(dat_b)
                                    ) {
    model_a <- dat_a %>% 
        getMutabilityModel(substitution_model=sub_mod_a)
    model_b <- dat_b %>% 
        getMutabilityModel(substitution_model=sub_mod_b)
    divergence <- getSumOfAbsoluteDifferences(model_a, model_b) 
    return(divergence)
}

#' Compare the substitution and mutability models of two datasets
#'
#' @param dat_a First dataset
#' @param dat_b Second dataset
#' @return List of two divergences: one for the substitution models and 
#'   another for the mutability models
compareSubstitutionAndMutabilityModels <- function(dat_a, dat_b) {
    sub_model_a <- dat_a %>% 
        getSubstitutionModel
    sub_model_b <- dat_b %>% 
        getSubstitutionModel
    sub_divergence <- getSumOfAbsoluteDifferences(sub_model_a, sub_model_b)

    mut_model_a <- dat_a %>% 
        getMutabilityModel(substitution_model=sub_model_a)
    mut_model_b <- dat_b %>% 
        getMutabilityModel(substitution_model=sub_model_b)
    mut_divergence <- getSumOfAbsoluteDifferences(mut_model_a, mut_model_b)

    divergences <- list(SubstitutionModel=sub_divergence,
                        MutabilityModel=mut_divergence)
    return(divergences)
}

getDeletionLengths <- function(dat, column) {
    lengths <- dat %>% 
        dplyr::select_(column) %>% 
        unlist(use.names=FALSE)
    return(lengths)
}

getVGene3PrimeDeletionLengths <- function(dat) {
    return(getDeletionLengths(dat, "v_3p_del"))
}

plotVGene3PrimeDeletionLengths <- function(dat_list,
                                           do_exact=FALSE,
                                           names=NULL
                                          ) { 
    p <- plotDistribution(dat_list,
                          getVGene3PrimeDeletionLengths,
                          do_exact=do_exact,
                          x_label="V gene 3' deletion length",
                          names=names
                         )
    return(p)
}

getVGene5PrimeDeletionLengths <- function(dat) {
    return(getDeletionLengths(dat, "v_5p_del"))
}

plotVGene5PrimeDeletionLengths <- function(dat_list,
                                           do_exact=FALSE,
                                           names=NULL
                                          ) { 
    p <- plotDistribution(dat_list,
                          getVGene5PrimeDeletionLengths,
                          do_exact=do_exact,
                          x_label="V gene 5' deletion length",
                          names=names
                         )
    return(p)
}

getDGene3PrimeDeletionLengths <- function(dat) {
    return(getDeletionLengths(dat, "d_3p_del"))
}

plotDGene3PrimeDeletionLengths <- function(dat_list,
                                           do_exact=FALSE,
                                           names=NULL
                                          ) { 
    p <- plotDistribution(dat_list,
                          getDGene3PrimeDeletionLengths,
                          do_exact=do_exact,
                          x_label="D gene 3' deletion length",
                          names=names
                         )
    return(p)
}

getDGene5PrimeDeletionLengths <- function(dat) {
    return(getDeletionLengths(dat, "d_5p_del"))
}

plotDGene5PrimeDeletionLengths <- function(dat_list,
                                           do_exact=FALSE,
                                           names=NULL
                                          ) { 
    p <- plotDistribution(dat_list,
                          getDGene5PrimeDeletionLengths,
                          do_exact=do_exact,
                          x_label="D gene 5' deletion length",
                          names=names
                         )
    return(p)
}

getJGene3PrimeDeletionLengths <- function(dat) {
    return(getDeletionLengths(dat, "j_3p_del"))
}

plotJGene3PrimeDeletionLengths <- function(dat_list,
                                           do_exact=FALSE,
                                           names=NULL
                                          ) { 
    p <- plotDistribution(dat_list,
                          getJGene3PrimeDeletionLengths,
                          do_exact=do_exact,
                          x_label="J gene 3' deletion length",
                          names=names
                         )
    return(p)
}

getJGene5PrimeDeletionLengths <- function(dat) {
    return(getDeletionLengths(dat, "j_5p_del"))
}

plotJGene5PrimeDeletionLengths <- function(dat_list,
                                           do_exact=FALSE,
                                           names=NULL
                                          ) { 
    p <- plotDistribution(dat_list,
                          getJGene5PrimeDeletionLengths,
                          do_exact=do_exact,
                          x_label="J gene 5' deletion length",
                          names=names
                         )
    return(p)
}

compareDeletionLengths <- function(dat_a, dat_b, gene, end) {
    deletion_length_function <- paste0("get", gene, end, "DeletionLengths") %>%
        get
    dist_a <- dat_a %>% 
        deletion_length_function
    dist_b <- dat_b %>% 
        deletion_length_function
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
    lengths <- dat %>% 
        dplyr::select_(column) %>% 
        unlist %>% 
        sapply(toString) %>% 
        sapply(nchar) %>% 
        unname
    return(lengths)
}

getVDInsertionLengths <- function(dat) {
    return(getInsertionLengths(dat, "vd_insertion"))
}

plotVDInsertionLengths <- function(dat_list,
                                   do_exact=FALSE,
                                   names=NULL
                                  ) { 
    p <- plotDistribution(dat_list,
                          getVDInsertionLengths,
                          do_exact=do_exact,
                          x_label="VD insertion length",
                          names=names
                         )
    return(p)
}

getDJInsertionLengths <- function(dat) {
    return(getInsertionLengths(dat, "dj_insertion"))
}

plotDJInsertionLengths <- function(dat_list,
                                   do_exact=FALSE,
                                   names=NULL
                                  ) { 
    p <- plotDistribution(dat_list,
                          getDJInsertionLengths,
                          do_exact=do_exact,
                          x_label="DJ insertion length",
                          names=names
                         )
    return(p)
}

compareInsertionLengths <- function(dat_a, dat_b, genes) {
    insertion_length_function <- paste0("get", genes, "InsertionLengths") %>% 
        get
    dist_a <- dat_a %>% 
        insertion_length_function
    dist_b <- dat_b %>% 
        insertion_length_function
    divergence <- getJSDivergence(dist_a, dist_b)
    return(divergence)
}

compareVDInsertionLengths <- function(dat_a, dat_b) {
    return(compareInsertionLengths(dat_a, dat_b, "VD"))
}

compareDJInsertionLengths <- function(dat_a, dat_b) {
    return(compareInsertionLengths(dat_a, dat_b, "DJ"))
}

#' Constructs a Markov matrix for a vector of DNA sequences
#'
#' Here, a transition is simply the next character in a sequence, so that
#' the sequence ACGT represents the transitions A -> C -> G -> T. The function
#' would record 1 A -> C transition, 1 C -> G transition, and 1 G -> T 
#' transition. Each sequence is assumed to follow the same transition matrix.
#' Any transition to or from a character other than a, c, g, or t is ignored.
#' @param seq_list Vector of DNA sequence strings
#' @return The empirical transition matrix for each (base, base) pair. 
getMarkovMatrix <- function(seq_list) {
    counts <- matrix(0, 4, 4)
    dna_chars <- c('a', 'c', 'g', 't')
    rownames(counts) <- dna_chars
    colnames(counts) <- dna_chars
    for(seq in seq_list) {
        char_vector <- seq %>% 
            toString %>% 
            strsplit('') %>% 
            unlist
        vec_length <- length(char_vector)
        for(i in 1:vec_length) {
            char <- char_vector[i]
            next_char <- char_vector[i + 1]
            if(char %in% dna_chars && next_char %in% dna_chars) {
                counts[char, next_char] <- counts[char, next_char] + 1
            }
        }
    }
    probs <- counts/sum(counts)
    return(probs)
}

#' Get the Markov transition matrix for either VD or DJ insertions
#'
#' @param dat Data.table with inferred insertion sequences
#' @param column The name of the column corresponding to the inserted sequences
#' @return The empirical transition matrix for each (base, base) pair.
getInsertionMatrix <- function(dat, column) {
    mat <- dat %>%
        dplyr::select_(column) %>%
        sapply(tolower) %>%
        getMarkovMatrix
    return(mat)
}

#' Get the Markov transition matrix for VD insertions
#'
#' @inheritParams getInsertionMatrix
#' @return The empirical transition matrix for each (base, base) pair.
getVDInsertionMatrix <- function(dat) {
    return(getInsertionMatrix(dat, "vd_insertion"))
}

#' Get the Markov transition matrix for DJ insertions
#'
#' @inheritParams getInsertionMatrix
#' @return The empirical transition matrix for each (base, base) pair.
getDJInsertionMatrix <- function(dat) {
    return(getInsertionMatrix(dat, "dj_insertion"))
}

#' Compare the transition matrices for VD insertions for two datasets
#'
#' @param dat_a First data.table with insertion annotations
#' @param dat_b First data.table with insertion annotations
#' @return The mean absolute difference of matrix entries, taken elementwise
compareVDInsertionMatrices <- function(dat_a, dat_b) {
    matrix_a <- dat_a %>% 
        getVDInsertionMatrix
    matrix_b <- dat_b %>% 
        getVDInsertionMatrix
    divergence <- getSumOfAbsoluteDifferences(matrix_a, matrix_b)
    return(divergence)
}

#' Compare the transition matrices for DJ insertions for two datasets
#'
#' @param dat_a First data.table with insertion annotations
#' @param dat_b First data.table with insertion annotations
#' @return The mean absolute difference of matrix entries, taken elementwise
compareDJInsertionMatrices <- function(dat_a, dat_b) {
    matrix_a <- dat_a %>% 
        getDJInsertionMatrix
    matrix_b <- dat_b %>% 
        getDJInsertionMatrix
    divergence <- getSumOfAbsoluteDifferences(matrix_a, matrix_b)
    return(divergence)
}

getClusterSizes <- function(dat) {
    print(dat$clone %>% length)
    sizes <- dat %$%
        clone %>% 
        table %>% 
        unname %>% 
        c %>% 
        sort
    return(sizes)
}

plotClusterSizeDistribution <- function(dat_list,
                                        do_exact=FALSE,
                                        names=NULL
                                       ) { 
    p <- plotDistribution(dat_list,
                          getClusterSizes,
                          do_exact=do_exact,
                          x_label="Cluster size",
                          names=names
                         )
    return(p)
}

compareClusterSizes <- function(dat_a, dat_b) {
    dist_a <- dat_a %>% 
        getClusterSizes
    dist_b <- dat_b %>% 
        getClusterSizes
    divergence <- getJSDivergence(dist_a, dist_b)
    return(divergence)
}

getHillNumbers <- function(dat, diversity_orders=c(0, 1, 2)) {
    counts <- dat %>% 
        getClusterSizes
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
    hill_numbers_a <- dat_a %>% 
        getHillNumbers(diversity_orders)
    hill_numbers_b <- dat_b %>% 
        getHillNumbers(diversity_orders)
    distance <- (hill_numbers_a - hill_numbers_b) %>% 
        abs %>% 
        mean
    return(distance)
}

#' Get the percentage of sequences in a dataset which, after VDJ recombination 
#'   and SHM, have the start of the CDR3 in line with the start of the germline 
#'   V sequence
#'
#' @param dat Annotated dataset
#' @return The percentage of in-frame sequences in \code{dat}
getInFramePercentage <- function(dat) {
    percentage <- 100*(dat %$%
                       in_frames %>% 
                       mean)
    return(percentage)
}

#' Compare the percentage of in-frame sequences of two datasets
#'
#' @param dat_a First dataset to compare
#' @param dat_b Second dataset to compare
#' @return Absolute differences of in-frame percentages
compareInFramePercentages <- function(dat_a, dat_b) {
    percent_a <- dat_a %>%
        getInFramePercentage
    percent_b <- dat_b %>%
        getInFramePercentage
    divergence <- abs(percent_a - percent_b)
    return(divergence)
}

#' Get selection estimates via the shazam Baseline models
#'
#' shazam::summarizeBaseline appends mutliple columns to \code{dat}, including
#'   BASELINE_SIGMA, which is the estimated selection strength of the 
#'   repertoire
#' @param dat Annotated dataset
#' @return Distribution of baseline values, which are measures of selection 
#'   strength. Thus, Sigma = 0 corresponds to no selection, Sigma > 0 
#'   corresponds to positive selection, and Sigma < 0 corresponds to negative
#'   selection
getSelectionEstimate <- function(dat) {
    baseline <- shazam::calcBaseline(dat,
                                     sequenceColumn="mature_seq",
                                     germlineColumn="naive_seq"
                                     ) %>%
        shazam::summarizeBaseline(returnType="df") %$%
        BASELINE_SIGMA 
    return(baseline)
}

#' Compare selection estimates for two datasets
#'
#' @param dat_a First dataset
#' @param dat_b Second datset
#' @return The estimated JS divergence of selection strength distributions
compareSelectionEstimates <- function(dat_a, dat_b) {
    baseline_a <- dat_a %>% 
        getSelectionEstimate
    baseline_b <- dat_b %>% 
        getSelectionEstimate
    divergence <- getJSDivergence(baseline_a, baseline_b, continuous=TRUE)
    return(divergence)
}

#' Compute Sackin's index for the given tree
#'
#' @param tree A phylo object corresponding to a phylogeny
#' @return Sackin's index of tree balance
getSackinIndex <- function(tree) {
    index <- tree %>% 
        CollessLike::sackin.index(norm=TRUE)
    return(index)
}

#' Compare Sackin's indices for two phylogenetic trees
#'
#' @param tree_1 The first phylo object
#' @param tree_2 The second phylo object
#' @return The absolute difference of the respective Sackin's indices
compareSackinIndices <- function(tree_1, tree_2) {
    index_1 <- tree_1 %>%
        getSackinIndex
    index_2 <- tree_2 %>%
        getSackinIndex
    difference <- abs(index_1 - index_2)
    return(difference)
}

#' Compute Colless-like index for the given tree
#'
#' @param tree A phylo object corresponding to a phylogeny
#' @return Colless-like index of tree balance
getCollessLikeIndex <- function(tree) {
    index <- tree %>%
        CollessLike::colless.like.index(norm=TRUE)
    return(index)
}

#' Compare Colless-like indices for two phylogenetic trees
#'
#' @param tree_1 The first phylo object
#' @param tree_2 The second phylo object
#' @return The absolute difference of the respective Colless-like indices
compareCollessLikeIndices <- function(tree_1, tree_2) {
    index_1 <- tree_1 %>%
        getCollessLikeIndex
    index_2 <- tree_2 %>%
        getCollessLikeIndex
    difference <- abs(index_1 - index_2)
    return(difference)
}

#' Compute cophenetic index for the given tree
#'
#' @param tree A phylo object corresponding to a phylogeny
#' @return cophenetic index of tree balance
getCopheneticIndex <- function(tree) {
    index <- tree %>%
        CollessLike::cophen.index(norm=TRUE) 
    return(index)
}

#' Compare cophenetic indices for two phylogenetic trees
#'
#' @param tree_1 The first phylo object
#' @param tree_2 The second phylo object
#' @return The absolute difference of the respective cophenetic indices
compareCopheneticIndices <- function(tree_1, tree_2) {
    index_1 <- tree_1 %>%
        getCopheneticIndex
    index_2 <- tree_2 %>%
        getCopheneticIndex
    difference <- abs(index_1 - index_2)
    return(difference)
}

#' Compute a table of counts for each amino acid in the list of sequences
#'
#' @param sequences Vector of amino acid sequences
#' @return table of counts for each amino acid in \code{sequences}
getAminoAcidDistribution <- function(dat,
                                     standardize=TRUE
                                    ) {
    sequences <- dat %$% cdr3_aa
    aa_dist <- paste(sequences[!is.na(sequences)], collapse="") %>%
        strsplit(split="") %>%
        unlist %>%
        table

    if(standardize) {
        aa_dist <- aa_dist/sum(aa_dist)
    }

    return(aa_dist)
}

compareCategoricalDistributions <- function(d1,
                                          d2) {
    full_names <- union(names(d1), names(d2))
    missing_1 <- setdiff(full_names, names(d1))
    d1[missing_1] <- 0
    missing_2 <- setdiff(full_names, names(d2))
    d2[missing_2] <- 0

    divergence <- abs(d1[sort(names(d1))] - d2[sort(names(d2))]) %>% sum
    return(divergence)
}

compareAminoAcidDistributions <- 
    function(dat_a,
             dat_b,
             aa_dist_1=getAminoAcidDistribution(dat_a),
             aa_dist_2=getAminoAcidDistribution(dat_b)) 
{
    divergence <- compareCategoricalDistributions(aa_dist_1, aa_dist_2)
    return(divergence)    
}

getSequence2mers <- function(sequence) {
    split_seq <- sequence %>% 
        strsplit(split="") %>%
        unlist
    
    seq_2mers <- {}
    seq_length <- length(split_seq)
    for(i in 1:(seq_length - 1)) {
        seq_2mers[i] <- paste(split_seq[i:(i + 1)], collapse="")
    }
    return(seq_2mers)
}

getAminoAcid2merDistribution <- function(dat,
                                         standardize=TRUE
                                        ) {
    sequences <- dat %$% cdr3_aa
    dist <- sequences[!is.na(sequences)] %>%
        lapply(getSequence2mers) %>%
        unlist %>%
        table

    if(standardize) {
        dist <- dist/sum(dist)
    }

    return(dist)
}

compareAminoAcid2merDistributions <- 
    function(dat_a,
             dat_b,
             aa_dist_1=getAminoAcid2merDistribution(dat_a),
             aa_dist_2=getAminoAcid2merDistribution(dat_b)
            )  {
    divergence <- compareCategoricalDistributions(aa_dist_1, aa_dist_2)
    return(divergence)    
}

#' @param dat List of annotation datasets for which to plot univariate summary 
#'   distributions
#' @param tall_plot Make the plot portrait-oriented rather than landscape
plotUnivariateDistributions <- function(dat_list,
                                        tall_plot=FALSE,
                                        do_exact=FALSE,
                                        names=NULL
                                       ) {
    plot_function_strings <- list("plotPairwiseDistanceDistribution",
                                  "plotNearestNeighborDistribution",
                                  "plotGCContentDistribution",
                                  "plotHotspotCountDistribution",
                                  "plotColdspotCountDistribution",
                                  "plotDistanceFromNaiveToMatureDistribution",
                                  "plotCDR3Lengths",
                                  "plotHydrophobicityDistribution",
                                  "plotAliphaticIndexDistribution",
                                  "plotGRAVYDistribution",
                                  "plotDistanceBetweenMutationsDistribution",
                                  "plotVGene3PrimeDeletionLengths",
                                  "plotDGene3PrimeDeletionLengths",
                                  "plotDGene5PrimeDeletionLengths",
                                  "plotJGene5PrimeDeletionLengths",
                                  "plotVDInsertionLengths",
                                  "plotDJInsertionLengths",
                                  "plotClusterSizeDistribution"
                                 )
    plots <- {}
    for(f_string in plot_function_strings) {
        plot_function <- eval(parse(text=f_string))
        plots[[f_string]] <- dat_list %>% 
            plot_function(do_exact=do_exact,
                          names=names
                          )
    }
    multiplot <- multiplot(plotlist=plots, tall_plot=tall_plot)
}
