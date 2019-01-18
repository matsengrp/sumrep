require(alakazam)
require(ape)
require(Biostrings)
require(CollessLike)
require(data.table)
require(dplyr)
require(HDMD)
require(jsonlite)
require(magrittr)
require(pegas)
require(Peptides)
require(RecordLinkage)
require(shazam)
require(seqinr)
require(stringdist)
require(textmineR)
require(yaml)

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
#' @param sequence_list vector of sequence strings (DNA, AA, etc.)
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
#' @param dat A \code{data.table} corresponding to repertoire annotations
#' @param column the column name of \code{dat} containing the strings on which
#'   the distribution should be computed
#' @param approximate if TRUE, approximate the distribution by subsampling
#'   and averaging
#' @return vector of integer-valued distances
getPairwiseDistanceDistribution <- function(dat,
                                            column="sequence",
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

#' Plot a summary distribution of one or more datasets
#'
#' @param dat_list A \code{list} of \code{data.table} objects corresponding to 
#'    repertoire annotations
#' @param summary_function A function that is applied to each dataset in 
#'   \code{dat_list} and whose output values are plotted
#' @param do_exact Display exact distribution plots rather than histograms
#' @param x_label The text label for the x-axis
#' @param names Strings to be displayed by the legend corresponding to the 
#'   elements of \code{dat_list}
plotDistribution <- function(dat_list,
                             summary_function,
                             plot_type,
                             x_label="Value",
                             names,
                             show_legend=TRUE,
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
    if(plot_type == "barplot") {
        # Actually plot the frequency of each value (which can be >= 200 bars)
        distance_table <- table(distribution_list)
        distribution <- distance_table/sum(distance_table) 
        d <- distribution %>%
            data.frame %>%
            setNames(c("Value", "Frequency"))
        p <- ggplot(d) +
            geom_bar(aes(x=as.numeric(Value), y=Frequency), stat="identity") 
    } else if(plot_type == "histogram") {
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
    } else if(plot_type == "freqpoly") {
        dat <- distribution_list %>%
            Map(function(x, y) {
                    data.frame(Value=x, Dataset=y)
                },
                .,
                names
               ) %>%
            do.call("rbind", .)

        p <- ggplot(dat,
                    aes(x=Value,
                        y=..density..,
                        group=Dataset,
                        colour=Dataset
                        )
                    ) +
            geom_freqpoly()
    }
    p <- p + 
        theme(panel.background=element_blank(),
                  panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),
                  axis.line=element_line(colour="grey")) +
        xlab(x_label) +
        ylab("Frequency")
    if(!show_legend) {
        p <- p + theme(legend.position="none")
    }
    return(p)
}
                                     

#' Plot the pairwise distance distribution of one or more datasets
#'
#' @inheritParams plotDistribution
plotPairwiseDistanceDistribution <- function(dat_list,
                                             plot_type,
                                             names=NULL,
                                             show_legend=TRUE,
                                             ...
                                            ) {
    p <- plotDistribution(dat_list=dat_list,
                          summary_function=getPairwiseDistanceDistribution,
                          plot_type=plot_type,
                          x_label="Pairwise distance",
                          names=names,
                          show_legend=show_legend,
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
#' @param dat_a,dat_b A \code{data.table} corresponding to repertoire annotations
#' @param column the column name of \code{dat} containing the strings on which
#'   the distribution should be computed
#' @param approximate If TRUE, uses approximate pairwise distance distributions
#' @return Estimated JS divergence of the distributions inferred from list_a
#'   and list_b
comparePairwiseDistanceDistributions <- function(dat_a, 
                                                 dat_b,
                                                 column="sequence",
                                                 approximate=TRUE,
                                                 ...
                                                ) {
    dist_a <- dat_a %>% 
        getPairwiseDistanceDistribution(approximate=approximate,
                                        column=column,
                                        ...)
    dist_b <- dat_b %>% 
        getPairwiseDistanceDistribution(approximate=approximate,
                                        column=column,
                                        ...)
    divergence <- getJSDivergence(dist_a, dist_b)
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

#' Get exact or approximate nearest neighbor distance distribution
#'
#' @param dat A \code{data.table} corresponding to repertoire annotations
#' @param column the column name of \code{dat} containing the strings on which
#'   the distribution should be computed
#' @param approximate If TRUE, approximate distribution by subsampling.
#'   Note: this inhibits interpretability since any subsample will have a 
#'   different definition of a "nearest neighbor"
#' @inheritParams getNearestNeighborDistances
#' @return vector of integer-value distances
getNearestNeighborDistribution <- function(dat,
                                           column="sequence", 
                                           k=1,
                                           approximate=TRUE,
                                           ...
                                           ) {
    sequence_list <- dat[[column]]
    if(approximate) {
        if(k == 1) {
            distribution <- getApproximateNearestNeighborDistribution(
                dat=dat,
                column=column,
                ...
            )
        } else {
            stop("k must be 1 to get approximate nearest neighbor distribution")
        }
    } else {
        distribution <- sequence_list %>%
            getNearestNeighborDistances(k=k)
    }
    return(distribution)
}

#' Plot the nearest neighbor distribution of one or more datasets
#'
#' @inheritParams plotDistribution
plotNearestNeighborDistribution <- function(dat_list,
                                            plot_type,
                                            names=NULL,
                                            show_legend=TRUE,
                                            ...
                                           ) {
    p <- plotDistribution(dat_list=dat_list,
                          summary_function=getNearestNeighborDistribution,
                          plot_type=plot_type,
                          x_label="Nearest neighbor distance",
                          names=names,
                          show_legend=show_legend,
                          ...
                         )
    return(p)
}

#' Compare kth nearest neighbor distance distributions of two lists of
#'   sequences
#' \code{compareNNDistanceDistribution} computes the JS divergence of
#'   the kth nearest neighbor distance distribution of two lists of
#'   DNA sequences
#' @param dat_a,dat_b A \code{data.table} corresponding to repertoire annotations
#' @param column the column name of \code{dat} containing the strings on which
#'   the distribution should be computed
#' @param k The separation depth for the nearest neighbor distances.
#'   k = 1 corresponds to the nearest neighbor, k = 2 corresponds to
#'   the second-nearest neighbor, etc.
#' @return Estimated JS divergence of the distributions inferred from list_a
#'   and list_b
compareNNDistanceDistributions <- function(dat_a, 
                                           dat_b, 
                                           column="sequence",
                                           k=1,
                                           approximate=TRUE,
                                           ...
                                          ) {
    dist_a <- getNearestNeighborDistribution(dat=dat_a,
                                             column=column,
                                             k=k,
                                             approximate=approximate,
                                             ...
                                            )
    dist_b <- getNearestNeighborDistribution(dat=dat_b,
                                             column=column,
                                             k=k,
                                             approximate=approximate,
                                             ...
                                            )
    divergence <- getJSDivergence(dist_a, dist_b)
    return(divergence)
}

#' Get the GC content of each sequence in a vector of strings
#' 
#' @param raw_sequences Vector of sequence strings
#' @return Vector of GC content values
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
#' @param dat A \code{data.table} corresponding to repertoire annotations
#' @param column the column name of \code{dat} containing the strings on which
#'   the distribution should be computed
#' @param approximate If TRUE, approximate distribution by subsampling.
#' @return A vector of GC content values
getGCContentDistribution <- function(dat,
                                     column="sequence",
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

#' Plot the GC content distribution of one or more datasets
#'
#' @inheritParams plotDistribution
plotGCContentDistribution <- function(dat_list,
                                      plot_type,
                                      names=NULL,
                                      show_legend=TRUE
                                     ) {
    p <- plotDistribution(dat_list,
                          getGCContentDistribution,
                          plot_type=plot_type,
                          x_label="GC content",
                          show_legend=show_legend,
                          names=names
                         )
    return(p)
}

#' Compare the GC distributions of two lists of DNA sequences
#'
#' @param dat_a,dat_b A \code{data.table} corresponding to repertoire annotations
#' @param column the column name of \code{dat} containing the strings on which
#'   the distribution should be computed
#' @return JS divergence of the GC content distributions inferred from list_a
#'   and list_b
compareGCContentDistributions <- function(dat_a, 
                                          dat_b, 
                                          column="sequence"
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
getMotifCount <- function(motif, 
                          dna_sequences
                         ) {
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
getSpotCount <- function(dna_sequences, 
                         spots
                        ) {
    spot_counts <- spots %>% 
        # Get a length(dna_sequences) x length(spots) matrix of counts
        sapply(getMotifCount, dna_sequences=dna_sequences)
        # Sum over each spot count for each sequence
    if(dim(spot_counts) %>% is.null) {
        total_count <- spot_counts %>% sum
    } else {
        total_count <- spot_counts %>% apply(1, sum)
    }
    return(total_count)
}

#' Get the number of occurrences of AID hotspots in a column of a dataset
#' 
#' @param dat A \code{data.table} corresponding to repertoire annotations
#' @param column the column name of \code{dat} containing the strings on which
#'   the distribution should be computed
#' @param hotspots Vector of hotspots of interest
#' @return The number of AID hotspot occurrences in \code{dna_sequences}
getHotspotCount <- function(dat,
                            column="sequence",
                            hotspots=c("WRC", "WA")
                           ) {
    return(getSpotCount(spots=hotspots,
                        dna_sequences=dat[[column]]
                       )
    )
}

#' Get full or approximate distribution of coldspot counts for a column of
#'   a given dataset
#'
#' @inheritParams getHotspotCount
#' @param approximate If TRUE, approximate distribution by subsampling.
getHotspotCountDistribution <- function(dat,
                                        column="sequence",
                                        hotspots=c("WRC", "WA"),
                                        approximate=FALSE,
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

#' Plot the hotspot count distribution of one or more datasets
#'
#' @inheritParams plotDistribution
plotHotspotCountDistribution <- function(dat_list,
                                         plot_type,
                                         names=NULL,
                                         show_legend=TRUE
                                        ) {
    p <- plotDistribution(dat_list,
                          getHotspotCountDistribution,
                          plot_type=plot_type,
                          x_label="Hotspot count",
                          show_legend=show_legend,
                          names=names
                         )
    return(p)
}

#' Get the number of occurrences of AID coldspots in a column of a dataset
#' 
#' @param dat A \code{data.table} corresponding to repertoire annotations
#' @param column the column name of \code{dat} containing the strings on which
#'   the distribution should be computed
#' @param coldspots Vector of coldspots of interest
#' @return The number of AID coldhotspot occurrences in \code{dna_sequences}
getColdspotCount <- function(dat,
                             column="sequence", 
                             coldspots="SYC"
                            ) {
    return(getSpotCount(dna_sequences=dat[[column]], 
                        spots=coldspots))
}

#' Get full or approximate distribution of coldspot counts for a column of
#'   a given dataset
#'
#' @inheritParams getColdspotCount
#' @param approximate If TRUE, approximate distribution by subsampling.
getColdspotCountDistribution <- function(dat,
                                         column="sequence",
                                         coldspots="SYC",
                                         approximate=FALSE,
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

#' Plot the coldspot count distribution of one or more datasets
#'
#' @inheritParams plotDistribution
plotColdspotCountDistribution <- function(dat_list,
                                          plot_type,
                                          names=NULL,
                                          show_legend=TRUE
                                         ) {
    p <- plotDistribution(dat_list,
                          getColdspotCountDistribution,
                          plot_type=plot_type,
                          x_label="Coldspot count",
                          show_legend=show_legend,
                          names=names,
                          binwidth=1
                         )
    return(p)
}

#' Compare hotspot count distributions of two sets of mature BCR sequences
#'
#' @inheritParams compareCounts
#' @param dat_a,dat_b A \code{data.table} corresponding to repertoire annotations
#' @param column the column name of \code{dat} containing the strings on which
#'   the distribution should be computed
#' @return The JS divergence of the hotspot count distributions inferred from
#'   \code{dat_a$sequence} and \code{dat_b$sequence}, respectively
compareHotspotCountDistributions <- function(dat_a, 
                                             dat_b,
                                             column="sequence",
                                             ...
                                            ) {
    dist_a <- dat_a %>% getHotspotCountDistribution(column=column,
                                                    ...
                                                   )
    dist_b <- dat_b %>% getHotspotCountDistribution(column=column,
                                                    ...
                                                   )
    divergence <- getJSDivergence(dist_a, dist_b)
    return(divergence)
}

#' Compare coldspot count distributions of two sets of mature BCR sequences
#'
#' @inheritParams compareCounts
#' @param dat_a,dat_b A \code{data.table} corresponding to repertoire annotations
#' @param column the column name of \code{dat} containing the strings on which
#'   the distribution should be computed
#' @return The JS divergence of the coldspot count distributions inferred from
#'   \code{dat_a$sequence} and \code{dat_b$sequence}, respectively
compareColdspotCountDistributions <- function(dat_a, 
                                              dat_b,
                                              column="sequence",
                                              ...
                                             ) {
    dist_a <- dat_a %>% getColdspotCountDistribution(column=column,
                                                    ...
                                                   )
    dist_b <- dat_b %>% getColdspotCountDistribution(column=column,
                                                    ...
                                                   )
    divergence <- getJSDivergence(dist_a, dist_b)
    return(divergence)
}

#' Get the exact distribution of Levenshtein distances from the inferred 
#'   naive sequences to the observed mature ones, sequence by sequence
#'
#' @inheritParams getDistanceFromGermlineToSequenceDistribution
#' @return Vector of Levenshtein distances from naive to mature
getDistancesFromGermlineToSequence <- function(dat,
                                          v_gene_only=FALSE,
                                          sequence_column=ifelse(
                                              v_gene_only,
                                              "v_qr_seqs",
                                              "sequence_alignment"
                                          ),
                                          germline_column=ifelse(
                                              v_gene_only,
                                              "v_gl_seq",
                                              "germline_alignment"
                                          )
                                         ) {
    distances <- dat[[sequence_column]] %>% 
        mapply(FUN=stringdist::stringdist, 
               b=dat[[germline_column]], 
               method="lv") %>% 
        sort %>% 
        unname
    return(distances)
}

#' Get the exact or approximate distribution of Levenshtein distances from 
#'   the inferred naive sequences to the observed mature ones
#'
#' @param dat A \code{data.table} corresponding to repertoire annotations
#' @param approximate If TRUE, approximate distribution by subsampling.
#' @param v_gene_only If TRUE, restrict sequences to the V gene only
#' @return Vector of Levenshtein distances from naive to mature
getDistanceFromGermlineToSequenceDistribution <- function(dat,
                                                          approximate=FALSE,
                                                          v_gene_only=FALSE,
                                                          ...
                                                         ) {
    if(approximate) {
        distribution <- dat %>%
            getApproximateDistribution(
                summary_function=getDistancesFromGermlineToSequence,
                divergence_function=getJSDivergence,
                v_gene_only=v_gene_only,
                ...
            )
    } else {
        distribution <- dat %>%
            getDistancesFromGermlineToSequence(v_gene_only=v_gene_only)
    }
}

#' Plot the distance from naive to mature distribution of one or more datasets
#'
#' @inheritParams plotDistribution
plotDistanceFromGermlineToSequenceDistribution <- function(dat_list,
                                                      plot_type,
                                                      names=NULL,
                                                      show_legend=TRUE
                                                     ) {
    p <- plotDistribution(dat_list,
                          getDistanceFromGermlineToSequenceDistribution,
                          plot_type=plot_type,
                          x_label="Dist. from germ. to seq.",
                          show_legend=show_legend,
                          names=names
                         )
    return(p)
}

#' Compare Levenshtein distance distributions from naive sequences to 
#  their corresponding mature ones for two repertoires
#'
#' @inheritParams getDistanceFromGermlineToSequenceDistribution
#' @param dat_a,dat_b A \code{data.table} corresponding to repertoire annotations
compareDistanceFromGermlineToSequenceDistributions <- function(
    dat_a, 
    dat_b, 
    approximate=FALSE,
    v_gene_only=FALSE,
    ...
) {
    dist_a <- dat_a %>% 
        getDistanceFromGermlineToSequenceDistribution(approximate=approximate,
                                                      v_gene_only=v_gene_only,
                                                      ...
                                                     )
    dist_b <- dat_b %>% 
        getDistanceFromGermlineToSequenceDistribution(approximate=approximate,
                                                      v_gene_only=v_gene_only,
                                                      ...
                                                     )
    divergence <- getJSDivergence(dist_a, dist_b)
    return(divergence)
}

#' Get the distribution of inferred length of CDR3 regions of each sequence.
#'   Requires either the \code{junction_length} or \code{junction} column for
#'   \code{by_amino_acid=FALSE}, or either the \code{junction} or 
#'   \code{junction_aa} column for \code{by_amnio_acid=TRUE}.
#'
#' @param dat A \code{data.table} corresponding to repertoire annotations
#' @param by_amino_acid If TRUE, the length is computed in terms of amino 
#'   acids; otherwise, the length is computed in terms of nucleotides.
#' @return Vector of CDR3 lengths (in nt units)
getCDR3LengthDistribution <- function(dat,
                                      by_amino_acid=FALSE
                                     ) {
    if(by_amino_acid) {
        if("junction_aa" %in% names(dat)) {
            CDR3_lengths <- dat$junction_aa %>%
                sapply(nchar) %>%
                unname
        } else {
            if("junction" %in% names(dat)) {
                CDR3_lengths <- dat$junction %>%
                    convertNucleobasesToAminoAcids %>%
                    sapply(nchar) %>%
                    unname
            } else {
                stop("Neither junction or junction_aa column are present in dat.")
            }
        }
    } else {
        if("junction_length" %in% names(dat)) {
            CDR3_lengths <- dat$junction_length %>% 
                na.omit
        } else {
            if("junction" %in% names(dat)) {
                CDR3_lengths <- dat$junction %>%
                    sapply(nchar) %>%
                    unname
            } else {
                stop("junction column not present in dat.")
            }
        }
    }
    return(CDR3_lengths)
}

#' Plot the CDR3 length distribution of one or more datasets
#'
#' @inheritParams plotDistribution
plotCDR3LengthDistribution <- function(dat_list,
                                       plot_type,
                                       names=NULL,
                                       show_legend=TRUE,
                                       ...
                                      ) {
    p <- plotDistribution(dat_list,
                          getCDR3LengthDistribution,
                          plot_type=plot_type,
                          x_label="CDR3 length",
                          show_legend=show_legend,
                          names=names,
                          ...
                         )
    return(p)
}

#' Compare the distribution of CDR3 lengths for two datasets
#'
#' @param dat_a,dat_b A \code{data.table} corresponding to repertoire annotations
#' @inheritParams getCDR3LengthDistribution
#' @return The JS divergence of the two CDR3 length distributions
compareCDR3LengthDistributions <- function(dat_a, 
                                           dat_b,
                                           by_amino_acid=FALSE
                                          ) {
    a_lengths <- dat_a %>% 
        getCDR3LengthDistribution(by_amino_acid=by_amino_acid)
    b_lengths <- dat_b %>% 
        getCDR3LengthDistribution(by_amino_acid=by_amino_acid)
    divergence <- getJSDivergence(a_lengths, b_lengths)
    return(divergence)
}

#' Get usage table from full list of categorical variables, such as
#'   genes or amino acids, setting zero to unused ones
#' 
#' @param factor_list List of factors
#' @param full_factor_list Full list of reference factors under consideration
#' @param standardize If TRUE, return relative frequencies rather than counts
#' @return A table of usage counts of \code{factor_list}
getUsageTableFromFullList <- function(factor_list, 
                                      full_factor_list,
                                      standardize=TRUE
                                     ) {
    usage_table <- factor_list %>% 
        table
    usage_table[full_factor_list %>% 
                setdiff(factor_list) %>% 
                unlist] <- 0
    if(standardize) {
        usage_table <- usage_table/sum(usage_table)
    }

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
compareGeneUsage <- function(gene_list_a, 
                             gene_list_b, 
                             collapse_alleles,
                             standardize
                            ) {
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

tabulateGenes <- function(dat,
                          gene_type,
                          collapse_alleles,
                          standardize
                         ) {
    gene_table <- dat[[gene_type]] %>%
        (function(x) { 
             if(collapse_alleles) { 
                 x %>% collapseAlleles 
             } else {
                 x
             }
         }) %>%
        table
    if(standardize) {
        gene_table <- gene_table/sum(gene_table)
    }

    return(gene_table)
}

#' Compare germline V, D, or J gene usage for two repertoires
#' 
#' @param dat_a,dat_b A \code{data.table} corresponding to repertoire annotations
#' @param gene_type String for gene type, taken as a column name of the 
#'   annotated datasets. Must be "v_call", "d_call", or "j_call"
#' @inheritParams compareGeneUsage
#' @return l1 divergence of gene counts between the two repertoires
compareGermlineGeneDistributions <- function(dat_a, 
                                             dat_b, 
                                             gene_type,
                                             collapse_alleles,
                                             standardize
                                            ) {
    if(gene_type %>% missing) {
        stop("gene_type needs to be supplied.")
    }

    gene_table_a <- tabulateGenes(dat=dat_a,
                                  gene_type=gene_type,
                                  collapse_alleles=collapse_alleles,
                                  standardize=standardize
                                 )
    gene_table_b <- tabulateGenes(dat=dat_b,
                                  gene_type=gene_type,
                                  collapse_alleles=collapse_alleles,
                                  standardize=standardize
                                 )
    divergence <- compareCategoricalDistributions(gene_table_a,
                                                  gene_table_b)
    return(divergence)
}

#' Compare V gene distribution between two repertoires
#'
#' @inheritParams compareGermlineGeneDistributions
#' @return l1 divergence of V gene counts between the two repertoires
compareVGeneDistributions <- function(dat_a, 
                                      dat_b,
                                      collapse_alleles=TRUE,
                                      standardize=TRUE
                                     ) {
    return(compareGermlineGeneDistributions(dat_a, 
                                            dat_b, 
                                            gene_type="v_call",
                                            collapse_alleles=collapse_alleles,
                                            standardize=standardize
                                           )
    )
}

#' Compare D gene distribution between two repertoires
#'
#' @inheritParams compareGermlineGeneDistributions
#' @return l1 divergence of D gene counts between the two repertoires
compareDGeneDistributions <- function(dat_a, 
                                      dat_b,
                                      collapse_alleles=TRUE,
                                      standardize=TRUE
                                     ) {
    return(compareGermlineGeneDistributions(dat_a, 
                                            dat_b, 
                                            gene_type="d_call",
                                            collapse_alleles=collapse_alleles,
                                            standardize=standardize
                                           )
    )
}

#' Compare J gene distribution between two repertoires
#'
#' @inheritParams compareGermlineGeneDistributions
#' @return l1 divergence of D gene counts between the two repertoires
compareJGeneDistributions <- function(dat_a, 
                                      dat_b,
                                      collapse_alleles=TRUE,
                                      standardize=TRUE
                                     ) {
    return(compareGermlineGeneDistributions(dat_a, 
                                            dat_b, 
                                            gene_type="j_call",
                                            collapse_alleles=collapse_alleles,
                                            standardize=standardize
                                           )
    )
}

#' Get table of combined V, D, and J usage frequencies
#'
#' @param dat A \code{data.table} corresponding to repertoire annotations
#' @inheritParams compareGermlineGeneDistributions
#' @param by_frequency If true, scale table entries as relative frequencies.
#'   Otherwise, entries will be counts.
#' @return A table of joint gene IDs and usage counts
getJointGeneTable <- function(dat, 
                              gene_calls=c("v_call", "j_call"),
                              collapse_alleles=TRUE,
                              by_frequency=TRUE
                             ) {
    tmp_dat <- dat
    if(collapse_alleles) {
        tmp_dat[, gene_calls] <- gene_calls %>%
            sapply(function(x) { 
                       tmp_dat[[x]] <- tmp_dat[[x]] %>% 
                           collapseAlleles
                   }
            )

    }

    # Concatenate each gene call into one string to define a joint call
    gene_table <- do.call(paste0, 
                                 tmp_dat[, gene_calls]
                         ) %>%
        table
    
    if(by_frequency) {
        gene_table <- gene_table/sum(gene_table)
    }

    return(gene_table)
}

#' Compare joint gene distributions of two datasets
#'
#' @inheritParams compareGermlineGeneDistributions
compareJointGeneDistributions <- function(dat_a,
                                          dat_b,
                                          gene_calls,
                                          collapse_alleles,
                                          by_frequency
                                         ) {
    table_a <- getJointGeneTable(dat_a, 
                                 gene_calls=gene_calls,
                                 collapse_alleles=collapse_alleles,
                                 by_frequency=by_frequency
                                )
    table_b <- getJointGeneTable(dat_b, 
                                 gene_calls=gene_calls,
                                 collapse_alleles=collapse_alleles,
                                 by_frequency=by_frequency
                                )
    divergence <- compareCategoricalDistributions(table_a,
                                                  table_b)
    return(divergence)
}

compareVJDistributions <- function(dat_a,
                                   dat_b,
                                   collapse_alleles=TRUE,
                                   by_frequency=TRUE
                                  ) {
    divergence <- compareJointGeneDistributions(
         dat_a=dat_a,
         dat_b=dat_b,
         gene_calls=c("v_call", "j_call"),
         collapse_alleles=collapse_alleles,
         by_frequency=by_frequency
        )
    return(divergence)
}

#' Compare joint V, D, and J gene usage between two annotated repertoires
#'
#' @inheritParams compareGermlineGeneDistributions
#' @return l1 divergence of counts of VDJ triples between the two repertoires
compareVDJDistributions <- function(dat_a, 
                                    dat_b, 
                                    collapse_alleles=TRUE,
                                    by_frequency=TRUE
                                   ) {
    divergence <- compareJointGeneDistributions(
         dat_a=dat_a,
         dat_b=dat_b,
         gene_calls=c("v_call", "d_call", "j_call"),
         collapse_alleles=collapse_alleles,
         by_frequency=by_frequency
        )
    return(divergence)
}

#' Get the kidera factors of a single sequence.
#'   See the \code{Peptides::kideraFactors} documentation for more information.
#'
#' @param aa_sequence An amino acid sequence string
#' @return A named vector of the 10 Kidera factors for \code{aa_sequence}
getKideraFactorsBySequence <- function(aa_sequence) {
    kidera_factors <- aa_sequence %>% 
        Peptides::kideraFactors() %>% 
        unlist
    return(kidera_factors)
}

#' Get distributions for each of the ten Kidera factors for a dataset 
#'
#' @param dat A \code{data.table} corresponding to repertoire annotations
#' @param column the column name of \code{dat} containing the strings on which
#'   the distribution should be computed
#' @param as_list If TRUE, return a list of the ten distributions. Otherwise,
#'   return an (n)x(10) data.table, where n is the number of valid amino acid 
#'   sequences in \code{column}
#' @return List or data.table of Kidera factor distributions
getKideraFactorDistributions <- function(dat,
                                         column="junction_aa",
                                         as_list=TRUE
                                        ) {
    factor_distributions <- dat[[column]] %>% 
        filterAminoAcidSequences %>%
        lapply(getKideraFactorsBySequence) %>%
        do.call("rbind", .)

    if(as_list) {
        factor_distributions <- 1:10 %>%
            lapply(function(x) { factor_distributions[, x] }) %>%
            setNames(paste0("KF", 1:10))
    }

    return(factor_distributions) 
}

#' Compare the distributions of Kidera factors for two datasets
#'
#' @param dat_a,dat_b A \code{data.table} corresponding to repertoire 
#'   annotations
#' @return Vector of divergences of each Kidera factor distribution
compareKideraFactorDistributions <- function(dat_a,
                                             dat_b,
                                             column="junction_aa"
                                            ) {
    dists_a <- dat_a %>% getKideraFactorDistributions(column=column)
    dists_b <- dat_b %>% getKideraFactorDistributions(column=column)
    divergences <- mapply(getContinuousJSDivergence,
                          dists_a,
                          dists_b
                         ) %>% 
        as.list %>%
        setNames(paste0("KideraFactor", 1:10, "Divergence"))

    return(divergences)
}

#' Get the distribution of a given Atchley factor for a vector of amino acid 
#'   sequences
#'
#' @param aa_sequences A vector of amino acid sequence string
#' @param factor_number The Atchley factor to be applied over 
#'  \code{aa_sequences}. Must be 1, 2, 3, 4, or 5.
#' @return A vector of Atchley factors of each amino acid in each sequence
getAtchleyFactorDistribution <- function(aa_sequences, 
                                         factor_number
                                        ) {
    factor_list <- aa_sequences %>% 
        HDMD::FactorTransform(Factor=factor_number) %>%
        unlist %>%
        unname
    return(factor_list)
}

#' Get a list of the distributions of each of the five Atchley factors for a 
#'   dataset.
#'
#' @param dat A \code{data.table} corresponding to repertoire annotations
#' @param column the column name of \code{dat} containing the strings on which
#'   the distribution should be computed
#' @return List of Atchley factor distributions
getAtchleyFactorDistributions <- function(dat,
                                          column="junction_aa"
                                         ) {
    collapsed_sequences <- dat[[column]] %>% 
        filterAminoAcidSequences
    factor_numbers <- 1:5
    atchley_factors <- factor_numbers %>% 
        lapply(function(x) {
            collapsed_sequences %>% 
                getAtchleyFactorDistribution(factor_number=x) %>%
                unlist
        })

    return(atchley_factors)
}

#' Compare the distributions of Atchley factors for two datasets
#'
#' @param dat_a,dat_b A \code{data.table} corresponding to repertoire 
#'   annotations
#' @return Vector of divergences of each Atchley factor distribution
compareAtchleyFactorDistributions <- function(dat_a, 
                                              dat_b,
                                              column="junction_aa"
                                             ) {
    dists_a <- dat_a %>% getAtchleyFactorDistributions(column=column)
    dists_b <- dat_b %>% getAtchleyFactorDistributions(column=column)
    divergences <- mapply(getContinuousJSDivergence,
                          dists_a,
                          dists_b
                         ) %>%
        as.list %>%
        setNames(paste0("AtchleyFactor", 1:5, "Divergence"))
    return(divergences)
}

#' Get the aliphatic index of a DNA sequence
#'
#' @param aa_sequence An amino acid sequence string
#' @return The aliphatic index of \code{dna_sequence}, if applicable
getAliphaticIndex <- function(aa_sequence) {
    aliphatic_index <- aa_sequence %>% 
               Peptides::aIndex() %>%
               ifelse(. < 125, ., NA) # Don't bother with outliers

    return(aliphatic_index)
}

#' Get the distribution of aliphatic indices of a list of sequences
#'
#' @param dat A \code{data.table} corresponding to repertoire annotations
#' @param column the column name of \code{dat} containing the strings on which
#'   the distribution should be computed
#' @return Vector of aliphatic indices
getAliphaticIndexDistribution <- function(dat,
                                          column="junction_aa"
                                         ) {
    a_indices <- dat[[column]] %>%
        filterAminoAcidSequences %>%
        sapply(function(x) {
                   ifelse(!is.na(x),
                          getAliphaticIndex(x),
                          NA)
        }
    )

    return(a_indices)
}

#' Plot the aliphatic index distribution of one or more datasets
#'
#' @inheritParams plotDistribution
plotAliphaticIndexDistribution <- function(dat_list,
                                           plot_type,
                                           names=NULL,
                                           show_legend=TRUE
                                          ) {
    p <- plotDistribution(dat_list,
                          getAliphaticIndexDistribution,
                          plot_type=plot_type,
                          x_label="Aliphatic index",
                          show_legend=show_legend,
                          names=names
                         )
    return(p)
}

#' Compare the distributions of aliphatic indices of two datasets
#'
#' @param dat_a,dat_b A \code{data.table} corresponding to repertoire annotations
#' @return The JS divergence of the two distributions
compareAliphaticIndexDistributions <- function(dat_a, 
                                               dat_b,
                                               ...
                                              ) {
    dist_a <- dat_a %>% getAliphaticIndexDistribution(...)
    dist_b <- dat_b %>% getAliphaticIndexDistribution(...)
    divergence <- getContinuousJSDivergence(dist_a, dist_b)
    return(divergence)
}

#' Get the distribution of GRAVY values from a list or vector of sequences
#'
#' @param dat A \code{data.table} corresponding to repertoire annotations
#' @param column the column name of \code{dat} containing the strings on which
#'   the distribution should be computed
#' @return Vector of GRAVY values for \code{sequence_list}
getGRAVYDistribution <- function(dat,
                                 column="junction_aa"
                                ) {
    dist <- dat[[column]] %>%
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

#' Plot the GRAVY index distribution of one or more datasets
#'
#' @inheritParams plotDistribution
plotGRAVYDistribution <- function(dat_list,
                                  plot_type,
                                  names=NULL,
                                  show_legend=TRUE
                                 ) { 
    p <- plotDistribution(dat_list,
                          getGRAVYDistribution,
                          plot_type=plot_type,
                          x_label="GRAVY index",
                          show_legend=show_legend,
                          names=names
                         )
    return(p)
}


#' Compare the GRAVY distributions of two datasets
#'
#' @param dat_a,dat_b A \code{data.table} corresponding to repertoire annotations
#' @return The JS divergence of GRAVY distributions
compareGRAVYDistributions <- function(dat_a, dat_b) {
    dist_a <- dat_a %>% 
        getGRAVYDistribution
    dist_b <- dat_b %>%
        getGRAVYDistribution
    divergence <- getJSDivergence(dist_a, dist_b, continuous=TRUE)
    return(divergence)
}

#' Get a particular distribution of amino acid properties using the
#'   \code{alakazam::aminoAcidProperties} function
#'
#' @param dat A \code{data.table} corresponding to repertoire annotations
#' @param column the column name of \code{dat} containing the strings on which
#'   the distribution should be computed
#' @param suffix A property-specific suffix used to extract the desired
#'   property from the \code{data.table} returned by \code{alakazam}. 
#'   This is usually of the form "AA_PROPERTY" (e.g. "AA_POLARITY"),
#'   but not always (e.g. "AA_BASIC" for basicity)
#' @return distribution of the given amino acid property
getAminoAcidProperties <- function(dat,
                                   column="junction_aa",
                                   suffix
                                  ) {
    properties <- dat %>%
        alakazam::aminoAcidProperties(seq=column,
                                      nt=FALSE)
    properties_column <- paste(column,
                               suffix,
                               sep="_"
                              )

    return(properties[[properties_column]])
}

#' Get a vector of polarity values for a set of an amino acid sequences
#'
#' @inheritParams getAminoAcidProperties
getPolarityDistribution <- function(dat,
                                    column="junction_aa"
                                   ) {
    polarities <- dat %>%
        getAminoAcidProperties(column=column,
                               suffix="AA_POLARITY"
                              )
    return(polarities)
}

#' Compare the polarity distributions of two datasets
#'
#' @param dat_a,dat_b A \code{data.table} corresponding to repertoire annotations
comparePolarityDistributions <- function(dat_a,
                                         dat_b
                                        ) {
    dist_a <- dat_a %>% getPolarityDistribution
    dist_b <- dat_b %>% getPolarityDistribution
    divergence <- getContinuousJSDivergence(dist_a, dist_b)
    return(divergence)
}

#' Plot the polarity distribution of one or more datasets
#'
#' @inheritParams plotDistribution
plotPolarityDistribution <- function(dat_list,
                                     plot_type,
                                     names=NULL,
                                     show_legend=TRUE,
                                     ...
                                    ) {
    p <- plotDistribution(dat_list,
                          getPolarityDistribution,
                          plot_type=plot_type,
                          x_label="Polarity",
                          show_legend=show_legend,
                          names=names
                         )

    return(p)
}

#' Get a vector of charge values for a set of an amino acid sequences
#'
#' @inheritParams getAminoAcidProperties
getChargeDistribution <- function(dat,
                                  column="junction_aa"
                                 ) {
    charges <- dat %>%
        getAminoAcidProperties(column=column,
                               suffix="AA_CHARGE"
                              )
    return(charges)
}

#' Compare the charge distributions of two datasets
#'
#' @param dat_a,dat_b A \code{data.table} corresponding to repertoire annotations
compareChargeDistributions <- function(dat_a,
                                       dat_b
                                      ) {
    dist_a <- dat_a %>% getChargeDistribution
    dist_b <- dat_b %>% getChargeDistribution
    divergence <- getContinuousJSDivergence(dist_a, dist_b)
    return(divergence)
}

#' Plot the charge distribution of one or more datasets
#'
#' @inheritParams plotDistribution
plotChargeDistribution <- function(dat_list,
                                     plot_type,
                                     names=NULL,
                                     show_legend=TRUE,
                                     ...
                                    ) {
    p <- plotDistribution(dat_list,
                          getChargeDistribution,
                          plot_type=plot_type,
                          x_label="Charge",
                          show_legend=show_legend,
                          names=names
                         )

    return(p)
}

#' Get a vector of basicity values for a set of an amino acid sequences
#'
#' @inheritParams getAminoAcidProperties
getBasicityDistribution <- function(dat,
                                    column="junction_aa"
                                   ) {
    basicities <- dat %>%
        getAminoAcidProperties(column=column,
                               suffix="AA_BASIC"
                              )

    return(basicities)
}

#' Compare the basicity distributions of two datasets
#'
#' @param dat_a,dat_b A \code{data.table} corresponding to repertoire annotations
compareBasicityDistributions <- function(dat_a,
                                         dat_b
                                        ) {
    dist_a <- dat_a %>% getBasicityDistribution
    dist_b <- dat_b %>% getBasicityDistribution
    divergence <- getContinuousJSDivergence(dist_a, dist_b)
    return(divergence)
}

#' Plot the basicity distribution of one or more datasets
#'
#' @inheritParams plotDistribution
plotBasicityDistribution <- function(dat_list,
                                     plot_type,
                                     names=NULL,
                                     show_legend=TRUE,
                                     ...
                                    ) {
    p <- plotDistribution(dat_list,
                          getBasicityDistribution,
                          plot_type=plot_type,
                          x_label="Basicity",
                          show_legend=show_legend,
                          names=names
                         )

    return(p)
}

#' Get a vector of acidity values for a set of an amino acid sequences
#'
#' @inheritParams getAminoAcidProperties
getAcidityDistribution <- function(dat,
                                   column="junction_aa"
                                  ) {
    acidities <- dat %>%
        getAminoAcidProperties(column=column,
                               suffix="AA_ACIDIC"
                              )
    return(acidities)
}

#' Compare the acidity distributions of two datasets
#'
#' @param dat_a,dat_b A \code{data.table} corresponding to repertoire annotations
compareAcidityDistributions <- function(dat_a,
                                         dat_b
                                        ) {
    dist_a <- dat_a %>% getAcidityDistribution
    dist_b <- dat_b %>% getAcidityDistribution
    divergence <- getContinuousJSDivergence(dist_a, dist_b)
    return(divergence)
}

#' Plot the acidity distribution of one or more datasets
#'
#' @inheritParams plotDistribution
plotAcidityDistribution <- function(dat_list,
                                     plot_type,
                                     names=NULL,
                                     show_legend=TRUE,
                                     ...
                                    ) {
    p <- plotDistribution(dat_list,
                          getAcidityDistribution,
                          plot_type=plot_type,
                          x_label="Acidity",
                          show_legend=show_legend,
                          names=names
                         )

    return(p)
}

#' Get a vector of aromaticity values for a set of an amino acid sequences
#'
#' @inheritParams getAminoAcidProperties
getAromaticityDistribution <- function(dat,
                                       column="junction_aa"
                                      ) {
    aromaticities <- dat %>%
        getAminoAcidProperties(column=column,
                               suffix="AA_AROMATIC"
                              )

    return(aromaticities)
}

#' Compare the aromaticity distributions of two datasets
#'
#' @param dat_a,dat_b A \code{data.table} corresponding to repertoire annotations
compareAromaticityDistributions <- function(dat_a,
                                         dat_b
                                        ) {
    dist_a <- dat_a %>% getAromaticityDistribution
    dist_b <- dat_b %>% getAromaticityDistribution
    divergence <- getContinuousJSDivergence(dist_a, dist_b)
    return(divergence)
}

#' Plot the aromaticity distribution of one or more datasets
#'
#' @inheritParams plotDistribution
plotAromaticityDistribution <- function(dat_list,
                                     plot_type,
                                     names=NULL,
                                     show_legend=TRUE,
                                     ...
                                    ) {
    p <- plotDistribution(dat_list,
                          getAromaticityDistribution,
                          plot_type=plot_type,
                          x_label="Aromaticity",
                          show_legend=show_legend,
                          names=names
                         )

    return(p)
}

#' Get a vector of bulkiness values for a set of an amino acid sequences
#'
#' @inheritParams getAminoAcidProperties
getBulkinessDistribution <- function(dat,
                                     column="junction_aa"
                                    ) {
    bulkinesses <- dat %>%
        getAminoAcidProperties(column=column,
                               suffix="AA_BULK"
                              )
    return(bulkinesses)
}

#' Compare the bulkiness distributions of two datasets
#'
#' @param dat_a,dat_b A \code{data.table} corresponding to repertoire annotations
compareBulkinessDistributions <- function(dat_a,
                                         dat_b
                                        ) {
    dist_a <- dat_a %>% getBulkinessDistribution
    dist_b <- dat_b %>% getBulkinessDistribution
    divergence <- getContinuousJSDivergence(dist_a, dist_b)
    return(divergence)
}

#' Plot the bulkiness distribution of one or more datasets
#'
#' @inheritParams plotDistribution
plotBulkinessDistribution <- function(dat_list,
                                     plot_type,
                                     names=NULL,
                                     show_legend=TRUE,
                                     ...
                                    ) {
    p <- plotDistribution(dat_list,
                          getBulkinessDistribution,
                          plot_type=plot_type,
                          x_label="Bulkiness",
                          show_legend=show_legend,
                          names=names
                         )

    return(p)
}

#' Extract CDR3 codon start positions from partis-returned dictionary strings
#'
#' This is a helper function to \link{getCDR3PairwiseDistanceDistribution}.
#' Note: We add one to \code{positions} because partis returns position values
#'   that are zero-based (an artifact) of Python, whereas R is one-based by
#'   default.
#' @param dictionary_list A vector of strings corresponding to a Python 
#'   dictionary, present in the column \code{codon_positions} of a data.table 
#'   returned by \code{getPartisAnnotations}.
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
#'   \code{sequence}.
#' @param dat A \code{data.table} corresponding to repertoire annotations
#' @return Vector of CDR3 strings
getCDR3PairwiseDistanceDistribution <- function(dat,
                                                by_amino_acid=TRUE,
                                                column=ifelse(by_amino_acid,
                                                              "junction_aa",
                                                              "junction"),
                                                approximate=TRUE,
                                                ...
                                               ) {
    distances <- dat %>% 
        getPairwiseDistanceDistribution(column=column,
                                        approximate=approximate,
                                        ...
                                       )
    return(distances)
}

#' Compare levenshtein distance distributions of two CDR3 repertoires
#'
#' @param dat_a,dat_b A \code{data.table} corresponding to repertoire annotations
#' @return The JS divergence of the levenshtein distance distributions of the
#'   CDR3s of the two repertoires
compareCDR3PairwiseDistanceDistributions <- function(dat_a, 
                                                     dat_b, 
                                                     by_amino_acid=TRUE,
                                                     column=ifelse(by_amino_acid,
                                                                   "junction_aa",
                                                                   "junction"),
                                                     approximate=TRUE,
                                                     ...
                                                    ) {
    distances_a <- dat_a %>% 
        getCDR3PairwiseDistanceDistribution(column=column,
                                            approximate=approximate,
                                            ...
                                           )
    distances_b <- dat_b %>% 
        getCDR3PairwiseDistanceDistribution(column=column,
                                            approximate=approximate,
                                            ...
                                           )
    divergence <- getJSDivergence(distances_a, distances_b)
    return(divergence)
}

#' Get the positional distances between mutations given a single mature sequence
#'   and corresponding naive sequence
#'
#' @param naive A string corresponding to a naive sequence
#' @param mature A string corresponding to a mature sequence
#' @return A vector of positional distances between the mutations
getPositionalDistancesBetweenMutationsBySequence <- function(naive, mature) {
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
        distances <- mut_indices_list %>% diff
    } else {
        # If there are less than two mutations, this statistic is undefined,
        # so return NA and eventually ignore
        distances <- NA
    }
    return(distances)
}

#' Get the distribution of positional distances between mutations of a given
#'   dataset.
#'   Note that this requires \code{nchar} to match for each naive/mature 
#'   sequence pair. (TODO: relax this assumption by discarding cases which
#'   don't match)
#'
#' @param dat 
#' @param naive_column Column containing the inferred naive/germline sequences
#' @param mature_column Column containing the mature SHM-experienced sequences
#' @param return Vector of positional distances
getPositionalDistanceBetweenMutationsDistribution <- function(
    dat,
    naive_column="germline_alignment",
    mature_column="sequence_alignment"
) {
    dists <- mapply(function(x, y) {
                        if(nchar(x) == nchar(y)) {
                            res <- getPositionalDistancesBetweenMutationsBySequence(x, y)
                        } else {
                            res <- NA
                        }
                    },
                    dat[[naive_column]],
                    dat[[mature_column]]
                   ) %>% 
        unname %>% 
        unlist %>% 
        subset(!is.na(.))
    return(dists)
}

#' Plot the distance between mutation distribution of one or more datasets
#'
#' @inheritParams plotDistribution
plotPositionalDistanceBetweenMutationsDistribution <- function(dat_list,
                                                     plot_type,
                                                     names=NULL,
                                                     show_legend=TRUE
                                                    ) { 
    p <- plotDistribution(dat_list,
                          getPositionalDistanceBetweenMutationsDistribution,
                          plot_type=plot_type,
                          x_label="Pos. dist. between mutations",
                          show_legend=show_legend,
                          names=names
                         )
    return(p)
}

#' Compare the positional distance between mutations distributions of two
#'   datasets
#'
#' @param dat_a,dat_b A \code{data.table} corresponding to repertoire annotations
comparePositionalDistanceBetweenMutationsDistributions <- function(dat_a, dat_b) {
    dists_a <- getPositionalDistanceBetweenMutationsDistribution(dat_a)
    dists_b <- getPositionalDistanceBetweenMutationsDistribution(dat_b)
    divergence <- getJSDivergence(dists_a, dists_b)
    return(divergence)
}

#' Get mutation rates of each gene from a \code{mutation_rates} object.
#'   This object is inferred from some XCR repertoire and should be an element 
#'   of a \code{list}, e.g. as returned by \code{getPartisAnnotations}.
#'
#' @param rate_dat A \code{mutation_rates} object containing, e.g. as obtained
#'   in \code{getPartisAnnotations}
#' return A named vector of mutation rates for each gene
getPerGeneMutationRates <- function(rate_dat) {
    rates <- rate_dat %>% 
        sapply( function(gene) { 
                    gene$overall_mut_rate 
                }
        )
    return(rates)
}

#' Compare per gene mutation rates of two datasets
#'
#' @param rate_dat_a,rate_dat_b A \code{mutation_rates} object
#' @return l1 divergence of per gene mutation rates
comparePerGeneMutationRates <- function(rate_dat_a, 
                                        rate_dat_b
                                       ) {
    rates_a <- rate_dat_a %>% 
        getPerGeneMutationRates
    rates_b <- rate_dat_b %>% 
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
        sum
    return(divergence) 
}

#' Get mutation rates by position of each gene from a \code{mutation_rates}
#'   object.
#'
#' @inheritParams getPerGeneMutationRates
getPerGenePerPositionMutationRates <- function(rate_dat) {
    rates <- rate_dat %>% 
        sapply( function(gene) {
                    gene$mut_rate_by_position 
                } 
        )
    return(rates)
}

#' Compare per gene per position mutation rates of two datasets
#'
#' @inheritParams comparePerGeneMutationRates
comparePerGenePerPositionMutationRates <- function(rate_dat_a, 
                                                   rate_dat_b
                                                  ) {
    rates_a <- rate_dat_a %>% 
        getPerGenePerPositionMutationRates
    rates_b <- rate_dat_b %>% 
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
        sum
    return(divergence) 
}

#' Get the inferred substitution model of a repertoire
#'
#' @param dat A \code{data.table} corresponding to repertoire annotations
#' @return The default shazam substitution model 
getSubstitutionModel <- function(dat) {
    sub_mat <- dat %>% 
        removeSequencesWithDifferentGermlineAndSequenceLengths %>%
        shazam::createSubstitutionMatrix(sequenceColumn="sequence_alignment",
                                         germlineColumn="germline_alignment",
                                         vCallColumn="v_call") 
    return(sub_mat)
}

#' Get the inferred mutability model of a repetoire
#'
#' This requires the inferred substitution model of \code{dat} as input
#' @param dat A \code{data.table} corresponding to repertoire annotations
#' @param substitution_model The inferred substitution model, which can be
#'   obtained via \code{getSubstitutionModel}
#' @return The default shazam mutability model
getMutabilityModel <- function(dat, 
                               substitution_model=getSubstitutionModel(dat)
                              ) {
    mut_mat <- dat %>% 
        removeSequencesWithDifferentGermlineAndSequenceLengths %>%
        shazam::createMutabilityMatrix(substitutionModel=substitution_model,
                                       sequenceColumn="sequence_alignment",
                                       germlineColumn="germline_alignment",
                                       vCallColumn="v_call")
    return(mut_mat)
}

#' Compare the mutability models of two datasets
#'
#' @param dat_a,dat_b A \code{data.table} corresponding to repertoire annotations
#' @param sub_mod_a,sub_mod_b A shazam substitution model
#' @return ell-1 divergence of mutability model matrices
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
#' @param dat_a,dat_b A \code{data.table} corresponding to repertoire annotations
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

    divergences <- list(SubstitutionModelDivergence=sub_divergence,
                        MutabilityModelDivergence=mut_divergence)
    return(divergences)
}

#' Get the requested deletion lengths of a dataset
#'
#' @param dat A \code{data.table} corresponding to repertoire annotations
#' @param column the column name of \code{dat} containing the desired deletion
#'   lengths.
#' @return Vector of deletion lengths
getDeletionLengths <- function(dat, column) {
    lengths <- dat %>% 
        dplyr::select_(column) %>% 
        unlist(use.names=FALSE)
    return(lengths)
}

#' Get the distribution of V gene 3' deletion lengths of a dataset
#'
#' @inheritParams getDeletionLengths
#' @return Distribution of V gene 3' deletion lengths
getVGene3PrimeDeletionLengthDistribution <- function(dat) {
    return(getDeletionLengths(dat, "v_3p_del"))
}

#' Plot the V gene 3' deletion length distribution of one or more datasets
#'
#' @inheritParams plotDistribution
plotVGene3PrimeDeletionLengthDistribution <- function(dat_list,
                                           plot_type,
                                           names=NULL,
                                           show_legend=TRUE
                                          ) { 
    p <- plotDistribution(dat_list,
                          getVGene3PrimeDeletionLengthDistribution,
                          plot_type=plot_type,
                          x_label="V gene 3' deletion length",
                          show_legend=show_legend,
                          names=names
                         )
    return(p)
}

#' Get the distribution of V gene 5' deletion lengths of a dataset
#'
#' @inheritParams getDeletionLengths
#' @return Distribution of V gene 5' deletion lengths
getVGene5PrimeDeletionLengthDistribution <- function(dat) {
    return(getDeletionLengths(dat, "v_5p_del"))
}

#' Plot the V gene 5' deletion length distribution of one or more datasets
#'
#' @inheritParams plotDistribution
plotVGene5PrimeDeletionLengthDistribution <- function(dat_list,
                                           plot_type,
                                           names=NULL,
                                           show_legend=TRUE
                                          ) { 
    p <- plotDistribution(dat_list,
                          getVGene5PrimeDeletionLengthDistribution,
                          plot_type=plot_type,
                          x_label="V gene 5' deletion length",
                          show_legend=show_legend,
                          names=names
                         )
    return(p)
}

#' Get the distribution of D gene 3' deletion lengths of a dataset
#'
#' @inheritParams getDeletionLengths
#' @return Distribution of D gene 3' deletion lengths
getDGene3PrimeDeletionLengthDistribution <- function(dat) {
    return(getDeletionLengths(dat, "d_3p_del"))
}

#' Plot the D gene 3' deletion length distribution of one or more datasets
#'
#' @inheritParams plotDistribution
plotDGene3PrimeDeletionLengthDistribution <- function(dat_list,
                                           plot_type,
                                           names=NULL,
                                           show_legend=TRUE
                                          ) { 
    p <- plotDistribution(dat_list,
                          getDGene3PrimeDeletionLengthDistribution,
                          plot_type=plot_type,
                          x_label="D gene 3' deletion length",
                          show_legend=show_legend,
                          names=names
                         )
    return(p)
}

#' Get the distribution of D gene 5' deletion lengths of a dataset
#'
#' @inheritParams getDeletionLengths
#' @return Distribution of D gene 5' deletion lengths
getDGene5PrimeDeletionLengthDistribution <- function(dat) {
    return(getDeletionLengths(dat, "d_5p_del"))
}

#' Plot the D gene 5' deletion length distribution of one or more datasets
#'
#' @inheritParams plotDistribution
plotDGene5PrimeDeletionLengthDistribution <- function(dat_list,
                                           plot_type,
                                           names=NULL,
                                           show_legend=TRUE
                                          ) { 
    p <- plotDistribution(dat_list,
                          getDGene5PrimeDeletionLengthDistribution,
                          plot_type=plot_type,
                          x_label="D gene 5' deletion length",
                          show_legend=show_legend,
                          names=names
                         )
    return(p)
}

#' Get the distribution of J gene 3' deletion lengths of a dataset
#'
#' @inheritParams getDeletionLengths
#' @return Distribution of J gene 3' deletion lengths
getJGene3PrimeDeletionLengthDistribution <- function(dat) {
    return(getDeletionLengths(dat, "j_3p_del"))
}

#' Plot the J gene 3' deletion length distribution of one or more datasets
#'
#' @inheritParams plotDistribution
plotJGene3PrimeDeletionLengthDistribution <- function(dat_list,
                                           plot_type,
                                           names=NULL,
                                           show_legend=TRUE
                                          ) { 
    p <- plotDistribution(dat_list,
                          getJGene3PrimeDeletionLengthDistribution,
                          plot_type=plot_type,
                          x_label="J gene 3' deletion length",
                          show_legend=show_legend,
                          names=names
                         )
    return(p)
}

#' Get the distribution of J gene 5' deletion lengths of a dataset
#'
#' @inheritParams getDeletionLengths
#' @return Distribution of J gene 5' deletion lengths
getJGene5PrimeDeletionLengthDistribution <- function(dat) {
    return(getDeletionLengths(dat, "j_5p_del"))
}

#' Plot the J gene 5' deletion length distribution of one or more datasets
#'
#' @inheritParams plotDistribution
plotJGene5PrimeDeletionLengthDistribution <- function(dat_list,
                                           plot_type,
                                           names=NULL,
                                           show_legend=TRUE
                                          ) { 
    p <- plotDistribution(dat_list,
                          getJGene5PrimeDeletionLengthDistribution,
                          plot_type=plot_type,
                          x_label="J gene 5' deletion length",
                          show_legend=show_legend,
                          names=names
                         )
    return(p)
}

#' Compare the deletion length distributions of two datasets
#'
#' @param dat_a,dat_b A \code{data.table} corresponding to repertoire annotations
#' @param gene Gene label Either "VGene", "DGene", or "JGene".
#' @param end The end of gene over which the deletion occured. Either "3Prime" 
#'   or "5Prime".
#' @return JS divergence of deletion lengths between \code{dat_a} and \code{dat_b}
compareDeletionLengths <- function(dat_a, dat_b, gene, end) {
    deletion_length_function <- paste0("get", gene, end, "DeletionLengthDistribution") %>%
        get
    dist_a <- dat_a %>% 
        deletion_length_function
    dist_b <- dat_b %>% 
        deletion_length_function
    divergence <- getJSDivergence(dist_a, dist_b)
    return(divergence)
}

#' Compare the V gene 3' deletion length distributions of two datasets
#'
#' @inheritParams compareDeletionLengths
compareVGene3PrimeDeletionLengthDistributions <- function(dat_a, dat_b) {
    return(compareDeletionLengths(dat_a, dat_b, "VGene", "3Prime"))
}

#' Compare the D gene 3' deletion length distributions of two datasets
#'
#' @inheritParams compareDeletionLengths
compareDGene3PrimeDeletionLengthDistributions <- function(dat_a, dat_b) {
    return(compareDeletionLengths(dat_a, dat_b, "DGene", "3Prime"))
}

#' Compare the D gene 5' deletion length distributions of two datasets
#'
#' @inheritParams compareDeletionLengths
compareDGene5PrimeDeletionLengthDistributions <- function(dat_a, dat_b) {
    return(compareDeletionLengths(dat_a, dat_b, "DGene", "5Prime"))
}

#' Compare the J gene 5' deletion length distributions of two datasets
#'
#' @inheritParams compareDeletionLengths
compareJGene5PrimeDeletionLengthDistributions <- function(dat_a, dat_b) {
    return(compareDeletionLengths(dat_a, dat_b, "JGene", "5Prime"))
}

#' Get the insertion lengths corresponding to a column of a given datset
#'
#' @param column the column name of \code{dat} containing the strings on which
#'   the distribution should be computed
getInsertionLengths <- function(dat, column) {
    lengths <- dat %>% 
        dplyr::select_(column) %>% 
        unlist(use.names=FALSE)
    return(lengths)
}

#' Get the distribution of VD insertion lengths of a dataset
#'
#' @inheritParams getInsertionLengths
getVDInsertionLengthDistribution <- function(dat) {
    return(getInsertionLengths(dat, "np1_length"))
}

#' Plot the VD insertion length distribution of one or more datasets
#'
#' @inheritParams plotDistribution
plotVDInsertionLengthDistribution <- function(dat_list,
                                   plot_type,
                                   names=NULL,
                                   show_legend=TRUE
                                  ) { 
    p <- plotDistribution(dat_list,
                          getVDInsertionLengthDistribution,
                          plot_type=plot_type,
                          x_label="VD insertion length",
                          show_legend=show_legend,
                          names=names
                         )
    return(p)
}


#' Get the distribution of DJ insertion lengths of a dataset
#'
#' @inheritParams plotDistribution
getDJInsertionLengthDistribution <- function(dat) {
    return(getInsertionLengths(dat, "np2_length"))
}

#' Plot the DJ insertion length distribution of one or more datasets
#'
#' @inheritParams plotDistribution
plotDJInsertionLengthDistribution <- function(dat_list,
                                   plot_type,
                                   names=NULL,
                                   show_legend=TRUE
                                  ) { 
    p <- plotDistribution(dat_list,
                          getDJInsertionLengthDistribution,
                          plot_type=plot_type,
                          x_label="DJ insertion length",
                          show_legend=show_legend,
                          names=names
                         )
    return(p)
}


#' This calls the function for VD insertions since the AIRR field
#'   for VD and VJ insertion lengths is the same
#'
#' @inheritParams getVDInsertionLengthDistribution
getVJInsertionLengthDistribution <- function(dat) {
    return(getVDInsertionLengthDistribution(dat))
}

#' This calls the function for VD insertions since the AIRR field
#'   for VD and VJ insertion lengths is the same
#'
#' @inheritParams plotVDInsertionLengthDistribution
plotVJInsertionLengthDistribution <- function(dat_list,
                                              plot_type,
                                              names=NULL,
                                              show_legend=TRUE
                                             ) {
    return(plotVDInsertionLengthDistribution(dat,
                                             plot_type=plot_type,
                                             names=names,
                                             show_legend=show_legend
                                            )
          )
}

compareInsertionLengths <- function(dat_a, dat_b, genes) {
    insertion_length_function <- paste0("get", 
                                        genes, 
                                        "InsertionLengthDistribution"
                                       ) %>% 
        get
    dist_a <- dat_a %>% 
        insertion_length_function
    dist_b <- dat_b %>% 
        insertion_length_function
    divergence <- getJSDivergence(dist_a, dist_b)
    return(divergence)
}

compareVDInsertionLengthDistributions <- function(dat_a, dat_b) {
    return(compareInsertionLengths(dat_a, dat_b, "VD"))
}

compareDJInsertionLengthDistributions <- function(dat_a, dat_b) {
    return(compareInsertionLengths(dat_a, dat_b, "DJ"))
}

compareVJInsertionLengthDistributions <- function(dat_a, dat_b) {
    return(compareInsertionLengths(dat_a, dat_b, "VJ"))
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
#' @param dat A \code{data.table} corresponding to repertoire annotations
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

#' Get the Markov transition matrix for DJ insertions
#'
#' @inheritParams getInsertionMatrix
#' @return The empirical transition matrix for each (base, base) pair.
getVJInsertionMatrix <- function(dat) {
    return(getInsertionMatrix(dat, "vj_insertion"))
}

#' Compare the transition matrices for VD insertions for two datasets
#'
#' @param dat_a,dat_b A \code{data.table} corresponding to repertoire annotations
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
#' @param dat_a,dat_b A \code{data.table} corresponding to repertoire annotations
#' @return The mean absolute difference of matrix entries, taken elementwise
compareDJInsertionMatrices <- function(dat_a, dat_b) {
    matrix_a <- dat_a %>% 
        getDJInsertionMatrix
    matrix_b <- dat_b %>% 
        getDJInsertionMatrix
    divergence <- getSumOfAbsoluteDifferences(matrix_a, matrix_b)
    return(divergence)
}

#' Compare the transition matrices for VJ insertions for two datasets
#'
#' @param dat_a,dat_b A \code{data.table} corresponding to repertoire annotations
#' @return The mean absolute difference of matrix entries, taken elementwise
compareVJInsertionMatrices <- function(dat_a, dat_b) {
    matrix_a <- dat_a %>% 
        getVJInsertionMatrix
    matrix_b <- dat_b %>% 
        getVJInsertionMatrix
    divergence <- getSumOfAbsoluteDifferences(matrix_a, matrix_b)
    return(divergence)
}

getClusterSizeDistribution <- function(dat,
                            column="clone_id"
                           ) {
    sizes <- dat[[column]] %>%
        table %>% 
        unname %>% 
        c %>% 
        sort
    return(sizes)
}

#' Plot the cluster size distribution of one or more datasets
#'
#' @inheritParams plotDistribution
plotClusterSizeDistribution <- function(dat_list,
                                        plot_type,
                                        names=NULL,
                                        show_legend=TRUE
                                       ) { 
    p <- plotDistribution(dat_list,
                          getClusterSizeDistribution,
                          plot_type=plot_type,
                          x_label="Cluster size",
                          show_legend=show_legend,
                          names=names,
                          binwidth=1
                         )
    return(p)
}

compareClusterSizeDistributions <- function(dat_a, dat_b) {
    dist_a <- dat_a %>% 
        getClusterSizeDistribution
    dist_b <- dat_b %>% 
        getClusterSizeDistribution
    divergence <- getJSDivergence(dist_a, dist_b)
    return(divergence)
}

getHillNumbers <- function(dat, diversity_orders=c(0, 1, 2)) {
    counts <- dat %>% 
        getClusterSizeDistribution
    diversity <- alakazam::calcDiversity(counts, diversity_orders)
    return(diversity)
}

#' Compare one or multiple Hill numbers of two datasets. The measure of distance
#' is a sum of absolute differences (or just the absolute difference if there is
#' only one diversity_order value
#' @param dat_a,dat_b A \code{data.table} corresponding to repertoire annotations
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
#' @param dat A \code{data.table} corresponding to repertoire annotations
#' @return The percentage of in-frame sequences in \code{dat}
getInFramePercentage <- function(dat,
                                 column="vj_in_frame") {
    percentage <- 100*(dat[[column]] %>%
                       mean)
    return(percentage)
}

#' Compare the percentage of in-frame sequences of two datasets
#'
#' @param dat_a,dat_b A \code{data.table} corresponding to repertoire annotations
#' @return Absolute differences of in-frame percentages
compareInFramePercentages <- function(dat_a, 
                                      dat_b,
                                      column="vj_in_frame"
                                     ) {
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
#' @param dat A \code{data.table} corresponding to repertoire annotations
#' @return Distribution of baseline values, which are measures of selection 
#'   strength. Thus, Sigma = 0 corresponds to no selection, Sigma > 0 
#'   corresponds to positive selection, and Sigma < 0 corresponds to negative
#'   selection
getSelectionEstimate <- function(dat) {
    baseline <- shazam::calcBaseline(dat,
                                     sequenceColumn="sequence_alignment",
                                     germlineColumn="germline_alignment"
                                     ) %>%
        shazam::summarizeBaseline(returnType="df") %$%
        BASELINE_SIGMA 
    return(baseline)
}

#' Compare selection estimates for two datasets
#'
#' @param dat_a,dat_b A \code{data.table} corresponding to repertoire annotations
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
                                     column="junction_aa",
                                     standardize=TRUE
                                    ) {
    sequences <- dat[[column]]
    aa_dist <- paste(sequences[!is.na(sequences)], collapse="") %>%
        strsplit(split="") %>%
        unlist %>%
        table

    if(standardize) {
        aa_dist <- aa_dist/sum(aa_dist)
    }

    return(aa_dist)
}

#' Get the l1 divergence of two categorical distributions, manifest as 
#'   \code{table} objects
#'
#' @param d1,d2 Table containing empirical frequencies 
compareCategoricalDistributions <- function(d1,
                                            d2
                                           ) {
    full_names <- union(names(d1), names(d2))
    missing_1 <- setdiff(full_names, names(d1))
    d1[missing_1] <- 0
    missing_2 <- setdiff(full_names, names(d2))
    d2[missing_2] <- 0

    # If some genes are unnamed, their name is converted to NA.
    # Right now we just ignore these in sum, but we should explicitly remove these
    # at some point
    divergence <- abs(d1[sort(names(d1))] - d2[sort(names(d2))]) %>% sum(na.rm=TRUE)
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
                                         column="junction_aa",
                                         standardize=TRUE
                                        ) {
    sequences <- dat[[column]]
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

#' @inheritParams plotUnivariateDistributions
getUnivariateDistributionDataTable <- function(dat_list,
                                               plot_type="freqpoly",
                                               names=NULL,
                                               plot_function_strings=NULL,
                                               do_all_plots=FALSE
                                              ) {
    if(plot_function_strings %>% is.null) {
        plot_function_strings <- list(
                                      "getPairwiseDistanceDistribution",
                                      "getGCContentDistribution",
                                      "getHotspotCountDistribution",
                                      "getColdspotCountDistribution",
                                      "getDistanceFromGermlineToSequenceDistribution",
                                      "getCDR3LengthDistribution",
                                      "getAliphaticIndexDistribution",
                                      "getGRAVYDistribution",
                                      "getPolarityDistribution",
                                      "getChargeDistribution",
                                      "getBasicityDistribution",
                                      "getAcidityDistribution",
                                      "getAromaticityDistribution",
                                      "getBulkinessDistribution",
                                      "getPositionalDistanceBetweenMutationsDistribution",
                                      "getVGene3PrimeDeletionLengthDistribution",
                                      "getDGene3PrimeDeletionLengthDistribution",
                                      "getDGene5PrimeDeletionLengthDistribution",
                                      "getJGene5PrimeDeletionLengthDistribution",
                                      "getVDInsertionLengthDistribution",
                                      "getDJInsertionLengthDistribution",
                                      "getClusterSizeDistribution"
                                     )
        if(do_all_plots) {
            plot_function_strings <- c(plot_function_strings,
                                       "getNearestNeighborDistribution"
                                      )
        }
    } 

    plot_names <- plot_function_strings %>%
        sapply(getNameFromFunctionString)

    distribution_df <- dat_list %>%
        Map(
            function(dat, dat_name) {
                dat_df <- Map(
                    function(f_string, name) {
                        f <- eval(parse(text=f_string))
                        distribution <- f(dat=dat)
                        dist_dat <- data.table(Value=distribution, 
                                               Dataset=dat_name,
                                               Name=name)
                        return(dist_dat)
                    },
                    plot_function_strings,
                    plot_names
                ) %>%
                    do.call("rbind", .)
               return(dat_df)
            },
            .,
            names
        ) %>%
            do.call("rbind", .)

}

#' Generate a gridded plot of each univariate distribution corresponding to
#'   one or more annotations datasets
#'
#' @inheritParams plotDistribution
#' @param names Strings to be displayed by the legend corresponding to the 
#'   elements of \code{dat_list}
plotUnivariateDistributions <- function(dat_list,
                                        plot_type,
                                        names=NULL,
                                        plot_function_strings=NULL,
                                        do_all_plots=FALSE
                                       ) {
    distribution_dat <- getUnivariateDistributionDataTable(
        dat_list,
        plot_type=plot_type,
        names=names,
        plot_function_strings=plot_function_strings,
        do_all_plots=do_all_plots
    )
    p <- distribution_dat %>% ggplot(aes(x=Value,
                                   group=Dataset,
                                   colour=Dataset
                                  )) +
        facet_wrap( ~ Name, scales="free")
    if(plot_type == "freqpoly") {
        p <- p + geom_freqpoly(aes(y=..density..))
    } else if(plot_type == "ecdf") {
        p <- p + stat_ecdf()
    }

    return(p)
}

getNameFromFunctionString <- function(function_string) {
    name_hash <- list(
        getNearestNeighborDistribution="Nearest neighbor distance",
        getNNDistanceDistribution="Nearest neighbor distance",
        getPairwiseDistanceDistribution="Pairwise distance",
        getCDR3PairwiseDistanceDistribution="CDR3 pairwise dist.",
        getGCContentDistribution="GC content",
        getHotspotCountDistribution="Hotspot count",
        getColdspotCountDistribution="Coldspot count",
        getDistanceFromGermlineToSequenceDistribution="Dist. from germ. to seq.",
        getCDR3LengthDistribution="CDR3 length",
        getAliphaticIndexDistribution="Aliphatic index",
        getGRAVYDistribution="GRAVY index",
        getPolarityDistribution="Polarity",
        getChargeDistribution="Charge",
        getBasicityDistribution="Basicity",
        getAcidityDistribution="Acidity",
        getAromaticityDistribution="Aromaticity",
        getBulkinessDistribution="Bulkiness",
        getInFramePercentage="In frame %",
        getAminoAcidDistribution="AA frequencies",
        getAminoAcid2merDistribution="AA 2mer frequencies",
        getPositionalDistanceBetweenMutationsDistribution="Pos. dist. b/w muts.",
        getVGeneDistribution="V usage",
        getDGeneDistribution="D usage",
        getJGeneDistribution="J usage",
        getVDJDistribution="VDJ usage",
        getVGene3PrimeDeletionLengthDistribution="V 3' del. length",
        getDGene3PrimeDeletionLengthDistribution="D 3' del. length",
        getDGene5PrimeDeletionLengthDistribution="D 5' del. length",
        getJGene5PrimeDeletionLengthDistribution="J 5' del. length",
        getVDInsertionLengthDistribution="VD insertion length",
        getDJInsertionLengthDistribution="DJ insertion length",
        getVDInsertionMatrice="VD insertion matrix",
        getDJInsertionMatrice="DJ insertion matrix",
        getClusterSizeDistribution="Cluster size",

        SubstitutionModelDivergence="Substitution model",
        MutabilityModelDivergence="Mutability model",
        getSelectionEstimate="Selection estimate",

        KideraFactor1Divergence="Kidera factor 1",
        KideraFactor2Divergence="Kidera factor 2",
        KideraFactor3Divergence="Kidera factor 3",
        KideraFactor4Divergence="Kidera factor 4",
        KideraFactor5Divergence="Kidera factor 5",
        KideraFactor6Divergence="Kidera factor 6",
        KideraFactor7Divergence="Kidera factor 7",
        KideraFactor8Divergence="Kidera factor 8",
        KideraFactor9Divergence="Kidera factor 9",
        KideraFactor10Divergence="Kidera factor 10",

        AtchleyFactor1Divergence="Atchley factor 1",
        AtchleyFactor2Divergence="Atchley factor 2",
        AtchleyFactor3Divergence="Atchley factor 3",
        AtchleyFactor4Divergence="Atchley factor 4",
        AtchleyFactor5Divergence="Atchley factor 5"
    )

    return(name_hash[[function_string]])
}
