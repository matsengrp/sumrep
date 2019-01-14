#' Subsample a vector
#' 
#' @inheritParams subsample
#' @param v Vector to subsample
#' @return A subsampled vector
#' 
#' @export
subsampleVector <- function(v, 
                            sample_count,
                            replace=TRUE
                           ) {
    new_vector <- {}
    if(length(v) > sample_count) {
        new_vector <- sample(v, 
                             sample_count,
                             replace
                            )
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
#'
#' @export
subsample <- function(dataset, 
                      sample_count,
                      replace=TRUE
                     ) {
    new_dataset <- {}
    if(is.null(dim(dataset))) {
        new_dataset <- subsampleVector(dataset, 
                                       sample_count,
                                       replace=replace
                                      )
    } else if(nrow(dataset) > sample_count) {
        new_dataset <- dataset[sample(1:nrow(dataset), 
                                      sample_count,
                                      replace=replace
                                     ), ]
    } else {
        new_dataset <- dataset
    }
    return(new_dataset)
}

#' Approximate the distribution of \code{summary_function} for \code{dat} by 
#    iteratively sampling from \code{dat}, computing the distribution of 
#'   \code{summary_function} for the subsample, and accumulating a representative
#'   distribution. The function stops when \code{divergence_function} between two
#'   distribution iterates is sufficiently small.
#'
#' @param dat The dataset for which the distribution is desired
#' @param summary_function A function which computes a summary statistic of \code{dat},
#'   which yields a distribution
#' @param sample_count The size of the subsample for each iteration
#' @param tol The threshold for \code{divergence_function} to be considered converged
#' @param divergence_function The divergence computed between successive iterates.
#'   This should probably be the JS divergence, specified for discrete or continuous data
#'   depending on the nature of \code{summary_function}
#'
#' @export
getApproximateDistribution <- function(dat,
                                       summary_function,
                                       sample_count=100,
                                       tol=0.001,
                                       divergence_function=getContinuousJSDivergence,
                                       ...
                                       ) {
    dist <- subsample(dat, sample_count) %>%
        summary_function(...)

    error <- Inf
    while(error > tol) {
        dist_prev <- dist
        sample_dat <- subsample(dat, sample_count)
        sample_dist <- sample_dat %>% 
            summary_function(...)
        dist <- c(dist_prev, sample_dist)
        error <- divergence_function(dist, dist_prev)
    }

    return(dist)
}

#' Runs a special approximation routine since nearest neighbor distances are
#'   not a regular distribution for subsampling. (Since the nearest neighbor
#'   distance of a sample is necessarily larger than that of the population)
#'
#' @param dat 
#' @param dat The dataset for which the distribution is desired
#' @param column The column over which the nearest neighbor distance
#'   should be computed. Must correspond to sequence strings, e.g. DNA or
#'   amino acid sequences.
#' @param tol The allowed tolerance for determining convergence.
#' @param batch_size The number of sequences sampled during each iteration.
#' @return Approximate nearest neighbor distance distribution of sequences
#'   given by \code{dat[[column]]}
#'
#' @export
getApproximateNearestNeighborDistribution <- function(dat,
                                                      column,
                                                      tol=10e-4,
                                                      batch_size=50
                                                     ) {
    d <- dat %>%
        doNNSubsamplingBatchStep(column=column)
    error <- Inf
    while(error > tol) {
        dist_prev <- d
        d <- c(dist_prev, 
               dat %>% 
                   doNNSubsamplingBatchStep(column=column)
               )
        error <- getJSDivergence(dist_prev, d)
    }
    
    return(d)
}

#' Compute the nearest neighbor distribution of a sample of \code{column} of 
#'   size \code{batch_size} to the full set of sequences in \code{column}.
#'
#' @inheritParams getApproximateNearestNeighborDistribution
#' @return Sample of \code{batch_size} nearest neighbor distances from 
#'   \code{dat[[column]]} 
#'
#' @export
doNNSubsamplingBatchStep <- function(dat,
                                     batch_size=10,
                                     column
                                    ) {
    nn_sample <- replicate(batch_size,
        {
            unique_id <- sample(nrow(dat), 1)
            seq_dat <- dat[unique_id, ]
            nn_dist <- stringdist::stringdist(dat[-unique_id, ][[column]],
                                              seq_dat[[column]],
                                              method="lv"
                                             ) %>%
                min
        }
    )

    return(nn_sample)
}
