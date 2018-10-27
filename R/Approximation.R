#' Subsample a vector
#' 
#' @inheritParams subsample
#' @param v Vector to subsample
#' @return A subsampled vector 
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
#'   So, we cluster the sequences using vsearch and sample sequences based
#'   on cluster frequency. This ensures no nearest neighbor distance is higher
#'   than the given cutoff.
getApproximateNearestNeighborDistribution <- function(dat,
                                                      column,
                                                      k=1,
                                                      pairwise_identity=0.5,
                                                      tol=0.01,
                                                      batch_size=100
                                                     ) {
    seq_fasta_filename <- "approx_nn_seq_tmp.fa"
    cluster_out_filename <- "approx_nn_clustered_tmp.fa"
    dat[[column]] %>%
        as.list %>%
        write.fasta(sequences=.,
                    names=1:nrow(dat),
                    file.out=seq_fasta_filename
                   )

    cluster_command <- paste("vsearch",
                             "--cluster_fast",
                             seq_fasta_filename,
                             "--clusterout_id",
                             "--consout",
                             cluster_out_filename,
                             "--id",
                             pairwise_identity,
                             "--usersort"
                            )

    print(cluster_command)

    cluster_command %>%
        system

    usearch_out_filename <- "approx_nn_usearch_tmp.uc"
    usearch_command <- paste("vsearch",
                             "--usearch_global",
                             seq_fasta_filename,
                             "--db",
                             cluster_out_filename,
                             "--id",
                             pairwise_identity,
                             "--uc",
                             usearch_out_filename,
                             "--uc_allhits"
                            )
    
    print(usearch_command)
    usearch_command %>% system

    cluster_df <- usearch_out_filename %>%
        data.table::fread()

    # vsearch doesn't include column names in the output tsv :(
    entry_type_index <- 1
    cluster_id_index <- 2
    sequence_id_index <- 9
    names(cluster_df)[entry_type_index] <- "entry_type"
    names(cluster_df)[cluster_id_index] <- "cluster_id" 
    names(cluster_df)[sequence_id_index] <- "sequence_id" 

    cluster_df <- cluster_df %>%
        filter(entry_type == "H") # Keep only sequences, not cluster centroids

    # Include sequences in clustered data.table based on sequence_ids from
    #   vsearch
    cluster_df$sequence <- dat[[column]][cluster_df$sequence_id]
    cluster_freqs <- (cluster_df$cluster_id %>% 
                          table)/sum(cluster_df$cluster_id %>% 
                                         table)
    cluster_df$cluster_freq <- cluster_df$cluster_id %>%
        sapply(function(x) {
                   cluster_freqs[x %>% toString]
               }
        )

    cluster_id_table <- cluster_df$cluster_id %>% 
        table
    singleton_cluster_ids <- names(cluster_id_table)[which(cluster_id_table == 1)] %>%
        as.numeric
    singleton_df <- cluster_df[cluster_df$cluster_id %in% 
                                         singleton_cluster_ids, ]
    singleton_nn_dists <- singleton_df$sequence %>%
        sapply(function(sequence) {
                   stringdist::stringdist(
                       a=singleton_df[singleton_df$sequence != sequence, ]$sequence,
                       b=sequence,
                       method="lv"
                   ) %>%
                   min
        }) %>% c
    singleton_df$nn_dist <- singleton_nn_dists

    # Now, sample sequences by cluster frequency and compute nearest neighbor 
    #   distance within its cluster. Do this for batch_size sequences at a 
    #   time, and stop when the KL-divergence of successive iterates is below
    #   tol.
    error <- Inf
    dist_current <- cluster_df %>% 
        doNNSubsamplingBatchStep(batch_size=batch_size,
                                 singleton_df=singleton_df
                                )
    while(error > tol) {
        dist_prev <- dist_current
        dist_current <- c(dist_current,
                          cluster_df %>%
                              doNNSubsamplingBatchStep(batch_size=batch_size,
                                                       singleton_df=singleton_df
                                                      )
                         )
        dist_current <- c(dist_prev, dist_current)
        error <- getJSDivergence(dist_prev, dist_current)
    }
    
    file.remove(seq_fasta_filename,
                cluster_out_filename,
                usearch_out_filename
               )

    return(dist_current)
}

doNNSubsamplingBatchStep <- function(dat,
                                     batch_size,
                                     singleton_df
                                    ) {
    nn_sample <- replicate(batch_size,
        {
            sequence_id <- sample(dat$sequence_id, 
                                  1, 
                                  prob=dat$cluster_freqs
                                 )
            if(sequence_id %in% singleton_df$sequence_id) {
                seq_df <- singleton_df[singleton_df$sequence_id == sequence_id, ]
                nn_dist <- singleton_df[singleton_df$sequence_id == sequence_id, ]$nn_dist
            } else {
                seq_dat <- dat[dat$sequence_id == sequence_id, ]
                cluster <- dat[dat$cluster_id == seq_dat$cluster_id &
                               dat$sequence_id != seq_dat$sequence_id, ]
                nn_dist <- stringdist::stringdist(cluster$sequence,
                                                   seq_dat$sequence,
                                                   method="lv"
                                                  ) %>%
                    min
            }
            nn_dist
        }
    )

    return(nn_sample)
}


