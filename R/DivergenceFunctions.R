library(entropy)

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
    tabulateData <- function(data_list) {
        result <- data_list %>%
            table %>%
            unname %>%
            as.vector
        return(result)
    }

    a_length <- list_a %>% 
        length
    b_length <- list_b %>% 
        length
    bin_count <- min(a_length, b_length) %>% 
        sqrt %>% 
        ceiling
    bin_count <- ifelse(bin_count < 2, 2, bin_count)
    bins <- c(list_a, list_b) %>% 
        cut(breaks=bin_count, labels=1:bin_count)
    table_a <- bins[1:a_length] %>%
        tabulateData
    table_b <- bins[-(1:a_length)] %>%
        tabulateData
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
getContinuousJSDivergenceByIntegration <- function(sample_1, sample_2) {
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

    p <- sample_1 %>% 
        density %>% 
        approxfun
    q <- sample_2 %>% 
        density %>% 
        approxfun
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
getJSDivergence <- function(list_a, list_b, continuous=FALSE, KL=FALSE) {
    tabulateDiscreteData <- function(data_list, factor_levels) {
        result <- data_list %>%
            factor(levels=factor_levels) %>%
            table %>%
            as.vector
        return(result)
    }

    if(continuous) {
        binned <- binContinuousListsAsDiscrete(list_a, list_b)
        divergence <- textmineR::CalcJSDivergence(binned[[1]], binned[[2]])
    } else {
        factor_levels <- c(list_a, list_b) %>%
            unique %>%
            as.factor
        table_a <- list_a %>%
            tabulateDiscreteData(factor_levels=factor_levels)
        table_b <- list_b %>%
            tabulateDiscreteData(factor_levels=factor_levels)
        if(KL) {
            divergence <- entropy::KL.empirical(table_a, table_b)
        } else {
            divergence <- textmineR::CalcJSDivergence(table_a, table_b)
        }
    }
    return(divergence)
}

getContinuousJSDivergence <- function(list_a, list_b) {
    return(getJSDivergence(list_a, list_b, continuous=TRUE))
}

#' Get the mean absolute difference of two sequences 
#' 
#' @param sequence_a First sequence, either a vector or matrix
#' @param sequence_b Second sequence
#' @param ignore_na Should we ignore na values when computing the mean?
getSumOfAbsoluteDifferences <- function(sequence_a, sequence_b, ignore_na=TRUE) {
    difference <- (sequence_a - sequence_b) %>% 
        abs %>% 
        sum(na.rm=ignore_na)
    return(difference)
}

#' Estimate a divergence by subsampling and averaging over the comparison
#'   function output
#'
#' Arguments to \code{func} may be supplied after the other required 
#'   parameters.
#' @param dataset_a First dataset (vector, matrix, or data.table) input to 
#'   \code{func}
#' @param dataset_b Second dataset input to \code{func}
#' @param summary_function Comparison function for which to calculate the JS 
#'   divergence
#' @param subsample_count Number of samples to extract on each trial
#' @param trial_count Number of trials to average over
#' @param divergenceFunction The divergence function to apply to the two 
#'   datasets. JS divergence by default, although this only makes sense for
#'   vectors corresponding to empirical distributions
#' @return The average computed divergence across trials, used as an
#'   estimate of the true JS divergence
#' @examples
#' seq_1 <- c("AAAAA", "AATTT", "TTATT", "ACATG", "GGGAA")
#' seq_2 <- c("AAAAC", "ACATG", "CGGGA", "ACATG", "GGACA")
#' getAverageDivergence(seq_1, seq_2, getNearestNeighborDistances, 
#'   subsample_count=3, trial_count=20, k=2)
getAverageDivergence <- function(dataset_a, 
                                 dataset_b, 
                                 summary_function, 
                                 subsample_count, 
                                 trial_count, 
                                 divergenceFunction=getJSDivergence,
                                 ...) {
    divergences <- {}
    for(trial in 1:trial_count) {
        summary_a <- dataset_a %>% 
            subsample(sample_count=subsample_count) %>%
            summary_function
        summary_b <- dataset_b %>% 
            subsample(sample_count=subsample_count) %>%
            summary_function
        divergences[trial] <- divergenceFunction(summary_a, summary_b)
    }
    return(divergences %>% mean)
}

#' Estimate the divergence by subsampling and averaging until a stable estimate
#'   is attained
#'
#' @inheritParams getAverageDivergence
#' @param tolerance The criterion to determine convergence based on relative
#'   difference between iterates
#' @return The average computed divergence across trials, used as an
#'   estimate of the true JS divergence
getAutomaticAverageDivergence <- function(dataset_a,
                                          dataset_b,
                                          summary_function,
                                          subsample_count,
                                          divergenceFunction=getJSDivergence,
                                          tolerance=1e-2,
                                          ...) {

    rollingAverage <- function(old_average, n, new_value) {
        new_average <- (old_average*(n - 1) + new_value)/n
        return(new_average)
    }

    n <- 0
    error <- Inf
    average_divergence <- 0
    divergences <- {}
    while(error >= tolerance) {
        old_average <- average_divergence
        n <- n + 1
        summary_a <- dataset_a %>%
            subsample(sample_count=subsample_count) %>%
            summary_function
        summary_b <- dataset_b %>%
            subsample(sample_count=subsample_count) %>%
            summary_function
        divergence <- divergenceFunction(summary_a, summary_b)
        average_divergence <- rollingAverage(old_average, n, divergence)
        error <- abs(average_divergence - old_average)/average_divergence
    }
    return(average_divergence)
}
