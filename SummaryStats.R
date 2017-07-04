library(alakazam)
library(dplyr)
library(RecordLinkage)
library(shazam)
library(stringdist)
library(textmineR)

get.distance.matrix <- function(sequence.list) {
    distance.matrix <- sequence.list %>% stringdistmatrix %>% as.matrix
    return(distance.matrix)
}

get.distance.vector <- function(sequence.list) {
    mat <- sequence.list %>% get.distance.matrix
    vec <- mat[mat %>% lower.tri] %>% as.vector
    return(vec)
}

compare.pairwise.distance.distribution <- function(list.a, list.b) {
    distance.vector.a <- list.a %>% get.distance.vector
    distance.vector.b <- list.b %>% get.distance.vector
    divergence <- CalcJSDivergence(distance.vector.a, distance.vector.b)
    return(divergence)
}

get.nearest.neighbor.distances <- function(sequence.list) {
    mat <- sequence.list %>% get.distance.matrix
    n <- sequence.list %>% length
    distances <- rep(NA, n)
    for(i in 1:n) {
        distances[i] <- mat[i, which.min(mat[i, -i])]
    }
    return(distances)
}

compare.NN.distance.distribution <- function(list.a, list.b) {
    distances.a <- get.nearest.neighbor.distances(list.a) 
    distances.b <- get.nearest.neighbor.distances(list.b)
    divergence <- CalcJSDivergence(distances.a, distances.b)
    return(divergence)
}

