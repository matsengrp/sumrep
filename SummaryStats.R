library(alakazam)
library(dplyr)
library(RecordLinkage)
library(shazam)
library(stringdist)
library(textmineR)

get.JS.divergence <- function(l1, l2) {
    p <- l1/(l1 %>% sum)
    q <- l2/(l2 %>% sum)
    m <- (p + q)/2
    KL.div.1 <- (p*log(p/m)) %>% sum
    KL.div.2 <- (q*log(q/m)) %>% sum
    divergence <- (KL.div.1 + KL.div.2)/2
    return(divergence)
}

get.distance.matrix <- function(sequence.list) {
    mat <- sequence.list %>% stringdistmatrix %>% as.matrix
    return(mat)
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
        distances[i] <- min(mat[i, -i]) 
    }
    return(distances)
}

compare.NN.distance.distribution <- function(list.a, list.b) {
    distances.a <- get.nearest.neighbor.distances(list.a) 
    distances.b <- get.nearest.neighbor.distances(list.b)
    divergence <- CalcJSDivergence(distances.a, distances.b)
    return(divergence)
}

