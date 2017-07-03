library(alakazam)
library(dplyr)
library(RecordLinkage)
library(shazam)
library(stringdist)
library(textmineR)

get.distance.vector <- function(distance.matrix) {
    vec <- distance.matrix[distance.matrix %>% lower.tri] %>% as.vector
    return(vec)
}

compare.pairwise.distance.distribution <- function(list.a, list.b) {
    distance.vector.a <- stringdistmatrix(list.a) %>% as.matrix %>% get.distance.vector
    distance.vector.b <- stringdistmatrix(list.b) %>% as.matrix %>% get.distance.vector
    divergence <- CalcJSDivergence(distance.vector.a, distance.vector.b)
    return(divergence)
}

compare.NN.distance.distribution <- function(list.a, list.b) {
     
}

