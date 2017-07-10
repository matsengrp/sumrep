library(alakazam)
library(ape)
library(dplyr)
library(flexmix)
library(RecordLinkage)
library(shazam)
library(stringdist)
library(textmineR)

get.JS.divergence <- function(l1, l2, continuous=FALSE) {
    p <- l1/(l1 %>% sum)
    q <- l2/(l2 %>% sum)
    m <- (p + q)/2
    KL.div.1 <- (p*log(p/m)) %>% sum
    KL.div.2 <- (q*log(q/m)) %>% sum
    divergence <- (KL.div.1 + KL.div.2)/2
    return(divergence)
}

get.distance.matrix <- function(sequence.list) {
    length.count <- sequence.list %>% sapply(nchar) %>% table %>% length 
    comparison.method <- ifelse(length.count > 1, "lv", "hamming")
    mat <- sequence.list %>% stringdistmatrix(method=comparison.method) %>% as.matrix
    return(mat)
}

get.distance.vector <- function(sequence.list) {
    mat <- sequence.list %>% get.distance.matrix
    vec <- mat[mat %>% lower.tri] %>% as.vector %>% sort
    return(vec)
}

compare.pairwise.distance.distribution <- function(list.a, list.b) {    
    distances.a <- list.a %>% get.distance.vector
    distances.b <- list.b %>% get.distance.vector
    max.val <- max(distances.a, distances.b)
    table.a <- factor(distances.a, levels=0:max.val) %>% table %>% as.vector
    table.b <- factor(distances.b, levels=0:max.val) %>% table %>% as.vector
    divergence <- CalcJSDivergence(table.a, table.b)
    return(divergence)
}

get.nearest.neighbor.distances <- function(sequence.list, k=1) {
    mat <- sequence.list %>% get.distance.matrix
    n <- sequence.list %>% length
    distances <- rep(NA, n)
    for(i in 1:n) {
        distances[i] <- sort(mat[i, -i], partial=k)[k]
    }
    return(distances)
}

compare.NN.distance.distribution <- function(list.a, list.b) {
    distances.a <- get.nearest.neighbor.distances(list.a) 
    distances.b <- get.nearest.neighbor.distances(list.b)
    divergence <- CalcJSDivergence(distances.a, distances.b)
    return(divergence)
}

get.GC.distribution <- function(sequence.list) {
    dna.list <- sequence.list %>% strsplit(split='') %>% lapply(as.DNAbin)
    return(sapply(dna.list, GC.content))
}

# In-progress
compare.GC.distributions <- function(list.a, list.b) {
    density.a <- get.GC.distribution(list.a) %>% density
    density.b <- get.GC.distribution(list.b) %>% density
    divergence <- CalcJSDivergence(density.a, density.b)
    return(divergence)
}

