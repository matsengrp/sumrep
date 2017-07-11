library(alakazam)
library(ape)
library(dplyr)
library(flexmix)
library(RecordLinkage)
library(shazam)
library(stringdist)
library(textmineR)

if(!exists("biocLite")) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("Biostrings")
}

# Hot/cold spot motif data from S5F model given by Yaari et al. 2013
motifs <- read.table("Data/Mutability.csv", header=TRUE)
motifs.2 <- read.table("Data/Substitution.csv", header=TRUE)

motifs$Rate <- motifs$Mutability/sum(motifs$Mutability)

get.continuous.JS.divergence <- function(sample.1, sample.2) {
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

    p <- sample.1 %>% density %>% approxfun
    q <- sample.2 %>% density %>% approxfun
    lower <- max(min(sample.1), min(sample.2))
    upper <- min(max(sample.1), max(sample.2))
    KL.div.1 <- integrate(integrand(p, m), lower, upper)$value
    KL.div.2 <- integrate(integrand(q, m), lower, upper)$value
    JS.divergence <- 0.5*(KL.div.1 + KL.div.2)
    return(JS.divergence)
}

get.JS.divergence <- function(l1, l2, continuous=FALSE) {
    if(continuous) {
        divergence <- get.continuous.JS.divergence(l1, l2)
    } else {
        p <- l1/(l1 %>% sum)
        q <- l2/(l2 %>% sum)
        m <- (p + q)/2
        KL.div.1 <- (p*log(p/m)) %>% sum
        KL.div.2 <- (q*log(q/m)) %>% sum
        divergence <- (KL.div.1 + KL.div.2)/2
    }
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
    table.a <- distances.a %>% factor(levels=0:max.val) %>% table %>% as.vector
    table.b <- distances.b %>% factor(levels=0:max.val) %>% table %>% as.vector
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

compare.NN.distance.distribution <- function(list.a, list.b, k=1) {
    distances.a <- list.a %>% get.nearest.neighbor.distances(k=k)
    distances.b <- list.b %>% get.nearest.neighbor.distances(k=k)
    divergence <- CalcJSDivergence(distances.a, distances.b)
    return(divergence)
}

get.GC.distribution <- function(sequence.list) {
    dna.list <- sequence.list %>% strsplit(split='') %>% lapply(as.DNAbin)
    gc.dist <- dna.list %>% sapply(GC.content)
    return(gc.dist)
}

# In-progress
compare.GC.distributions <- function(list.a, list.b) {
    density.a <- list.a %>% get.GC.distribution
    density.b <- list.b %>% get.GC.distribution
    divergence <- get.continuous.JS.divergence(density.a, density.b)
    return(divergence)
}

