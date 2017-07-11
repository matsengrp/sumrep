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

# Hot/cold spot motif data from S5F model given my Yaari et al. 2013
motifs <- read.table("Data/Mutability.csv", header=TRUE)
motifs.2 <- read.table("Data/Substitution.csv", header=TRUE)

motifs$Rate <- motifs$Mutability/sum(motifs$Mutability)

get.continuous.JS.divergence <- function(s1, s2) {
    p <- s1 %>% density %>% approxfun()
    q <- s2 %>% density %>% approxfun()
    m <- function(x) {
        res <- 0.5*(p(x) + q(x))
        return(res)
    }

    integrand <- function(f, g) {
        func <- function(x) {
            res <- ifelse(f(x) == 0, 0, ifelse(g(x) == 0, Inf, f(x)*(log(f(x)) - log(g(x)))))
            return(res)
        }
        return(func)
    }
    lower <- max(min(s1), min(s2))
    upper <- min(max(s1), max(s2))
    KL.div.1 <- integrate(integrand(p, m), lower, upper)$value
    KL.div.2 <- integrate(integrand(q, m), lower, upper)$value
    divergence <- 0.5*(KL.div.1 + KL.div.2)
    return(divergence)
}

get.JS.divergence <- function(l1, l2, continuous=FALSE) {
    if(continuous) {
        return(get.continuous.JS.divergence(l1, l2))
    }
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
    density.a <- get.GC.distribution(list.a)
    density.b <- get.GC.distribution(list.b)
    divergence <- get.continuous.JS.divergence(density.a, density.b)
    return(divergence)
}

