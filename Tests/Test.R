source("SummaryStats.R")
library(dplyr)

seq.1 <- c(A="AT", B="AT")
seq.2 <- list(A="AA", B="AC")
seq.3 <- c(A="", B="C-TA")
seq.4 <- c(A="", B="")
seq.5 <- c("AAA", "AAT", "TTT")
seq.6 <- c("A-T-", "GCAC", "CGGC", "ACGT")
seq.7 <- list(c("G", "T"), c("T", "A"))
seq.8 <- c("AAAAAA", "TTTTTT", "AAAAGGGG", "AAAAACCCC", "CCCCCCCCC", "GGGCCCCCGG")
seq.9 <- list("AAAAAC", c("T", "T", "T", "C", "C", "T"), "AAAAGGGG", "AAAAACCCC", "CCCCTTTC", 
              list("G", "G", "G", "C", "C", "C", "C", "C", "G", "G"))

apply.row.column.names <- function(m) {
    y <- m
    rownames(y) <- colnames(y) <- 1:dim(y)[1]
    return(y)
}

test.get.distance.vector <- function() {
    checkEquals(get.distance.vector(seq.1), 0)
    checkEquals(get.distance.vector(seq.2), 1)
    checkEquals(get.distance.vector(seq.3), 4)
    checkEquals(get.distance.vector(seq.4), 0)
    checkEquals(get.distance.vector(seq.7), 2)
}

test.get.distance.matrix <- function() {
    m <- matrix(0, 2, 2) %>% apply.row.column.names
    checkEquals(get.distance.matrix(seq.1), m)

    m2 <- matrix(c(0, 1, 1, 0), 2, 2) %>% apply.row.column.names
    checkEquals(get.distance.matrix(seq.2), m2)

    m3 <- matrix(c(0, 4, 4, 0), 2, 2) %>% apply.row.column.names
    checkEquals(get.distance.matrix(seq.3), m3)

    m4 <- matrix(0, 2, 2) %>% apply.row.column.names
    checkEquals(get.distance.matrix(seq.4), m4)
}

test.get.nearest.neighbor.distances <- function() {
    d1 <- c(0, 0)
    checkEquals(get.nearest.neighbor.distances(seq.1), d1)

    d2 <- c(1, 1, 2)    
    checkEquals(get.nearest.neighbor.distances(seq.5) %>% sort, d2)

    checkEquals(get.nearest.neighbor.distances(seq.6, k=1), c(3, 3, 3, 3))
    checkEquals(get.nearest.neighbor.distances(seq.6, k=2), c(4, 3, 3, 3))
    checkEquals(get.nearest.neighbor.distances(seq.6, k=3), c(4, 4, 4, 3))
}

test.compare.pairwise.distance.distribution <- function() {
    s1 <- c("AAA", "AAT", "ATT")
    s2 <- c("ATT", "AAA", "AAT")
    s3 <- c("AAA", "AAT", "TTT")
    s4 <- c("AAA", "ATC", "GGG")
    checkEquals(compare.pairwise.distance.distribution(s1, s2), 0)    
    
    c1 <- compare.pairwise.distance.distribution(s1, s3)
    c2 <- compare.pairwise.distance.distribution(s1, s4)
    checkTrue(c1 < c2)

    c3 <- compare.pairwise.distance.distribution(s4, s1)
    checkEquals(c2, c3)
}

test.get.GC.distribution <- function() {
    d1 <- rep(0, 3)
    checkEquals(get.GC.distribution(seq.5), d1)

    d2 <- c(0, 0.75, 1.0, 0.5)
    checkEquals(get.GC.distribution(seq.6), d2)
}

test.compare.GC.distributions <- function() {
    checkEquals(compare.GC.distributions(seq.2, seq.7), 0)    
    c1 <- compare.GC.distributions(seq.8, seq.9)
    c2 <- compare.GC.distributions(seq.9, seq.8)
    checkTrue(c1 > 0)
    checkTrue(c1 == c2)
}

test.get.hotspot.count <- function() {
    seq.a <- "TTTTT"
    seq.b <- list("AAC", "TTTT")
    seq.c <- "NWRC"
    seq.d <- c("WRC", "WRC", "WGC")
    checkEquals(get.hotspot.count(seq.a), 0)
    checkEquals(get.hotspot.count(seq.b), 2)
    checkEquals(get.hotspot.count(seq.c), 3)
    checkEquals(get.hotspot.count(seq.d), 5)
}

test.get.coldspot.count <- function() {
    seq.a <- "TTTTT"
    seq.b <- list("GCC", "AAA", "AAAA")
    seq.c <- c("AAA", "CTC")
    seq.d <- "SYCSYC"
    checkEquals(get.coldspot.count(seq.a), 0)
    checkEquals(get.coldspot.count(seq.b), 1)
    checkEquals(get.coldspot.count(seq.c), 1)
    checkEquals(get.coldspot.count(seq.d), 4)
}
