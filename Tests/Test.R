source("SummaryStats.R")
library(dplyr)

seq.1 <- c(A="AT", B="AT")
seq.2 <- list(A="AA", B="AC")
seq.3 <- c(A="", B="C-TA")
seq.4 <- c(A="", B="")
seq.5 <- c("AAA", "AAT", "TTT")
seq.6 <- c("A-T-", "GCAC", "CGGC", "ACGT")
seq.7 <- list(c("G", "T"), c("T", "A"))
seq.8 <- c("AAAAAA", "TTTTTT", "AAAA", "AAAAACC", "AAAAAAATTTT", "GGGCCCCCGG")
seq.9 <- list("AAAAAC", c("T", "T", "T", "C", "C", "T"), "AAAAGGGG", "AAAAACCCC", "CCCCTTTC", 
              list("G", "G", "G", "C", "C", "C", "C", "C", "G", "G"))

apply.row.column.names <- function(m) {
    y <- m
    rownames(y) <- colnames(y) <- 1:dim(y)[1]
    return(y)
}

test.bin.continuous.lists.as.discrete <- function() {
    l1 <- 1:20
    l2 <- rep(1, 20)
    binned <- bin.continuous.lists.as.discrete(l1, l2)
    bin.a <- binned[[1]]
    bin.b <- binned[[2]]
    checkEquals(bin.a, rep(4, 5))
    checkEquals(bin.b, c(20, rep(0, 4)))
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

test.get.nucleotide.diversity <- function() {
    d1 <- get.nucleotide.diversity(list("AT", "AT"))
    d2 <- get.nucleotide.diversity(list("AT", "AC"))
    d3 <- get.nucleotide.diversity(list("AT", "GC"))
    d4 <- get.nucleotide.diversity(list("ATATATATAT", "AAAAAAAAAA"))
    d5 <- get.nucleotide.diversity(list("ATATATATAT", "CGTACGTAAT"))
    checkEquals(d1, 0)
    checkEquals(d2, 0.5)
    checkEquals(d3, 1)
    checkTrue(d4 < d5)
}

test.compare.nucleotide.diversities <- function() {
    r1 <- list("AT", "AT")
    r2 <- list("AT", "AC")
    r3 <- list("AT", "GC")
    r4 <- list("ATATATATAT", "AAAAAAAAAA")
    r5 <- list("ATATATATAT", "CGTACGTAAT")
    checkEquals(compare.nucleotide.diversities(r1, r1), 0)
    checkEquals(compare.nucleotide.diversities(r1, r2), 0.5)
    checkEquals(compare.nucleotide.diversities(r2, r3), 0.5)
}

test.get.distances.from.naive.to.mature <- function() {
    naive <- c("AAAAAA")
    m1 <- c("AAAAAT")
    m2 <- c("CGCAAA")
    m3 <- c("GGGGGG")
    m4 <- c("AAAAAA")

    checkEquals(get.distances.from.naive.to.mature(naive, list(m1, m2, m3, m4)), c(0, 1, 3, 6))
}

test.compare.distances.from.naive.to.mature <- function() {
    naive.a <- c("AAAAAA")
    naive.b <- c("GGGGGG")
    m1 <- c("AAAAAA")
    m2 <- c("AAAAAG")
    m3 <- c("GAAAAA")
    m4 <- c("AAAGGG")
    m5 <- c("GGGG")
    m.list.1 <- list(m1, m2, m3)
    m.list.2 <- list(m2, m3, m4) 
    c1 <- compare.distances.from.naive.to.mature(naive.a, m.list.1, naive.a, m.list.1)
    c2 <- compare.distances.from.naive.to.mature(naive.a, m.list.1, naive.b, m.list.2)
    c3 <- compare.distances.from.naive.to.mature(naive.b, m.list.2, naive.a, m.list.1)
    checkEquals(c1, 0)
    checkEquals(c2, c3)
    checkTrue(c2 > 0)
}
