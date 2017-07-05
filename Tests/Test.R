source("SummaryStats.R")
library(dplyr)

db <- read.csv("Data/test_data.csv")

seqs <- db$SEQUENCE_VDJ %>% sapply(toString)

seq.1 <- c(A="AT", B="AT")
seq.2 <- c(A="AA", B="AC")
seq.3 <- c(A="", B="C-TA")
seq.4 <- c(A="", B="")
seq.5 <- c("AAA", "AAT", "TTT")

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
}
