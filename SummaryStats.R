library(alakazam)
library(shazam)
library(textmineR)

compare.pairwise.distance.distribution <- function(list.a, list.b) {
    distance.matrix.a <- pairwiseDist(list.a)
    distance.vector.a <- as.vector(distance.matrix.a[lower.tri(distance.matrix.a)])
    distance.matrix.b <- pairwiseDist(list.b)
    distance.vector.b <- as.vector(distance.matrix.b[lower.tri(distance.matrix.b)])
    divergence <- CalcJSDivergence(distance.vector.a, distance.vector.b)
}

sequence.1 <- c(A="ATGGC", B="ATGGG", C="ATGGG", D="AT--C")
sequence.2 <- c(A="TTGGC", B="ATGGG", C="ATGGG", D="ATGAC")

sequence.3 <- c(A="CGAAC", B="GGAAT", C="TACCC", D="ATAGC")

c1 <- compare.pairwise.distance.distribution(sequence.1, sequence.2)
c2 <- compare.pairwise.distance.distribution(sequence.1, sequence.3)

cat(c1, c2, '\n')
