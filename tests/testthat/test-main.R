library(sumrep)
context("Summary statistics")

test_that("test.binContinuousListsAsDiscrete", {
    l1 <- 1:20
    l2 <- rep(1, 20)
    binned <- binContinuousListsAsDiscrete(l1, l2)
    bin_a <- binned[[1]]
    bin_b <- binned[[2]]
    expect_equal(rep(4, 5), bin_a)
    expect_equal(c(20, rep(0, 4)), bin_b)
})

test_that("test.compareDistancesFromNaiveToMature", {
    naive_a <- c("AAAAAA")
    naive_b <- c("GGGGGG")
    m1 <- c("AAAAAA")
    m2 <- c("AAAAAG")
    m3 <- c("GAAAAA")
    m4 <- c("AAAGGG")
    m5 <- c("GGGG")
    dat_a <- data.table(naive_seq=naive_a, mature_seq=list(m1, m2, m3))
    dat_b <- data.table(naive_seq=naive_b, mature_seq=list(m2, m3, m4))
    c1 <- compareDistancesFromNaiveToMature(dat_a, dat_a)
    c2 <- compareDistancesFromNaiveToMature(dat_a, dat_b)
    c3 <- compareDistancesFromNaiveToMature(dat_b, dat_a)
    expect_equal(0, c1)
    expect_equal(c3, c2)
    expect_true(c2 > 0)
})

test_that("test.compareGCDistributions", {
    expect_equal(0, compareGCDistributions(seq_2, seq_7))
    c1 <- compareGCDistributions(seq_8, seq_9)
    c2 <- compareGCDistributions(seq_9, seq_8)
    expect_true(c1 > 0)
    expect_true(c1 == c2)
})

test_that("test.compareNucleotideDiversities", {
    r1 <- list("AT", "AT")
    r2 <- list("AT", "AC")
    r3 <- list("AT", "GC")
    r4 <- list("ATATATATAT", "AAAAAAAAAA")
    r5 <- list("ATATATATAT", "CGTACGTAAT")
    expect_equal(0, compareNucleotideDiversities(r1, r1))
    expect_equal(0.5, compareNucleotideDiversities(r1, r2))
    expect_equal(0.5, compareNucleotideDiversities(r2, r3))
})

test_that("test.comparePairwiseDistanceDistributions", {
    s1 <- c("AAA", "AAT", "ATT")
    s2 <- c("ATT", "AAA", "AAT")
    s3 <- c("AAA", "AAT", "TTT")
    s4 <- c("AAA", "ATC", "GGG")
    expect_equal(0, comparePairwiseDistanceDistributions(s1, 
        s2))
    c1 <- comparePairwiseDistanceDistributions(s1, s3)
    c2 <- comparePairwiseDistanceDistributions(s1, s4)
    expect_true(c1 < c2)
    c3 <- comparePairwiseDistanceDistributions(s4, s1)
    expect_equal(c3, c2)
})

test_that("test.getColdspotCount", {
    seq_a <- "TTTTT"
    seq_b <- list("GCC", "AAA", "AAAA")
    seq_c <- c("AAA", "CTC")
    seq_d <- "SYCSYC"
    expect_equal(0, getColdspotCount(seq_a))
    expect_equal(1, getColdspotCount(seq_b))
    expect_equal(1, getColdspotCount(seq_c))
    expect_equal(4, getColdspotCount(seq_d))
})

test_that("test.getDistanceMatrix", {
    m <- matrix(0, 2, 2) %>% applyRowAndColumnNames
    expect_equal(m, getDistanceMatrix(seq_1))
    m2 <- matrix(c(0, 1, 1, 0), 2, 2) %>% applyRowAndColumnNames
    expect_equal(m2, getDistanceMatrix(seq_2))
    m3 <- matrix(c(0, 4, 4, 0), 2, 2) %>% applyRowAndColumnNames
    expect_equal(m3, getDistanceMatrix(seq_3))
    m4 <- matrix(0, 2, 2) %>% applyRowAndColumnNames
    expect_equal(m4, getDistanceMatrix(seq_4))
})

test_that("test.getDistancesFromNaiveToMature", {
    naives <- c("AAAAAC", rep(c("AAAAAA"), 3))
    m1 <- c("AAAAAT")
    m2 <- c("CGCAAA")
    m3 <- c("GGGGGG")
    m4 <- c("AAAAAA")
    dat <- data.table(naive_seq=naives, mature_seq=c(m1, m2, m3, m4))
    expect_equal(c(0, 1, 3, 6), getDistancesFromNaiveToMature(dat))
})

test_that("test.getDistanceVector", {
    expect_equal(0, getDistanceVector(seq_1))
    expect_equal(1, getDistanceVector(seq_2))
    expect_equal(4, getDistanceVector(seq_3))
    expect_equal(0, getDistanceVector(seq_4))
    expect_equal(2, getDistanceVector(seq_7))
})

test_that("test.getGCDistribution", {
    d1 <- rep(0, 3)
    expect_equal(d1, getGCDistribution(seq_5))
    d2 <- c(0, 0.75, 1, 0.5)
    expect_equal(d2, getGCDistribution(seq_6))
})

test_that("test.getGRAVYDistribution", {
    s <- c("AAA", "GGG", "ACGTACGTACGT")
    expect_equal(c(-0.4, 0.8, 1.8), getGRAVYDistribution(s) %>% 
        sort)
})

test_that("test.getHotspotCount", {
    seq_a <- "TTTTT"
    seq_b <- list("AAC", "TTTT")
    seq_c <- "NWRC"
    seq_d <- c("WRC", "WRC", "WGC")
    expect_equal(0, getHotspotCount(seq_a))
    expect_equal(2, getHotspotCount(seq_b))
    expect_equal(3, getHotspotCount(seq_c))
    expect_equal(5, getHotspotCount(seq_d))
})

test_that("test.getNearestNeighborDistances", {
    d1 <- c(0, 0)
    expect_equal(d1, getNearestNeighborDistances(seq_1))
    d2 <- c(1, 1, 2)
    expect_equal(d2, getNearestNeighborDistances(seq_5) %>% sort)
    expect_equal(c(3, 3, 3, 3), getNearestNeighborDistances(seq_6, 
        k = 1))
    expect_equal(c(4, 3, 3, 3), getNearestNeighborDistances(seq_6, 
        k = 2))
    expect_equal(c(4, 4, 4, 3), getNearestNeighborDistances(seq_6, 
        k = 3))
})

test_that("test.getNucleotideDiversity", {
    d1 <- getNucleotideDiversity(list("AT", "AT"))
    d2 <- getNucleotideDiversity(list("AT", "AC"))
    d3 <- getNucleotideDiversity(list("AT", "GC"))
    d4 <- getNucleotideDiversity(list("ATATATATAT", "AAAAAAAAAA"))
    d5 <- getNucleotideDiversity(list("ATATATATAT", "CGTACGTAAT"))
    expect_equal(0, d1)
    expect_equal(0.5, d2)
    expect_equal(1, d3)
    expect_true(d4 < d5)
})

test_that("test.getPositionalDistancesBetweenMutations", {
    naive_seq <- c("AAAA", "AAAA", "TTTT", "GGGGGGG")
    mature_seq <- c("AAAA", "ATAG", "TTAT", "CGGCGCG")

    expect_equal(getDistancesBetweenMutations(naive_seq[1], mature_seq[1]), logical(0))
    expect_equal(getDistancesBetweenMutations(naive_seq[2], mature_seq[2]), 1)
    expect_equal(getDistancesBetweenMutations(naive_seq[3], mature_seq[3]), logical(0))
    expect_equal(getDistancesBetweenMutations(naive_seq[4], mature_seq[4]) %>% sort,
                 c(1, 2))
    expect_equal(getDistancesBetweenMutations(naive_seq, mature_seq) %>% sort, 
                 c(1, 1, 2))
})

test_that("getClusterSizes returns the correct cluster size distribution", {
    clones <- "1;2:3:4;5;6;7:8:9:10"
    dat <- data.frame(partition=clones)
    expect_equal(getClusterSizes(dat) %>% sort, c(1, 1, 1, 3, 4))

    expect_equal(getHillNumbers(dat, c(0, 2)), c(5, 3.571429), tolerance=0.001)
})
