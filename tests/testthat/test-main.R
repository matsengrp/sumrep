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

test_that("test.compareDistanceFromGermlineToSequenceDistributions", {
    naive_a <- c("AAAAAA")
    naive_b <- c("GGGGGG")
    m1 <- c("AAAAAA")
    m2 <- c("AAAAAG")
    m3 <- c("GAAAAA")
    m4 <- c("AAAGGG")
    m5 <- c("GGGG")
    dat_a <- data.table(germline_alignment=naive_a, sequence_alignment=list(m1, m2, m3))
    dat_b <- data.table(germline_alignment=naive_b, sequence_alignment=list(m2, m3, m4))
    c1 <- compareDistanceFromGermlineToSequenceDistributions(dat_a, 
                                            dat_a, 
                                            approximate=FALSE,
                                            v_gene_only=FALSE
                                            )
    c2 <- compareDistanceFromGermlineToSequenceDistributions(dat_a, 
                                            dat_b, 
                                            approximate=FALSE,
                                            v_gene_only=FALSE
                                            )
    c3 <- compareDistanceFromGermlineToSequenceDistributions(dat_b, 
                                            dat_a, 
                                            approximate=FALSE,
                                            v_gene_only=FALSE)
    expect_equal(0, c1)
    expect_equal(c3, c2)
    expect_true(c2 > 0)
})

test_that("test.compareGCContentDistributions", {
    dat_a <- data.table(sequence=seq_2)
    dat_b <- data.table(sequence=seq_7)
    dat_c <- data.table(sequence=seq_8)
    dat_d <- data.table(sequence=seq_9)
    expect_equal(0, compareGCContentDistributions(dat_a, dat_b))
    c1 <- compareGCContentDistributions(dat_c, dat_d)
    c2 <- compareGCContentDistributions(dat_d, dat_c)
    expect_true(c1 > 0)
    expect_true(c1 == c2)
})

test_that("test.comparePairwiseDistanceDistributions", {
    s1 <- data.table(sequence=c("AAA", "AAT", "ATT"))
    s2 <- data.table(sequence=c("ATT", "AAA", "AAT"))
    s3 <- data.table(sequence=c("AAA", "AAT", "TTT"))
    s4 <- data.table(sequence=c("AAA", "ATC", "GGG"))
    expect_equal(0, 
                 comparePairwiseDistanceDistributions(s1,
                                                      s2,
                                                      approximate=FALSE
                                                      ))
    c1 <- comparePairwiseDistanceDistributions(s1, 
                                               s3, 
                                               approximate=FALSE)
    c2 <- comparePairwiseDistanceDistributions(s1, 
                                               s4, 
                                               approximate=FALSE)
    expect_true(c1 < c2)
    c3 <- comparePairwiseDistanceDistributions(s4, 
                                               s1, 
                                               approximate=FALSE)
    expect_equal(c3, c2)
})

test_that("test.getColdspotCountDistribution", {
    seq_a <- "TTTTT"
    seq_b <- c("GCC", "AAA", "AAAA")
    seq_c <- c("AAA", "CTC")
    seq_d <- "SYCSYC"
    expect_equal(0, getColdspotCountDistribution(data.table(junction=seq_a), 
                                                 column="junction"))
    expect_equal(c(1, 0, 0), getColdspotCountDistribution(data.table(junction=seq_b), 
                                                 column="junction"))
    expect_equal(c(0, 1), getColdspotCountDistribution(data.table(junction=seq_c), 
                                                 column="junction"))
    expect_equal(4, getColdspotCountDistribution(data.table(junction=seq_d), 
                                                 column="junction"))
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

test_that("test.getDistanceFromGermlineToSequenceDistribution", {
    naives <- c("AAAAAC", rep(c("AAAAAA"), 3))
    m1 <- c("AAAAAT")
    m2 <- c("CGCAAA")
    m3 <- c("GGGGGG")
    m4 <- c("AAAAAA")
    dat <- data.table(germline_alignment=naives, sequence_alignment=c(m1, m2, m3, m4))
    expect_equal(c(0, 1, 3, 6), getDistanceFromGermlineToSequenceDistribution(dat,
                                                              v_gene_only=FALSE
                                                              ))
    dat$v_gl_seq <- c("AAA", "CGC", "GGG", "AAA")
    dat$v_qr_seqs <- rep("AAA", 4)
    expect_equal(c(0, 0, 3, 3), 
                 getDistanceFromGermlineToSequenceDistribution(dat,
                                                          v_gene_only=TRUE
                                                         )
                 )
})

test_that("test.getDistanceVector", {
    expect_equal(0, getDistanceVector(seq_1))
    expect_equal(1, getDistanceVector(seq_2))
    expect_equal(4, getDistanceVector(seq_3))
    expect_equal(0, getDistanceVector(seq_4))
    expect_equal(2, getDistanceVector(seq_7))
})

test_that("CDR3 pairwise distances functions", {
    dat <- data.table(junction=c("AAACCC", "AAACCC", "ACTCAT"), 
                        junction_aa=c("KP", "KP", "TH"))
    dist_1 <- dat %>% 
        getCDR3PairwiseDistanceDistribution(
            approximate=FALSE,
            by_amino_acid=FALSE
        )
    expect_equal(c(0, 4, 4), dist_1)
    dist_2 <- dat %>% 
        getCDR3PairwiseDistanceDistribution(
            approximate=FALSE,
            by_amino_acid=TRUE
        )
    expect_equal(c(0, 2, 2), dist_2)
})

test_that("test.getGCContentDistribution", {
    d1 <- rep(0, 3)
    expect_equal(d1, getGCContentDistribution(dat_1, approximate=FALSE))
    d2 <- c(0, 0.75, 1, 0.5)
    expect_equal(d2, getGCContentDistribution(dat_2, approximate=FALSE))
})

test_that("test.getGRAVYDistribution", {
    # Sample GRAVY values for a few particular amino acids obtained from 
    # https://github.com/PRIDE-Utilities/pride-utilities/wiki/1.2-GRAVY-Calculations 
    alanine <- data.table(junction_aa=convertNucleobasesToAminoAcids("GCC"))
    glutamine <- data.table(junction_aa=convertNucleobasesToAminoAcids("CAA"))
    isoleucine <- data.table(junction_aa=convertNucleobasesToAminoAcids("ATT"))
    expect_equal(1.8, getGRAVYDistribution(alanine))
    expect_equal(-3.5, getGRAVYDistribution(glutamine))
    expect_equal(4.5, getGRAVYDistribution(isoleucine))
    expect_equal(c(1.8, -3.5, 4.5), 
                 getGRAVYDistribution(rbind(alanine, glutamine, isoleucine)))
})

test_that("test.getHotspotCountDistribution", {
    seq_a <- "TTTTT"
    seq_b <- c("AAC", "TTTT")
    seq_c <- "NWRC"
    seq_d <- c("WRC", "WRC", "WGC")
    expect_equal(0, getHotspotCountDistribution(data.table(junction=seq_a), 
                                                column="junction"))
    expect_equal(c(2, 0), getHotspotCountDistribution(data.table(junction=seq_b), 
                                                column="junction"))
    expect_equal(3, getHotspotCountDistribution(data.table(junction=seq_c), 
                                                column="junction"))
    expect_equal(c(2, 2, 1), getHotspotCountDistribution(data.table(junction=seq_d), 
                                                column="junction"))
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

test_that("test.getPositionalPositionalDistanceBetweenMutations", {
    germline_alignment <- c("AAAA", "AAAA", "TTTT", "GGGGGGG")
    sequence_alignment <- c("AAAA", "ATAG", "TTAT", "CGGCGCG")

    expect_equal(getPositionalDistancesBetweenMutationsBySequence(germline_alignment[1], sequence_alignment[1]), NA)
    expect_equal(getPositionalDistancesBetweenMutationsBySequence(germline_alignment[2], sequence_alignment[2]), 2)
    expect_equal(getPositionalDistancesBetweenMutationsBySequence(germline_alignment[3], sequence_alignment[3]), NA)
    expect_equal(getPositionalDistancesBetweenMutationsBySequence(germline_alignment[4], sequence_alignment[4]) %>% 
                 sort,
                 c(2, 3))
    expect_equal(getPositionalDistanceBetweenMutationsDistribution(
                     data.table(germline_alignment=germline_alignment, 
                                sequence_alignment=sequence_alignment
                               )
                 ) %>% sort, 
                 c(2, 2, 3))
})

test_that("getClusterSizes returns the correct cluster size distribution", {
    dat <- data.frame(clone=c(13, 14, 15, 13, 13, 14, 14, 13, 16, 17))
    cluster_sizes <- dat %>% getClusterSizes %>% sort
    expect_equal(c(1, 1, 1, 3, 4), cluster_sizes)

    expect_equal(getHillNumbers(dat, c(0, 2)), c(5, 3.571429), tolerance=0.001)
})

test_that("getInsertionMatrix functions return correct transition matrices", {
    dat <- data.frame(vd_insertion=c("aa", "CC", "xttx", "gga", "xyz"),
                      dj_insertion=c("aAa", "cc", "hello", "yo", "sup"))

    vd_truth <- matrix(0, 4, 4)
    rownames(vd_truth) <- c('a', 'c', 'g', 't')
    colnames(vd_truth) <- c('a', 'c', 'g', 't')
    vd_truth['a', 'a'] <- 0.2
    vd_truth['c', 'c'] <- 0.2
    vd_truth['t', 't'] <- 0.2
    vd_truth['g', 'g'] <- 0.2
    vd_truth['g', 'a'] <- 0.2
    vd_mat <- dat %>% getVDInsertionMatrix
    expect_equal(vd_truth, vd_mat)

    dj_truth <- matrix(0, 4, 4)
    rownames(dj_truth) <- c('a', 'c', 'g', 't')
    colnames(dj_truth) <- c('a', 'c', 'g', 't')
    dj_truth['a', 'a'] <- 2/3
    dj_truth['c', 'c'] <- 1/3
    dj_mat <- dat %>% getDJInsertionMatrix
    expect_equal(dj_truth, dj_mat)
})

test_that("getInFramePercentage returns the correct percentage of in-frame sequences", {
    dat_a <- data.frame(vj_in_frame=c(TRUE, FALSE, TRUE, FALSE))
    dat_b <- data.frame(vj_in_frame=c(TRUE, FALSE, TRUE, TRUE))
    expect_equal(getInFramePercentage(dat_a), 50)
    expect_equal(getInFramePercentage(dat_b), 75)
    expect_equal(compareInFramePercentages(dat_a, dat_b), 25)
})

test_that("Gene usage comparison", {
    dat_a <- data.frame(v_call=c("IGHV3-30*18",
                                 "IGHV4-31*03",
                                 "IGHV1-3*01"),
                        d_call=c("IGHD6-19*01",
                                 "IGHD2-15*01",
                                 "IGHD1-1*01"),
                        j_call=c("IGHJ6*02",
                                 "IGHJ5*02",
                                 "IGHJ4*02")                        
                        )
    dat_b <- data.frame(v_call=c("IGHV3-30*18",
                                 "IGHV2-01*03",
                                 "IGHV1-3*02"),
                        d_call=c("IGHD6-19*02",
                                 "IGHD2-11*01",
                                 "IGHD1-1*02"),
                        j_call=c("IGHJ6*02",
                                 "IGHJ5*02",
                                 "IGHJ4*02")                        
                        )
    expect_equal(4, compareVGeneDistributions(dat_a, 
                                              dat_b,
                                              collapse_alleles=FALSE,
                                              standardize=FALSE
                                             )
    )
    expect_equal(4/3, compareVGeneDistributions(dat_a, 
                                              dat_b,
                                              collapse_alleles=FALSE,
                                              standardize=TRUE
                                             )
    )
    expect_equal(2, compareVGeneDistributions(dat_a, 
                                              dat_b,
                                              collapse_alleles=TRUE,
                                              standardize=FALSE
                                             )
    )
    expect_equal(2/3, compareVGeneDistributions(dat_a, 
                                              dat_b,
                                              collapse_alleles=TRUE,
                                              standardize=TRUE
                                             )
    )
    expect_equal(2/3, compareVGeneDistributions(dat_a, 
                                              dat_b,
                                              collapse_alleles=TRUE
                                             )
    )
    expect_equal(6/3, compareDGeneDistributions(dat_a, 
                                              dat_b,
                                              collapse_alleles=FALSE
                                             )
    )
    expect_equal(2/3, compareDGeneDistributions(dat_a, 
                                              dat_b,
                                              collapse_alleles=TRUE
                                             )
    )
    expect_equal(0, compareJGeneDistributions(dat_a, 
                                              dat_b,
                                              collapse_alleles=FALSE
                                             )
    )
    expect_equal(6, compareVDJDistributions(dat_a, 
                                            dat_b,
                                            collapse_alleles=FALSE,
                                            by_frequency=FALSE
                                           )
    )
    expect_equal(2, compareVDJDistributions(dat_a, 
                                            dat_b,
                                            collapse_alleles=FALSE,
                                            by_frequency=TRUE
                                           )
    )
    expect_equal(2/3, compareVDJDistributions(dat_a, 
                                            dat_b,
                                            collapse_alleles=TRUE,
                                            by_frequency=TRUE
                                           )
    )
    expect_equal(2, compareVDJDistributions(dat_a, 
                                            dat_b,
                                            collapse_alleles=TRUE,
                                            by_frequency=FALSE
                                           )
    )
    expect_equal(2/3, compareVJDistributions(dat_a,
                                           dat_b,
                                           collapse_alleles=TRUE,
                                           by_frequency=TRUE
                                          )
   ) 
    expect_equal(4/3, compareVJDistributions(dat_a,
                                           dat_b,
                                           collapse_alleles=FALSE,
                                           by_frequency=TRUE
                                          )
   ) 
    expect_equal(4, compareVJDistributions(dat_a,
                                           dat_b,
                                           collapse_alleles=FALSE,
                                           by_frequency=FALSE
                                          )
   ) 
})
