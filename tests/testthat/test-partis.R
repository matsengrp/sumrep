library(sumrep)
context("Partis annotations and partitions")

test_data_path <- system.file("data/test_data.fa", package="sumrep")

test_that("sumrep correctly calls partis annotate", {
    dat_a <- annotateSequences(test_data_path, 
                               do_full_annotation=FALSE)
    dat_b <- annotateSequences(test_data_path, 
                               output_filename="blah.csv",
                               output_path="_output_arbitrary", 
                               num_procs=8)

    expect_equal(dat_a %>% names, c("annotations", "mutation_rates"))
    expect_equal(dat_a$mutation_rates[[1]] %>% names, 
                 c("overall_mut_rate", "mut_rate_by_position"))

    expect_equal(ncol(dat_a$annotations), 36)
    expect_equal(nrow(dat_a$annotations), 17)
    expect_equal(ncol(dat_b$annotations), 45)
    expect_equal(nrow(dat_b$annotations), 17)

    dir.create("_output")
    file.create("_output/blah.txt")
    expect_error(annotateSequences(test_data_path))
    "_output" %>% unlink(recursive=TRUE)
})

test_that("sumrep correctly calls partis partition", {
    dat <- partitionSequences(test_data_path, 
                              output_filename="blah2.csv",
                              output_path="_output_arbitrary_again", 
                              num_procs=8)

    expect_equal(dat %>% names, c("logprob", "n_clusters", "n_procs", 
                                  "partition"))

    dir.create("_output")
    file.create("_output/blah2.txt")
    expect_error(partitionSequences(test_data_path))
    "_output" %>% unlink(recursive=TRUE)
})
