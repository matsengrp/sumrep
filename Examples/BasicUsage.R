# The path of the sumrep directory. This will need to be changed if the sumrep
# directory is not the working directory.
path_to_sumrep <- getwd()

# Make sure sumrep is loaded
devtools::load_all(path_to_sumrep)

# Read in an annotations dataset from a csv file
dat_a <- data.table::fread("data/test_annotations.csv")

# Compute the distribution of pairwise distances for dat_a.
# By default, this uses the sequence_alignment column
pairwise_distances <- getPairwiseDistanceDistribution(dat_a)

# Compute the distribution of pairwise distances of the junction column for dat_a
junction_pairwise_distances <- getPairwiseDistanceDistribution(dat_a, column="junction")

# Read in a second annotations dataset (this time from an RDS file)
dat_b <- readRDS("data/test_dat_boot.rds")

# Compute the JS-divergence of the pairwise distance distributions of dat_a and dat_b
pairwise_distance_divergence <- comparePairwiseDistanceDistributions(dat_a, dat_b)

# Run a full repertoire comparison of dat_a and dat_b
repertoire_a <- list(annotations=dat_a)
repertoire_b <- list(annotations=dat_b)
compareRepertoires(repertoire_a, repertoire_b, locus="igh")

# Plot all relevant univariate distributions of dat_a and dat_b
plot_1 <- plotUnivariateDistributions(list(dat_a, dat_b), 
                                      locus="igh"
                                     )

# Plot only the pairwise distance and aromaticity distributions of dat_a and dat_b
plot_2 <- plotUnivariateDistributions(list(dat_a, dat_b), 
                                      locus="igh", 
                                      plot_types="freqpoly", 
                                      plot_function_strings=c("getPairwiseDistanceDistribution", 
                                                              "getAromaticityDistribution"
                                                             )
                                     )


