# The path of the sumrep directory. This will need to be changed if the sumrep
#   directory is not the working directory.
path_to_sumrep <- getwd()

# Make sure sumrep is loaded
devtools::load_all(path_to_sumrep)

# Load annnotation and clonal partition dataset
test_dat <- readRDS(file.path(path_to_sumrep,
                              "data-raw/test_dat.rds")) 

# Load a simulated dataset based on the observed DNA sequences using partis
test_simu <- readRDS(file.path(path_to_sumrep, 
                               "data-raw/test_simu.rds"))

# Run repertoire comparison for the observed and simulated data,
# and also compare the observed to bootstrapped data for refernce
comparison <- compareRepertoires(test_dat, 
                                 test_simu, 
                                 locus="igh")
