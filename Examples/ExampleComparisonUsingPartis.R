# The path of the sumrep directory. This will need to be changed if the sumrep
#   directory is not the working directory.
path_to_sumrep <- getwd()

# Make sure sumrep is loaded
devtools::load_all(path_to_sumrep)

# Get annnotation and clonal partition dataset for IgH sequences
test_dat <- getPartisAnnotations(file.path(path_to_sumrep,
                                           "data/test_dat.fa"),
                                 locus="igh",
                                 cleanup=FALSE,
                                 output_path="tmp_output"
                                )

# Simulate a dataset based on the observed DNA sequences
test_simu <- getPartisSimulation(parameter_dir="tmp_output",
                                 num_events=nrow(test_dat[["annotations"]])
                                )

# Delete partis output folder
"tmp_output" %>% unlink(recursive=TRUE)

# Resample original fasta sequences and save to fasta file
test_dat_boot <- bootstrapFasta(file.path(path_to_sumrep,
                                          "data/compare_data.fa"),
                                output_filename="boot.fa"
                               )

# Get annotation and clonal partition dataset from resampled sequences
test_dat_boot <- getPartisAnnotations("boot.fa",
                                      locus="igh",
                                      output_path="tmp_output_boot"
                                     )

# Delete resampled fasta file
"boot.fa" %>% file.remove

# Run repertoire comparison for the observed and simulated data,
# and also compare the observed to bootstrapped data for refernce
comparison <- compareRepertoires(test_dat, 
                                 test_simu, 
                                 test_dat_boot,
                                 locus="igh"
                                )
