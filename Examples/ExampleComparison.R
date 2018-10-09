# Get annnotation and clonal partition dataset
test_dat <- getPartisAnnotations("data/compare_data.fa", 
                                 cleanup=FALSE,
                                 output_path="tmp_output"
                                )

# Simulate a dataset based on the observed DNA sequences
test_simu <- simulateDataset(parameter_dir="tmp_output",
                             num_events=round(nrow(test_dat$annotations)/4)
                            )

# Delete partis output folder
"tmp_output" %>% unlink(recursive=TRUE)

# Resample original fasta sequences and save to fasta file
test_dat_boot <- bootstrapFasta("data/compare_data.fa",
                                output_filename="boot.fa"
                               )

# Get annotation and clonal partition dataset from resampled sequences
test_dat_boot <- getPartisAnnotations("boot.fa",
                                      output_path="tmp_output_boot"
                                     )

# Delete resampled fasta file
"boot.fa" %>% file.remove

# Run repertoire comparison for the observed and simulated data,
# and also compare the observed to bootstrapped data for refernce
comparison <- compareRepertoires(test_dat, test_simu, test_dat_boot)
