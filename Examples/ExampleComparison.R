# Get annnotation and clonal partition dataset
test_dat <- annotateAndPartitionSequences("data/compare_data.fa", 
                                          cleanup=FALSE)

# Simulate a dataset based on the observed DNA sequences
test_simu <- simulateDataset(parameter_dir="_output/params/")

# Delete partis output folder
"_output" %>% unlink(recursive=TRUE)

# Resample original fasta sequences and save to fasta file
test_dat_boot <- bootstrapFasta("data/compare_data.fa",
                                output_filename="boot.fa")

# Get annotation and clonal partition dataset from resampled sequences
test_dat_boot <- annotateAndPartitionSequences("boot.fa")

# Delete resampled fasta file
"boot.fa" %>% file.remove

# Run repertoire comparison for the observed and simulated data,
# and also compare the observed to bootstrapped data for refernce
compareRepertoires(test_dat, test_simu, test_dat_boot)
