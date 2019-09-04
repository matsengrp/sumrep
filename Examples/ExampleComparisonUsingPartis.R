# The path of the sumrep directory. This will need to be changed if the sumrep
#   directory is not the working directory.
path_to_sumrep <- getwd()

# Make sure sumrep is loaded
devtools::load_all(path_to_sumrep)

# Get annnotation and clonal partition dataset for IgH sequences
test_dat <- getPartisAnnotations(file.path(path_to_sumrep,
                                           "data/test_data.fa"),
                                 locus="igh",
                                 cleanup=FALSE,
                                 output_path="tmp_output"
                                )

# Simulate a dataset based on the observed DNA sequences
test_simu <- getPartisSimulation(parameter_dir="tmp_output",
                                 num_events=nrow(test_dat[["annotations"]]),
                                 extra_command_args="--min-observations-per-gene 15"
                                 # ^ partis has troubles simulating from small
                                 # datasets, so give it some wiggle room
                                )

# Delete partis output folder
"tmp_output" %>% unlink(recursive=TRUE)

# Run repertoire comparison for the observed and simulated data,
# and also compare the observed to bootstrapped data for refernce
comparison <- compareRepertoires(test_dat, 
                                 test_simu, 
                                 locus="igh"
                                )
