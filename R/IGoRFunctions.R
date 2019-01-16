#' Run igor shell script from R
#'
#' @param dir_name Absolute path to desired output directory
runIgor <- function(input_filename,
                    dir_name,
                    num_scenarios=num_scenarios,
                    num_gen_sequences=num_gen_sequences,
                    eval_batch_name="eval_batch",
                    gen_batch_name="gen_batch",
                    species=species,
                    chain=chain
                    ) {
    igor_command <- paste("sh inst/run_igor.sh",
                          "-w", dir_name,
                          "-i", input_filename,
                          "-n", num_scenarios,
                          "-g", num_gen_sequences,
                          "-e", eval_batch_name,
                          "-b", gen_batch_name,
                          "-c", chain,
                          "-s", species
                          )
    cat(igor_command, '\n')
    igor_command %>% system
}

getIgorAnnotations <- function(input_filename,
                               igor_wd_name,
                               chain,
                               output_filename="annotations.csv"
                              ) {
    python_command <- paste("python",
                            "inst/run_igor.py",
                            input_filename,
                            igor_wd_name,
                            chain,
                            output_filename
                           )
    python_command %>% system
    annotations <- data.table::fread(file.path(igor_wd_name,
                                               output_filename
                                              )
    )

    if(chain == "beta") {
        annotations[["vd_insertion"]] <- annotations[["vd_insertion"]] %>%
            gsub(pattern="[^A-Z]", replacement="") %>%
            tolower
        annotations[["dj_insertion"]] <- annotations[["dj_insertion"]] %>%
            gsub(pattern="[^A-Z]", replacement="") %>%
            tolower
    } else {
        annotations[["vj_insertion"]] <- annotations[["vj_insertion"]] %>%
            gsub(pattern="[^A-Z]", replacement="") %>%
            tolower
    }

    indexed_seqs <- fread(file.path(igor_wd_name,
                                    "aligns",
                                    "igor_indexed_sequences.csv")
    )

    # Get query sequences from an igor alignment file. Add one to indices 
    # since annotations$seq_index contains zero-based positions
    annotations$sequence <- 
        indexed_seqs[["sequence"]][1 + annotations[["seq_index"]]] %>%
        tolower
    annotations$sequence_alignment <- annotations$sequence
    return(annotations)
}
