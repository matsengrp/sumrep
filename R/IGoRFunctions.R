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
    python_command <- paste("python3",
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

    annotations <- processIgorAnnotations(annotations,
                                         chain=chain)

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

    simulations <- fread(file.path(igor_wd_name,
                                   "sim.csv"
                                   )
    )

    simulations <- processIgorAnnotations(simulations,
                                          chain=chain)

    sim_indexed_seqs <- fread(file.path(igor_wd_name,
                                        "igor_generated",
                                        "generated_seqs_werr.csv"
                                       )
    )

    sim_aux_df <- fread(file.path(igor_wd_name,
                                  "igor_generated",
                                  "generated_seqs_werr_CDR3_info.csv"
                                 )
    )

    sim_annotations <- merge(sim_aux_df, simulations)
    sim_annotations$sequence <-
        sim_indexed_seqs[["nt_sequence"]][1 + sim_annotations[["seq_index"]]] %>%
        tolower

    names(sim_annotations)[which(names(sim_annotations) == "nt_CDR3")] <- 
        "junction"
    sim_annotations$junction <- sim_annotations$junction %>% tolower
    return(list(annotations=annotations,
                simulations=sim_annotations
               )
    )
}

processIgorAnnotations <- function(annotations,
                                   chain
                                  ) {
    if(chain == "beta") {
        annotations[["vd_insertion"]] <- annotations[["vd_insertion"]] %>%
            gsub(pattern="[^A-Z]", replace="") %>%
            tolower
        annotations[["dj_insertion"]] <- annotations[["dj_insertion"]] %>%
            gsub(pattern="[^A-Z]", replace="") %>%
            tolower
    } else {
        annotations[["vj_insertion"]] <- annotations[["vj_insertion"]] %>%
            gsub(pattern="[^A-Z]", replace="") %>%
            tolower
    }

    if("is_inframe" %in% names(annotations)) {
        annotations$vj_in_frame <- annotations$is_inframe %>%
            as.logical
    }
    return(annotations)
}
