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
                               dir_name="tmp",
                               num_scenarios=1,
                               species="human",
                               chain="beta",
                               cleanup=TRUE
                               ) {
    if(!(dir_name %in% list.files())) {
        dir.create(dir_name) 
    }

    input_seqs <- input_filename %>%
        data.table::fread()
    num_gen_sequences <- input_seqs %>%
        nrow
    runIgor(input_filename,
            file.path(getwd(), dir_name),
            num_scenarios=num_scenarios,
            num_gen_sequences=num_gen_sequences,
            species=species,
            chain=chain
           )
}
