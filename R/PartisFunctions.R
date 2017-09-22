#' Call partis
#'
#' \code{callPartis} calls partis from within the sumrep package.
#' This is done through the bash script run_partis.sh.
#' It is assumed that a "SHELL" environmental variable is set to run the above 
#' script.
#' @param action The desired partis command
#' @param input_filename The input (fasta) file, i.e. --infname argument
#' @param output_filename Desired output (csv) filename, i.e. --outfname 
#'   argument
#' @param output_path Desired output directory, which will contain 
#'   \code{output_filename}
#' @param partis_path The full path to the partis executable
#' @param num_procs The number of processors to use in parallel
#' @param cleanup Flag to delete all interim files created by partis
#' @return A data.table object containing the output of the partis call
callPartis <- function(action, input_filename, output_filename, output_path, 
                       partis_path, num_procs, cleanup) {
    shell <- Sys.getenv("SHELL")
    script.file <- system.file("run_partis.sh", package="sumrep")
    command <- paste(shell, script.file, 
                     "-p", partis_path, 
                     "-a", action, 
                     "-i", input_filename, 
                     "-o", output_filename, 
                     "-n", num_procs,
                     "-h", file.path(output_path, "params"))
    command %>% system
    partis_dataset <- output_filename %>% 
        data.table::fread(stringsAsFactors=TRUE)
    return(partis_dataset)
} 

#' Get per-gene and per-gene-per-position mutation information from partis
#'
#' \code{getMutationInfo} is run if \code{doFullAnnotation} is true in the 
#'   \link{annotateSequences} call. This function extracts mutation information
#'   from the given YAML file in the parameter directory output by partis.
#' @param filename YAML file from partis output folder corresponding to a gene
#' @return Object containing the gene's overall mutation rate as well as a 
#'   vector of positional mutation rates
getMutationInfo <- function(filename) {
    getSubstitutionRate <- function(state) {
        position <- state$name %>% strsplit("_") %>% unlist %>% last
        if(!grepl("[0-9]+", gsub("[^0-9]+", "", position))) {
            # The state in question does not represent an allele-position entry
            # (partis returns various metadata as 'states')
            result <- NA
        } else {
            # Partis returns zero-based position values, so add one
            position <- as.numeric(position) + 1
            germline_base <- state$extras$germline
            pos_emissions <- state$emissions
            highest_base_sub <- pos_emissions$probs %>% which.max %>% names

            # If the germline is not the most frequent base, the rate is 
            # unreliable, so treat as missing info
            result <- ifelse(highest_base_sub == germline_base,
                             1 - get(highest_base_sub, pos_emissions$probs),
                             NA)
        }
        names(result) <- position
        return(result)
    }

    result <- list()
    yaml_object <- yaml.load_file(filename)
    gene_info <- yaml_object$name %>% strsplit("_star_") %>% unlist
    result$gene <- gene_info[1]
    result$allele <- gene_info[2]
    states <- yaml_object$states[-(1:2)]

    overall_substitution_rate <- yaml_object$extras$per_gene_mute_freq
    position_substitution_rates <- states %>% sapply(getSubstitutionRate) %>% 
        subset(!is.na(.))
    substitution_rates <- {}
    substitution_rates$overall_mut_rate <- overall_substitution_rate
    substitution_rates$mut_rate_by_position <- position_substitution_rates
    return(substitution_rates)
}

#' Ensure sumrep does not accidentally delete any previous _output folders 
#' created by partis
#' @param output_path Desired partis output directory, which is \code{_output}
#'   by default
#' @param cleanup Flag to delete all interim files created by partis
preventOutputOverwrite <- function(output_path, cleanup) {
    if(length(list.files("_output")) > 0 && 
       output_path == "_output" && cleanup) {
        stop(paste("_output path already exists.",
                   "Please remove it, use another path name,",
                   "or set 'cleanup' to 'FALSE'."))
    }
}

#' Do extended partis processing of annotations
#'
#' @inheritParams annotateSequences
doFullAnnotation <- function(output_path, output_file, partis_path) {
    extended_output_filename <- "new_output.csv"
    extended_output_filepath <- file.path(output_path, extended_output_filename)
    script_file <- system.file("process_output.py", package="sumrep")
    script_command <-  paste("python", script_file, output_file, 
                             extended_output_filepath, 
                             partis_path %>% dirname %>% dirname,
                             sep=' ')
    script_command %>% system
    annotated_data <- extended_output_filepath %>% 
        data.table::fread(stringsAsFactors=TRUE) %>% 
        subset(select=which(!duplicated(names(.))))
    annotated_data$naive_seq <- annotated_data$naive_seq %>% 
        sapply(toString) %>% tolower
    extended_output_filepath %>% file.remove
    return(annotated_data)
}

#' Perform sequence annotation with partis
#'
#' @param input_filename The input (fasta) file, i.e. --infname argument
#' @param output_filename Desired output (csv) filename, i.e. --outfname 
#'   argument
#' @param partis_path The full path to the partis executable
#' @param num_procs The number of processors to use in parallel
#' @param cleanup Flag to delete all interim files created by partis
#' @param do_full_annotation Include per-gene and per-gene-per-position mutation
#'   rate information
#' @param output_path Desired output directory, which will contain 
#'   output_filename
#' @return A data.table object containing the output of the partis call
annotateSequences <- function(input_filename, output_filename="partis_output.csv", 
                              partis_path=Sys.getenv("PARTIS_PATH"), num_procs=4, 
                              cleanup=TRUE, 
                              do_full_annotation=TRUE, output_path="_output") {
    preventOutputOverwrite(output_path, cleanup)

    output_file <- file.path(output_path, output_filename)
    annotated_data <- callPartis("annotate", input_filename, output_file, 
                                 output_path, partis_path, num_procs, cleanup)


    hmm_yaml_filepath <- file.path(output_path, "params/hmm/hmms")
    yaml_files <- hmm_yaml_filepath %>% list.files
    yaml_filepath_and_files <- sapply(yaml_files, 
                                      function(x) { file.path(hmm_yaml_filepath,
                                                              x) }) 

    mutation_rates <- {}
    for(yaml_file in yaml_files) {
        allele <- yaml_file %>% gsub(pattern=".yaml", replace='') %>%
            gsub(pattern="_star_", replace="\\*")
        mutation_rates[[allele]] <- hmm_yaml_filepath %>% 
            file.path(yaml_file) %>% getMutationInfo
    }

    if(do_full_annotation) {
        annotated_data <- doFullAnnotation(output_path, output_file, partis_path)
    }

    if(cleanup) {
        output_path %>% unlink(recursive=TRUE)
    }

    raw_sequences <- input_filename %>% getSequenceListFromFasta
    annotated_data$mature_seq <- raw_sequences[annotated_data$unique_ids]

    annotation_object <- {}
    annotation_object$annotations <- annotated_data
    annotation_object$mutation_rates <- mutation_rates
    return(annotation_object)
}

#' Perform clonal partitioning with partis
#'
#' @param input_filename The input (fasta) file, i.e. --infname argument
#' @param output_filename Desired output (csv) filename, i.e. --outfname 
#'   argument
#' @param partis_path The full path to the partis executable
#' @param num_procs The number of processors to use in parallel
#' @param cleanup Flag to delete all interim files created by partis
#' @param output_path Desired output directory, which will contain 
#'   output_filename
#' @return A data.table object containing the output of the partis call
partitionSequences <- function(input_filename, 
                               output_filename="partis_output.csv", 
                               partis_path=Sys.getenv("PARTIS_PATH"), 
                               num_procs=4, 
                               cleanup=TRUE, output_path="_output") {
    preventOutputOverwrite(output_path, cleanup)

    output_file <- file.path(output_path, output_filename)
    partitioned_data <- callPartis("partition", input_filename, output_file, 
                                    output_path, partis_path, num_procs, 
                                    cleanup)

    if(cleanup) {
        output_path %>% unlink(recursive=TRUE)
    }

    return(partitioned_data)
}

annotateAndPartitionSequences <- function(input_filename, 
                                          output_filename="partis_output.csv", 
                                          partis_path=Sys.getenv("PARTIS_PATH"),
                                          num_procs=4, 
                                          cleanup=TRUE, 
                                          do_full_annotation=TRUE, 
                                          output_path="_output") {
    annotation_object <- annotateSequences(input_filename, output_filename, 
                                           partis_path, num_procs, cleanup, 
                                           do_full_annotation,
                                     output_path)
    partitions <- partitionSequences(input_filename, output_filename, 
                                     partis_path, num_procs, cleanup, 
                                     output_path)
    annotation_object$annotations <- 
        includeClonalMemberships(annotation_object$annotations, partitions)
    return(annotation_object)
}

simulateDataset <- function(parameter_dir,
                            partis_path=Sys.getenv("PARTIS_PATH"),
                            output_file="simu.csv",
                            num_events=100,
                            cleanup=TRUE,
                            do_full_annotation=TRUE) {
    partis_command <- paste(partis_path, "simulate", 
                            "--parameter-dir", parameter_dir,
                            "--outfname", output_file,
                            "--n-sim-events", num_events)
    partis_command %>% system
    sim_annotations <- output_file %>% fread(stringsAsFactors=FALSE)
    if(do_full_annotation) {
        sim_annotations <- doFullAnnotation(output_file,
                                            output_path=".",
                                            partis_path)
    }

    if(cleanup) {
        output_file %>% unlink
    }

    names(sim_annotations)[which(names(sim_annotations) == "input_seqs")] <- 
        "mature_seq"
    sim_annotations$mature_seq <- sim_annotations$mature_seq %>% 
        sapply(toString) %>% sapply(tolower)
    sim_data <- list(annotations=sim_annotations)
    return(sim_data)
}

