#' Call partis
#'
#' \code{callPartis} calls partis from within the sumrep package.
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
#' @param locus String denoting the locus for the receptors given in 
#'   \code{input_filename}. 
#'   Either "tra", "trb", "trd", or "trg" for TCRs,
#'   or "igl", "igk", or "igh" for BCRs
#' @param cleanup Flag to delete all interim files created by partis
#' @return A data.table object containing the output of the partis call
#' 
#' @export
callPartis <- function(action, 
                       input_filename, 
                       output_filename, 
                       output_path, 
                       partis_path, 
                       num_procs, 
                       locus,
                       germline_dir,
                       extra_columns,
                       seed=NULL
                       ) {
    command <- paste(partis_path,
                     action,
                     "--infname",
                     input_filename,
                     "--outfname",
                     output_filename,
                     "--n-procs",
                     num_procs,
                     "--locus",
                     locus,
                     "--parameter-dir",
                     file.path(output_path, "params")
                    )
    
    if(!missing(germline_dir) && !is.null(germline_dir)) {
        command <- paste(command,
                         "--initial-germline-dir",
                         germline_dir
                        )
    }

    if(!missing(extra_columns) && !is.null(extra_columns)) {
        command <- paste(command,
                         "--extra-annotation-columns",
                         extra_columns
                        )
    }

    if(!is.null(seed)) {
        command <- paste(command,
                         "--seed",
                         seed
                        )
    }
    
    if(action == "partition") {
        command <- paste(command,
                         "--count-parameters"
                        )
    }

    command %>% cat
    command %>% 
        system

    # Manually add query sequences to dataset from fasta file.
    # partis by default alters these, and it's not straightforward to
    # get both the raw and altered sequences (we need both).
    annotation_filename <- output_filename %>%
        ifelse(action == "partition",
               gsub(., 
                    pattern='.csv',
                    replacement='-cluster-annotations.csv'
                   ),
               .
              )
    partis_dataset <- appendQuerySequencesToPartisAnnotationsFile(
        input_filename,
        annotation_filename
    )

    return(partis_dataset)
} 

appendQuerySequencesToPartisAnnotationsFile <- function(input_filename,
                                                        annotation_filename
                                                       ) {
    partis_dataset <- annotation_filename %>% 
        data.table::fread(stringsAsFactors=TRUE)
    query_seqs <- input_filename %>%
        read.fasta(as.string=TRUE, seqonly=TRUE) %>%
        unlist %>%
        tolower
    partis_dataset[["sequence"]] <- query_seqs[partis_dataset[["unique_ids"]]]
    partis_dataset %>%
        data.table::fwrite(file=annotation_filename)
    return(partis_dataset)
}


#' Get per-gene and per-gene-per-position mutation information from partis
#'
#' \code{getMutationInfo} is run if \code{getFullPartisAnnotation} is true in the 
#'   \link{getPartisAnnotations} call. This function extracts mutation information
#'   from the given YAML file in the parameter directory output by partis.
#' @param filename YAML file from partis output folder corresponding to a gene
#' @return Object containing the gene's overall mutation rate as well as a 
#'   vector of positional mutation rates
#'   
#' @export
getMutationInfo <- function(filename) {
    getSubstitutionRate <- function(state) {
        position <- state[["name"]] %>%
            strsplit("_") %>% 
            unlist %>% 
            last
        if(!grepl("[0-9]+", gsub("[^0-9]+", "", position))) {
            # The state in question does not represent an allele-position entry
            # (partis returns various metadata as 'states')
            result <- NA
        } else {
            # Partis returns zero-based position values, so add one
            position <- as.numeric(position) + 1
            germline_base <- state[["extras"]][["germline"]]
            pos_emissions <- state[["emissions"]]
            highest_base_sub <- pos_emissions[["probs"]] %>%
                which.max %>% 
                names

            # If the germline is not the most frequent base, the rate is 
            # unreliable, so treat as missing info
            result <- ifelse(highest_base_sub == germline_base,
                             1 - get(highest_base_sub, pos_emissions[["probs"]]),
                             NA)
        }
        names(result) <- position
        return(result)
    }

    result <- list()
    yaml_object <- yaml.load_file(filename)
    gene_info <- yaml_object[["name"]] %>%
        strsplit("_star_") %>% 
        unlist
    result[["gene"]] <- gene_info[1]
    result[["allele"]] <- gene_info[2]
    states <- yaml_object[["states"]][-(1:2)]

    overall_substitution_rate <- yaml_object[["extras"]][["per_gene_mute_freq"]]
    position_substitution_rates <- states %>% 
        sapply(getSubstitutionRate) %>% 
        subset(!is.na(.))
    substitution_rates <- list()
    substitution_rates[["overall_mut_rate"]] <- overall_substitution_rate
    substitution_rates[["mut_rate_by_position"]] <- position_substitution_rates
    return(substitution_rates)
}

#' Ensure sumrep does not accidentally delete any previous _output folders 
#' created by partis
#'
#' @inheritParams callPartis
#' 
#' @export
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
#' @inheritParams getPartisAnnotations
#' 
#' @export
getFullPartisAnnotation <- function(output_path, 
                                    output_file, 
                                    partis_path
                                   ) {
    extended_output_filename <- "new_output.csv"
    extended_output_filepath <- file.path(output_path, extended_output_filename)
    script_file <- system.file("exec", "process_output.py", package="sumrep")
    script_command <-  paste("python", script_file, output_file, 
                             extended_output_filepath, 
                             partis_path %>% 
                                 dirname %>% 
                                 dirname,
                             sep=' ')
    script_command %>% system
    annotated_data <- extended_output_filepath %>% 
        data.table::fread(stringsAsFactors=TRUE) %>% 
        subset(select=which(!duplicated(names(.))))
    extended_output_filepath %>% 
        file.remove
    return(annotated_data)
}

# partis' output contains an indel_reversed_seqs column which best corresponds
# to the sequence_alignment field in the AIRR schema. However, this column is 
# the empty string if no indels were detected. Thus, we must manually construct 
# the sequence_alignment column to contain the indel_reversed_seqs string if 
# there was an indel, and the input_seqs string if there was not an indel.
processPartisMatureSequences <- function(dat) {
    indel_seqs <- which(dat[["indel_reversed_seqs"]] != "")
    dat[["sequence_alignment"]] <- dat[["input_seqs"]] %>%
        sapply(toString) %>%
        tolower
    dat[indel_seqs, ][["sequence_alignment"]] <- 
        dat[["indel_reversed_seqs"]][indel_seqs] %>%
        tolower
    return(dat)
}

# Here we use stringr::str_split to get empty strings after a ':', if necessary.
# E.g., with strsplit, '::' would yield c('', ''), whereas
# stringr::str_split would yield c('', '', '')
collapseColonedList <- function(coloned_list,
                                type_conversion=as.numeric
                               ) {
    collapsed_vector <- coloned_list %>% 
        toString %>%
        stringr::str_split(':') %>%
        unlist %>%
        unname

    return(collapsed_vector)
}

collapseClones <- function(partition_dataset) {
    coloned_columns <- c("unique_ids",
                         "mut_freqs",
                         "n_mutations",
                         "sequence",
                         "input_seqs",
                         "indel_reversed_seqs",
                         "mutated_invariants",
                         "in_frames",
                         "stops",
                         "aligned_j_seqs",
                         "aligned_v_seqs",
                         "v_qr_seqs",
                         "cdr3_seqs"
                         )
    partition_dataset[["clone_id"]] <- 0
    all_columns <- partition_dataset %>% names

    clone_count <- partition_dataset %>% nrow
    column_count <- partition_dataset %>% ncol
    collapsed_dat <- matrix(nrow=0, ncol=column_count) %>%
        data.table::data.table() %>%
        setNames(all_columns)

    for(clone in 1:clone_count) {
        clone_annotations <- partition_dataset[clone, ]
        unique_ids <- clone_annotations[["unique_ids"]] %>% 
            collapseColonedList

        if(length(unique_ids) > 1) {
            clone_list <- {} 
            for(column in all_columns) {
                if(column %in% coloned_columns) {
                    clone_list[[column]] <- clone_annotations[[column]] %>%
                        collapseColonedList
                } else {
                    clone_list[[column]] <- clone_annotations[[column]]
                }
            }

            clone_dat <- clone_list %>% as.data.table
            clone_dat[["clone_id"]] <- clone
        }  else {
            clone_dat <- clone_annotations
            clone_dat[["clone_id"]] <- clone
        }
        collapsed_dat <- rbind(collapsed_dat, clone_dat)
    }

    return(collapsed_dat)
}

#' Get partis annotations from the partis output directory.
#' Some column names are changed to match the AIRR standard, including:
#' This function does NOT run partis from R. 
#'
#' @param filepath The output directory from the partis call
#' @param output_filename The --infname argument to partis
#' @param partis_path The path to the partis binary
#' @param do_full_annotation Does further processing than the default by partis
#' @param collapse_clones Convert a row of colon-separated clonal families to
#'   separate rows for each member
#'   
#' @export
readPartisAnnotations <- function(output_path,
                                  annotation_filename,
                                  partis_path=Sys.getenv("PARTIS_PATH"),
                                  do_full_annotation=FALSE,
                                  collapse_clones=TRUE
                                ) {
    annotation_file <- file.path(output_path,
                                 annotation_filename)

    if(do_full_annotation) {
        annotated_data <- getFullPartisAnnotation(output_path, 
                                                  annotation_file, 
                                                  partis_path
                                                 )
        # annotated_data[["junction"]] <- annotated_data %>% getCDR3s
    } else {
        annotated_data <- annotation_file %>%
            data.table::fread()
    }

    if(collapse_clones) {
        annotated_data <- annotated_data %>%
            collapseClones

        collapsed_filename <- annotation_filename %>%
            gsub(pattern='.csv',
                 replacement='-collapsed.csv')
        collapsed_file=file.path(output_path,
                                 collapsed_filename)
        write.csv(annotated_data,
                  file=collapsed_file)
    }

    hmm_yaml_filepath <- file.path(output_path, "params/hmm/hmms")
    yaml_files <- hmm_yaml_filepath %>% 
        list.files
    yaml_filepath_and_files <- sapply(yaml_files, 
                                      function(x) { 
                                          file.path(hmm_yaml_filepath,
                                                    x)
                                      }
                               ) 

    mutation_rates <- {}
    for(yaml_file in yaml_files) {
        allele <- yaml_file %>% gsub(pattern=".yaml", replacement='') %>%
            gsub(pattern="_star_", replacement="\\*")
        mutation_rates[[allele]] <- hmm_yaml_filepath %>% 
            file.path(yaml_file) %>% 
            getMutationInfo
    }

    annotated_data <- annotated_data %>%
        processPartisSequences

    annotation_object <- {}
    annotation_object[["annotations"]] <- annotated_data
    annotation_object[["mutation_rates"]] <- mutation_rates
    return(annotation_object)
}

processPartisSequences <- function(annotated_data) {
    if("sequence" %in% names(annotated_data)) {
        annotated_data[["sequence"]] <- annotated_data[["sequence"]] %>%
            sapply(toString)
    }

    annotated_data[["germline_alignment"]] <- 
        annotated_data[["naive_seq"]] %>% 
        sapply(toString) %>% 
        tolower

    annotated_data[["v_qr_seqs"]] <- annotated_data[["v_qr_seqs"]] %>%
        sapply(toString) %>% 
        tolower

    annotated_data[["v_gl_seq"]] <- annotated_data[["v_gl_seq"]] %>%
        sapply(toString) %>% 
        tolower

    annotated_data[["junction"]] <- annotated_data[["cdr3_seqs"]] %>%
        sapply(toString) %>%
        tolower

    annotated_data[["junction_aa"]] <- annotated_data[["junction"]] %>%
        convertNucleobasesToAminoAcids

    annotated_data <- annotated_data %>% 
        processPartisMatureSequences

    annotated_data[["vj_in_frame"]] <- annotated_data[["in_frames"]] %>% 
        as.logical

    annotated_data[["np1_length"]] <- annotated_data[["vd_insertion"]] %>%
        sapply(nchar)
    annotated_data[["np2_length"]] <- annotated_data[["dj_insertion"]] %>%
        sapply(nchar)

    names(annotated_data)[which(names(annotated_data) == "cdr3_length")] <- 
        "junction_length"
    names(annotated_data)[which(names(annotated_data) == "v_gene")] <- "v_call"
    names(annotated_data)[which(names(annotated_data) == "d_gene")] <- "d_call"
    names(annotated_data)[which(names(annotated_data) == "j_gene")] <- "j_call"
    names(annotated_data)[which(names(annotated_data) == "stops")] <- "stop_codon" 
    annotated_data[["stop_codon"]] <- annotated_data[["stop_codon"]] %>%
        sapply(as.logical)

    return(annotated_data)
}

checkForValidLocus <- function(locus) {
    valid_loci <- c("tra", "trb", "trd", "trg", "igl", "igk", "igh")
    if(!(locus %in% valid_loci)) {
        stop(paste("locus must be one of",
                   paste(valid_loci %>%
                            sapply(function(x) { paste0("'", x, "'") }) %>%
                            unname, 
                         collapse=", ")
                  )
        )
    }
}

#' Perform sequence annotation with partis.
#' This function does call partis.
#'
#' @inheritParams callPartis
#' @return A data.table object containing the output of the partis annotate call
#' 
#' @export
getPartisAnnotations <- function(input_filename, 
                                 output_filename="partis_output.csv", 
                                 partis_path=Sys.getenv("PARTIS_PATH"), 
                                 num_procs=16, 
                                 locus,
                                 partition=ifelse(
                                     stringr::str_sub(locus, 0, 2) == "ig",
                                     TRUE,
                                     FALSE
                                 ),
                                 collapse_clones=TRUE,
                                 cleanup=TRUE, 
                                 do_full_annotation=FALSE, 
                                 output_path="_output",
                                 germline_dir=NULL,
                                 extra_columns="v_gl_seq:v_qr_seqs:cdr3_seqs"
                                ) {
    locus <- locus %>% tolower
    locus %>% checkForValidLocus

    preventOutputOverwrite(output_path, cleanup)

    if(partition) {
        partition_data <- getPartisPartitions(input_filename=input_filename,
                                              output_filename=output_filename,
                                              partis_path=partis_path,
                                              num_procs=num_procs,
                                              locus=locus,
                                              cleanup=FALSE,
                                              output_path=output_path,
                                              germline_dir=germline_dir,
                                              extra_columns=extra_columns
                                             )
        annotation_filename <- output_filename %>%
            gsub(pattern='.csv',
                 replacement='-cluster-annotations.csv')
    } else {
        output_file <- file.path(output_path, output_filename)
        annotated_data <- callPartis(action="annotate", 
                                     input_filename=input_filename, 
                                     output_file=output_file, 
                                     output_path=output_path, 
                                     partis_path=partis_path, 
                                     num_procs=num_procs,
                                     locus=locus,
                                     extra_columns=extra_columns
                                     )
        annotation_filename <- output_filename
    }

    annotation_object <- readPartisAnnotations(output_path,
                                               annotation_filename=annotation_filename,
                                               partis_path=partis_path,
                                               do_full_annotation=do_full_annotation,
                                               collapse_clones=collapse_clones
                                              )

    if(cleanup) {
        output_path %>% 
            unlink(recursive=TRUE)
    }

    return(annotation_object)
}

#' Perform clonal partitioning with partis
#'
#' @inheritParams callPartis
#' @return A data.table object containing the output of the partis partition 
#'   call
#'   
#' @export
getPartisPartitions <- function(input_filename, 
                                output_filename="partition.csv", 
                                partis_path=Sys.getenv("PARTIS_PATH"), 
                                num_procs=4, 
                                locus,
                                cleanup=TRUE, 
                                output_path="_output",
                                germline_dir=NULL,
                                extra_columns
                               ) {
    preventOutputOverwrite(output_path, cleanup)

    output_file <- file.path(output_path, output_filename)
    partitioned_data <- callPartis("partition", 
                                   input_filename=input_filename, 
                                   output_file=output_file, 
                                   output_path=output_path ,
                                   partis_path=partis_path, 
                                   num_procs=num_procs,
                                   locus=locus,
                                   germline_dir=germline_dir,
                                   extra_columns=extra_columns
                                  )

    if(cleanup) {
        output_path %>% 
            unlink(recursive=TRUE)
    }

    return(partitioned_data)
}

#' Simulate a dataset based on parameters from partis annotations
#'
#' @inheritParams callPartis
#' @param parameter_dir The parent output folder for partis (which is '_output'
#'    by default. The function cd's into the params directory
#' @param num_events The desired number of VDJ rearragement events for the
#'    simulation
#' @param num_leaves The exponent a of the zipf p.m.f, which is of the form
#'    x^(-a). Thus, a higher value of a leads to fewer leaves in a given 
#'    clonal family, since the probability mass will be highly distributed near
#'    zero. a = 2 leads to larger clonal families, and this can be very, very
#'    slow if num_events is high.
#' @param seed The random generator seed to be supplied to 
#'    \code{partis simulate}
#' @param subsample_to_unique_clones If TRUE, one member of each clonal family
#'    is subsampled, and the function returns the resultant dataset of unique
#'    clones
#'    
#' @export
getPartisSimulation <- function(parameter_dir,
                                partis_path=Sys.getenv("PARTIS_PATH"),
                                output_file=file.path(parameter_dir,
                                                      "simu.csv"),
                                num_events=NULL,
                                num_leaves=NULL,
                                cleanup=TRUE,
                                do_full_annotation=FALSE,
                                extra_columns="v_gl_seq:v_qr_seqs:cdr3_seqs:naive_seq",
                                seed=NULL,
                                subsample_to_unique_clones=FALSE,
                                do_multi_hmm=FALSE,
                                extra_command_args=""
                               ) {
    partis_command <- paste(partis_path, 
                            "simulate", 
                            "--parameter-dir", 
                            file.path(parameter_dir, "params"),
                            "--outfname", 
                            output_file,
                            "--n-leaf-distribution",
                            "zipf"
                           )

    if(!missing(num_events)) {
        partis_command <- paste(partis_command,
                                "--n-sim-events", num_events)
    }
    
    if(!missing(num_leaves)) {
        partis_command <- paste(partis_command,
                                "--n-leaves", num_leaves)
    }

    if(!is.null(extra_columns)) {
        partis_command <- paste(partis_command,
                                "--extra-annotation-columns", extra_columns)
    }

    if(!is.null(seed)) {
        partis_command <- paste(partis_command,
                                "--seed",
                                seed
                               )
    }

    if(do_multi_hmm) {
        partis_command <- paste(partis_command,
                                "--parameter-type",
                                "multi-hmm"
                               )
    }

    partis_command <- paste(partis_command,
                            extra_command_args
                           )

    print(partis_command)

    partis_command %>% 
        system

    sim_annotations <- output_file %>% 
        fread(stringsAsFactors=FALSE)
    if(do_full_annotation) {
        sim_annotations <- getFullPartisAnnotation(output_file,
                                                   output_path=".",
                                                   partis_path
                                                  )
        # sim_annotations[["junction"]] <- sim_annotations %>% getCDR3s
    }

    if(cleanup) {
        output_file %>% 
            unlink
    }

    sim_annotations <- sim_annotations %>% 
        processPartisSequences

    sim_annotations[["clone_id"]] <- sim_annotations[["reco_id"]] %>% 
        sapply(as.numeric) %>%
        rank
    sim_annotations[["sequence"]] <- sim_annotations[["sequence_alignment"]]

    names(sim_annotations)[which(names(sim_annotations) == "cdr3_length")] <- 
        "junction_length"
    
    if(subsample_to_unique_clones) {
        sim_annotations <- sim_annotations %>%
            subsampleToUniqueClones
    }

    sim_data <- list(annotations=sim_annotations)
    return(sim_data)
}

#' Run partis annotate on vector or list of sequences directly (i.e. without
#'   saving as a fasta file)
#'
#' @param sequences List or vector of query BCR sequences for annotation
#' 
#' @export
getPartisAnnotationsFromStrings <- function(sequences, 
                                            ...
                                           ) {
    annotations <- tryCatch({
        filename <- Sys.time() %>%
            gsub(pattern=' |:|-', replacement='') %>%
            paste0('tmp', ., '.fasta')
        write.fasta(sequences, names="tmp", file.out=filename)
        getPartisAnnotations(filename, 
                             ...
                            )
    }, warning = function(warning_condition) {
        message(warning_condition)
        return(NA)
    }, error = function(error_condition) {
        message(error_condition)
        return(NULL)
    }, finally={
        filename %>% unlink
    })
    return(annotations)
}
