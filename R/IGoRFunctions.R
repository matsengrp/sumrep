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
    unlink(file.path(dir.name, "aligns"), recursive=TRUE)
    unlink(file.path(dir.name, 
                     gen_batch_name), 
           recursive=TRUE
          )
    unlink(file.path(dir.name, 
                     eval_batch_name), 
           recursive=TRUE
          )

    igor_prefix <- paste("igor",
                         "-set_wd",
                         dir_name
                        )

    read_seqs_command <- paste(igor_prefix,
                               "-batch",
                               eval_batch_name,
                               "-read_seqs",
                               input_filename
                              )
    read_seqs_command %>% system

    igor_prefix <- paste(igor_prefix,
                         "-species",
                         species,
                         "-chain",
                         chain
                        )

    align_command <- paste(igor_prefix,
                           "-batch",
                           eval_batch_name,
                           "-align",
                           "--all"
                          )

    align_command %>% system

    evaluate_command <- paste(igor_prefix,
                              "-batch",
                              eval_batch_name,
                              "-evaulate",
                              "-output",
                              "--scenarios",
                              num_scenarios,
                              "--CDR3"
                             )

    evaluate_command %>% system

    generate_command <- paste(igor_prefix,
                              "-batch",
                              gen_batch_name,
                              "-generate",
                              num_gen_sequences,
                              "--CDR3"
                             )

    generate_command %>% system
}

getIgorAnnotations <- function(input_filename,
                               locus,
                               igor_wd_name="igor_wd",
                               organism="human",
                               output_filename="annotations.csv",
                               nproc=4,
                               cleanup=TRUE
                              ) {
    chain_hash <- list("tra"="alpha",
                       "trb"="beta"
                      )
    chain <- chain_hash[locus]
    igor_input_filename <- "tmp.txt"
    convertFastaToTxt(input_filename,
                      output_filename=igor_input_filename)
    igor_script <- system.file("run_igor.py",
                               package="sumrep")
    python_command <- paste("python3",
                            igor_script,
                            igor_input_filename,
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
    # since annotations[["seq_index"]] contains zero-based positions
    annotations[["sequence"]] <- 
        indexed_seqs[["sequence"]][1 + annotations[["seq_index"]]] %>%
        tolower
    annotations[["sequence_alignment"]] <- annotations[["sequence"]]

    seq_align_file <- file.path(igor_wd_name, "seq_align.fa")
    write.fasta(annotations[["sequence_alignment"]] %>% as.list,
                names=1:nrow(annotations),
                file.out=seq_align_file,
                as.string=TRUE
               )

    igblast_outfile <- file.path(igor_wd_name, "igblast_out.tsv")
    igblast_dir <- Sys.getenv("IGBLAST_DIR")
    igblast_command <- paste(file.path(Sys.getenv("CHANGEO_DIR"),
                                       "AssignGenes.py"
                                      ),
                             "igblast",
                             "-s",
                             seq_align_file,
                             "-b",
                             igblast_dir,
                             "--organism",
                             "human",
                             "--loci",
                             stringr::str_sub(locus, 0, 2),
                             "--format",
                             "airr",
                             "--exec",
                             file.path(igblast_dir,
                                       "bin",
                                       "igblastn"
                                      ),
                             "-o",
                             igblast_outfile
                            )
    
    igblast_command %>% cat('\n')
    igblast_command %>% system

    igblast_dat <- data.table::fread(igblast_outfile)
    annotations[["junction"]] <- igblast_dat[["junction"]] %>% 
        tolower
    annotations[["junction_aa"]] <- igblast_dat[["junction_aa"]]


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
    sim_annotations[["sequence"]] <-
        sim_indexed_seqs[["nt_sequence"]][1 + sim_annotations[["seq_index"]]] %>%
        tolower
        
    sim_annotations[["sequence_alignment"]] <- sim_annotations[["sequence"]]

    names(sim_annotations)[which(names(sim_annotations) == "nt_CDR3")] <- 
        "junction"
    sim_annotations[["junction"]] <- sim_annotations[["junction"]] %>% 
        tolower
    sim_annotations[["junction_aa"]] <- sim_annotations[["junction"]] %>%
        convertNucleobasesToAminoAcids

    if(cleanup) {
        igor_wd_name %>% 
            unlink(recursive=TRUE)
    }

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
        annotations[["vj_in_frame"]] <- annotations[["is_inframe"]] %>%
            as.logical
    }
    return(annotations)
}

convertFastaToTxt <- function(input_filename,
                              output_filename
                             ) {
    tmp_dat <- input_filename %>%
        seqinr::read.fasta(seqonly=TRUE, as.string=TRUE) %>%
        data.table::data.table() %>%
        data.table::fwrite(file=output_filename,
                           col.names=FALSE)
}
