library(data.table)
library(shazam)

#' Get IgBLAST annotations from an input fasta file
#'
#' Currently kind of hacky, as it changes directories to the location of
#'   igblast and its germline databases, runs some commands, and changes back
#'   to the initial directory.
#'   Change-O is used to parse the raw output of igblastn.
#'   TODO: Delete interim files, try-catch the success of changeo
#'
#' @param input_filename The path to the input fasta file
#' @param output_filename The desired output filename WITHOUT an extension
#' @param locus The locus from which the sequences are derived
#' @param organism The organism/species to which the sequences belong
#' @param domain_system (input to IgBLAST)
#' @param ig_seqtype (input to IgBLAST
getIgBlastAnnotations <- function(input_filename,
                                  output_filename="igblast_out",
                                  locus="igh",
                                  organism="human",
                                  domain_system="imgt",
                                  ig_seqtype="Ig",
                                  num_threads=8
                                  ) {

    input_filename <- input_filename %>% 
        normalizePath
    igblast_dir <- "~/Software/igblast/bin"
    changeo_dir <- "~/.local/bin"

    # IgBlast needs to be in its own directory to find the germline databases
    # easily, so let's remember where we are currently and chnage directories
    initial_wd <- getwd()
    igblast_dir %>% setwd
    igblast_exec <- file.path(igblast_dir, "igblastn")
    arguments <- igblast_exec
    regions <- c("V", "D", "J")
    for(region in regions) {
        arguments <- c(arguments, 
                     paste0("-germline_db_", region),
                     file.path(
                               locus, 
                               paste0(locus,
                                   region %>% tolower,
                                   "-unaligned.fasta")
                               )
                     )
    }

    full_output_filename <- 
                   paste(output_filename,
                         "fmt7",
                         sep='.')

    arguments <- c(arguments,
                   "-auxiliary_data",
                   file.path(
                             "optional_file",
                             paste(organism, 
                                   "gl.aux", 
                                   sep="_")
                             ),
                   "-domain_system",
                   domain_system,
                   "-ig_seqtype",
                   ig_seqtype,
                   "-organism",
                   organism,
                   "-outfmt",
                   "'7 std qseq sseq btop'", # Suggested by Change-O docs
                   "-query",
                   # Give absolute path for 
                   # input file, since we have changed directories
                   input_filename %>% normalizePath, 
                   "-out",
                   full_output_filename,
                   "-num_threads",
                   num_threads
                   )
    
    command <- arguments %>%
        paste(collapse=' ')

    command %>% system

     germline_sequences <- regions %>%
                               sapply(tolower) %>%
                               sapply(function(x) {
                                          paste0(locus,
                                                 '/',
                                                 locus,
                                                 x,
                                                 '.fasta')
                                   }
                               ) %>%
                               paste(collapse=' ')

    changeo_command <- paste(file.path(changeo_dir,
                                       "MakeDb.py"),
                             "igblast",
                             "-i",
                             full_output_filename,
                             "-s",
                             input_filename,
                             "-r",
                             germline_sequences
                             )

    changeo_command %>% system
    changeo_filename <- output_filename %>%
        paste0("_db-pass.tab")

    annotations_without_naive <- changeo_filename %>%
        data.table::fread(stringsAsFactors=TRUE)

    # Need to do clustering with some help from shazam and changeo
    #   before getting naive sequences
    cluster_threshold <- annotations_without_naive %>%
        shazam::distToNearest() %$%
        DIST_NEAREST %>%
        shazam::findThreshold() %>%
        slot(name="threshold")

    cluster_command <- paste(file.path(changeo_dir,
                                       "DefineClones.py"),
                             "bygroup",
                             "-d",
                             changeo_filename,
                             "--act",
                             "set",
                             "--model",
                             "ham",
                             "--sym",
                             "min",
                             "--norm",
                             "len",
                             "--dist",
                             cluster_threshold
                             )

    cluster_command %>% system

    changeo_clone_filename <- changeo_filename %>%
        gsub(pattern=".tab", replace="_clone-pass.tab")

    germline_command <- paste(file.path(changeo_dir,
                                        "CreateGermlines.py"),
                              "-d",
                              changeo_clone_filename,
                              "-r",
                              germline_sequences,
                              "-g",
                              "full",
                              "--cloned"
                              )

    germline_command %>% system

    germline_filename <- changeo_clone_filename %>%
        gsub(pattern=".tab", replace="_germ-pass.tab")

    annotations <- germline_filename %>%
        data.table::fread(stringsAsFactors=TRUE) 

    names(annotations)[which(names(annotations) == "V_CALL")] <- "v_gene"
    names(annotations)[which(names(annotations) == "D_CALL")] <- "d_gene"
    names(annotations)[which(names(annotations) == "J_CALL")] <- "j_gene"
    names(annotations)[which(names(annotations) == "SEQUENCE_INPUT")] <- "mature_seq"
    names(annotations)[which(names(annotations) == "SEQUENCE_VDJ")] <- "naive_seq"
    names(annotations)[which(names(annotations) == "JUNCTION")] <- "cdr3s"
    names(annotations)[which(names(annotations) == "JUNCTION_LENGTH")] <- "cdr3_length"

    annotations$mature_seq <- annotations$mature_seq %>%
        sapply(toString)
    annotations$naive_seq <- annotations$naive_seq %>%
        sapply(toString)

    # Go back to initial working directory
    initial_wd %>% setwd

    annotation_object <- list(annotations=annotations)
    return(annotation_object)
}
