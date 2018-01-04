library(data.table)

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
                                  ig_seqtype="Ig"
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
                   full_output_filename
                   )
    
    command <- arguments %>%
        paste(collapse=' ')

    command %>% system

    changeo_command <- paste(file.path(changeo_dir,
                                       "MakeDb.py"),
                             "igblast",
                             "-i",
                             full_output_filename,
                             "-s",
                             input_filename,
                             "-r",
                             regions %>%
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
                             )

    changeo_command %>% system
    changeo_filename <- output_filename %>%
        paste0("_db-pass.tab")
    annotations <- changeo_filename %>%
        data.table::fread(stringsAsFactors=TRUE)

    # Go back to initial working directory
    initial_wd %>% setwd

    return(annotations)
}
