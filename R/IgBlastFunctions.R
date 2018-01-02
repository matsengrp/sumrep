#' Run IgBLAST on an input fasta file
#'
#' Currently kind of hacky, as it changes directories to the location of
#'   igblast and its germline databases, runs the command, and changes back
#'   to the initial directory
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

    # IgBlast needs to be in its own directory to find the germline databases
    # easily, so let's remember where we are currently and chnage directories
    initial_wd <- getwd()
    igblast_dir <- "~/Software/igblast/bin"
    igblast_dir %>% setwd
    igblast_exec <- file.path(igblast_dir, "igblastn")
    arguments <- igblast_exec
    for(region in c("V", "D", "J")) {
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
                   paste(output_filename,
                         "fmt7",
                         sep='.')
                   )
    
    command <- arguments %>%
        paste(collapse=' ')

    command %>% system

    initial_wd %>% setwd

    return(command)
}
