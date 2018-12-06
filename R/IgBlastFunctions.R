require(data.table)
require(shazam)

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
                                  output_filename="igblast_out.tsv",
                                  locus="IGH",
                                  organism="human",
                                  domain_system="imgt",
                                  loci="ig",
                                  num_threads=8,
                                  igblast_dir,
                                  changeo_dir,
                                  airr_format=TRUE
                                 ) {

    input_filename <- input_filename %>% 
        normalizePath

    igblast_bin_dir <- file.path(igblast_dir,
                                  "bin")

    # IgBlast needs to be in its own directory to find the germline databases
    # easily, so let's remember where we are currently and chnage directories
    initial_wd <- getwd()
    tryCatch({
        igblast_dir %>% setwd
        igblast_exec <- file.path(igblast_bin_dir, "igblastn")
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

        full_output_filename <- file.path(initial_wd,
                                          output_filename
                                         )

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
                       loci,
                       "-organism",
                       organism,
                       "-outfmt",
                       ifelse(airr_format, 
                              19, # AIRR tab-delimited file
                              "'7 std qseq sseq btop'" # blast-style tabular output
                              ),
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

        command <- paste(file.path(changeo_dir,
                                   "AssignGenes.py"),
                         "igblast",
                         "-s",
                         input_filename %>% normalizePath,
                         "-b",
                         igblast_dir,
                         "--organism",
                         organism,
                         "--loci",
                         loci,
                         "--format",
                         "airr",
                         "--exec",
                         igblast_exec,
                         "-o",
                         full_output_filename,
                         sep=" "
                        )



        cat("\n", command, "\n")
        command %>% system

         germline_sequences <- regions %>%
                                   sapply(function(x) {
                                              file.path(
                                                     igblast_dir,
                                                     "germlines/imgt",
                                                     organism,
                                                     "vdj",
                                                     paste0(
                                                         paste(
                                                           "imgt",
                                                           organism,
                                                           paste0(locus, x),
                                                           sep='_'),
                                                         '.fasta')
                                                     )
                                       }
                                   ) %>%
                                   paste(collapse=' ')

        changeo_filename <- full_output_filename

        annotations_without_naive <- changeo_filename %>%
            data.table::fread(stringsAsFactors=TRUE)

        # Need to do clustering with some help from shazam and changeo
        #   before getting naive sequences
        try({
        cluster_threshold <- annotations_without_naive %>%
            shazam::distToNearest(sequenceColumn="junction",
                                  vCallColumn="v_call",
                                  jCallColumn="j_call"
                                 ) %$%
            DIST_NEAREST %>%
            shazam::findThreshold() %>%
            slot(name="threshold")
        })
        if(!exists("cluster_threshold")) {
            cluster_threshold <- 0.5
        }

        cluster_filename <-  file.path(
                                      gsub(changeo_filename,
                                      pattern=".tsv",
                                      replace="_clustered.tsv"
                                     )
                                      )
        cluster_command <- paste(file.path(changeo_dir,
                                           "DefineClones.py"),
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
                                 cluster_threshold,
                                 "-o",
                                 cluster_filename
                                )

        cluster_command %>% system

        annotations <- cluster_filename %>%
            data.table::fread() 

        names(annotations) <- names(annotations) %>%
            sapply(tolower)
    }, error = function(e) {
        print(e)
    }, finally = {
        setwd(initial_wd)
    })

    # Go back to initial working directory
    initial_wd %>% setwd

    annotation_object <- list(annotations=annotations)
    return(annotation_object)
}

#' Add n's to positions for which the input (mature) sequence did not contain
#'   a base, but for which the inferred naive (germline) sequence does contain
#'   a base. IgBlast/Change-O by default place '.'s in these positions, which
#'   is also the symbol to denote IMGT gaps. Thus, manual processing to account
#'   for this is required.
getMatureSequences <- function(germline_imgt,
                               sequence_imgt
                              ) {
    naive_sequence <- germline_imgt %>%
        toString %>%
        tolower

    mature_sequence <- sequence_imgt %>%
        toString %>%
        tolower

    missing_indices <- getMissingDataIndices(naive_sequence,
                                             mature_sequence)

    mature_seq_split <- mature_sequence %>%
        strsplit('') %>%
        first 

    mature_seq_split[missing_indices] <- "n"

    mature_sequence <- mature_seq_split %>%
        paste(collapse='') %>%
        gsub(pattern='\\.', replace='')

    return(mature_sequence) 
}

getMissingDataIndices <- function(germline, mature) {
    germline_dots <- germline %>%
        getDotIndices

    mature_dots <- mature %>%
        getDotIndices

    missing_indices <- setdiff(mature_dots,
                               germline_dots)
    return(missing_indices)
}

getDotIndices <- function(seq) {
    indices <- seq %>%
        strsplit('') %>%
        first %>%
        grep(pattern='\\.')
    
    return(indices)
}

processGenes <- function(genes) {
    genes <- genes %>%
        sapply(toString) %>%
        # Remove all but first gene
        sapply(gsub, pattern=',.*', replace='') %>%
        # Remove gene after OR, e.g. IGHV3/OR12-1*01 should be IGHV3-1*01
        sapply(gsub, pattern='\\/OR.*-', replace='-') %>% 
        as.factor %>%
        unname
    return(genes)
}
