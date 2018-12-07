require(data.table)
require(dplyr)
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
#' @param output_filename The desired output filename (should be a tsv file)
#' @param organism The organism/species to which the sequences belong
#' @param domain_system (input to IgBLAST -- alter with caution)
#' @param receptor_type String denoting the type of receptors given in 
#'   \code{input_filename}. Either "BCR" or "TCR".
#' @param igblast_dir Path to parent directory of the igblast bin folder
#'   (e.g. "path/to/parent", NOT "path/to/parent/bin").
#' @param changeo_dir The path of the changeo bin folder (e.g.
#'   "path/to/bin").
#' @param cleanup If TRUE, remove all interim files created
getIgBlastAnnotations <- function(input_filename,
                                  output_filename="igblast_out.tsv",
                                  organism="human",
                                  domain_system="imgt",
                                  receptor_type,
                                  num_threads=8,
                                  igblast_dir,
                                  changeo_dir,
                                  cleanup=TRUE
                                 ) {
    if(receptor_type %>% missing) {
        stop("receptor_type argument must be specified.")
    }
    if(igblast_dir %>% missing) {
        stop("igblast_dir argument must be specified.")
    }
    if(changeo_dir %>% missing) {
        stop("changeo_dir argument must be specified.")
    }

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
        full_output_filename <- file.path(initial_wd,
                                          output_filename
                                         )

        loci_hash <- list(BCR="ig",
                          TCR="tr"
                         )
        
        igblast_command <- paste(file.path(changeo_dir,
                                           "AssignGenes.py"),
                                 "igblast",
                                 "-s",
                                 input_filename %>% normalizePath,
                                 "-b",
                                 igblast_dir,
                                 "--organism",
                                 organism,
                                 "--loci",
                                 loci_hash[[receptor_type]],
                                 "--format",
                                 "airr",
                                 "--exec",
                                 igblast_exec,
                                 "-o",
                                 full_output_filename,
                                 sep=" "
                                )



        cat("\n", igblast_command, "\n")
        igblast_command %>% system

        changeo_filename <- full_output_filename

        annotations <- changeo_filename %>%
            data.table::fread(stringsAsFactors=TRUE)

        if(receptor_type == "BCR") {
            try({
                cluster_threshold <- annotations %>%
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

            names(annotations)[which(names(annotations) == "clone")] <- 
                "clone_id"
        }
    }, error = function(e) {
        print(e)
        setwd(initial_wd)
    }, finally = {
        if(cleanup) {
           c(
             "changeo_filename",
             "cluster_filename"
            ) %>%
                sapply(function(x) { 
                           if(exists(x)) { 
                               eval(parse(text=x)) %>%
                                   file.remove
                           } 
                       } 
                )
        }
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
