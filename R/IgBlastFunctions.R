require(data.table)
require(dplyr)
require(shazam)

#' Get IgBLAST annotations from an input fasta file
#'
#' The routine changes directories to the provided igblast directory,
#'   runs some Change-O commands, reads in the resultant tsv files,
#'   changes back to the initial directory, and returns the annotations
#'   list object.
#'
#' @param input_filename The path to the input fasta file
#' @param output_filename The desired output filename (should be a tsv file)
#' @param organism The organism/species to which the sequences belong
#' @param domain_system (input to IgBLAST -- alter with caution)
#' @param locus String denoting the locus for the receptors given in 
#'   \code{input_filename}. 
#'   Either "tra", "trb", "trd", or "trg" for TCRs,
#'   or "igl", "igk", or "igh" for BCRs
#' @param igblast_dir Path to parent directory of the igblast bin folder
#'   (e.g. "path/to/parent", NOT "path/to/parent/bin").
#' @param changeo_dir The path of the changeo bin folder (e.g.
#'   "path/to/bin").
#' @param cleanup If TRUE, remove all interim files created
#' @param nproc The number of processors for use within the igblast and
#'   shazam::distToNearest function calls. Defaults to 4.
#'   Setting nproc=alakazam::cpuCount() attempts to automatically detect the
#'   number of available cores -- see the alakazam documtation for details.
#' @return A \code{list} of a \code{data.table} containing the annotations.
getIgBlastAnnotations <- function(input_filename,
                                  output_filename="igblast_out.tsv",
                                  organism="human",
                                  domain_system="imgt",
                                  locus,
                                  nproc=4,
                                  igblast_dir=Sys.getenv("IGBLAST_DIR"),
                                  changeo_dir=Sys.getenv("CHANGEO_DIR"),
                                  cleanup=TRUE
                                 ) {
    checkForValidLocus(locus)
    if(igblast_dir == "") {
        stop(paste("Empty string given as igblast_dir argument.",
                   "Please specify the correct igblast directory,",
                   "or make sure IGBLAST_DIR environment variable is set."
                  )
        )
    }
    if(changeo_dir == "") {
        stop(paste("Empty string given as changeo_dir argument.",
                   "Please specify the correct Change-O directory,",
                   "or make sure CHANGEO_DIR environment variable is set."
                  )
        )
    }

    input_filename <- input_filename %>% 
        normalizePath

    igblast_bin_dir <- file.path(igblast_dir,
                                  "bin")

    initial_wd <- getwd()
    tryCatch({
        igblast_dir %>% setwd
        igblast_exec <- file.path(igblast_bin_dir, "igblastn")
        full_output_filename <- file.path(initial_wd,
                                          output_filename
                                         )

        receptor_type <- stringr::str_sub(locus, 0, 2)

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
                                 receptor_type,
                                 "--format",
                                 "blast",
                                 "--exec",
                                 igblast_exec,
                                 "-o",
                                 full_output_filename,
                                 "--nproc",
                                 nproc,
                                 sep=" "
                                )



        cat("\n", igblast_command, "\n")
        igblast_command %>% system

        gene_segments <- c("v", "j")
        if(locus %in% c("igh", "trb", "trd")) {
            gene_segments <- c(gene_segments, "d")
        }

        database_files <- gene_segments %>%
            sapply(function(x) { 
                       file.path(igblast_dir,
                                 "germlines",
                                 domain_system,
                                 organism,
                                 "vdj",
                                 paste(domain_system,
                                     organism,
                                     paste0(locus %>% toupper,
                                            x %>% toupper,
                                            ".fasta"
                                           )
                                     ,
                                     sep="_"
                                 )
                                 )
                   }
                ) %>%
            paste(collapse=" ")
                         
        make_db_command <- paste(file.path(changeo_dir,
                                           "MakeDb.py"),
                                 "igblast",
                                 "-i",
                                 full_output_filename,
                                 "-s",
                                 input_filename %>% normalizePath,
                                 "-r",
                                 database_files,
                                 "--regions",
                                 "--format",
                                 "airr"
                                ) 

        cat("\n", make_db_command, "\n")
        make_db_command %>% system 

        changeo_filename <- full_output_filename %>%
            gsub(pattern=".tsv",
                 replacement="_db-pass.tsv"
                )

        annotations <- changeo_filename %>%
            data.table::fread()


        if(receptor_type == "ig") {
            tryCatch({
                cluster_threshold <- annotations %>%
                    shazam::distToNearest(sequenceColumn="junction",
                                          vCallColumn="v_call",
                                          jCallColumn="j_call",
                                          nproc=nproc
                                         ) %$%
                    DIST_NEAREST %>%
                    shazam::findThreshold() %>%
                    slot(name="threshold")
            }, error = function(e) {
                stop(e)
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

    annotations[["vj_in_frame"]] <- annotations[["vj_in_frame"]] %>%
        as.logical
    annotations[["stop_codon"]] <- annotations[["stop_codon"]] %>%
        as.logical

    initial_wd %>% setwd

    annotation_object <- list(annotations=annotations)
    return(annotation_object)
}
