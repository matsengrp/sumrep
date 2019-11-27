#' Remove empty strings from a list or vector
#'
#' @param l List or vector of strings
#' @return The given list or vector with all empty strings removed
#' 
#' @export
removeEmptyStrings <- function(l) {
    return(l[l != ""])
}

maskEmptyStringWithNA <- function(s) {
    return(ifelse(s != "", s, NA))
}

#' Convert lists of factors or vectors of characters into a vector of strings
#' 
#' @param List or vector of either strings or char vectors of DNA sequences
#' @return Vector of strings of DNA sequences
#' 
#' @export
standardizeList <- function(l) {
    new_list <- l %>% 
        sapply(toString) %>%
        gsub(pattern=", *", replacement="") %>%
        sapply(paste, collapse='') %>% 
        unname
    return(new_list)
}

#' Import a vector of DNA sequence strings from a fasta file
#' 
#' @param filename Name of fasta file including the sequences
#' @return A vector of DNA sequence strings
#' 
#' @export
getSequenceListFromFasta <- function(filename) {
    sequences <- filename %>% 
        seqinr::read.fasta() %>% 
        lapply(paste, collapse="") %>% 
        unlist %>% 
        unname
    return(sequences)
}

correspondsToProperAASequence <- function(x) {
    condition <- (nchar(x) >= 3) && (nchar(x) %% 3 == 0)
    return(condition) 
}

# Convert each sequence with < 3 bases into NA
filterStringsForAAFunctions <- function(sequence_list) {
    filtered_list <- sequence_list %>%
        sapply(function(x) {
                            ifelse(correspondsToProperAASequence(x),
                                   x,
                                   NA)
                            }
                            )
    return(filtered_list)
}

#' Convert a single DNA sequence string to an amino acid sequence string
#' 
#' @param sequence String of DNA bases
#' @return String of single-letter amino acid codes
#' 
#' @export
convertNucleobasesToAminoAcidsBySequence <- function(sequence) {
    aa <- ifelse(!is.na(sequence),
                 sequence %>%
        strsplit(split='') %>%
        unlist %>%
        seqinr::translate() %>%
        unname %>%
        paste0(collapse=''),
        NA)
    
    return(aa)
}

#' Convert a vector of DNA strings to a strings of single-letter amino acid codes.
#' If the string does not correspond to a valid amino acid, convert it to NA.
#'
#' @param sequence_list Vector of string of DNA bases
#' @return Vector of strings of single-letter amino acid codes
#' 
#' @export
convertNucleobasesToAminoAcids <- function(sequence_list) {
    aa_sequences <- sequence_list %>% 
        filterStringsForAAFunctions %>%
        sapply(convertNucleobasesToAminoAcidsBySequence)
    return(aa_sequences)
}

#' Check if an amino acid sequence has unrecognized codons
#' 
#' @param aa_sequence String of single-letter amino acid codons
#' @return A boolean stating whether any unrecognized AA codons are present 
#'   in \code{aa_sequence}
#'   
#' @export
hasUnrecognizedAminoAcids <- function(aa_sequence) {
    return(grepl("X", aa_sequence))  
}

#' Check if an amino acid sequence has a stop codon
#' 
#' @param aa_sequence String of single-letter amino acid codes
#' @return A boolean stating whether any stop codons are present 
#'   in \code{aa_sequence}
#'   
#' @export
hasStopCodon <- function(aa_sequence) {
    return(grepl("\\*", aa_sequence))
}

#' Filter amino acid sequences by removing any 'X' (unrecognized) values,
#' as well as converting sequences with stop codons to NA
#'
#' @param aa_sequences Vector or list of amino acid sequences
#' @return Filtered list with "bad" sequences removed
#' 
#' @export
filterAminoAcidSequences <- function(aa_sequences) {
    filtered_seqs <- aa_sequences %>%
        sapply(function(x) {
        ifelse(hasStopCodon(x) || hasUnrecognizedAminoAcids(x),
                            NA,
                            x
                            )
        })

    return(filtered_seqs)
}

#' Parse a python dictionary string into an R list
#'
#' @param dictionary String that represents a dictionary in Python's syntax
#' @return An R list containing the same information as \code{dictionary}
#' 
#' @export
parsePythonDictionary <- function(dictionary) {
    parsed <- dictionary %>% gsub(pattern="'", replacement='"') %>% 
        gsub(pattern="\\(", replacement="\\[") %>% 
        gsub(pattern="\\)", replacement="\\]") %>% 
        jsonlite::fromJSON()
    return(parsed)
}

removeSequencesWithDifferentGermlineAndSequenceLengths <- 
    function(dat,
             germline_column="germline_alignment",
             sequence_column="sequence_alignment"
            ) {
    return(dat %>% 
               subset(nchar(getColumnValues(dat, germline_column)) == 
                      nchar(getColumnValues(dat, sequence_column))
                     )
          )
}

#' Load datasets that are not already in the workspace.
#'   The variable names are taken to be exactly the name of the corresponding
#'   .rds files
#' 
#' @param data_dir The directory which contains .rds files to be loaded
#' 
#' @export
loadNewDatasets <- function(data_dir,
                            pattern=""
                           ) {
    for(data_file in list.files(data_dir,
                                pattern=pattern
                               )) {
        var_name <- data_file %>%
            gsub(pattern="-", replacement="_") %>%
            gsub(pattern=".rds", replacement="")
        if(!exists(var_name)) {
            assign(var_name, 
                readRDS(file.path(data_dir, data_file)),
                envir=.GlobalEnv)
        }
    }
}

#' Handy automatic dimension retrieval given number of grid cells 
#'
#' @export
getGridDims <- function(n) {
    cols <- n %>% sqrt %>% floor 
    rows <- ceiling(n/cols) %>% floor 
    return(c(cols, rows))
}

#' Display an array of ggplots within one figure
#' @param plotlist List of plots created with ggplot
#' @param cols Number of columns for the figure
#' @param rows Number of rows for the figure
#' 
#' @export
multiplot <- function(plotlist=NULL, 
                      cols, 
                      rows, 
                      layout=NULL, 
                      tall_plot=FALSE, 
                      ...
                     ) {
    plots <- c(list(...), plotlist)
    numPlots = length(plots)

    if(missing(cols) || missing(rows)) {
        grid_dims <- numPlots %>% 
            getGridDims
        if(tall_plot) {
            cols <- grid_dims[1]
            rows <- grid_dims[2]
        } else {
            rows <- grid_dims[1]
            cols <- grid_dims[2]
        }
        
    }

    if(is.null(layout)) {
        layout <- matrix(seq(1, cols * rows),
                         ncol = cols, nrow = rows)
    }
    if (numPlots==1) {
       print(plots[[1]])
    } else {
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        for(i in 1:numPlots) {
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx[["row"]],
                                            layout.pos.col = matchidx[["col"]]))
        }   
    }     
} 

#' A wrapper of \code{alakazam::checkColumns} that stops execution when an
#'   error is thrown.
#' Checks whether \code{column} is a column name present in \code{dat},
#'   and whether this column contains at least some non-missing values.
#'
#' @param dat A \code{data.table} corresponding to repertoire annotations 
#' @param column A string naming the expected column of \code{dat}
#' 
#' @export
checkColumn <- function(dat, column) {
    column_check <- alakazam::checkColumns(dat, column)
    if(column_check != TRUE) {
        stop(column_check)
    }
}

#' A wrapper of \code{alakazam::checkColumns} that stops execution when an
#'   error is thrown, and returns the vector of column values if available.
#'
#' @param dat A \code{data.table} corresponding to repertoire annotations 
#' @param column A string naming the expected column of \code{dat}
#' @return A vector of values given by \code{dat[[column]]}, if available
#' 
#' @export
getColumnValues <- function(dat, column) {
    checkColumn(dat, column)
    column_values <- dat[[column]]
    return(column_values)
}

#' This function attempts to access the provided column from the provided
#'   dat, perfoming any optional filters that were specified, and returns
#'   a vector of the relevant sequence strings.
#'
#' @param dat A \code{data.table} corresponding to repertoire annotations 
#' @param column A string naming the expected column of \code{dat}
#' @param drop_gaps If TRUE, removes '.' and '-' gap characters from each
#'   sequence
#' @param remove_stop_codon_seqs If TRUE, removes any sequence with a 
#'   stop codon. This assumes the \code{stop_codon} column is present in
#'   \code{dat}, as defined by the AIRR rearrangement schema.
#' @param in_frame_only If TRUE, removes any sequence with an
#'   out-of-frame V or J segment.
#' @return A vector of filtered strings from by \code{dat[[column]]}, 
#'   if available
#'   
#' @export
getColumnSequences <- function(dat,
                               column,
                               drop_gaps,
                               remove_stop_codon_seqs=TRUE,
                               in_frame_only=TRUE
                              ) {
    if(remove_stop_codon_seqs) {
        tryCatch({
            checkColumn(dat, "stop_codon")
            dat <- dat[!dat[["stop_codon"]], ]
        }, error=function(e) {
            cat("\nWarning: stop_codon column not present.\n")
            cat("  Unable to remove sequences with stop codons.\n\n")
        })
    }

    if(in_frame_only) {
        tryCatch({
            checkColumn(dat, "vj_in_frame")
            dat <- dat[dat[["vj_in_frame"]], ]
        }, error=function(e) {
            cat("Warning: vj_in_frame column not present.\n")
            cat("  Unable to remove sequences with out-of-frame V or J segments.\n\n")
        })
    }

    sequences <- dat %>%
        getColumnValues(column=column)

    if(drop_gaps) {
        sequences <- sequences %>%
            sapply(function(x) { 
                       gsub(x,
                            pattern="\\.|\\-",
                            replacement=""
                           )
                   }
            ) %>%
            unname
    }

    return(sequences)
}

subsampleToUniqueClones <- function(dat,
                                    clone_id_column="clone_id"
                                   ) {
    clone_list <- dat$clone_id %>% unique
    clone_sample_ids <- {}
    for(clone in clone_list) {
        clone_sample_ids <-
            c(clone_sample_ids,
              which(getColumnValues(dat, clone_id_column) == clone) %>%
              subsample(1)
             )
    }
    unique_clone_dat <- dat[clone_sample_ids, ]
    return(unique_clone_dat)
}
