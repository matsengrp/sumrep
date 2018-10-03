#' Remove empty strings from a list or vector
#'
#' @param l List or vector of strings
#' @return The given list or vector with all empty strings removed
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
standardizeList <- function(l) {
    new_list <- l %>% 
        sapply(toString) %>%
        gsub(pattern=", *", replace="") %>%
        sapply(paste, collapse='') %>% 
        unname
    return(new_list)
}
#' Subsample a vector
#' 
#' @inheritParams subsample
#' @param v Vector to subsample
#' @return A subsampled vector 
subsampleVector <- function(v, 
                            sample_count,
                            replace=TRUE
                           ) {
    new_vector <- {}
    if(length(v) > sample_count) {
        new_vector <- sample(v, 
                             sample_count,
                             replace
                            )
    } else {
        new_vector <- v
    }
    return(new_vector)
}

#' Subsample a dataset
#'
#' @param dataset A data.table, data.frame, or vector object
#' @param sample_count Number of samples to retain in the subsampled data.
#'   Samples refers to elements in a vector or rows in a data.table/data.frame
#' @return A subsampled dataset of the same type given by \code{dataset}
subsample <- function(dataset, 
                      sample_count,
                      replace=TRUE
                     ) {
    new_dataset <- {}
    if(is.null(dim(dataset))) {
        new_dataset <- subsampleVector(dataset, 
                                       sample_count,
                                       replace=replace
                                      )
    } else if(nrow(dataset) > sample_count) {
        new_dataset <- dataset[sample(1:nrow(dataset), 
                                      sample_count,
                                      replace=replace
                                     ), ]
    } else {
        new_dataset <- dataset
    }
    return(new_dataset)
}

#' Import a vector of DNA sequence strings from a fasta file
#' 
#' @param filename Name of fasta file including the sequences
#' @return A vector of DNA sequence strings
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
hasUnrecognizedAminoAcids <- function(aa_sequence) {
    return(grepl("X", aa_sequence))  
}

#' Check if an amino acid sequence has a stop codon
#' 
#' @param aa_sequence String of single-letter amino acid codes
#' @return A boolean stating whether any stop codons are present 
#'   in \code{aa_sequence}
hasStopCodon <- function(aa_sequence) {
    return(grepl("\\*", aa_sequence))
}

#' Filter amino acid sequences by removing any 'X' (unrecognized) values,
#' as well as converting sequences with stop codons to NA
#'
#' @param aa_sequences Vector or list of amino acid sequences
#' @return Filtered list with "bad" sequences removed
filterAminoAcidSequences <- function(aa_sequences) {
    filtered_seqs <- ifelse(hasStopCodon(aa_sequences),
                            NA,
                            aa_sequences
                            ) %>%
        gsub(pattern="X", replace="")
    return(filtered_seqs)
}

#' Parse a python dictionary string into an R list
#'
#' @param dictionary String that represents a dictionary in Python's syntax
#' @return An R list containing the same information as \code{dictionary}
parsePythonDictionary <- function(dictionary) {
    parsed <- dictionary %>% gsub(pattern="'", replacement='"') %>% 
        gsub(pattern="\\(", replacement="\\[") %>% 
        gsub(pattern="\\)", replacement="\\]") %>% 
        jsonlite::fromJSON()
    return(parsed)
}

removeSequencesWithDifferentNaiveAndMatureLengths <- function(dat) {
    return(dat %>% subset(nchar(dat$mature_seq) == nchar(dat$naive_seq)))
}

#' Load datasets that are not already in the workspace.
#'   The variable names are taken to be exactly the name of the corresponding
#'   .rds files
#' 
#' @param data_dir The directory which contains .rds files to be loaded
loadNewDatasets <- function(data_dir) {
    for(data_file in list.files(data_dir)) {
        var_name <- data_file %>%
            gsub(pattern="-", replace="_") %>%
            gsub(pattern=".rds", replace="")
        if(!exists(var_name)) {
            assign(var_name, 
                readRDS(file.path(data_dir, data_file)),
                envir=.GlobalEnv)
        }
    }
}

#' Approximate the distribution of \code{summary_function} for \code{dat} by 
#    iteratively sampling from \code{dat}, computing the distribution of 
#'   \code{summary_function} for the subsample, and accumulating a representative
#'   distribution. The function stops when \code{divergence_function} between two
#'   distribution iterates is sufficiently small.
#'
#' @param dat The dataset for which the distribution is desired
#' @param summary_function A function which computes a summary statistic of \code{dat},
#'   which yields a distribution
#' @param sample_count The size of the subsample for each iteration
#' @param tol The threshold for \code{divergence_function} to be considered converged
#' @param divergence_function The divergence computed between successive iterates.
#'   This should probably be the JS divergence, specified for discrete or continuous data
#'   depending on the nature of \code{summary_function}
getApproximateDistribution <- function(dat,
                                       summary_function,
                                       sample_count=100,
                                       tol=0.001,
                                       divergence_function=getContinuousJSDivergence,
                                       ...
                                       ) {
    dist <- subsample(dat, sample_count) %>%
        summary_function(...)

    error <- Inf
    while(error > tol) {
        dist_prev <- dist
        sample_dat <- subsample(dat, sample_count)
        sample_dist <- sample_dat %>% 
            summary_function(...)
        dist <- c(dist_prev, sample_dist)
        error <- divergence_function(dist, dist_prev)
    }

    return(dist)
}

#' Handy automatic dimension retrieval given number of grid cells 
#'
getGridDims <- function(n) {
    cols <- n %>% sqrt %>% floor 
    rows <- ceiling(n/cols) %>% floor 
    return(c(cols, rows))
}

#' Display an array of ggplots within one figure
#' @param plotlist List of plots created with ggplot
#' @param cols Number of columns for the figure
#' @param rows Number of rows for the figure
multiplot <- function(plotlist=NULL, 
                      cols, 
                      rows, 
                      layout=NULL, 
                      tall_plot=FALSE, 
                      ...
                     ) {
    library(grid)
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
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }   
    }     
} 
