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
subsampleVector <- function(v, sample_count) {
    new_vector <- {}
    if(length(v) > sample_count) {
        new_vector <- sample(v, sample_count)
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
subsample <- function(dataset, sample_count) {
    new_dataset <- {}
    if(is.null(dim(dataset))) {
        new_dataset <- subsampleVector(dataset, sample_count)
    } else if(nrow(dataset) > sample_count) {
        new_dataset <- dataset[sample(1:nrow(dataset), sample_count), ]
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
