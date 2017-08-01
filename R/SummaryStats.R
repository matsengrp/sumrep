library(alakazam)
library(ape)
library(Biostrings)
library(data.table)
library(dplyr)
library(flexmix)
library(jsonlite)
library(magrittr)
library(pegas)
library(RecordLinkage)
library(shazam)
library(seqinr)
library(stringdist)
library(textmineR)
library(yaml)

compareCategoricalDistributions <- function(table_a, table_b) {
    table_a_dims <- table_a %>% dim
    table_b_dims <- table_b %>% dim
    dim_count_a <- table_a_dims %>% length
    dim_count_b <- table_b_dims %>% length
    if(dim_count_a != dim_count_b || !all(table_a_dims == table_b_dims)) {
        stop("Table dimensions must agree")
    }
    dims <- table_a_dims
    entrywise_sum_of_absolute_differences <- (table_a - table_b) %>% abs %>% sum
    entrywise_comparison <- entrywise_sum_of_absolute_differences/prod(dims)
    dimension_count <- dims %>% length
    dimension_comparisons <- rep(NA, dimension_count)
    for(i in 1:dimension_count) {
        dimension_sums_a <- apply(table_a, i, sum)
        dimension_sums_b <- apply(table_b, i, sum)
        dimension_sum_of_absolute_differences <- 
            (dimension_sums_a - dimension_sums_b) %>% abs %>% sum
        other_dims <- dims[-i]
        dimension_comparisons[[i]] <- 
            dimension_sum_of_absolute_differences/prod(other_dims)
    }
    total_comparison <- entrywise_comparison + sum(dimension_comparisons)
    return(total_comparison)
}

binContinuousListsAsDiscrete <- function(list_a, list_b) {
    a_length <- list_a %>% length
    b_length <- list_b %>% length
    bin_count <- min(a_length, b_length) %>% sqrt %>% ceiling
    bin_count <- ifelse(bin_count  < 2, 2, bin_count)
    bins <- c(list_a, list_b) %>% cut(breaks=bin_count, labels=1:bin_count)
    table_a <- bins[1:a_length] %>% table %>% unname %>% as.vector
    table_b <- bins[-(1:a_length)] %>% table %>% unname %>% as.vector
    return(list(table_a, table_b))
}

getContinuousJSDivergence <- function(sample_1, sample_2) {
    m <- function(x) {
        result <- 0.5*(p(x) + q(x))
        return(result)
    }

    integrand <- function(f, g) {
        func <- function(x) {
            result <- ifelse(f(x) == 0, 0, ifelse(g(x) == 0, Inf, 
                                                  f(x)*(log(f(x)) - log(g(x)))))
            return(result)
        }
        return(func)
    }

    p <- sample_1 %>% density %>% approxfun
    q <- sample_2 %>% density %>% approxfun
    lower <- max(min(sample_1), min(sample_2))
    upper <- min(max(sample_1), max(sample_2))
    KL_div_1 <- integrate(integrand(p, m), lower, upper)$value
    KL_div_2 <- integrate(integrand(q, m), lower, upper)$value
    JS_divergence <- 0.5*(KL_div_1 + KL_div_2)
    return(JS_divergence)
}

getJSDivergence <- function(list_a, list_b, continuous=FALSE) {
    if(continuous) {
        binned <- binContinuousListsAsDiscrete(list_a, list_b)
        divergence <- CalcJSDivergence(binned[[1]], binned[[2]])
    } else {
        max_val <- max(list_a, list_b)
        table_a <- list_a %>% factor(levels=0:max_val) %>% table %>% as.vector
        table_b <- list_b %>% factor(levels=0:max_val) %>% table %>% as.vector
        divergence <- CalcJSDivergence(table_a, table_b)
    }
    return(divergence)
}

removeEmptyStrings <- function(l) {
    return(l[l != ""])
}

standardizeList <- function(l) {
    new_list <- l %>% 
                sapply(toString) %>%
                sapply(paste, collapse='') %>% unname
    return(new_list)
}

determineComparisonMethod <- function(sequence_list) {
    length_count <- sequence_list %>% standardizeList %>% sapply(nchar) %>% 
        table %>% length 
    comparison_method <- ifelse(length_count > 1, "lv", "hamming")
    return(comparison_method)
}

getDistanceMatrix <- function(raw_sequences) {
    sequence_list <- raw_sequences %>% standardizeList
    comparison_method <- sequence_list %>% determineComparisonMethod
    mat <- sequence_list %>% stringdistmatrix(method=comparison_method) %>% as.matrix
    return(mat)
}

getDistanceVector <- function(sequence_list) {
    mat <- sequence_list %>% getDistanceMatrix
    vec <- mat[mat %>% lower.tri] %>% as.vector %>% sort
    return(vec)
}

comparePairwiseDistanceDistributions <- function(list_a, list_b) {    
    distances_a <- list_a %>% getDistanceVector
    distances_b <- list_b %>% getDistanceVector
    divergence <- getJSDivergence(distances_a, distances_b)
    return(divergence)
}

getNearestNeighborDistances <- function(sequence_list, k=1) {
    mat <- sequence_list %>% getDistanceMatrix
    n <- sequence_list %>% length
    distances <- rep(NA, n)
    for(i in 1:n) {
        distances[i] <- sort(mat[i, -i], partial=k)[k]
    }
    return(distances)
}

compareNNDistanceDistribution <- function(list_a, list_b, k=1) {
    distances_a <- list_a %>% getNearestNeighborDistances(k=k)
    distances_b <- list_b %>% getNearestNeighborDistances(k=k)
    divergence <- getJSDivergence(distances_a, distances_b)
    return(divergence)
}

getGCDistribution <- function(raw_sequences) {
    sequence_list <- raw_sequences %>% sapply(paste, collapse='') %>% unname
    dna_list <- sequence_list %>% strsplit(split='') %>% lapply(as.DNAbin)
    gc_dist <- dna_list %>% sapply(GC.content)
    return(gc_dist)
}

compareGCDistributions <- function(list_a, list_b) {
    density_a <- list_a %>% getGCDistribution
    density_b <- list_b %>% getGCDistribution
    divergence <- getJSDivergence(density_a, density_b, continuous=TRUE)
    return(divergence)
}

getMotifCount <- function(motif, subject) {
    dna_strings <- subject %>% unlist %>% DNAStringSet
    count <- motif %>% vcountPattern(dna_strings, fixed=FALSE) %>% sum
    return(count)
}

getHotspotCount <- function(dna_sequence) {
    hotspots <- c("WRC", "WA")
    count <- hotspots %>% sapply(getMotifCount, subject=dna_sequence) %>% sum
    return(count)
}

getColdspotCount <- function(dna_sequence) {
    coldspot <- "SYC"
    count <- coldspot %>% getMotifCount(dna_sequence)
    return(count)
}

getNucleotideDiversity <- function(repertoire) {
    diversity <- repertoire %>% sapply(strsplit, split='') %>% as.DNAbin %>% nuc.div
    return(diversity)
}

compareNucleotideDiversities <- function(repertoire_a, repertoire_b) {
    nuc_div_a <- getNucleotideDiversity(repertoire_a)
    nuc_div_b <- getNucleotideDiversity(repertoire_b)    
    distance <- (nuc_div_b - nuc_div_a) %>% abs
    return(distance)
}

getDistancesFromNaiveToMature <- function(dt) {
    distances <- dt$mature_seq %>% 
        mapply(FUN=stringdist, b=dt$naive_seq, method="lv") %>% 
        sort %>% unname
    return(distances)
}

compareDistancesFromNaiveToMature <- function(dt_a, dt_b) {
    distances_a <- getDistancesFromNaiveToMature(dt_a)
    distances_b <- getDistancesFromNaiveToMature(dt_b)
    divergence <- getJSDivergence(distances_a, distances_b)
    return(divergence)
}

getSequenceListFromFasta <- function(filename) {
    sequences <- filename %>% read.fasta %>% lapply(paste, collapse="") %>% 
        unlist %>% unname
    return(sequences)
}

callPartis <- function(action, input_filename, output_filename, output_path, 
                       partis_path, num_procs, cleanup) {
    shell <- Sys.getenv("SHELL")
    script.file <- system.file("run_partis.sh", package="sumrep")
    command <- paste(shell, script.file, 
                     "-p", partis_path, 
                     "-a", action, 
                     "-i", input_filename, 
                     "-o", output_filename, 
                     "-n", num_procs,
                     "-h", file.path(output_path, "params"))
    command %>% system
    partis_dataset <- output_filename %>% fread(stringsAsFactors=TRUE)
    return(partis_dataset)
} 

getMutationInfo <- function(filename) {
    getSubstitutionRate <- function(state) {
        position <- state$name %>% strsplit("_") %>% unlist %>% last
        if(!grepl("[0-9]+", gsub("[^0-9]+", "", position))) {
            # The state in question does not represent an allele-position entry
            # (partis returns various metadata as 'states')
            result <- NA
        } else {
            # Partis returns zero-based position values, so add one
            position <- as.numeric(position) + 1
            germline_base <- state$extras$germline
            pos_emissions <- state$emissions
            highest_base_sub <- pos_emissions$probs %>% which.max %>% names

            # If the germline is not the most frequent base, the rate is unreliable,
            # so treat as missing info
            result <- ifelse(highest_base_sub == germline_base,
                             1 - get(highest_base_sub, pos_emissions$probs),
                             NA)
        }
        names(result) <- position
        return(result)
    }

    result <- list()
    yaml_object <- yaml.load_file(filename)
    gene_info <- yaml_object$name %>% strsplit("_star_") %>% unlist
    result$gene <- gene_info[1]
    result$allele <- gene_info[2]
    states <- yaml_object$states[-(1:2)]

    overall_substitution_rate <- yaml_object$extras$overall_mute_freq
    position_substitution_rates <- states %>% sapply(getSubstitutionRate) %>% subset(!is.na(.))
    substitution_rates <- {}
    substitution_rates$overall_mut_rate <- overall_substitution_rate
    substitution_rates$mut_rate_by_position <- position_substitution_rates
    return(substitution_rates)
}

annotateSequences <- function(input_filename, output_filename="partis_output.csv", 
                              partis_path='partis', num_procs=4, cleanup=TRUE, 
                              do_full_annotation=TRUE, output_path="_output") {
    if(length(list.files("_output")) > 0 && output_path == "_output" && cleanup) {
        stop("_output path already exists. Please remove it, use another path name, or set 'cleanup' to 'FALSE'.")
    }
    output_file <- file.path(output_path, output_filename)
    annotated_data <- callPartis("annotate", input_filename, output_file, output_path,
                                  partis_path, num_procs, cleanup)


    hmm_yaml_filepath <- file.path(output_path, "params/hmm/hmms")
    yaml_files <- hmm_yaml_filepath %>% list.files
    yaml_filepath_and_files <- sapply(yaml_files, 
                                      function(x) { file.path(hmm_yaml_filepath, x) }) 

    mutation_rates <- {}
    for(yaml_file in yaml_files) {
        allele <- yaml_file %>% gsub(pattern=".yaml", replace='') %>%
            gsub(pattern="_star_", replace="\\*")
        mutation_rates[[allele]] <- hmm_yaml_filepath %>% file.path(yaml_file) %>%
            getMutationInfo
    }

    if(do_full_annotation) {
        extended_output_filename <- file.path(output_path, "new_output.csv")
        script_file <- system.file("process_output.py", package="sumrep")
        system(paste("python", script_file, output_file, 
                     extended_output_filename, sep=' '))
        annotated_data <- extended_output_filename %>% fread(stringsAsFactors=TRUE) %>% 
            subset(select=which(!duplicated(names(.))))
    }

    if(cleanup) {
        output_path %>% unlink(recursive=TRUE)
    }

    raw_sequences <- input_filename %>% getSequenceListFromFasta
    annotated_data$mature_seq <- raw_sequences[annotated_data$unique_ids]

    annotation_object <- {}
    annotation_object$annotations <- annotated_data
    annotation_object$mutation_rates <- mutation_rates
    return(annotation_object)
}

partitionSequences <- function(input_filename, output_filename="partis_output.csv", 
                                partis_path='partis', num_procs=4, cleanup=TRUE) {
    partitioned_data <- callPartis("partition", input_filename, output_filename, 
                                    partis_path, num_procs, cleanup)
    return(partitioned_data)
}

getCDR3Lengths <- function(dt) {
    CDR3_lengths <- dt$cdr3_length %>% na.omit
    return(CDR3_lengths)
}

compareCDR3Lengths <- function(dt_a, dt_b) {
    a_lengths <- getCDR3Lengths(dt_a)
    b_lengths <- getCDR3Lengths(dt_b)
    divergence <- getJSDivergence(a_lengths, b_lengths)
    return(divergence)
}

compareGermlineGeneDistributions <- function(data_table_a, data_table_b, gene_type) {
    column_name <- switch(gene_type,
                          V="v_gene",
                          D="d_gene",
                          J="j_gene")
    levels_a <- data_table_a[, column_name, with=FALSE] %>% sapply(as.numeric) 
    levels_b <- data_table_b[, column_name, with=FALSE] %>% sapply(as.numeric) 
    divergence <- getJSDivergence(levels_a, levels_b)
    return(divergence)
}

compareVDJDistributions <- function(dt_a, dt_b) {
    table_a <- table(dt_a$v_gene, dt_a$d_gene, dt_a$j_gene)
    table_b <- table(dt_b$v_gene, dt_b$d_gene, dt_b$j_gene)
    divergence <- compareCategoricalDistributions(table_a, table_b)
    return(divergence)
}

getGRAVYDistribution <- function(sequence_list) {
    dist <- sequence_list %>% removeEmptyStrings %>% 
            standardizeList %>%
            sapply(gravy) %>% unname
    return(dist)
}

compareGRAVYDistributions <- function(list_a, list_b) {
    dist_a <- getGRAVYDistribution(list_a)
    dist_b <- getGRAVYDistribution(list_b)
    divergence <- getJSDivergence(list_a, list_b,
                                    continuous=TRUE)
    return(divergence)
}

parsePythonDictionary <- function(dictionary) {
    parsed <- dictionary %>% gsub(pattern="'", replacement='"') %>% fromJSON
    return(parsed)
}

extractCDR3CodonStartPositions <- function(dictionary) {
    dict_length <- length(dictionary)
    positions <- dictionary %>% sapply(toString) %>% lapply(parsePythonDictionary) %>% 
        sapply(extract, "v") %>% unlist %>% as.numeric

    # partis returns zero-based positions, so add one. 
    positions <- positions + 1
    return(positions)
}

getCDR3s <- function(dt) {
    codon_starts <- dt$codon_positions %>% extractCDR3CodonStartPositions
    codon_ends <- codon_starts + dt$cdr3_length
    collapsed_seqs <- dt$mature_seq %>% sapply(paste, collapse='')
    ordered_seqs <- collapsed_seqs[dt$unique_ids]
    cdr3s <- ordered_seqs %>% substr(codon_starts, codon_ends - 1)
    return(cdr3s)
}
