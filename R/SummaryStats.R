library(alakazam)
library(ape)
library(Biostrings)
library(data.table)
library(dplyr)
library(flexmix)
library(HDMD)
library(jsonlite)
library(magrittr)
library(pegas)
library(Peptides)
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
    bin_count <- ifelse(bin_count < 2, 2, bin_count)
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

getDistancesFromNaiveToMature <- function(dat) {
    distances <- dat$mature_seq %>% 
        mapply(FUN=stringdist, b=dat$naive_seq, method="lv") %>% 
        sort %>% unname
    return(distances)
}

compareDistancesFromNaiveToMature <- function(dat_a, dat_b) {
    distances_a <- getDistancesFromNaiveToMature(dat_a)
    distances_b <- getDistancesFromNaiveToMature(dat_b)
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

    overall_substitution_rate <- yaml_object$extras$per_gene_mute_freq
    position_substitution_rates <- states %>% sapply(getSubstitutionRate) %>% 
        subset(!is.na(.))
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
        annotated_data$naive_seq <- annotated_data$naive_seq %>% sapply(toString) %>%
            tolower
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

getCDR3Lengths <- function(dat) {
    CDR3_lengths <- dat$cdr3_length %>% na.omit
    return(CDR3_lengths)
}

compareCDR3Lengths <- function(dat_a, dat_b) {
    a_lengths <- getCDR3Lengths(dat_a)
    b_lengths <- getCDR3Lengths(dat_b)
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

compareVDJDistributions <- function(dat_a, dat_b) {
    table_a <- table(dat_a$v_gene, dat_a$d_gene, dat_a$j_gene)
    table_b <- table(dat_b$v_gene, dat_b$d_gene, dat_b$j_gene)
    divergence <- compareCategoricalDistributions(table_a, table_b)
    return(divergence)
}

convertDNAToAminoAcids <- function(sequence) {
    aa_list <- sequence %>% sapply(strsplit, '') %>% sapply(translate) %>%
        paste0(collapse='')
    return(aa_list)
}

getKideraFactorsBySequence <- function(sequence) {
    kideraFactors <- sequence %>% convertDNAtoAminoAcids %>% kideraFactors %>% 
        unlist %>% as.list
    return(kideraFactors)
}

getKideraFactors <- function(sequence_list) {
    kidera_factors <- sequence_list %>% sapply(getKideraFactorsBySequence) %>% t
    return(kidera_factors)
}

getHydrophobicityDistribution <- function(sequence_list) {
    hydrophobicity_list <- sequence_list %>% getKideraFactors %>%
        data.table %>% select(KF4) %>% unlist
    return(hydrophobicity_list)
}

compareHydrophobicityDistributions <- function(dat_a, dat_b) {
    dist_a <- dat_a$naive_seq %>% getHydrophobicityDistribution
    dist_b <- dat_b$naive_seq %>% getHydrophobicityDistribution
    divergence <- getJSDivergence(dist_a, dist_b, continuous=TRUE)
    return(divergence)
}

getMeanAtchleyFactorDistribution <- function(sequence_list) {
    atchley_factors <- sequence_list %>% sapply(convertDNAToAminoAcids) %>% 
        FactorTransform %>% sapply(mean, na.rm=TRUE) %>% unname
    return(atchley_factors)
}

compareAtchleyFactorDistributions <- function(dat_a, dat_b) {
    dist_a <- dat_a$naive_seq %>% getMeanAtchleyFactorDistribution
    dist_b <- dat_b$naive_seq %>% getMeanAtchleyFactorDistribution
    divergence <- getJSDivergence(dist_a, dist_b, continuous=TRUE)
    return(divergence)
}

getAliphaticIndexDistribution <- function(sequence_list) {
    a_indices <- sequence_list %>% sapply(convertDNAToAminoAcids) %>% sapply(aIndex)
    return(a_indices)
}

compareAliphaticIndexDistributions <- function(dat_a, dat_b) {
    dist_a <- dat_a$naive_seq %>% getAliphaticIndexDistribution
    dist_b <- dat_b$naive_seq %>% getAliphaticIndexDistribution
    divergence <- getJSDivergence(dist_a, dist_b, continuous=TRUE)
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
    parsed <- dictionary %>% gsub(pattern="'", replacement='"') %>% 
        gsub(pattern="\\(", replacement="\\[") %>% 
        gsub(pattern="\\)", replacement="\\]") %>% fromJSON
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

getCDR3s <- function(dat) {
    codon_starts <- dat$codon_positions %>% extractCDR3CodonStartPositions
    codon_ends <- codon_starts + dat$cdr3_length
    collapsed_seqs <- dat$mature_seq %>% sapply(paste, collapse='')
    cdr3s <- collapsed_seqs %>% substr(codon_starts, codon_ends - 1)
    return(cdr3s)
}

compareCDR3Distributions <- function(dat_a, dat_b, subsample=TRUE, 
                                     subsample_count=10000) {
    subsample_count <- min(nrow(dat_a), nrow(dat_b), subsample_count)
    dist_a <- dat_a[sample(nrow(dat_a), subsample_count), ] %>% getCDR3s %>% 
        getDistanceVector
    dist_b <- dat_b[sample(nrow(dat_b), subsample_count), ] %>% getCDR3s %>%
        getDistanceVector
    divergence <- getJSDivergence(dist_a, dist_b)
}

getDistancesBetweenMutationsBySequence <- function(naive, mature) {
    if(nchar(naive) != nchar(mature)) {
        stop(paste0("nchar(naive) [", nchar(naive), "] != ", "nchar(mature) [", 
                    nchar(mature), "]"))
    }
    naive_char_list <- naive %>% strsplit(split='') %>% unlist
    mature_char_list <- mature %>% strsplit(split='') %>% unlist
    mut_indices_list <- which(naive_char_list != mature_char_list)
    if(length(mut_indices_list) > 1) {
        distances <- (mut_indices_list %>% diff - 1) %>% list
    } else {
        # If there are less than two mutations, this statistic is undefined,
        # so return NA and eventually ignore
        distances <- NA
    }
    return(distances)
}

getDistancesBetweenMutations <- function(naive_list, mature_list) {
    dists <- mapply(function(x, y) {
                        if(nchar(x) == nchar(y)) {
                            res <- getDistancesBetweenMutationsBySequence(x, y)
                        } else {
                            res <- NA
                        }
                    },
                    naive_list,
                    mature_list) %>% unname %>% unlist %>% subset(!is.na(.))
    return(dists)
}

compareDistanceBetweenMutationsDistributions <- function(dat_a, dat_b) {
    dists_a <- getDistancesBetweenMutations(dat_a$naive_seq, dat_a$mature_seq)
    dists_b <- getDistancesBetweenMutations(dat_b$naive_seq, dat_b$mature_seq)
    divergence <- getJSDivergence(dists_a, dists_b)
    return(divergence)
}

getPerGeneMutationRates <- function(dat) {
    rates <- dat$mutation_rates %>% sapply( function(gene) { 
                                                gene$overall_mut_rate } )
    return(rates)
}

comparePerGeneMutationRates <- function(dat_a, dat_b) {
    rates_a <- dat_a %>% getPerGeneMutationRates
    rates_b <- dat_b %>% getPerGeneMutationRates
    common_genes <- intersect(rates_a %>% names, rates_b %>% names)
    rates_a_common <- rates_a[names(rates_a) %in% common_genes][common_genes] %>% unname
    rates_b_common <- rates_b[names(rates_b) %in% common_genes][common_genes] %>% unname
    divergence <- (rates_a_common - rates_b_common) %>% abs %>% sum
    return(divergence/length(common_genes)) 
}

getPerGenePerPositionMutationRates <- function(dat) {
    rates <- dat$mutation_rates %>% sapply( function(gene) {
                                                gene$mut_rate_by_position } )
    return(rates)
}

comparePerGenePerPositionMutationRates <- function(dat_a, dat_b) {
    rates_a <- dat_a %>% getPerGenePerPositionMutationRates
    rates_b <- dat_b %>% getPerGenePerPositionMutationRates
    common_genes <- intersect(rates_a %>% names, rates_b %>% names)
    rates_a_common <- rates_a[names(rates_a) %in% common_genes][common_genes] %>% unname
    rates_b_common <- rates_b[names(rates_b) %in% common_genes][common_genes] %>% unname
    divergence <- mapply(function(positions_a, positions_b) { 
                            common_positions <- intersect(positions_a %>% names, 
                                                          positions_b %>% names)
                            a_common <- positions_a[names(positions_a) %in% 
                                                    common_positions][common_positions]
                            b_common <- positions_b[names(positions_b) %in% 
                                                    common_positions][common_positions]
                            abs(a_common - b_common)/length(common_positions)
                         }, 
                         rates_a_common, rates_b_common) %>% unlist %>% sum
    return(divergence/length(common_genes)) 
}

