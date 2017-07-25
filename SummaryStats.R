library(alakazam)
library(ape)
library(Biostrings)
library(data.table)
library(dplyr)
library(flexmix)
library(pegas)
library(RecordLinkage)
library(shazam)
library(stringdist)
library(textmineR)

compareCategoricalDistributions <- function(a, b) {
    dims <- dim(a)
    entry.compare <- (a - b) %>% abs %>% sum
    entry.compare <- entry.compare/prod(dims)
    dimension.count <- dims %>% length
    dimension.diffs <- rep(NA, dimension.count)
    for(i in 1:dimension.count) {
        dimension.sums.a <- apply(a, i, sum)
        dimension.sums.b <- apply(b, i, sum)
        dimension.diff <- (dimension.sums.a - dimension.sums.b) %>% abs %>% sum
        dimension.lengths <- dims[-i]
        dimension.diffs[[i]] <- dimension.diff/prod(dimension.lengths)
    }
    total.compare <- entry.compare + sum(dimension.diffs)
    return(total.compare)
}

binContinuousListsAsDiscrete <- function(list.a, list.b) {
    a.length <- list.a %>% length
    b.length <- list.b %>% length
    bin.count <- min(a.length, b.length) %>% sqrt %>% ceiling
    bin.count <- ifelse(bin.count  < 2, 2, bin.count)
    bins <- c(list.a, list.b) %>% cut(breaks=bin.count, labels=1:bin.count)
    table.a <- bins[1:a.length] %>% table %>% unname %>% as.vector
    table.b <- bins[-(1:a.length)] %>% table %>% unname %>% as.vector
    return(list(table.a, table.b))
}

getContinuousJSDivergence <- function(sample.1, sample.2) {
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

    p <- sample.1 %>% density %>% approxfun
    q <- sample.2 %>% density %>% approxfun
    lower <- max(min(sample.1), min(sample.2))
    upper <- min(max(sample.1), max(sample.2))
    KL.div.1 <- integrate(integrand(p, m), lower, upper)$value
    KL.div.2 <- integrate(integrand(q, m), lower, upper)$value
    JS.divergence <- 0.5*(KL.div.1 + KL.div.2)
    return(JS.divergence)
}

getJSDivergence <- function(list.a, list.b, continuous=FALSE) {
    if(continuous) {
        binned <- binContinuousListsAsDiscrete(list.a, list.b)
        divergence <- CalcJSDivergence(binned[[1]], binned[[2]])
    } else {
        max.val <- max(list.a, list.b)
        table.a <- list.a %>% factor(levels=0:max.val) %>% table %>% as.vector
        table.b <- list.b %>% factor(levels=0:max.val) %>% table %>% as.vector
        divergence <- CalcJSDivergence(table.a, table.b)
    }
    return(divergence)
}

removeEmptyStrings <- function(l) {
    return(l[l != ""])
}

standardizeList <- function(l) {
    new.list <- l %>% 
                sapply(toString) %>%
                sapply(paste, collapse='') %>% unname
    return(new.list)
}

determineComparisonMethod <- function(sequences) {
    length.count <- sequences %>% standardizeList %>% sapply(nchar) %>% table %>% length 
    comparison.method <- ifelse(length.count > 1, "lv", "hamming")
    return(comparison.method)
}

getDistanceMatrix <- function(raw.sequences) {
    sequence.list <- raw.sequences %>% standardizeList
    comparison.method <- sequence.list %>% determineComparisonMethod
    mat <- sequence.list %>% stringdistmatrix(method=comparison.method) %>% as.matrix
    return(mat)
}

getDistanceVector <- function(sequence.list) {
    mat <- sequence.list %>% getDistanceMatrix
    vec <- mat[mat %>% lower.tri] %>% as.vector %>% sort
    return(vec)
}

comparePairwiseDistanceDistributions <- function(list.a, list.b) {    
    distances.a <- list.a %>% getDistanceVector
    distances.b <- list.b %>% getDistanceVector
    divergence <- getJSDivergence(distances.a, distances.b)
    return(divergence)
}

getNearestNeighborDistances <- function(sequence.list, k=1) {
    mat <- sequence.list %>% getDistanceMatrix
    n <- sequence.list %>% length
    distances <- rep(NA, n)
    for(i in 1:n) {
        distances[i] <- sort(mat[i, -i], partial=k)[k]
    }
    return(distances)
}

compareNNDistanceDistribution <- function(list.a, list.b, k=1) {
    distances.a <- list.a %>% getNearestNeighborDistances(k=k)
    distances.b <- list.b %>% getNearestNeighborDistances(k=k)
    divergence <- getJSDivergence(distances.a, distances.b)
    return(divergence)
}

getGCDistribution <- function(raw.sequences) {
    sequence.list <- raw.sequences %>% sapply(paste, collapse='') %>% unname
    dna.list <- sequence.list %>% strsplit(split='') %>% lapply(as.DNAbin)
    gc.dist <- dna.list %>% sapply(GC.content)
    return(gc.dist)
}

compareGCDistributions <- function(list.a, list.b) {
    density.a <- list.a %>% getGCDistribution
    density.b <- list.b %>% getGCDistribution
    divergence <- getJSDivergence(density.a, density.b, continuous=TRUE)
    return(divergence)
}

getMotifCount <- function(motif, subject) {
    dna.strings <- subject %>% unlist %>% DNAStringSet
    count <- motif %>% vcountPattern(dna.strings, fixed=FALSE) %>% sum
    return(count)
}

getHotspotCount <- function(dna.sequence) {
    hotspots <- c("WRC", "WA")
    count <- hotspots %>% sapply(getMotifCount, subject=dna.sequence) %>% sum
    return(count)
}

getColdspotCount <- function(dna.sequence) {
    coldspot <- "SYC"
    count <- coldspot %>% getMotifCount(dna.sequence)
    return(count)
}

getNucleotideDiversity <- function(repertoire) {
    diversity <- repertoire %>% sapply(strsplit, split='') %>% as.DNAbin %>% nuc.div
    return(diversity)
}

compareNucleotideDiversities <- function(rep.a, rep.b) {
    nuc.div.a <- getNucleotideDiversity(rep.a)
    nuc.div.b <- getNucleotideDiversity(rep.b)    
    distance <- (nuc.div.b - nuc.div.a) %>% abs
    return(distance)
}

getDistancesFromNaiveToMature <- function(naive, mature.list) {
    comparison.method <- mature.list %>% standardizeList %>% determineComparisonMethod
    distances <- mature.list %>% 
                 sapply(stringdist, b=naive, method=comparison.method) %>% sort
    return(distances)
}

compareDistancesFromNaiveToMature <- function(naive.a, mature.list.a, 
                                                   naive.b, mature.list.b) {
    distances.a <- getDistancesFromNaiveToMature(naive.a, mature.list.a)
    distances.b <- getDistancesFromNaiveToMature(naive.b, mature.list.b)
    divergence <- getJSDivergence(distances.a, distances.b)
    return(divergence)
}

callPartis <- function(action, input.filename, output.filename, partis.path, num.procs, 
                        cleanup) {
    shell <- Sys.getenv("SHELL")
    command <- paste(shell, "run_partis.sh", 
                     "-p", partis.path, 
                     "-a", action, 
                     "-i", input.filename, 
                     "-o", output.filename, 
                     "-n", num.procs)
    command %>% system
    partis.dataset <- output.filename %>% fread(stringsAsFactors=TRUE)
    return(partis.dataset)
} 

annotateSequences <- function(input.filename, output.filename="partis_output.csv", 
                              partis.path='partis', num.procs=4, cleanup=TRUE, 
                              do.full.annotation=TRUE) {
    annotated.data <- callPartis("annotate", input.filename, output.filename, 
                                  partis.path, num.procs, cleanup)
    if(do.full.annotation) {
        extended.output.filename <- "new_output.csv"
        system(paste("python process_output.py", output.filename, 
                     extended.output.filename, sep=' '))
        annotated.data <- extended.output.filename %>% fread(stringsAsFactors=TRUE)
        
        if(cleanup) {
            extended.output.filename %>% file.remove
        }
    }

    if(cleanup) {
        paste("rm", output.filename, sep=' ') %>% system

        if("_output" %in% list.files()) {
            if(length(list.files("_output")) == 0) {
                "_output" %>% unlink(recursive=TRUE)
            }
        }

        files.to.remove <- dir(path=".", pattern="*cluster-annotations.csv")
        files.to.remove %>% file.remove
    }

    return(annotated.data)
}

partitionSequences <- function(input.filename, output.filename="partis_output.csv", 
                                partis.path='partis', num.procs=4, cleanup=TRUE) {
    partitioned.data <- callPartis("partition", input.filename, output.filename, 
                                    partis.path, num.procs, cleanup)
    return(partitioned.data)
}

getCDR3Lengths <- function(dt) {
    CDR3.lengths <- dt$cdr3_length %>% na.omit
    return(CDR3.lengths)
}

compareCDR3Lengths <- function(dt.a, dt.b) {
    a.lengths <- getCDR3Lengths(dt.a)
    b.lengths <- getCDR3Lengths(dt.b)
    divergence <- getJSDivergence(a.lengths, b.lengths)
    return(divergence)
}

compareGermlineGeneDistributions <- function(data.table.a, data.table.b, gene.type) {
    column.name <- switch(gene.type,
                          V="v_gene",
                          D="d_gene",
                          J="j_gene")
    levels.a <- data.table.a[, column.name, with=FALSE] %>% sapply(as.numeric) 
    levels.b <- data.table.b[, column.name, with=FALSE] %>% sapply(as.numeric) 
    divergence <- getJSDivergence(levels.a, levels.b)
    return(divergence)
}

compareVDJDistributions <- function(dt.a, dt.b) {
    table.a <- table(dt.a$v_gene, dt.a$d_gene, dt.a$j_gene)
    table.b <- table(dt.b$v_gene, dt.b$d_gene, dt.b$j_gene)
    divergence <- compareCategoricalDistributions(table.a, table.b)
    return(divergence)
}

getGRAVYDistribution <- function(sequence.list) {
    dist <- sequence.list %>% removeEmptyStrings %>% 
            standardizeList %>%
            sapply(gravy) %>% unname
    return(dist)
}


compareGRAVYDistributions <- function(list.a, list.b) {
    dist.a <- getGRAVYDistribution(list.a)
    dist.b <- getGRAVYDistribution(list.b)
    divergence <- getJSDivergence(list.a, list.b,
                                    continuous=TRUE)
    return(divergence)
}

extractCDR3CodonStartPositions <- function(dictionary) {
    parsePythonKeyValue <- function(x) {
        splits <- x %>% strsplit(split=",") 
        pair <- gsub("[^0-9]", " ", splits) %>% strsplit("\\s+") %>% unlist %>% 
            removeEmptyStrings
        start <- pair[2]
        return(start)
    }

    dict.length <- length(dictionary)
    positions <- dictionary %>% sapply(toString) %>% sapply(parsePythonKeyValue) %>% 
        unlist %>% as.numeric

    # partis returns zero-based positions, so add one. 
    positions <- positions + 1
    return(positions)
}

getCDR3s <- function(dt, raw.seqs) {
    codon_starts <- dt$codon_positions %>% extractCDR3CodonStartPositions
    codon_ends <- codon_starts + dt$cdr3_length
    collapsed.seqs <- raw.seqs %>% sapply(paste, collapse='')
    ordered.seqs <- collapsed.seqs[dt$unique_ids]
    cdr3s <- ordered.seqs %>% substr(codon_starts, codon_ends - 1)
    return(cdr3s)
}
