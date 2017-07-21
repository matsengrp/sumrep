library(alakazam)
library(ape)
library(Biostrings)
library(dplyr)
library(flexmix)
library(data.table)
library(pegas)
library(RecordLinkage)
library(shazam)
library(stringdist)
library(textmineR)

compare.categorical.distributions <- function(a, b) {
    dims <- dim(a)
    entry.compare <- (a - b) %>% abs %>% sum
    entry.compare <- entry.compare/(dims %>% prod)
    dimension <- dims %>% length
    dimension.diffs <- rep(NA, dimension)
    for(i in 1:dimension) {
        dimension.sums.a <- apply(a, i, sum)
        dimension.sums.b <- apply(b, i, sum)
        dimension.diff <- (dimension.sums.a - dimension.sums.b) %>% abs %>% sum
        dimension.length <- dims[i]
        dimension.diffs[[i]] <- dimension.diff/dimension.length
    }
    total.compare <- entry.compare + sum(dimension.diffs)
    return(total.compare)
}

bin.continuous.lists.as.discrete <- function(list.a, list.b) {
    a.length <- list.a %>% length
    b.length <- list.b %>% length
    bin.count <- min(a.length, b.length) %>% sqrt %>% ceiling
    bin.count <- ifelse(bin.count  < 2, 2, bin.count)
    bins <- c(list.a, list.b) %>% cut(breaks=bin.count, labels=1:bin.count)
    table.a <- bins[1:a.length] %>% table %>% unname %>% as.vector
    table.b <- bins[-(1:a.length)] %>% table %>% unname %>% as.vector
    return(list(table.a, table.b))
}

get.continuous.JS.divergence <- function(sample.1, sample.2) {
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

get.JS.divergence <- function(list.a, list.b, continuous=FALSE) {
    if(continuous) {
        binned <- bin.continuous.lists.as.discrete(list.a, list.b)
        divergence <- get.continuous.JS.divergence(binned[[1]], binned[[2]])
    } else {
        max.val <- max(list.a, list.b)
        table.a <- list.a %>% factor(levels=0:max.val) %>% table %>% as.vector
        table.b <- list.b %>% factor(levels=0:max.val) %>% table %>% as.vector
        divergence <- CalcJSDivergence(table.a, table.b)
    }
    return(divergence)
}

standardize.list <- function(l) {
    new.list <- l %>% sapply(paste, collapse='') %>% unname
    return(new.list)
}

determine.comparison.method <- function(sequences) {
    length.count <- sequences %>% standardize.list %>% sapply(nchar) %>% table %>% length 
    comparison.method <- ifelse(length.count > 1, "lv", "hamming")
    return(comparison.method)
}

get.distance.matrix <- function(raw.sequences) {
    sequence.list <- raw.sequences %>% standardize.list
    comparison.method <- sequence.list %>% determine.comparison.method
    mat <- sequence.list %>% stringdistmatrix(method=comparison.method) %>% as.matrix
    return(mat)
}

get.distance.vector <- function(sequence.list) {
    mat <- sequence.list %>% get.distance.matrix
    vec <- mat[mat %>% lower.tri] %>% as.vector %>% sort
    return(vec)
}

compare.pairwise.distance.distribution <- function(list.a, list.b) {    
    distances.a <- list.a %>% get.distance.vector
    distances.b <- list.b %>% get.distance.vector
    divergence <- get.JS.divergence(distances.a, distances.b)
    return(divergence)
}

get.nearest.neighbor.distances <- function(sequence.list, k=1) {
    mat <- sequence.list %>% get.distance.matrix
    n <- sequence.list %>% length
    distances <- rep(NA, n)
    for(i in 1:n) {
        distances[i] <- sort(mat[i, -i], partial=k)[k]
    }
    return(distances)
}

compare.NN.distance.distribution <- function(list.a, list.b, k=1) {
    distances.a <- list.a %>% get.nearest.neighbor.distances(k=k)
    distances.b <- list.b %>% get.nearest.neighbor.distances(k=k)
    divergence <- get.JS.divergence(distances.a, distances.b)
    return(divergence)
}

get.GC.distribution <- function(raw.sequences) {
    sequence.list <- raw.sequences %>% sapply(paste, collapse='') %>% unname
    dna.list <- sequence.list %>% strsplit(split='') %>% lapply(as.DNAbin)
    gc.dist <- dna.list %>% sapply(GC.content)
    return(gc.dist)
}

compare.GC.distributions <- function(list.a, list.b) {
    density.a <- list.a %>% get.GC.distribution
    density.b <- list.b %>% get.GC.distribution
    divergence <- get.JS.divergence(density.a, density.b, continuous=TRUE)
    return(divergence)
}

get.motif.count <- function(motif, subject) {
    dna.strings <- subject %>% unlist %>% DNAStringSet
    count <- motif %>% vcountPattern(dna.strings, fixed=FALSE) %>% sum
    return(count)
}

get.hotspot.count <- function(dna.sequence) {
    hotspots <- c("WRC", "WA")
    count <- hotspots %>% sapply(get.motif.count, subject=dna.sequence) %>% sum
    return(count)
}

get.coldspot.count <- function(dna.sequence) {
    coldspot <- "SYC"
    count <- coldspot %>% get.motif.count(dna.sequence)
    return(count)
}

get.nucleotide.diversity <- function(repertoire) {
    diversity <- repertoire %>% sapply(strsplit, split='') %>% as.DNAbin %>% nuc.div
    return(diversity)
}

compare.nucleotide.diversities <- function(rep.a, rep.b) {
    nuc.div.a <- get.nucleotide.diversity(rep.a)
    nuc.div.b <- get.nucleotide.diversity(rep.b)    
    distance <- (nuc.div.b - nuc.div.a) %>% abs
    return(distance)
}

get.distances.from.naive.to.mature <- function(naive, mature.list) {
    comparison.method <- mature.list %>% standardize.list %>% determine.comparison.method
    distances <- mature.list %>% 
                 sapply(stringdist, b=naive, method=comparison.method) %>% sort
    return(distances)
}

compare.distances.from.naive.to.mature <- function(naive.a, mature.list.a, 
                                                   naive.b, mature.list.b) {
    distances.a <- get.distances.from.naive.to.mature(naive.a, mature.list.a)
    distances.b <- get.distances.from.naive.to.mature(naive.b, mature.list.b)
    divergence <- get.JS.divergence(distances.a, distances.b)
    return(divergence)
}

call.partis <- function(action, input.filename, output.filename, partis.path, num.procs, 
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
    return(partis.dataset)
} 

annotate.sequences <- function(input.filename, output.filename="partis_output.csv", 
                               partis.path='partis', num.procs=4, cleanup=TRUE) {
    annotated.data <- call.partis("annotate", input.filename, output.filename, 
                                  partis.path, num.procs, cleanup)
    return(annotated.data)
}

partition.sequences <- function(input.filename, output.filename="partis_output.csv", 
                                partis.path='partis', num.procs=4, cleanup=TRUE) {
    partitioned.data <- call.partis("partition", input.filename, output.filename, 
                                    partis.path, num.procs, cleanup)
    return(partitioned.data)
}

get.CDR3.lengths <- function(dt) {
    CDR3.lengths <- dt$cdr3_length %>% na.omit
    return(CDR3.lengths)
}

compare.CDR3.lengths <- function(dt.a, dt.b) {
    a.lengths <- get.CDR3.lengths(dt.a)
    b.lengths <- get.CDR3.lengths(dt.b)
    divergence <- get.JS.divergence(a.lengths, b.lengths)
    return(divergence)
}

compare.germline.gene.distributions <- function(data.table.a, data.table.b, gene.type) {
    column.name <- switch(gene.type,
                          V="v_gene",
                          D="d_gene",
                          J="j_gene")
    levels.a <- data.table.a[, column.name, with=FALSE] %>% sapply(as.numeric) 
    levels.b <- data.table.b[, column.name, with=FALSE] %>% sapply(as.numeric) 
    divergence <- get.JS.divergence(levels.a, levels.b)
    return(divergence)
}

compare.vdj.distributions <- function(dt.a, dt.b) {
    table.a <- table(dt.a$v_gene, dt.a$d_gene, dt.a$j_gene)
    table.b <- table(dt.b$v_gene, dt.b$d_gene, dt.b$j_gene)
    divergence <- compare.categorical.distributions(table.a, table.b)
    return(divergence)
}
