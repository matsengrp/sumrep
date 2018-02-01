#' Get and write annotations from fasta as RDS objects
#'
#' Currently supports partis and igblast annotations
#' For partis, a simulation is also written by default
#'
#' @param filename Full path of fasta file
#' @param outname Desired path for rds object
#' @param method Annotation method, currently either 'partis' or 'igblast'
#' @param num_procs Number of cores for annotation
writeAnnotations <- function(filename, 
                             outname,
                             method,
                             num_procs=100) {
    if(method == "partis") {
        annotations <- annotateSequences(filename, 
                                         num_procs=num_procs,
                                         output_path="tmp_output",
                                         cleanup=FALSE)
        saveRDS(annotations, outname)
        simulation <- simulateDataset("tmp_output")
        saveRDS(simulation, outname %>% gsub(pattern='.rds',
                                             replace='-sim.rds'))
        "tmp_output" %>% unlink
    } else if(method == "igblast") {
        annotations <- getIgBlastAnnotations(filename, num_threads=num_procs)
        saveRDS(annotations, outname)
    }
}

if(FALSE) {
writeAnnotations("~/Data/FV-igh-m1h.fa", "data/Annotations/igb_fv1.rds", "igblast")
writeAnnotations("~/Data/FV-igh-m8d.fa", "data/Annotations/igb_fv2.rds", "igblast")
writeAnnotations("~/Data/GMC-igh-m1h.fa", "data/Annotations/igb_gmc1.rds", "igblast")
}

writeAnnotations("~/Data/FV-igh-m1h.fa", "data/Annotations/partis_fv1.rds", "partis")
writeAnnotations("~/Data/FV-igh-m8d.fa", "data/Annotations/partis_fv2.rds", "partis")
writeAnnotations("~/Data/GMC-igh-m1h.fa", "data/Annotations/partis_gmc1.rds", "partis")
