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
                             germline_dir=NULL,
                             num_procs=100) {
    if(method == "partis") {
        annotations <- getPartisAnnotations(input_filename=filename, 
                                            num_procs=num_procs,
                                            output_path="tmp_output",
                                            cleanup=FALSE,
                                            germline_dir=germline_dir
                                           )
        saveRDS(annotations, outname)

        num_clones <- annotations$annotations$clone %>% 
            unique %>% 
            length

        num_leaves <- nrow(annotations$annotations)/num_clones

        simulation <- getPartisSimulation("tmp_output",
                                          num_events=num_clones,
                                          num_leaves=num_leaves
                                         )
        saveRDS(simulation, outname %>% gsub(pattern='.rds',
                                             replace='-sim.rds'))
        "tmp_output" %>% unlink(recursive=TRUE)
    } else if(method == "igblast") {
        annotations <- getIgBlastAnnotations(filename, 
                                             num_threads=num_procs,
                                             igblast_dir="~/Software/igblast",
                                             changeo_dir="~/.local/bin")
        saveRDS(annotations, outname)
    }
}

igb_germline_dir <- "~/Software/igblast/partis_friendly_bin"

write_partis_annotations <- FALSE
if(write_partis_annotations) {
    writeAnnotations("~/Data/GMC-igh-m1h.fa", 
                     "data/Annotations/p_g1.rds", 
                     "partis")
    writeAnnotations("~/Data/FV-igh-m1h.fa", 
                     "data/Annotations/p_f1.rds", 
                     "partis")
    writeAnnotations("~/Data/FV-igh-m8d.fa", 
                     "data/Annotations/p_f2.rds", 
                     "partis")
}

write_partis_igb_annotations <- FALSE
if(write_partis_igb_annotations) {
    writeAnnotations("~/Data/FV-igh-m1h.fa", 
                     "data/Annotations/pi_f1.rds", 
                     "partis",
                     germline_dir=igb_germline_dir)
    writeAnnotations("~/Data/FV-igh-m8d.fa", 
                     "data/Annotations/pi_f2.rds", 
                     "partis",
                     germline_dir=igb_germline_dir)
    writeAnnotations("~/Data/GMC-igh-m1h.fa", 
                     "data/Annotations/pi_g1.rds", 
                     "partis",
                     germline_dir=igb_germline_dir)
}

write_igb_annotations <- TRUE
if(write_igb_annotations) {
    writeAnnotations("~/Data/FV-igh-m1h.fa", 
                     "data/Annotations/i_f1.rds", 
                     "igblast")
    writeAnnotations("~/Data/FV-igh-m8d.fa", 
                     "data/Annotations/i_f2.rds", 
                     "igblast")
    writeAnnotations("~/Data/GMC-igh-m1h.fa", 
                     "data/Annotations/i_g1.rds", 
                     "igblast")
}
