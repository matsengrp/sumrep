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
                             dat_name,
                             method,
                             outname=paste0("data/Annotations/",
                                            dat_name,
                                            ".rds"),
                             germline_dir=NULL,
                             num_procs=100,
                             locus=NULL
                            ) {
    if(method == "partis") {
        output_path <- paste0("_output_", dat_name)
        annotations <- getPartisAnnotations(input_filename=filename, 
                                            num_procs=num_procs,
                                            output_path=output_path,
                                            cleanup=FALSE,
                                            germline_dir=germline_dir,
                                            locus="igh"
                                           )
        saveRDS(annotations, outname)

        num_clones <- annotations$annotations$clone %>% 
            unique %>% 
            length

        num_leaves <- nrow(annotations$annotations)/num_clones


        # Get yaml simulation output for cft
        yaml_command <- paste(Sys.getenv("PARTIS_PATH"),
              "simulate",
              "--parameter-dir",
              file.path(output_path, "params"),
              "--n-sim-events",
              num_clones,
              "--n-leaves",
              num_leaves,
              "--outfname",
              file.path(output_path, "simu.yaml"),
              "--seed",
              13
             )
        yaml_command %>% system

        simulation <- getPartisSimulation(parameter_dir=output_path,
                                          num_events=num_clones,
                                          num_leaves=num_leaves,
                                          cleanup=F,
                                          seed=13
                                         )

        saveRDS(simulation, outname %>% gsub(pattern='.rds',
                                            replace='-sim.rds'))
        "tmp_output" %>% unlink(recursive=TRUE)
    } else if(method == "igblast") {
        annotations <- getIgBlastAnnotations(input_filename=filename, 
                                             num_threads=num_procs,
                                             igblast_dir="~/Software/igblast",
                                             changeo_dir="~/.local/bin",
                                             receptor_type="BCR"
                                            )
        saveRDS(annotations, outname)
    } else if(method == "igor") {
        ann_sim <- getIgorAnnotations(input_filename=filename,
                                      output_filename=paste0(dat_name, '.csv'),
                                      igor_wd=dat_name,
                                      locus=locus
                                     )
        annotations <- ann_sim$annotations
        saveRDS(list(annotations=annotations), outname)
        simulation <- ann_sim$simulation
        saveRDS(list(annotations=simulation), outname %>% gsub(pattern='.rds',
                                             replacement='-sim.rds')
        )
        
    }
}

igb_germline_dir <- "~/Software/igblast/partis_friendly_bin"

write_partis_annotations <- FALSE
if(write_partis_annotations) {
    writeAnnotations("~/Data/GMC-igh-m1h.fa", 
                      "p_g1",
                     "partis")
    writeAnnotations("~/Data/GMC-igh-m8d.fa", 
                     "p_g2",
                     "partis")
    writeAnnotations("~/Data/IB-igh-m1h.fa", 
                     "p_i1",
                     "partis")
    writeAnnotations("~/Data/IB-igh-m8d.fa", 
                     "p_i2",
                     "partis")
    writeAnnotations("~/Data/FV-igh-m1h.fa", 
                     "p_f1",
                     "partis")
    writeAnnotations("~/Data/FV-igh-m8d.fa", 
                     "p_f2",
                     "partis")
}

write_partis_igb_annotations <- FALSE
if(write_partis_igb_annotations) {
    writeAnnotations("~/Data/FV-igh-m1h.fa", 
                     "pi_f1",
                     "partis",
                     germline_dir=igb_germline_dir)
    writeAnnotations("~/Data/FV-igh-m8d.fa", 
                     "pi_f2", 
                     "partis",
                     germline_dir=igb_germline_dir)
    writeAnnotations("~/Data/GMC-igh-m1h.fa", 
                     "pi_g1", 
                     "partis",
                     germline_dir=igb_germline_dir)
    writeAnnotations("~/Data/GMC-igh-m8d.fa", 
                     "pi_g2", 
                     "partis",
                     germline_dir=igb_germline_dir)
    writeAnnotations("~/Data/IB-igh-m1h.fa", 
                     "pi_i1", 
                     "partis",
                     germline_dir=igb_germline_dir)
    writeAnnotations("~/Data/IB-igh-m8d.fa", 
                     "pi_i2", 
                     "partis",
                     germline_dir=igb_germline_dir)
}

write_igb_annotations <- FALSE
if(write_igb_annotations) {
    writeAnnotations("~/Data/FV-igh-m1h.fa", 
                     "i_f1", 
                     "igblast")
    writeAnnotations("~/Data/FV-igh-m8d.fa", 
                     "i_f2", 
                     "igblast")
    writeAnnotations("~/Data/GMC-igh-m1h.fa", 
                     "i_g1", 
                     "igblast")
    writeAnnotations("~/Data/GMC-igh-m8d.fa", 
                     "i_g2", 
                     "igblast")
    writeAnnotations("~/Data/IB-igh-m1h.fa", 
                     "i_i1", 
                     "igblast")
    writeAnnotations("~/Data/IB-igh-m8d.fa", 
                     "i_i2", 
                     "igblast")
}

write_igor_annotations <- TRUE
if(write_igor_annotations) {
    tcr_dir <- "/fh/fast/matsen_e/bolson2/mike_aging/converted"
    writeAnnotations(file.path(tcr_dir, "A5-S22_R2.fa"),
                     "A5_S22",
                     "igor",
                     locus="trb"
                    )
    writeAnnotations(file.path(tcr_dir, "A5-S9_R2.fa"),
                     "A5_S9",
                     "igor",
                     locus="trb"
                    )
    writeAnnotations(file.path(tcr_dir, "A5-S10_R2.fa"),
                     "A5_S10",
                     "igor",
                     locus="trb"
                    )
    writeAnnotations(file.path(tcr_dir, "A5-S15_R2.fa"),
                     "A5_S15",
                     "igor",
                     locus="trb"
                    )
    writeAnnotations(file.path(tcr_dir, "A4-i194_R2.fa"),
                     "A4_i194",
                     "igor",
                     locus="trb"
                    )
    writeAnnotations(file.path(tcr_dir, "A4-i107_R2.fa"),
                     "A4_i107",
                     "igor",
                     locus="trb"
                    )
}
