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

        num_leaves <- 4

        simulation <- getPartisSimulation(parameter_dir=output_path,
                                          num_events=num_clones,
                                          num_leaves=num_leaves,
                                          cleanup=F,
                                          seed=13,
                                          subsample_to_unique_clones=TRUE
                                         )

        saveRDS(simulation, outname %>% gsub(pattern='.rds',
                                            replace='-sim.rds'))
    } else if(method == "igblast") {
        annotations <- getIgBlastAnnotations(input_filename=filename, 
                                             nproc=16,
                                             locus="igh"
                                            )
        saveRDS(annotations, outname)
    } else if(method == "igor") {
        ann_sim <- getIgorAnnotations(input_filename=filename,
                                      output_filename=paste0(dat_name, '.csv'),
                                      igor_wd=dat_name,
                                      locus=locus,
                                      cleanup=FALSE
                                     )
        annotations <- ann_sim$annotations
        saveRDS(list(annotations=annotations), outname)
        simulation <- ann_sim$simulation
        saveRDS(list(annotations=simulation), outname %>% gsub(pattern='.rds',
                                             replacement='-sim.rds')
        )
        
    }
}
