loadNewDatasets("data/Annotations")

do_partis_plots <- FALSE
if(do_partis_plots) {
    partis_dats <- list(p_f1$annotations,
                        p_f1_sim$annotations,
                        p_g1$annotations,
                        p_g1_sim$annotations
                       )
    partis_names <- c("p_f1", "p_f1_sim", "p_g1", "p_g1_sim")
    
    p <- plotUnivariateDistributions(partis_dats,
                                     plot_type="freqpoly",
                                     locus="igh",
                                     names=partis_names
                                    )
    ggsave("~/Manuscripts/sumrep-ms/Figures/partis_freqpoly.pdf", width=14, height=14)
    
    p <- plotUnivariateDistributions(partis_dats,
                                     plot_type="ecdf",
                                     locus="igh",
                                     names=partis_names
                                    )
    ggsave("~/Manuscripts/sumrep-ms/Figures/partis_ecdf.pdf", width=14, height=14)
}

do_igor_plots <- TRUE
if(do_igor_plots) {
    igor_dats <- list(A5_S9$annotations,
                      A5_S9_sim$annotations,
                      A5_S22$annotations,
                      A5_S22_sim$annotations
                     )
    igor_names <- c("A5_S9", "A5_S9_sim",
                    "A5_S22", "A5_S22_sim")
    p <- plotUnivariateDistributions(igor_dats,
                                     plot_type="freqpoly",
                                     locus="trb",
                                     names=igor_names
                                    )
    ggsave("~/Manuscripts/sumrep-ms/Figures/igor_freqpoly.pdf", width=14, height=14)
    
    p <- plotUnivariateDistributions(igor_dats,
                                     plot_type="ecdf",
                                     locus="trb",
                                     names=igor_names
                                    )
    ggsave("~/Manuscripts/sumrep-ms/Figures/igor_ecdf.pdf", width=14, height=14)
}
