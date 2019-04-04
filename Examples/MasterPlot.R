loadNewDatasets("data/Annotations")

dataset_list <- list(p_f1$annotations,
                     p_f1_sim$annotations,
                     p_g1$annotations,
                     p_g1_sim$annotations
                    )
plot_names <- c("p_f1", "p_f1_sim", "p_f2", "p_f2_sim")

p <- plotUnivariateDistributions(dataset_list,
                            plot_type="freqpoly",
                            names=plot_names
                           )
ggsave("~/Manuscripts/sumrep-ms/Figures/master_plot_freqpoly.pdf", width=14, height=14)

p <- plotUnivariateDistributions(dataset_list,
                                 plot_type="ecdf",
                                 names=plot_names
                                )
ggsave("~/Manuscripts/sumrep-ms/Figures/master_plot_ecdf.pdf", width=14, height=14)
