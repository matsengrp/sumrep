library(viridis)

plotDistribution <- function(function_string,
                             datasets,
                             column_name,
                             dat_names,
                             group_names,
                             filename=paste0(function_string, ".pdf")) {

    summary_function <- get(function_string)
    d_full <- {}
    num_datasets <- length(datasets)
    for(i in 1:num_datasets) {
        density_i <- datasets[[i]][[column_name]] %>%
            summary_function
        d_full <- rbind(d_full,
                        data.table(Value=density_i,
                                   Dataset=dat_names[[i]],
                                   Group=group_names[[i]])
                        )
    }

    ggplot(d_full, aes(x=Value, colour=Dataset, lty=Group)) + 
        geom_density() +
        scale_color_viridis(discrete=TRUE) +
        ggtitle(paste("Densities of", function_string))

    ggsave(filename, width=10, height=6)
}

dats <- list(p_f1$annotations,
             p_f2$annotations,
             p_g1$annotations,
             p_f1_sim$annotations,
             p_f2_sim$annotations,
             p_g1_sim$annotations)
groups <- c(rep("Observed", 3), rep("Simulated", 3))

dat_names=c("Individual 1, time 1",
            "Individual 1, time 2",
            "Individual 2, time 1",
            "Individual 1, time 1",
            "Individual 1, time 2",
            "Individual 2, time 1")

getUnivariateSummaries <- function(dat) {
    gc_content <- dat %$%
        sequence %>%
        getGCContentDistribution

    cdr3_lengths <- dat %>%
        getCDR3Lengths

    hydrophobicities <- dat %$%
        naive_seq %>%
        getHydrophobicityDistribution(include_NA=TRUE)

    summary_dat <- data.table(
                              GCContent=gc_content,
                              CDR3Length=cdr3_lengths,
                              Hydrophobicity=hydrophobicities
                             )

    return(summary_dat)
}

loadNewDatasets("data/Annotations")
summary_dat_f1 <- p_f1$annotations %>% getUnivariateSummaries

stop()

plotDistribution(function_string="getGCContentDistribution",
                 datasets=dats,
                 column_name="sequence",
                 dat_names=c("Individual 1, time 1",
                             "Individual 1, time 2",
                             "Individual 2, time 1",
                             "Individual 1, time 1",
                             "Individual 1, time 2",
                             "Individual 2, time 1"),
                 group_names=groups
                )

                 
plotDistribution(function_string="getHydrophobicityDistribution",
                 datasets=dats,
                 column_name="cdr3s",
                 dat_names=c("Individual 1, time 1",
                             "Individual 1, time 2",
                             "Individual 2, time 1",
                             "Individual 1, time 1",
                             "Individual 1, time 2",
                             "Individual 2, time 1"),
                 group_names=groups
                )

plotDistribution(function_string="getGRAVYDistribution",
                 datasets=dats,
                 column_name="cdr3s",
                 dat_names=c("Individual 1, time 1",
                             "Individual 1, time 2",
                             "Individual 2, time 1",
                             "Individual 1, time 1",
                             "Individual 1, time 2",
                             "Individual 2, time 1"),
                 group_names=groups
                )
