library(viridis)
ggplot2::theme_set(theme_gray(base_size = 14))

getComparisonValue <- function(dat, comparison) {
    value <- dat %>% 
        filter(Comparison == comparison) %$%
        Divergence
    return(value)
}

scoreStatistics <- function(sim_dats, obs_dats) {
    comparison_types <- sim_dats[[1]]$Comparison %>% unique
    score_dat <- matrix(NA, nrow=0, ncol=2) %>% 
        data.table %>%
        setNames(c("Comparison", "Score"))
    for(c_type in comparison_types) {
            sim_score <- sim_dats %>% 
                sapply(getComparisonValue, c_type) %>% 
                unlist %>%
                mean
            obs_score <- obs_dats %>%
                sapply(getComparisonValue, c_type) %>%
                unlist %>%
                mean
            c_score <- sim_score/obs_score
            comparison_name <- c_type %>%
                gsub(pattern="compare", replace="get") %>%
                ifelse(str_sub(., -1) == "s", 
                       substr(., start=1, stop=nchar(.) - 1), 
                       .
                      ) %>%
                getNameFromFunctionString
            
            score_dat <- rbind(score_dat,
                               data.table(
                                          Comparison=comparison_name,
                                          Score=c_score
                                         )
                              )
    }
    return(score_dat)
}

shortenName <- function(string) {
    shortened_name <- string %>% 
        gsub(pattern="compare", replace="get") %>%
        getNameFromFunctionString
    return(shortened_name)
}


plotComparisons <- function(dat, filename, cols=1, rows=1) {
    comparison_types <- dat$Comparison %>% unique
    plot_list <- list()
    for(c_type in comparison_types) {
        dat_sub <- dat[dat$Comparison == c_type, ]
        plot_list[[c_type]] <- ggplot(dat_sub, aes(Type1, Type2)) +
            geom_tile(aes(fill=Divergence)) +
            ggtitle(c_type %>% shortenName) +
            theme(
                  legend.key.size=unit(0.8, "cm"),
                  legend.title=element_text(size=16),
                  legend.text=element_text(size=14),
                  axis.text.x=element_text(size=14),
                  axis.text.y=element_text(size=14),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank(),
                  plot.title=element_text(size=18)
                  ) +
            scale_fill_viridis(direction=-1)
    }

    grid_dims <- comparison_types %>% 
        length %>%
        getGridDims
    
    pdf(filename, width=24, height=18)
    multiplot(plotlist=plot_list, cols=grid_dims[1], rows=grid_dims[2])
    dev.off()

}

loadNewDatasets("data/Comparisons")

part_igb_dat <- rbind(compare_pi_f1_pi_f1_sim,
                      compare_pi_f2_pi_f2_sim,
                      compare_pi_g1_pi_g1_sim,
                      compare_i_f1_pi_f1_sim,
                      compare_i_f2_pi_f2_sim,
                      compare_i_g1_pi_g1_sim,
                      compare_pi_f1_pi_f2,
                      compare_pi_f1_pi_g1,
                      compare_pi_f2_pi_g1,
                      compare_i_f1_i_f2,
                      compare_i_f1_i_g1,
                      compare_i_f2_i_g1
                     )

part_igb_dat <- part_igb_dat[!(part_igb_dat$Comparison %in% 
                               c(
                                 "comparePerGeneMutationRates",
                                 "comparePerGenePerPositionMutationRates"
                                 )), ]

part_igb_dat$Type2 <- factor(part_igb_dat$Type2,
                             levels=c(
                                      "pi_f1_sim",
                                      "pi_f2_sim",
                                      "pi_g1_sim",
                                      "i_f1",
                                      "i_f2",
                                      "i_g1",
                                      "pi_f1",
                                      "pi_f2",
                                      "pi_g1"
                                      ))

obs_sim_dat <- rbind(compare_p_f1_p_f1_sim,
                     compare_p_f2_p_f2_sim,
                     compare_p_g1_p_g1_sim,
                     compare_p_f1_p_f2,
                     compare_p_f1_p_g1,
                     compare_p_f2_p_g1
                    )

obs_sim_dat$Type2 <- factor(obs_sim_dat$Type2, 
                      levels=c("p_f1_sim",
                               "p_f2_sim",
                               "p_g1_sim",
                               "p_f1",
                               "p_f2",
                               "p_g1"))

obs_sim_dats <- list(
                 compare_p_f1_p_f1_sim,
                 compare_p_f2_p_f2_sim,
                 compare_p_g1_p_g1_sim,
                 compare_p_g2_p_g2_sim,
                 compare_p_i1_p_i1_sim,
                 compare_p_i2_p_i2_sim
                )

obs_obs_dats <- list(
                 compare_p_f1_p_f2,
                 compare_p_f1_p_g1,
                 compare_p_f1_p_g2,
                 compare_p_f1_p_i1,
                 compare_p_f1_p_i2,

                 compare_p_f2_p_g1,
                 compare_p_f2_p_g2,
                 compare_p_f2_p_i1,
                 compare_p_f2_p_i2,

                 compare_p_g1_p_g2,
                 compare_p_g1_p_i1,
                 compare_p_g1_p_i2,

                 compare_p_g2_p_i1,
                 compare_p_g2_p_i2,

                 compare_p_i1_p_i2
                )

sim_sim_dats <- list(
                 compare_p_f1_sim_p_f2_sim,
                 compare_p_f1_sim_p_g1_sim,
                 compare_p_f1_sim_p_g2_sim,
                 compare_p_f1_sim_p_i1_sim,
                 compare_p_f1_sim_p_i2_sim,

                 compare_p_f2_sim_p_g1_sim,
                 compare_p_f2_sim_p_g2_sim,
                 compare_p_f2_sim_p_i1_sim,
                 compare_p_f2_sim_p_i2_sim,

                 compare_p_g1_sim_p_g2_sim,
                 compare_p_g1_sim_p_i1_sim,
                 compare_p_g1_sim_p_i2_sim,

                 compare_p_g2_sim_p_i1_sim,
                 compare_p_g2_sim_p_i2_sim,

                 compare_p_i1_sim_p_i2_sim
                    )

obs_score_dat <- scoreStatistics(obs_sim_dats, obs_obs_dats)
obs_score_dat <- obs_score_dat[order(obs_score_dat$Score)]
obs_score_plot <- obs_score_dat %>% ggplot(aes(x=reorder(Comparison, Score), 
                                       y=log(Score))) +
                                       # fill=reorder(Comparison, Score))) +
    geom_bar(stat="identity") +
    xlab("Statistic") + ylab("log(Relative deviance)") +
    theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1),
          plot.margin=unit(c(1, 1, 1, 1.5), "cm"))
ggsave("Images/obs_score_plot.pdf", width=20, height=12)

sim_score_dat <- scoreStatistics(obs_sim_dats, sim_sim_dats)
sim_score_dat <- sim_score_dat[order(sim_score_dat$Score)]
sim_score_plot <- sim_score_dat %>% ggplot(aes(x=reorder(Comparison, Score), 
                                       y=log(Score))) +
                                       # fill=reorder(Comparison, Score))) +
    geom_bar(stat="identity") +
    xlab("Statistic") + ylab("log(Relative deviance)") +
    theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1),
          plot.margin=unit(c(1, 1, 1, 1.5), "cm"))
ggsave("Images/sim_score_plot.pdf", width=20, height=12)

plotComparisons(obs_sim_dat, "Images/sim_obs.pdf")
plotComparisons(part_igb_dat, "Images/partis_igb.pdf")

# Save plots to sumrep ms
sumrep_ms_dir <- "/home/bolson2/Manuscripts/sumrep-ms/Figures"
ggsave(filename=file.path(sumrep_ms_dir, "obs_score_plot.pdf"), 
       plot=obs_score_plot,
       width=12, height=6)
ggsave(filename=file.path(sumrep_ms_dir, "sim_score_plot.pdf"), 
       plot=sim_score_plot,
       width=12, height=6)
plotComparisons(obs_sim_dat, file.path(sumrep_ms_dir, "sim_obs.pdf"))
plotComparisons(part_igb_dat, file.path(sumrep_ms_dir, "partis_igb.pdf"))

