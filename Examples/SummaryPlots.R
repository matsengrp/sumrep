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
            score_dat <- rbind(score_dat,
                               data.table(Comparison=shortenName(c_type),
                                          Score=c_score))
    }
    return(score_dat)
}

shortenName <- function(string) {
    shortened_name <- string %>% 
        gsub(pattern="compare", replace="") %>%
        gsub(pattern="Distribution", replace="")
    return(shortened_name)
}

#' Handy automatic dimension retrieval given number of grid cells 
#'
getGridDims <- function(n) {
    cols <- n %>% sqrt %>% floor 
    rows <- ceiling(n/cols) %>% floor 
    return(c(cols, rows))
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

sim_dats <- list(
                 compare_p_f1_p_f1_sim,
                 compare_p_f2_p_f2_sim,
                 compare_p_g1_p_g1_sim
                )

obs_dats <- list(
                 compare_p_f1_p_f2,
                 compare_p_f1_p_g1,
                 compare_p_f2_p_g1
                )

score_dat <- scoreStatistics(sim_dats, obs_dats)
score_dat <- score_dat[order(score_dat$Score)]
score_plot <- score_dat %>% ggplot(aes(x=reorder(Comparison, Score), 
                                       y=log(Score))) +
                                       # fill=reorder(Comparison, Score))) +
    geom_bar(stat="identity") +
    xlab("Comparison") + ylab("log(Relative deviance)") +
    ggtitle(paste0("Comparison scores of observations-to-simulations with ",
                   "respect to comparisons of observations-to-observations")) +
    theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1),
          plot.margin=unit(c(1, 1, 1, 1.5), "cm"))
ggsave("Images/score_plot.pdf", width=20, height=12)

plotComparisons(obs_sim_dat, "Images/sim_obs.pdf")
plotComparisons(part_igb_dat, "Images/partis_igb.pdf")

# Save plots to sumrep ms
sumrep_ms_dir <- "/home/bolson2/Manuscripts/sumrep-ms/Figures"
ggsave(filename=file.path(sumrep_ms_dir, "score_plot.pdf"), 
       plot=score_plot,
       width=10, height=6)
plotComparisons(obs_sim_dat, file.path(sumrep_ms_dir, "sim_obs.pdf"))
plotComparisons(part_igb_dat, file.path(sumrep_ms_dir, "partis_igb.pdf"))

