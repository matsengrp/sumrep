library(viridis)
ggplot2::theme_set(theme_gray(base_size = 14))

source("Examples/ScoreStatistics.R")

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

obs_obs_igb_dats <- list(
                 compare_i_f1_i_f2,
                 compare_i_f1_i_g1,
                 compare_i_f1_i_g2,
                 compare_i_f1_i_i1,
                 compare_i_f1_i_i2,

                 compare_i_f2_i_g1,
                 compare_i_f2_i_g2,
                 compare_i_f2_i_i1,
                 compare_i_f2_i_i2,

                 compare_i_g1_i_g2,
                 compare_i_g1_i_i1,
                 compare_i_g1_i_i2,

                 compare_i_g2_i_i1,
                 compare_i_g2_i_i2,

                 compare_i_i1_i_i2
                )

obs_sim_igb_dats <- list(
                         compare_i_f1_pi_f1_sim,
                         compare_i_f2_pi_f2_sim,
                         compare_i_g1_pi_g1_sim,
                         compare_i_g2_pi_g2_sim,
                         compare_i_i1_pi_i1_sim,
                         compare_i_i2_pi_i2_sim
                        )


part_igb_dat <- part_igb_dat[!(part_igb_dat[["Comparison"]] %in% 
                               c(
                                 "comparePerGeneMutationRates",
                                 "comparePerGenePerPositionMutationRates"
                                 )), ]

part_igb_dat[["Type2"]] <- factor(part_igb_dat[["Type2"]],
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

obs_sim_partis_dat <- rbind(compare_p_f1_p_f1_sim,
                     compare_p_f2_p_f2_sim,
                     compare_p_g1_p_g1_sim,
                     compare_p_f1_p_f2,
                     compare_p_f1_p_g1,
                     compare_p_f2_p_g1
                    )

obs_sim_partis_dat[["Type2"]] <- factor(obs_sim_partis_dat[["Type2"]], 
                      levels=c("p_f1_sim",
                               "p_f2_sim",
                               "p_g1_sim",
                               "p_f1",
                               "p_f2",
                               "p_g1"))

obs_sim_partis_dats <- list(
                 compare_p_f1_p_f1_sim,
                 compare_p_f2_p_f2_sim,
                 compare_p_g1_p_g1_sim,
                 compare_p_g2_p_g2_sim,
                 compare_p_i1_p_i1_sim,
                 compare_p_i2_p_i2_sim
                )

obs_obs_partis_dats <- list(
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

sim_sim_partis_dats <- list(
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

sumrep_ms_partis_dir <- "/home/bolson2/Manuscripts/sumrep-ms/Figures/PartisScores"
sumrep_ms_partis_dir %>% dir.create

plotSummaryScores(dats_1=obs_sim_partis_dats,
                  dats_2=obs_obs_partis_dats,
                  filename=file.path(sumrep_ms_partis_dir,
                                     "obs_score_plot.pdf"
                                    )
                 )

plotSummaryScores(dats_1=obs_sim_partis_dats,
                  dats_2=sim_sim_partis_dats,
                  filename=file.path(sumrep_ms_partis_dir,
                                     "sim_score_plot.pdf"
                                    )
                 )

plotSummaryScores(dats_1=obs_sim_igb_dats,
                  dats_2=obs_obs_igb_dats,
                  filename=file.path(sumrep_ms_partis_dir,
                                     "obs_score_plot_igb.pdf"
                                    )
                 )

plotComparisons(obs_sim_partis_dat, "Images/sim_obs.pdf")
plotComparisons(part_igb_dat, "Images/partis_igb.pdf")

# Save plots to sumrep ms
plotComparisons(obs_sim_partis_dat, file.path(sumrep_ms_dir, "sim_obs.pdf"))
plotComparisons(part_igb_dat, file.path(sumrep_ms_dir, "partis_igb.pdf"))

