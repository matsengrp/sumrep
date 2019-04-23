library(viridis)
ggplot2::theme_set(theme_gray(base_size = 14))

source("Examples/ScoreStatistics.R")

loadNewDatasets("data/Comparisons")

obs_sim_igor_dats <- list(
                     compare_A4_i107_A4_i107_sim,
                     compare_A4_i194_A4_i194_sim,
                     compare_A5_S10_A5_S10_sim,
                     compare_A5_S15_A5_S15_sim,
                     compare_A5_S22_A5_S22_sim,
                     compare_A5_S9_A5_S9_sim
                )


obs_obs_igor_dats <- list(
                     compare_A4_i194_A4_i107,
                     compare_A5_S10_A4_i107,
                     compare_A5_S10_A4_i194,
                     compare_A5_S10_A5_S15,
                     compare_A5_S15_A4_i107,
                     compare_A5_S15_A4_i194,
                     compare_A5_S22_A4_i107,
                     compare_A5_S22_A4_i194,
                     compare_A5_S22_A5_S10,
                     compare_A5_S22_A5_S15,
                     compare_A5_S22_A5_S9,
                     compare_A5_S9_A4_i107,
                     compare_A5_S9_A4_i194,
                     compare_A5_S9_A5_S10,
                     compare_A5_S9_A5_S15
                )

sim_sim_igor_dats <- list(
                     compare_A4_i194_sim_A4_i107_sim,
                     compare_A5_S10_sim_A4_i107_sim,
                     compare_A5_S10_sim_A4_i194_sim,
                     compare_A5_S10_sim_A5_S15_sim,
                     compare_A5_S15_sim_A4_i107_sim,
                     compare_A5_S15_sim_A4_i194_sim,
                     compare_A5_S22_sim_A4_i107_sim,
                     compare_A5_S22_sim_A4_i194_sim,
                     compare_A5_S22_sim_A5_S10_sim,
                     compare_A5_S22_sim_A5_S15_sim,
                     compare_A5_S22_sim_A5_S9_sim,
                     compare_A5_S9_sim_A4_i107_sim,
                     compare_A5_S9_sim_A4_i194_sim,
                     compare_A5_S9_sim_A5_S10_sim,
                     compare_A5_S9_sim_A5_S15_sim
                    )

sumrep_ms_igor_dir <- "/home/bolson2/Manuscripts/sumrep-ms/Figures/IgorScores"
plotSummaryScores(dats_1=obs_sim_igor_dats,
                  dats_2=obs_obs_igor_dats,
                  filename=file.path(sumrep_ms_igor_dir,
                                     "obs_score_plot.pdf"
                                    )
                 )
plotSummaryScores(dats_1=obs_sim_igor_dats,
                  dats_2=sim_sim_igor_dats,
                  filename=file.path(sumrep_ms_igor_dir,
                                     "sim_score_plot.pdf"
                                    )
                 )
