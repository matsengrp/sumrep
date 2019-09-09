library(viridis)
ggplot2::theme_set(theme_gray(base_size = 14))

source("Examples/ScoreStatistics.R")

loadNewDatasets("data/Annotations", pattern="A")
loadNewDatasets("data/Comparisons", pattern="compare_A")

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
sumrep_ms_igor_dir %>% dir.create

# Since IGoR does not output a "sequence_alignment" column (or any column
#   for the full variable-region sequence), or a "vj_in_frame" column
#   for annotations, let's omit the corresponding summaries from the
#   analysis
summaries_to_omit <- c("getPairwiseDistanceDistribution",
                       "getGCContentDistribution",
                       "getInFramePercentage"
                      )
comparisons_to_omit <- c("comparePairwiseDistanceDistributions",
                         "compareGCContentDistributions",
                         "compareInFramePercentages"
                        )

plotSummaryScores(dats_1=obs_sim_igor_dats,
                  dats_2=obs_obs_igor_dats,
                  filename=file.path(sumrep_ms_igor_dir,
                                     "obs_score_plot.pdf"
                                    ),
                  comparisons_to_omit=comparisons_to_omit
                 )
plotSummaryScores(dats_1=obs_sim_igor_dats,
                  dats_2=sim_sim_igor_dats,
                  filename=file.path(sumrep_ms_igor_dir,
                                     "sim_score_plot.pdf"
                                    ),
                  comparisons_to_omit=comparisons_to_omit
                 )

igor_plots <- plotUnivariateDistributions(
    list(
         A4_i194[["annotations"]],
         (A4_i194_sim %>% subsetToInFrame)[["annotations"]],
         A5_S10[["annotations"]],
         (A5_S10_sim %>% subsetToInFrame)[["annotations"]],
         A5_S9[["annotations"]],
         (A5_S9_sim %>% subsetToInFrame)[["annotations"]]
    ),
    locus="trb",
    color=c(rep("A4_i194", 2),
            rep("A5_S10", 2),
            rep("A5_S9", 2)
           ),
    lty=rep(c("Observed", "Simulated"), 3),
    functions_to_omit=summaries_to_omit
)

ggsave(file.path(sumrep_ms_igor_dir,
                 "igor_freqpoly.pdf"
                ),
       plot=igor_plots[["freqpoly"]],
       width=14,
       height=14
)

ggsave(file.path(sumrep_ms_igor_dir,
                 "igor_ecdf.pdf"
                ),
       plot=igor_plots[["ecdf"]],
       width=14,
       height=14
)

