getComparisonValue <- function(dat, comparison) {
    value <- dat %>% 
        filter(Comparison == comparison) %$%
        Divergence
    return(value)
}

scoreStatistics <- function(dats_1, 
                            dats_2,
                            comparisons_to_omit=NULL
                           ) {
    comparison_types <- dats_1[[1]][["Comparison"]] %>% 
        unique
    comparison_types <- comparison_types[!(comparison_types %in% comparisons_to_omit)]
        
    score_dat <- matrix(NA, nrow=0, ncol=3) %>% 
        data.table %>%
        setNames(c("Comparison", "Score", "Group"))
    for(c_type in comparison_types) {
            score_1 <- dats_1 %>% 
                sapply(getComparisonValue, c_type) %>% 
                unlist %>%
                mean
            score_2 <- dats_2 %>%
                sapply(getComparisonValue, c_type) %>%
                unlist %>%
                mean
            c_score <- log(score_2/score_1)
            comparison_name <- c_type %>%
                gsub(pattern="compare", replacement="get") %>%
                ifelse(str_sub(., -1) == "s", 
                       substr(., start=1, stop=nchar(.) - 1), 
                       .
                      ) %>%
                getNameFromFunctionString
            summary_group <- c_type %>% getSummaryGroup
            
            score_dat <- rbind(score_dat,
                               data.table(
                                          Comparison=comparison_name,
                                          Score=c_score,
                                          Group=summary_group
                                         )
                              )
    }
    return(score_dat)
}

shortenName <- function(string) {
    shortened_name <- string %>% 
        gsub(pattern="compare", replacement="get") %>%
        getNameFromFunctionString
    return(shortened_name)
}


plotComparisons <- function(dat, filename, cols=1, rows=1) {
    comparison_types <- dat[["Comparison"]] %>% unique
    plot_list <- list()
    for(c_type in comparison_types) {
        dat_sub <- dat[dat[["Comparison"]] == c_type, ]
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

filterScores <- function(dat,
                         comparisons_to_omit
                        ) {
    return(dat[!(dat[["Comparison"]] %in% comparisons_to_omit), ])
}

plotSummaryScores <- function(
                              dats_1, 
                              dats_2,
                              filename,
                              plot_width=12,
                              plot_height=6,
                              comparisons_to_omit=NULL
                             ) {
    score_dat <- scoreStatistics(dats_1, dats_2, comparisons_to_omit=comparisons_to_omit)
    score_dat <- score_dat[order(score_dat[["Score"]])]

    score_dat[["Group"]] <- factor(score_dat[["Group"]],
                                   levels=c("Sequence",
                                            "Rearrangement",
                                            "Physiochemical (CDR3aa)"
                                           )
                                  )

    score_plot <- score_dat %>%
        ggplot(aes(x=reorder(Comparison, -Score), 
                   y=Score
                  )
              ) +
        geom_bar(stat="identity") +
        facet_grid( ~ Group, scales="free_x", space="free_x") +
        xlab("Summary statistic") +
        ylab("log(Relative average divergence)") +
        theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1),
              plot.margin=unit(c(1, 1, 1, 1.5), "cm"),
              panel.background=element_blank(),
              panel.grid.major=element_line(colour="lightgray")
             )
        
    ggsave(filename,
           width=plot_width,
           height=plot_height
          )
}

plotSummaryScoreDifferences <- function(
                                       partis_obs_sim_dats,
                                       partis_obs_obs_dats,
                                       igblast_obs_sim_dats,
                                       igblast_obs_obs_dats,
                                       filename,
                                       plot_width=12,
                                       plot_height=6
                                      ) { 
    partis_score_dat <- scoreStatistics(partis_obs_sim_dats,
                                        partis_obs_obs_dats)
    names(partis_score_dat)[which(names(partis_score_dat) == "Score")] <-
        "PartisScore"
    igblast_score_dat <- scoreStatistics(igblast_obs_sim_dats,
                                         igblast_obs_obs_dats
                                        )
    names(igblast_score_dat)[which(names(igblast_score_dat) == "Score")] <-
        "IgBlastScore"
    score_dat <- merge(partis_score_dat, igblast_score_dat)

    score_dat[["Difference"]] <- 
        score_dat[["PartisScore"]] - score_dat[["IgBlastScore"]]
    score_dat <- score_dat[order(score_dat[["Difference"]])]

    score_dat[["Group"]] <- factor(score_dat[["Group"]],
                                   levels=c("Sequence",
                                            "Rearrangement",
                                            "Physiochemical (CDR3aa)"
                                           )
                                  )
    score_plot <- score_dat %>%
        ggplot(aes(x=reorder(Comparison, -Difference), 
                   y=Difference
                  )
              ) +
        geom_bar(stat="identity") +
        xlab("Summary statistic") +
        ylab("Difference of scores") +
        facet_grid( ~ Group, scales="free_x", space="free_x") +
        theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1),
              plot.margin=unit(c(1, 1, 1, 1.5), "cm"),
              panel.background=element_blank(),
              panel.grid.major=element_line(colour="lightgray")
             )

    ggsave(filename,
           width=plot_width,
           height=plot_height
          )
}

getSummaryGroup <- function(comparison) {
    group <- {}
    if(comparison %in%
       c("compareCDR3LengthDistributions",
         "compareVGeneDistributions",
         "compareJGeneDistributions",
         "compareVGene3PrimeDeletionLengthDistributions",
         "compareJGene5PrimeDeletionLengthDistributions",
         "compareDGene3PrimeDeletionLengthDistributions",
         "compareDGene5PrimeDeletionLengthDistributions",
         "compareDGeneDistributions",
         "compareVDJDistributions",
         "compareVDInsertionLengthDistributions",
         "compareDJInsertionLengthDistributions",
         "compareVDInsertionMatrices",
         "compareDJInsertionMatrices",
         "compareInFramePercentages"
        ) ) {
        group <- "Rearrangement"
    } else if (comparison %in%
               c("comparePairwiseDistanceDistributions",
                 "compareGCContentDistributions",
                 "compareHotspotCountDistributions",
                 "compareColdspotCountDistributions",
                 "compareCDR3PairwiseDistanceDistributions",
                 "compareDistanceFromGermlineToSequenceDistributions",
                 "comparePositionalDistanceBetweenMutationsDistributions",
                 "comparePerGeneMutationRates",
                 "comparePerGenePerPositionMutationRates",
                 "compareAminoAcidDistributions",
                 "compareAminoAcid2merDistributions"
                ) )
    {
        group <- "Sequence"
    } else if (comparison %in%
               c("KideraFactor1Divergence",
                 "KideraFactor2Divergence",
                 "KideraFactor3Divergence",
                 "KideraFactor4Divergence",
                 "KideraFactor5Divergence",
                 "KideraFactor6Divergence",
                 "KideraFactor7Divergence",
                 "KideraFactor8Divergence",
                 "KideraFactor9Divergence",
                 "KideraFactor10Divergence",
                 "AtchleyFactor1Divergence",
                 "AtchleyFactor2Divergence",
                 "AtchleyFactor3Divergence",
                 "AtchleyFactor4Divergence",
                 "AtchleyFactor5Divergence",
                 "compareAliphaticIndexDistributions",
                 "compareGRAVYDistributions",
                 "comparePolarityDistributions",
                 "compareChargeDistributions",
                 "compareBasicityDistributions",
                 "compareAcidityDistributions",
                 "compareAromaticityDistributions",
                 "compareBulkinessDistributions",
                 "compareAminoAcidDistributions",
                 "compareAminoAcid2merDistributions",
                 "compareAliphaticIndexDistributions",
                 "compareGRAVYDistributions",
                 "comparePolarityDistributions",
                 "compareChargeDistributions",
                 "compareBasicityDistributions",
                 "compareAcidityDistributions",
                 "compareAromaticityDistributions",
                 "compareBulkinessDistributions",
                 "compareInFramePercentages"
               ) ) {
        group <- "Physiochemical (CDR3aa)"
    } else {
        stop(paste("Invalid comparison string",
                   comparison
                  )
        )
    } 

    return(group)
}
