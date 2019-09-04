library(glmnet)
library(plotmo)

devtools::load_all()
loadNewDatasets("data/Annotations", pattern="p_")
loadNewDatasets("data/Annotations", pattern="pi_")
loadNewDatasets("data/Annotations", pattern="A")

set.seed(1)

getSummaryDataTable <- function(dat) {
    summary_df <- data.frame(matrix(NA, nrow=nrow(dat), ncol=0))
    summary_df['cdr3_length'] <- dat %>% getCDR3LengthDistribution
    summary_df['v_3p_del'] <- dat %>% getVGene3PrimeDeletionLengthDistribution
    summary_df['d_5p_del'] <- dat %>% getDGene5PrimeDeletionLengthDistribution
    summary_df['d_3p_del'] <- dat %>% getDGene3PrimeDeletionLengthDistribution
    summary_df['j_5p_del'] <- dat %>% getJGene5PrimeDeletionLengthDistribution
    summary_df['vd_ins'] <- dat %>% getVDInsertionLengthDistribution
    summary_df['dj_ins'] <- dat %>% getDJInsertionLengthDistribution
    summary_df['charge'] <- dat %>% getChargeDistribution
    summary_df['polarity'] <- dat %>% getPolarityDistribution
    summary_df['basicity'] <- dat %>% getBasicityDistribution
    summary_df['bulkiness'] <- dat %>% getBulkinessDistribution
    summary_df['aromaticity'] <- dat %>% getAromaticityDistribution
    summary_df['acidity'] <- dat %>% getAcidityDistribution
    summary_df['gravy'] <- dat %>% getGRAVYDistribution
    summary_df['aliphatic_index'] <- dat %>% getAliphaticIndexDistribution

    return(summary_df)
}

getSummaryDesignMatrix <- function(dat_list) {
    summary_dat <- dat_list %>% 
        lapply(getSummaryDataTable) %>%
        rbindlist %>%
        as.matrix

    return(summary_dat)
}

getSubsampledBCRDatasets <- function(dat_list) {
    sub_dat_list <- dat_list %>%
        lapply( function(x) { 
                    x[['annotations']] %>% subsampleToUniqueClones
                }
        ) 

    return(sub_dat_list)
}

getResponseVector <- function(dat_list) {
    ids <- dat_list %>%
        sapply(nrow) %>%
        Map(function(x, y) { rep(x, y) },
            1:length(dat_list),
            .
           ) %>%
        unlist %>%
        as.factor
    
    return(ids)
}

getLassoModel <- function(dat_list) {
    summary_dat <- dat_list %>%
        getSummaryDesignMatrix

    rows_with_na <- which(is.na(rowSums(summary_dat)))
    summary_dat <- summary_dat[-rows_with_na, ]
    
    # Standardize each covariate column so that lasso coefficients are 
    #   interpretable on a common scale
    summary_dat <- summary_dat %>% scale
    
    ids <- dat_list %>%
        getResponseVector
    ids <- ids[-rows_with_na]
    
    model <- glmnet::glmnet(summary_dat,
                            ids,
                            family="multinomial",
                            type.multinomial="grouped",
                            intercept=FALSE # No need for intercept due to scaling
                           )

    return(model)
}

plotAllBCRLassoPaths <- function(dat_lists) {
    par(mfrow=c(4, 6))
    for(dat in names(dat_lists)) {
        model <- dat_lists[[dat]] %>%
            getLassoModel
    
        for(i in 1:6) {
            lasso_plot <- plot_glmnet(model,
                                      nresponse=i,
                                      label=T,
                                      main=dat
                                     )
        }
    }
}

getSummaryRanks <- function(model) {
    ranks <- list()
    for(i in 1:length(coef(model))) {
        coef_mat <- coef(model)[[i]]
        first_nonzero_indices <- list()
        # First row is for intercept which we don't care about
        for(s in rownames(coef_mat)[-1]) {
            nonzero_index <- which(coef_mat[s, ] != 0)[1] %>%
                ifelse(is.na(.), Inf, .) # If coefficient is always zero, 
                                         # set its branching time to Inf
            first_nonzero_indices[s] <- ifelse(nonzero_index < Inf,
                                               -log(model$lambda[ nonzero_index]),
                                               -Inf
                                              )
        }
        
        index_dat <- first_nonzero_indices %>% as.data.table
        order_dat <- index_dat
        order_dat <- rank(index_dat[1, ],
                               ties.method="random"
                              )
        ranks[[i]] <- first_nonzero_indices %>% unlist
    }
    
    rank_df <- ranks %>%
        do.call("rbind", .)
    
    melted_ranks <- rank_df %>% 
        melt %>%
        setNames(c("Response",
                   "Summary",
                   "Rank"
                  )
        )
    return(melted_ranks)
}

if(F) {
    p_dat_list <- list(p_f1, p_f2, p_g1, p_g2, p_i1, p_i2) %>%
        getSubsampledBCRDatasets
    p_sim_dat_list <- list(p_f1_sim, p_f2_sim, p_g1_sim, p_g2_sim, p_i1_sim, p_i2_sim) %>%
        getSubsampledBCRDatasets
    pi_dat_list <- list(pi_f1, pi_f2, pi_g1, pi_g2, pi_i1, pi_i2) %>%
        getSubsampledBCRDatasets
    pi_sim_dat_list <- list(pi_f1_sim, pi_f2_sim, pi_g1_sim, pi_g2_sim, pi_i1_sim, 
                            pi_i2_sim) %>%
        getSubsampledBCRDatasets
    dat_lists <- list("partis"=p_dat_list,
                      "partis sim"=p_sim_dat_list,
                      "imgt"=pi_dat_list,
                      "imgt sim"=pi_sim_dat_list
                     )
    
    igor_dat_list <- list(A5_S9[['annotations']],
                          A4_i107[['annotations']],
                          A4_i194[['annotations']],
                          A5_S10[['annotations']],
                          A5_S15[['annotations']],
                          A5_S22[['annotations']]
                         )
}


if(F) {
}

    igor_model <- igor_dat_list %>% 
        getLassoModel
    partis_model <- pi_dat_list %>% 
        getLassoModel

plotSummaryRanks <- function(model,
                             filepath
                            ) {
    ranks <- model %>% getSummaryRanks
    p <- ggplot(ranks, 
                aes(x=reorder(Summary, Rank, median), y=Rank)) + 
        geom_boxplot() +
        xlab("Summary") +
        theme_minimal() +
        theme(axis.text.x=element_text(angle=45, hjust=1))
    p
    ggsave(filepath, width=8, height=4,
           plot=p)
}

sumrep_ms_dir <- "/home/bolson2/Manuscripts/sumrep-ms/Figures/Lasso"
plotSummaryRanks(partis_model,
                 file.path(sumrep_ms_dir, "partis_lasso_scores.pdf")
                )
plotSummaryRanks(igor_model,
                 file.path(sumrep_ms_dir, "igor_lasso_scores.pdf")
                )

plotLassoPaths <- function(model, filepath) {
    pdf(filepath, width=10, height=6)
    par(mfrow=c(2, 3), cex=0.9)
    for(i in 1:6) {
        lasso_plot <- plot_glmnet(model,
                                  nresponse=i,
                                  label=T
                                 )
    }
    dev.off()
}

plotLassoPaths(partis_model, 
               file.path(sumrep_ms_dir, "partis_lasso_paths.pdf")
              )
plotLassoPaths(igor_model, 
               file.path(sumrep_ms_dir, "igor_lasso_paths.pdf")
              )
