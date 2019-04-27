devtools::load_all()

library(xtable)

printTable <- function(dat, filename, digits=4, hline_pos=0, label=NULL) {
    sink(filename)
    cat("\\begin{figure}", '\n')
    dat %>%
        xtable::xtable(digits=digits
                      ) %>%
        print(
              floating=FALSE,
              include.rownames=FALSE,
              latex.environments="center",
              hline.after=c(0, hline_pos),
              sanitize.text.function=function(x){x}
             )
    cat(paste0("\\label{", label, "}\n"))
    cat("\\end{figure}", '\n')
    sink()
}

getTrueDistributionDatasetInfo <- function(dat,
                                           summary_function,
                                           ...
                                           ) {
    pt <- proc.time()
    true_dist <- summary_function(dat,
                                  approximate=FALSE,
                                  ...
                                 ) 

    true_time <- (proc.time() - pt)["elapsed"]
    true_dat <- cbind.data.frame(true_dist, "True") %>%
        setNames(c("Value", "Setting"))
    return(list(true_dat=true_dat, true_time=true_time))
}

getNaiveNNDistribution <- function(dat,
                                column="sequence_alignment",
                                approximate=TRUE,
                                ...
                               ) {
    sequence_list <- dat[[column]]
    if(approximate) {
        distribution <- sequence_list %>%
            getApproximateDistribution(summary_function=getNearestNeighborDistances,
                                       divergence_function=getJSDivergence,
                                       ...
                                      )
    } else {
        distribution <- sequence_list %>%
            getNearestNeighborDistances
    }

    return(distribution) 
}

loadNewDatasets("data/Annotations", pattern="p_f1")

runPerformanceAnalysis <- function(dat,
                                   distribution_function,
                                   tols,
                                   trial_count,
                                   continuous,
                                   ...
                                  ) {
    print(dat %>% nrow)

    # Get true distribution
    pt <- proc.time()
    true_info <- getTrueDistributionDatasetInfo(dat,
                                                distribution_function,
                                                ...
                                               ) 
    true_dat <- true_info[["true_dat"]]
    true_time <- true_info[["true_time"]]

    # Get subsampled distributions for various tolerances
    times <- {}
    metric_names <- c("Tolerance", "Time", "Divergence", "Trial")
    metric_dat <- matrix(nrow=0, ncol=4) %>%
        data.frame %>%
        setNames(metric_names)
    divergences <- {}
    for(trial in 1:trial_count) {
        distributions <- list()
        
        for(i in 1:length(tols)) { 
            
            print(tols[i])
            pt <- proc.time()

            distributions[[i]] <- distribution_function(dat,
                                                        approximate=TRUE,
                                                        tol=tols[i],
                                                        ...
                                                       )
            times[i] <- (proc.time() - pt)["elapsed"]
            divergences[i] <- getJSDivergence(distributions[[i]], 
                                              true_dat[["Value"]],
                                              continuous=continuous,
                                              KL=TRUE
                                             )
            metric_dat <- rbind(metric_dat,
                                data.frame(Tolerance=tols[i],
                                           Time=times[i],
                                           Divergence=divergences[i],
                                           Trial=trial)
                                )
        }
    }
    
    dist_dat <- distributions %>%
        lapply(data.frame) %>%
        Map(cbind, ., paste("Tol = ", tols)) %>%
        lapply(setNames, c("Value", "Setting")) %>%
        do.call(rbind, .) %>%
        rbind(., true_dat)

    dist_dat[["Value"]] <- as.numeric(dist_dat[["Value"]])

    summary_table <- dist_dat[["Setting"]] %>% 
        unique %>% 
        lapply(function(x) { 
                  dist_dat[dist_dat[["Setting"]] == x, ][["Value"]] %>% summary 
               }) %>% 
        do.call(rbind, .) %>%
        as.data.table %>%
        cbind(Distribution=dist_dat[["Setting"]] %>% unique, .)

    return(list(Distributions=dist_dat,
                Metrics=metric_dat,
                TrueTime=true_time,
                Summary=summary_table
               )
          )
}

runSingleSummaryAnalysis <- function(dat,
                               distribution_function,
                               tols,
                               trial_count,
                               continuous,
                               column,
                               out_dir,
                               xlab="Value"
                              ) {
    perf_list <- runPerformanceAnalysis(
                           dat=dat,
                           distribution_function=distribution_function,
                           tols=tols,
                           trial_count=trial_count,
                           continuous=continuous,
                           column=column
                          )
    dist_dat <- perf_list[["Distributions"]]
    metric_dat <- perf_list[["Metrics"]]
    true_time <- perf_list[["TrueTime"]]
    summary_table <- perf_list[["Summary"]]

    dir.create(out_dir)
    
    freq_plot <- ggplot(dist_dat,
           aes(x=as.numeric(Value), group=Setting, color=Setting)) +
        geom_freqpoly(aes(y=..density..), binwidth=1) +
        xlim(quantile(dist_dat[["Value"]], c(0.01, 0.99), na.rm=T))  +
        xlab(xlab) +
        ylab("Density") +
        theme(panel.background=element_blank(),
              panel.grid.major=element_line(color="lightgray"),
              panel.grid.minor=element_line(color="lightgray")
              )
    ggsave(file.path(out_dir, "freqpoly_by_tol.pdf"), width=6, height=4)
    
    density_plot <- ggplot(dist_dat,
           aes(x=as.numeric(Value), group=Setting, color=Setting)) +
        geom_density(adjust=4) +
        xlim(quantile(dist_dat[["Value"]], c(0.01, 0.99), na.rm=T))  +
        xlab(xlab) +
        ylab("Density") +
        theme(panel.background=element_blank(),
              panel.grid.major=element_line(color="lightgray"),
              panel.grid.minor=element_line(color="lightgray")
             )
    ggsave(file.path(out_dir, "density_by_tol.pdf"), width=6, height=4)
    
    ecdf_plot <- ggplot(dist_dat,
           aes(x=as.numeric(Value), group=Setting, color=Setting)) +
        stat_ecdf() +
        xlab(xlab) +
        ylab("Density") +
        # We don't want outliers to stretch the x-axis too much
        xlim(quantile(dist_dat[["Value"]], c(0.01, 0.99), na.rm=T))  +
        theme(panel.background=element_blank(),
              panel.grid.major=element_line(color="lightgray"),
              panel.grid.minor=element_line(color="lightgray")
             )
    ggsave(file.path(out_dir, "ecdf_by_tol.pdf"), width=6, height=4)
    
    time_plot <- ggplot(d=metric_dat, 
                        aes(x=log10(Tolerance), y=Time, group=log10(Tolerance))) +
        geom_boxplot() +
        geom_hline(aes(yintercept=true_time, color="Full distribution")) +
        xlab("Log_10(tolerance)") +
        ylab("Time (seconds)") +
        theme(panel.background=element_blank(),
              panel.grid.major=element_line(color="lightgray"),
              panel.grid.minor=element_line(color="lightgray")
             )
    ggsave(file.path(out_dir, "time_by_tol.pdf"), width=6, height=4)
    
    log_time_plot <- ggplot(d=metric_dat, 
                            aes(x=log10(Tolerance), y=log(Time), group=log10(Tolerance))) +
        geom_boxplot() +
        geom_hline(aes(yintercept=log(true_time), color="Full distribution")) +
        xlab("Log_10(Tolerance)") +
        ylab("Time (log-seconds)") +
        theme(panel.background=element_blank(),
              panel.grid.major=element_line(color="lightgray"),
              panel.grid.minor=element_line(color="lightgray")
             )
    ggsave(file.path(out_dir, "log_time_by_tol.pdf"), width=6, height=4)
    
    div_plot <- ggplot(d=metric_dat, 
                       aes(x=log10(Tolerance), y=Divergence, group=log10(Tolerance))) +
        geom_boxplot() +
        xlab("Log_10(tolerance)") +
        ylab("KL-divergence") +
        ylim(0, max(metric_dat[["Divergence"]])) +
        theme(panel.background=element_blank(),
              panel.grid.major=element_line(color="lightgray"),
              panel.grid.minor=element_line(color="lightgray")
             )
    ggsave(file.path(out_dir, "div_by_tol.pdf"), width=6, height=4)
    
    summary_table %>%
        printTable(filename=file.path(out_dir, "summary_by_tol.tex"),
                   digits=c(0, 0, 0, 0, 0, 2, 0, 0),
                   label="tab:SummaryTable"
                  )
}

test_dat <- p_f1[["annotations"]] %>%
    subsample(10000, replace=FALSE)

rbind(test_dat[["junction_length"]] %>% summary,
                     test_dat[["sequence_alignment"]] %>% sapply(nchar) %>% summary
                    ) %>%
    printTable("~/Manuscripts/sumrep-ms/Tables/sequence_lengths.tex",
               digits=c(0,0,0, 0, 1, 0, 0)
              )

pairwise_dist_analysis <- runSingleSummaryAnalysis(
    dat=test_dat,
    distribution_function=getPairwiseDistanceDistribution,
    tols=10^seq(-1, -7),
    trial_count=10,
    continuous=FALSE,
    column="sequence_alignment",
    out_dir="~/Manuscripts/sumrep-ms/Figures/PairwiseDistance",
    xlab="Pairwise distance"
)

nn_dist_analysis_cdr3 <- runSingleSummaryAnalysis(
    dat=test_dat,
    distribution_function=getNearestNeighborDistribution,
    tols=10^seq(-1, -7),
    trial_count=10,
    continuous=FALSE,
    column="junction",
    out_dir="~/Manuscripts/sumrep-ms/Figures/NearestNeighbor/CDR3",
    xlab="NN distance"
)

nn_dist_analysis_sequence <- runSingleSummaryAnalysis(
    dat=test_dat,
    distribution_function=getNearestNeighborDistribution,
    tols=10^seq(-1, -7),
    trial_count=5,
    continuous=FALSE,
    column="sequence_alignment",
    out_dir="~/Manuscripts/sumrep-ms/Figures/NearestNeighbor/Sequence",
    xlab="NN distance"
)

naive_nn_dist_analysis <- runSingleSummaryAnalysis(
    dat=test_dat,
    distribution_function=getNaiveNNDistribution,
    tols=10^seq(-1, -7),
    trial_count=10,
    continuous=FALSE,
    column="sequence_alignment",
    out_dir="~/Manuscripts/sumrep-ms/Figures/NaiveNearestNeighbor",
    xlab="NN distance"
)

runAnalysisBySampleSize <- function(
    dat,
    sample_sizes,
    distribution_function,
    tols,
    trial_count,
    continuous,
    column,
    out_dir
) {
    dats <- sample_sizes %>%
        lapply(function(x) {
                   dat %>% subsample(x, replace=FALSE)
               })

    metric_dat <- dats %>%
        lapply(.,
               FUN=function(dat) {
                perf <- runPerformanceAnalysis(
                            dat=dat,
                            distribution_function=distribution_function,
                            tols=tols,
                            trial_count=trial_count,
                            continuous=continuous,
                            column=column
                )

                dist_dat <- cbind(perf[["Metrics"]],
                                  Size=nrow(dat),
                                  TrueTime=perf[["TrueTime"]]
                                 )
                return(dist_dat)
               }) %>%
        do.call(rbind, .)

    metric_dat[["Tolerance"]] <- metric_dat[["Tolerance"]] %>% as.factor
    metric_dat[["Efficiency"]] <- metric_dat[["TrueTime"]]/metric_dat[["Time"]]

    size_plot <- metric_dat %>%
        ggplot(aes(x=log(Size), 
                   y=Divergence, 
                   group=interaction(log(Size), Tolerance),
                   fill=Tolerance
                  )
              ) +
        geom_boxplot() +
        xlab("log(Size)") +
        ylab("KL-divergence") +
        ylim(0, max(metric_dat[["Divergence"]])) +
        theme(panel.background=element_blank(),
              panel.grid.major=element_line(color="lightgray"),
              panel.grid.minor=element_line(color="lightgray")
             )

    ggsave(file.path(out_dir,
                     "div_by_size_and_tol.pdf"), width=6, height=4)

    time_plot <- ggplot(d=metric_dat, 
                        aes(x=log(Size), 
                            y=log(Time), 
                            group=interaction(log(Size), Tolerance),
                            fill=Tolerance
                            )
                        )  +
        geom_boxplot() +
        xlab("log(Size)") +
        ylab("Time in log(seconds)") +
        theme(panel.background=element_blank(),
              panel.grid.major=element_line(color="lightgray"),
              panel.grid.minor=element_line(color="lightgray")
             )
    ggsave(file.path(out_dir, "time_by_size_and_tol.pdf"), width=6, height=4)

    efficiency_plot <- ggplot(d=metric_dat, 
                        aes(x=log(Size), 
                            y=log(Efficiency), 
                            group=interaction(log(Size), Tolerance),
                            fill=Tolerance
                            )
                        )  +
        geom_boxplot() +
        xlab("log(Size)") +
        ylab("log(Efficiency)") +
        theme(panel.background=element_blank(),
              panel.grid.major=element_line(color="lightgray"),
              panel.grid.minor=element_line(color="lightgray")
             )
    ggsave(file.path(out_dir, "efficiency_by_size_and_tol.pdf"), width=6, height=4)
}

pairwise_dist_size_analysis <- runAnalysisBySampleSize(
    dat=p_f1[["annotations"]],
    sample_sizes=c(6, 7, 8, 9, 10) %>% exp,
    distribution_function=getPairwiseDistanceDistribution,
    tols=10^seq(-1, -3),
    trial_count=10,
    continuous=FALSE,
    column="junction",
    out_dir="~/Manuscripts/sumrep-ms/Figures/PairwiseDistance"
)

nn_size_analysis <- runAnalysisBySampleSize(
    dat=p_f1[["annotations"]],
    sample_sizes=c(6, 7, 8, 9, 10) %>% exp,
    distribution_function=getNearestNeighborDistribution,
    tols=10^seq(-1, -5),
    trial_count=5,
    continuous=FALSE,
    column="sequence_alignment",
    out_dir="~/Manuscripts/sumrep-ms/Figures/NearestNeighbor"
)

runMultipleSummaryAnalysis <- function(
    dat,
    summaries,
    summary_names,
    continuous,
    tols,
    trial_count,
    column,
    out_dir
) {
    distribution_dat <- mapply(
            FUN=function(summary, continuous, summary_name) {
                    perf <- runPerformanceAnalysis(dat=dat,
                                           distribution_function=summary,
                                           tols=tols,
                                           trial_count=trial_count,
                                           continuous=continuous,
                                           column=column
                                          )
                    dist_dat <- cbind(perf[["Metrics"]], Summary=summary_name)
                    return(dist_dat)
                },
                summaries,
                continuous,
                summary_names ,
                SIMPLIFY=FALSE
            ) %>%
        do.call(rbind, .)

    distribution_dat[["Tolerance"]] <- distribution_dat[["Tolerance"]] %>% as.factor

    dir.create(out_dir)

    summ_plot <- distribution_dat %>%
        ggplot(aes(x=Summary,
                   y=Divergence,
                   group=interaction(Summary, Tolerance),
                   fill=Tolerance
                   )) +
               geom_boxplot() +
               theme(panel.background=element_blank(),
                     panel.grid.major=element_line(color="lightgray"),
                     panel.grid.minor=element_line(color="lightgray")
                    )
    ggsave(file.path(out_dir,
                     "div_by_summary_and_tol.pdf"),
           width=6, height=4)
}

multiple_summary_analysis <- runMultipleSummaryAnalysis(
    dat=p_f1[["annotations"]] %>% subsample(5000, replace=FALSE),
    summaries=c(getPairwiseDistanceDistribution,
                getGCContentDistribution,
                getHotspotCountDistribution,
                getColdspotCountDistribution
               ),
    summary_names=c(
                    "Pairwise distances", 
                    "GC content", 
                    "Hotspot count", 
                    "Coldspot count"
                    ),
    continuous=c(FALSE, TRUE, FALSE, FALSE),
    tols=10^seq(-1, -7),
    trial_count=10,
    column="junction",
    out_dir="~/Manuscripts/sumrep-ms/Figures/Multiple"
)
