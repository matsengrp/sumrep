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
    print(label)
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

loadNewDatasets("data/Annotations")

runPerformanceAnalysis <- function(dat,
                                   distribution_function,
                                   tols,
                                   trial_count,
                                   continuous,
                                   ...
                                  ) {
    # Get true distribution
    pt <- proc.time()
    true_info <- getTrueDistributionDatasetInfo(dat,
                                                distribution_function,
                                                ...
                                               ) 
    true_dat <- true_info$true_dat
    true_time <- true_info$true_time
    
    # Get subsampled distributions for various tolerances
    times <- {}
    metric_names <- c("Tolerance", "Time", "Divergence", "Trial")
    metric_dat <- matrix(nrow=0, ncol=4) %>%
        data.frame %>%
        setNames(metric_names)
    divergences <- {}
    for(trial in 1:trial_count) {
        distributions <- list()
        
        for(i in 1:length(tols)) { pt <- proc.time()
            distributions[[i]] <- distribution_function(dat,
                                                        approximate=TRUE,
                                                        tol=tols[i],
                                                        ...
                                                       )
            times[i] <- (proc.time() - pt)["elapsed"]
            divergences[i] <- getJSDivergence(distributions[[i]], 
                                              true_dat$Value,
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

    dist_dat$Value <- as.numeric(dist_dat$Value)

    summary_table <- dist_dat$Setting %>% 
        unique %>% 
        lapply(function(x) { 
                  dist_dat[dist_dat$Setting == x, ]$Value %>% summary 
               }) %>% 
        do.call(rbind, .) %>%
        as.data.table %>%
        cbind(Distribution=dist_dat$Setting %>% unique, .)

    return(list(Distributions=dist_dat,
                Metrics=metric_dat,
                TrueTime=true_time,
                Summary=summary_table
               )
          )
}

run <- TRUE
if(run) {
    perf_list <- runPerformanceAnalysis(
                           dat=p_f1$annotations %>% 
                               subsample(10000, replace=FALSE),
                           distribution_function=getPairwiseDistanceDistribution,
                           tols=10^seq(-1, -7),
                           trial_count=10,
                           continuous=FALSE,
                           column="cdr3s"
                          )
    dist_dat <- perf_list$Distributions
    metric_dat <- perf_list$Metrics
    true_time <- perf_list$TrueTime
    summary_table <- perf_list$Summary
}
    
freq_plot <- ggplot(dist_dat,
       aes(x=as.numeric(Value), group=Setting, colour=Setting)) +
    geom_freqpoly(aes(y=..density..), binwidth=1) +
    xlim(quantile(dist_dat$Value, c(0.01, 0.99)))  +
    xlab("Pairwise Distance") +
    ylab("Density")
ggsave("~/Manuscripts/sumrep-ms/Figures/freqpoly_by_tol.pdf", width=10)

density_plot <- ggplot(dist_dat,
       aes(x=as.numeric(Value), group=Setting, colour=Setting)) +
    geom_density(adjust=4) +
    xlim(quantile(dist_dat$Value, c(0.01, 0.99)))  +
    xlab("Pairwise Distance") +
    ylab("Density")
ggsave("~/Manuscripts/sumrep-ms/Figures/density_by_tol.pdf", width=10)

ecdf_plot <- ggplot(dist_dat,
       aes(x=as.numeric(Value), group=Setting, colour=Setting)) +
    stat_ecdf() +
    xlab("GC content") +
    ylab("Density") +
    # We don't want outliers to stretch the x-axis
    xlim(quantile(dist_dat$Value, c(0.01, 0.99))) 
ggsave("~/Manuscripts/sumrep-ms/Figures/ecdf_by_tol.pdf", width=10)

time_plot <- ggplot(d=metric_dat, 
                    aes(x=log10(Tolerance), y=Time, group=log10(Tolerance))) +
    geom_boxplot() +
    geom_hline(yintercept=true_time, colour="red") +
    geom_text(aes(0, true_time, label="Time for full dataset", hjust=1, 
                  vjust=-1)) +
    xlab("Log_10(tolerance)") +
    ylab("Time (seconds)") +
    ggtitle("Time complexity of distribution subsampling by tolerance")
ggsave("~/Manuscripts/sumrep-ms/Figures/time_by_tol.pdf", width=10)

log_time_plot <- ggplot(d=metric_dat, 
                        aes(x=log10(Tolerance), y=log(Time), group=log10(Tolerance))) +
    geom_boxplot() +
    geom_hline(yintercept=log(true_time), colour="red") +
    geom_text(aes(0, log(true_time), label="Log(Time) for full dataset", 
                  hjust=1, vjust=-1)) +
    xlab("Log_10(Tolerance)") +
    ylab("Log(time) (log-seconds)") +
    ggtitle(
     "Time complexity (in log-seconds) of distribution subsampling by tolerance"
    )
ggsave("~/Manuscripts/sumrep-ms/Figures/log_time_by_tol.pdf", width=10)

div_plot <- ggplot(d=metric_dat, 
                   aes(x=log10(Tolerance), y=Divergence, group=log10(Tolerance))) +
    geom_boxplot() +
    xlab("Log_10(tolerance)") +
    ylab("KL-divergence") +
    ggtitle("KL-divergence to true distribution by tolerance")
ggsave("~/Manuscripts/sumrep-ms/Figures/div_by_tol.pdf", width=10)

summary_table %>%
    printTable(filename="~/Manuscripts/sumrep-ms/Tables/summary_by_tol.tex",
               digits=c(0, 0, 0, 0, 0, 2, 0, 0)
              )
