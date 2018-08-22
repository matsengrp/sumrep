devtools::load_all()

getTrueDistributionDatasetInfo <- function(dat,
                                           summary_function,
                                           ...
                                           ) {
    pt <- proc.time()
    true_dist <- summary_function(dat,
                                  ...) 
    true_time <- (proc.time() - pt)["elapsed"]
    true_dat <- cbind(true_dist, "True") %>%
        data.frame %>%
        setNames(c("Value", "Setting"))
    return(list(true_dat=true_dat, true_time=true_time))
}

loadNewDatasets("data/Annotations")

do_dists <- TRUE
if(do_dists) {
    dat <- p_f1$annotations %>% subsample(10000)
    distribution_function <- getPairwiseDistanceDistribution

    # Get true distribution
    pt <- proc.time()
    true_info <- getTrueDistributionDatasetInfo(dat,
                                               distribution_function,
                                               column="cdr3s",
                                               approximate=FALSE) 
    true_dat <- true_info$true_dat
    true_time <- true_info$true_time
    
    # Get subsampled distributions for various tolerances
    times <- {}
    distributions <- list()
    tols <- 10^seq(-1, -7)
    
    for(i in 1:length(tols)) { pt <- proc.time()
        distributions[[i]] <- getPairwiseDistanceDistribution(dat,
                                                              column="cdr3s",
                                                              approximate=TRUE,
                                                              tol=tols[i])
        times[i] <- (proc.time() - pt)["elapsed"]
    }
    
    dist_dat <- distributions %>%
        lapply(data.frame) %>%
        Map(cbind, ., paste("Tol = ", tols)) %>%
        lapply(setNames, c("Value", "Setting")) %>%
        do.call(rbind, .) %>%
        rbind(., true_dat)
    
    print("Hi")
    dist_dat$Value <- as.numeric(dist_dat$Value)
}
    
dist_plot <- ggplot(dist_dat,
       aes(x=as.numeric(Value), group=Setting, colour=Setting)) +
    geom_density(adjust=3) +
    xlim(0, 60) +
    xlab("Pairwise Distance") +
    ylab("Density")
ggsave("~/Manuscripts/sumrep-ms/Figures/dists_by_tol.pdf", width=10)

time_dat <- data.frame(Time=times, Tolerance=tols)
time_plot <- ggplot(d=time_dat, aes(x=log10(Tolerance), y=Time)) +
    geom_point() +
    geom_hline(yintercept=true_time, colour="red") +
    geom_text(aes(0, true_time, label="Time for full dataset", hjust=1, 
                  vjust=-1)) +
    xlab("Log_10(tolerance)") +
    ylab("Time (seconds)") +
    ggtitle("Time complexity of distribution subsampling by tolerance")

log_time_plot <- ggplot(d=time_dat, aes(x=log10(Tolerance), y=log(Time))) +
    geom_point() +
    geom_hline(yintercept=log(true_time), colour="red") +
    geom_text(aes(0, log(true_time), label="Log(Time) for full dataset", 
                  hjust=1, vjust=-1)) +
    xlab("Log_10(Tolerance)") +
    ylab("Log(time) (log-seconds)") +
    ggtitle(
     "Time complexity (in log-seconds) of distribution subsampling by tolerance"
    )

m <- multiplot(plotlist=list(time_plot, log_time_plot), cols=2)
ggsave("~/Manuscripts/sumrep-ms/Figures/time_by_tol.pdf", width=10)
