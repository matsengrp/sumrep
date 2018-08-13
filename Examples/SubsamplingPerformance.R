devtools::load_all()

loadNewDatasets("data/Annotations")

do_dists <- TRUE
if(do_dists) {
    dat <- p_f1$annotations %>% subsample(5000)
    
    times <- {}
    distributions <- list()
    tols <- 10^seq(-1, -7)
    
    for(i in 1:length(tols)) { pt <- proc.time()
        distributions[[i]] <- getPairwiseDistanceDistribution(dat$cdr3s,
                                                              approximate=TRUE,
                                                              tol=tols[i])
        times[i] <- (proc.time() - pt)["elapsed"]
    }
    
    pt <- proc.time()
    true_dist <- getPairwiseDistanceDistribution(dat$cdr3s,
                                                 approximate=FALSE) 
    true_time <- (proc.time() - pt)["elapsed"]
    true_dat <- cbind(true_dist, "True", 500) %>%
        data.frame %>%
        setNames(c("Value", "Setting", "Smoothness"))
    
    dist_dat <- distributions %>%
        lapply(data.frame) %>%
        Map(cbind, ., paste("Tol = ", tols), smooth) %>%
        lapply(setNames, c("Value", "Setting", "Smoothness")) %>%
        do.call(rbind, .) %>%
        rbind(., true_dat)
    
    dist_dat$Value <- as.numeric(dist_dat$Value)
    dist_dat$Smoothness <- as.numeric(dist_dat$Smoothness)
}
    
p1 <- ggplot(dist_dat,
       aes(x=as.numeric(Value), group=Setting, colour=Setting)) +
    geom_density(adjust=3) +
    xlim(0, 60) +
    xlab("Pairwise Distance") +
    ylab("Density")

p2 <- ggplot(dist_dat,
       aes(x=as.numeric(Value), group=Setting, colour=Setting, adjust=Smoothness)) +
    stat_ecdf() +
    xlim(0, 60)

time_dat <- data.frame(Time=times, Tolerance=tols)
p3 <- ggplot(d=time_dat, aes(x=log10(Tolerance), y=Time)) +
    geom_point() +
    geom_hline(yintercept=true_time, colour="red") +
    geom_text(aes(0, true_time, label="Time for full dataset", hjust=1, vjust=-1)) +
    xlab("Log_10(tolerance)") +
    ylab("Time (seconds)") +
    ggtitle("Time complexity of distribution subsampling by tolerance")


p3_log_log <- ggplot(d=time_dat, aes(x=log10(Tolerance), y=log(Time))) +
    geom_point() +
    geom_hline(yintercept=log(true_time), colour="red") +
    geom_text(aes(0, log(true_time), label="Log(Time) for full dataset", hjust=1, vjust=-1)) +
    xlab("Log_10(Tolerance)") +
    ylab("Log(time) (log-seconds)") +
    ggtitle("Time complexity (in log-seconds) of distribution subsampling by tolerance")

multiplot(plotlist=list(p3, p3_log_log), cols=2)
