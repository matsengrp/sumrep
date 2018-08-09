devtools::load_all()

loadNewDatasets("data/Annotations")

dat <- p_f1$annotations %>% subsample(5000)

times <- {}
distributions <- list()
tols <- c(1, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001)
smooth <- c(1, 1, 1, 1, 2, 3, 4, 5)

for(i in 1:length(tols)) {
    pt <- proc.time()
    distributions[[i]] <- getPairwiseDistanceDistribution(dat$cdr3s,
                                                          approximate=TRUE,
                                                          tol=tols[i])
    times[i] <- (proc.time() - pt)["elapsed"]
}

pt <- proc.time()
true_dist <- getPairwiseDistanceDistribution(dat$cdr3s,
                                             approximate=FALSE) 
true_time <- (proc.time() - pt)["elapsed"]
true_dat <- cbind(true_dist, "True", smooth) %>%
    data.frame %>%
    setNames(c("Value", "Setting", "Smoothness"))

dist_dat <- distributions %>%
    lapply(data.frame) %>%
    Map(cbind, ., paste("Tol = ", tols), 10) %>%
    lapply(setNames, c("Value", "Setting", "Smoothness")) %>%
    do.call(rbind, .) %>%
    rbind(., true_dat)

p1 <- ggplot(dist_dat,
       aes(x=as.numeric(Value), group=Setting, colour=Setting, adjust=Smoothness)) +
    geom_density() +
    xlim(0, 60)

p2 <- ggplot(dist_dat,
       aes(x=as.numeric(Value), group=Setting, colour=Setting, adjust=Smoothness)) +
    stat_ecdf() +
    xlim(0, 60)

time_dat <- data.frame(Time=times, Tolerance=tols)
p3 <- ggplot(d=time_dat, aes(x=log(Tolerance), y=Time)) +
    geom_point() +
    geom_hline(yintercept=true_time)
