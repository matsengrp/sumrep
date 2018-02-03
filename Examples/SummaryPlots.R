library(viridis)

multiplot <- function(..., plotlist=NULL, file, cols=1, rows=1, layout=NULL) {
  library(grid)

  plots <- c(list(...), plotlist)
  numPlots = length(plots)

  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * rows),
                    ncol = cols, nrow = rows)
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

if(!exists("c1")) {
    c1 <- compare_partis_fv1_partis_fv1_sim
    c2 <- compare_partis_fv2_partis_fv2_sim
    c3 <- compare_partis_gmc1_partis_gmc1_sim
    c4 <- compare_igb_fv1_partis_fv1_sim
    c5 <- compare_igb_fv2_partis_fv2_sim
    c6 <- compare_igb_gmc1_partis_gmc1_sim
}

if(!("Type1" %in% names(c1))) {
    c1$Type1 <- "F1h-p"
    c1$Type2 <- "F1h-p-sim"

    c2$Type1 <- "F8d-p"
    c2$Type2 <- "F8d-p-sim"

    c3$Type1 <- "G1h-p"
    c3$Type2 <- "G1h-p-sim"

    c4$Type1 <- "F1h-i"
    c4$Type2 <- "F1h-p-sim"

    c5$Type1 <- "F8d-i"
    c5$Type2 <- "F8d-p-sim"

    c6$Type1 <- "G1h-i"
    c6$Type2 <- "G1h-p-sim"
}

cfull <- rbind(c1, c2, c3, c4, c5, c6)
cfull$Divergence <- cfull$Divergence %>% unlist


# cfull$Comparison <- factor(cfull$Comparison, levels=ordered_strings)

comparison_types <- c1$Comparison


getComparisonValue <- function(dat, comparison) {
    value <- dat %>% 
        filter(Comparison == comparison) %$%
        Value
    return(value)
}

scoreStatistics <- function(sim_dats, obs_dats) {
    score_dat <- matrix(NA, nrow=0, ncol=2) %>% 
        data.table %>%
        setNames(c("Comparison", "Score"))
    for(c_type in comparison_types) {
            sim_score <- sim_dats %>% 
                sapply(getComparisonValue, c_type) %>% 
                unlist %>%
                mean
            obs_score <- obs_dats %>%
                sapply(getComparisonValue, c_type) %>%
                unlist %>%
                mean
            cat(sim_score, obs_score, '\n')
            c_score <- sim_score/obs_score
            score_dat <- rbind(score_dat,
                               data.table(Comparison=shortenName(c_type),
                                          Score=c_score))
    }
    return(score_dat)
}

shortenName <- function(string) {
    shortened_name <- string %>% 
        gsub(pattern="compare", replace="") %>%
        gsub(pattern="Distribution", replace="")
    return(shortened_name)
}

partis_dats <- list(c1, c2, c3)
igb_dats <- list(c4, c5, c6)

#score_dat <- scoreStatistics(partis_dats, igb_dats)
score_dat <- score_dat[order(score_dat$Score)]
score_plot <- score_dat %>% ggplot(aes(x=reorder(Comparison, Score), 
                                       y=Score)) +
                                       # fill=reorder(Comparison, Score))) +
    geom_bar(stat="identity") +
    xlab("Comparison") + ylab("Relative deviance") +
    ggtitle(paste0("Comparison scores of observations-to-simulations with ",
                   "respect to comparisons of observations-to-observations")) +
    theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1),
          plot.margin=unit(c(1, 1, 1, 1.5), "cm"))

# Multi-plots
do.multiplot <- TRUE
if(do.multiplot) {
    plot_list <- list()
    for(c_type in comparison_types) {
        csub <- cfull[cfull$Comparison == c_type, ]
        plot_list[[c_type]] <- ggplot(csub, aes(Type1, Type2)) +
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
    
    pdf("multiplot.pdf", width=34, height=16)
    multiplot(plotlist=plot_list, cols=5, rows=6)
    dev.off()
}
