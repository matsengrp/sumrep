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
    c1 <- compareRepertoires(obs_1, sim_1)
    c2 <- compareRepertoires(obs_2, sim_2)
    c3 <- compareRepertoires(obs_3, sim_3)
    c4 <- compareRepertoires(obs_1, obs_2)
    c5 <- compareRepertoires(obs_1, obs_3)
    c6 <- compareRepertoires(obs_2, obs_3)
}

if(!("Type1" %in% names(c1))) {
    c1$Type1 <- "FV-1h"
    c1$Type2 <- "FV-1h-sim"

    c2$Type1 <- "FV-8d"
    c2$Type2 <- "FV-8d-sim"

    c3$Type1 <- "GMC-1h"
    c3$Type2 <- "GMC-1h-sim"

    c4$Type1 <- "FV-1h"
    c4$Type2 <- "FV-8d"

    c5$Type1 <- "FV-1h"
    c5$Type2 <- "GMC-1h"

    c6$Type1 <- "FV-8d"
    c6$Type2 <- "GMC-1h"
}

cfull <- rbind(c1, c2, c3, c4, c5, c6)

# Reorder the levels for plotting 
cfull$Type2 <- factor(cfull$Type2, 
                      levels=c("FV-1h-sim",
                               "FV-8d-sim",
                               "GMC-1h-sim",
                               "FV-1h",
                               "FV-8d",
                               "GMC-1h"))

cfull$Comparison <- factor(cfull$Comparison, levels=ordered_strings)

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
                mean
            obs_score <- obs_dats %>%
                sapply(getComparisonValue, c_type) %>%
                mean
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

sim_dats <- list(c1, c2, c3)
obs_dats <- list(c4, c5, c6)

score_dat <- scoreStatistics(sim_dats, obs_dats)
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
do.multiplot <- FALSE
if(do.multiplot) {
    plot_list <- list()
    for(c_type in comparison_types) {
        csub <- cfull[cfull$Comparison == c_type, ]
        plot_list[[c_type]] <- ggplot(csub, aes(Type1, Type2)) +
            geom_tile(aes(fill=Value)) +
            ggtitle(c_type %>% shortenName) +
            theme(
                  legend.key.size=unit(0.8, "cm"),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank(),
                  plot.title=element_text(size=10)
                  ) +
            scale_fill_viridis(direction=-1)
    }
    
    multiplot(plotlist=plot_list, cols=7, rows=4)
}
