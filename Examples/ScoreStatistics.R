getComparisonValue <- function(dat, comparison) {
    value <- dat %>% 
        filter(Comparison == comparison) %$%
        Divergence
    return(value)
}

scoreStatistics <- function(dats_1, dats_2) {
    comparison_types <- dats_1[[1]]$Comparison %>% unique
    score_dat <- matrix(NA, nrow=0, ncol=2) %>% 
        data.table %>%
        setNames(c("Comparison", "Score"))
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
                gsub(pattern="compare", replace="get") %>%
                ifelse(str_sub(., -1) == "s", 
                       substr(., start=1, stop=nchar(.) - 1), 
                       .
                      ) %>%
                getNameFromFunctionString
            
            score_dat <- rbind(score_dat,
                               data.table(
                                          Comparison=comparison_name,
                                          Score=c_score
                                         )
                              )
    }
    return(score_dat)
}

shortenName <- function(string) {
    shortened_name <- string %>% 
        gsub(pattern="compare", replace="get") %>%
        getNameFromFunctionString
    return(shortened_name)
}


plotComparisons <- function(dat, filename, cols=1, rows=1) {
    comparison_types <- dat$Comparison %>% unique
    plot_list <- list()
    for(c_type in comparison_types) {
        dat_sub <- dat[dat$Comparison == c_type, ]
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

plotSummaryScores <- function(
                              dats_1, 
                              dats_2,
                              filename,
                              plot_width=12,
                              plot_height=6
                             ) {
    score_dat <- scoreStatistics(dats_1, dats_2)
    score_dat <- score_dat[order(score_dat[["Score"]])]
    score_plot <- score_dat %>%
        ggplot(aes(x=reorder(Comparison, -Score), 
                   y=Score
                  )
              ) +
        geom_bar(stat="identity") +
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