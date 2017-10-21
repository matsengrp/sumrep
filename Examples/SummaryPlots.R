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
    c1 <- compareRepertoires(dat, simu)
    c2 <- compareRepertoires(dat, boot)
    c3 <- compareRepertoires(simu, boot)
}

if(!("Type1" %in% names(c1))) {
    c1$Type1 <- "Obs"
    c1$Type2 <- "Sim"

    c2$Type1 <- "Obs"
    c2$Type2 <- "Boot"

    c3$Type1 <- "Sim"
    c3$Type2 <- "Boot"
}

cfull <- rbind(c1, c2, c3)

comparison_types <- c1$Comparison

plot_list <- list()
for(c_type in comparison_types) {
    csub <- cfull[cfull$Comparison == c_type, ]
    plot_list[[c_type]] <- ggplot(csub, aes(Type1, Type2)) +
        geom_tile(aes(fill=Value)) +
        ggtitle(c_type %>% gsub(pattern="compare", replace="")
            %>% gsub(pattern="Distribution", replace="")) +
        theme(
              legend.key.size=unit(0.3, "cm"),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              plot.title=element_text(size=10)
              )
}

multiplot(plotlist=plot_list, cols=4, rows=6)
