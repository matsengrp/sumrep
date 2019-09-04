devtools::load_all()

# Need to pull the data from zenodo, e.g. via the following:
# wget https://zenodo.org/record/3381680/files/flu_rds.tar
# tar -C /path/to/sumrep/data/flu -xvf flu_rds.tar 
flu_dir <- "data/flu"
p_f1 <- readRDS(file.path(flu_dir, "p_fv_igh_m1h.rds"))
p_f2 <- readRDS(file.path(flu_dir, "p_fv_igh_m8d.rds"))
p_f7 <- readRDS(file.path(flu_dir, "p_fv_igh_p7d.rds"))
p_f28 <- readRDS(file.path(flu_dir, "p_fv_igh_p28d.rds"))
p_g1 <- readRDS(file.path(flu_dir, "p_gmc_igh_m1h.rds"))
p_g2 <- readRDS(file.path(flu_dir, "p_gmc_igh_m8d.rds"))
p_g7 <- readRDS(file.path(flu_dir, "p_gmc_igh_p7d.rds"))
p_g28 <- readRDS(file.path(flu_dir, "p_gmc_igh_p28d.rds"))
p_i1 <- readRDS(file.path(flu_dir, "p_ib_igh_m1h.rds"))
p_i2 <- readRDS(file.path(flu_dir, "p_ib_igh_m8d.rds"))
p_i7 <- readRDS(file.path(flu_dir, "p_ib_igh_p7d.rds"))
p_i28 <- readRDS(file.path(flu_dir, "p_ib_igh_p28d.rds"))

dat_list <- list(p_f1, p_f2, p_f7, p_f28,
                 p_g1, p_g2, p_g7, p_g28,
                 p_i1, p_i2, p_i7, p_i28)
num_dats <- length(dat_list)
colors <- c(rep("purple", 4),
            rep("green", 4),
            rep("orange", 4))
ltys <- rep(c("-8d", "-1hr", "+7d", "+28d"), 4) %>% as.factor

div_matrix <- matrix(0, nrow=num_dats, ncol=num_dats)
for(i in 1:num_dats) {
    for(j in 1:num_dats) {
        div_matrix[i, j] = 
            sumrep::compareCDR3LengthDistributions(dat_list[[i]][["annotations"]],
                                                   dat_list[[j]][["annotations"]]
                                                  )
    }
}

pdf("mds.pdf")
div_matrix %>%
    cmdscale %>%
    plot(col=colors, 
         pch=as.integer(ltys), 
         xlab="Coordinate 1", 
         ylab="Coordinate 2", 
         main="MDS of CDR3 length divergences"
        )
legend(x=-0.022, 
       y=0.018, 
       legend=c("FV", "GMC", "IB"), 
       col=c("purple", "green", "orange"), 
       lty=1
      )
legend(x=-0.03, 
       y=0.018, 
       legend=c("1h", "2h", "7d", "28d"), 
       pch=c(1, 2, 3, 4)
      )
dev.off()
