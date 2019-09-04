devtools::load_all()

# Need to pull the data from zenodo, e.g. via the following:
# wget https://zenodo.org/record/3381680/files/flu_rds.tar
# tar -C /path/to/sumrep/data/flu -xvf flu_rds.tar 
flu_dir <- "data/flu"
p_f_m1h <- readRDS(file.path(flu_dir, "p_fv_igh_m1h.rds"))
p_f_m8d <- readRDS(file.path(flu_dir, "p_fv_igh_m8d.rds"))
p_f_p7d <- readRDS(file.path(flu_dir, "p_fv_igh_p7d.rds"))
p_f_p28d <- readRDS(file.path(flu_dir, "p_fv_igh_p28d.rds"))
p_g_m1h <- readRDS(file.path(flu_dir, "p_gmc_igh_m1h.rds"))
p_g_m8d <- readRDS(file.path(flu_dir, "p_gmc_igh_m8d.rds"))
p_g_p7d <- readRDS(file.path(flu_dir, "p_gmc_igh_p7d.rds"))
p_g_p28d <- readRDS(file.path(flu_dir, "p_gmc_igh_p28d.rds"))
p_i_m1h <- readRDS(file.path(flu_dir, "p_ib_igh_m1h.rds"))
p_i_m8d <- readRDS(file.path(flu_dir, "p_ib_igh_m8d.rds"))
p_i_p7d <- readRDS(file.path(flu_dir, "p_ib_igh_p7d.rds"))
p_i_p28d <- readRDS(file.path(flu_dir, "p_ib_igh_p28d.rds"))

# Combine and label datasets
dat_list <- list(p_f_m1h, p_f_m8d, p_f_p7d, p_f_p28d,
                 p_g_m1h, p_g_m8d, p_g_p7d, p_g_p28d,
                 p_i_m1h, p_i_m8d, p_i_p7d, p_i_p28d)
num_dats <- length(dat_list)
colors <- c(rep("purple", 4),
            rep("green", 4),
            rep("orange", 4))
ltys <- rep(c("-8d", "-1hr", "+7d", "+28d"), 4) %>% as.factor

# Compare CDR3 length distribution of each pairwise dataset,
# and store in matrix of divergences (distances)
div_matrix <- matrix(0, nrow=num_dats, ncol=num_dats)
for(i in 1:num_dats) {
    for(j in 1:num_dats) {
        div_matrix[i, j] = 
            sumrep::compareCDR3LengthDistributions(dat_list[[i]][["annotations"]],
                                                   dat_list[[j]][["annotations"]]
                                                  )
    }
}

# Plot MDS coordinates 1 and 2 of div_matrix
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
