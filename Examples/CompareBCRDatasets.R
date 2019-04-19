loadNewDatasets("data/Annotations")

do_bcr <- TRUE
if(do_bcr) {
    data_to_compare <- list(
                            # Compare partis annotations to
                            # simulations using imgt database
                            c("pi_f1", "pi_f1_sim"),
                            c("pi_f2", "pi_f2_sim"),
                            c("pi_g1", "pi_g1_sim"),
                            c("pi_g2", "pi_g2_sim"),
                            c("pi_i1", "pi_i1_sim"),
                            c("pi_i2", "pi_i2_sim"),
    
                            # Compare partis annotations to simulations
                            c("p_f1", "p_f1_sim"),
                            c("p_f2", "p_f2_sim"),
                            c("p_g1", "p_g1_sim"),
                            c("p_g2", "p_g2_sim"),
                            c("p_i1", "p_i1_sim"),
                            c("p_i2", "p_i2_sim"),
    
                            # Compare partis annotations to annotations
                            c("p_f1", "p_f2"),
                            c("p_f1", "p_g1"),
                            c("p_f1", "p_g2"),
                            c("p_f1", "p_i1"),
                            c("p_f1", "p_i2"),
    
                            c("p_f2", "p_g1"),
                            c("p_f2", "p_g2"),
                            c("p_f2", "p_i1"),
                            c("p_f2", "p_i2"),
    
                            c("p_g1", "p_g2"),
                            c("p_g1", "p_i1"),
                            c("p_g1", "p_i2"),
    
                            c("p_g2", "p_i1"),
                            c("p_g2", "p_i2"),
                            
                            c("p_i1", "p_i2"),
    
                            # Compare partis simulations to simulations
                            c("p_f1_sim", "p_f2_sim"),
                            c("p_f1_sim", "p_g1_sim"),
                            c("p_f1_sim", "p_g2_sim"),
                            c("p_f1_sim", "p_i1_sim"),
                            c("p_f1_sim", "p_i2_sim"),
    
                            c("p_f2_sim", "p_g1_sim"),
                            c("p_f2_sim", "p_g2_sim"),
                            c("p_f2_sim", "p_i1_sim"),
                            c("p_f2_sim", "p_i2_sim"),
    
                            c("p_g1_sim", "p_g2_sim"),
                            c("p_g1_sim", "p_i1_sim"),
                            c("p_g1_sim", "p_i2_sim"),
    
                            c("p_g2_sim", "p_i1_sim"),
                            c("p_g2_sim", "p_i2_sim"),
                            
                            c("p_i1_sim", "p_i2_sim"),
    
                            # Compare partis annotations using imgt database
                            # to other partis annotations using imgt database
                            c("pi_f1", "pi_f2"),
                            c("pi_f1", "pi_g1"),
                            c("pi_f2", "pi_g1"),
    
                            # Compare igblast annotations to annotations
                            c("i_f1", "i_f2"),
                            c("i_f1", "i_g1"),
                            c("i_f2", "i_g1"),
    
                            # Compare igblast annotations to partis
                            # simulations using imgt database
                            c("i_f1", "pi_f1_sim"),
                            c("i_f2", "pi_f2_sim"),
                            c("i_g1", "pi_g1_sim")
                           )
}


for(dats in data_to_compare) {
    comparison_name <- dats %>% 
        paste(collapse="_") %>%
        paste0("compare_", .)
    if(!exists(comparison_name)) {
        dat_1 <- eval(parse(text=dats[1]))
        dat_2 <- eval(parse(text=dats[2]))
        comparison <-
               compareRepertoires(dat_1, 
                                  dat_2,
                                  locus="igh",
                                  do_full_comparison=FALSE
                                 ) %>%
               cbind(Type1=dats[1],
                     Type2=dats[2])
        saveRDS(comparison, 
                file.path("data/Comparisons", 
                          comparison_name %>% paste0('.rds')))
    }
}

