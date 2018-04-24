data_to_compare <- list(
                        # Compare partis annotations to
                        # simulations using imgt database
                        c("pi_f1", "pi_f1_sim"),
                        c("pi_f2", "pi_f2_sim"),
                        c("pi_g1", "pi_g1_sim"),

                        # Compare igblast annotations to partis
                        # simulations using imgt database
                        c("i_f1", "pi_f1_sim"),
                        c("i_f2", "pi_f2_sim"),
                        c("i_g1", "pi_g1_sim"),

                        # Compare partis annotations to simulations
                        c("p_f1", "p_f1_sim"),
                        c("p_f2", "p_f2_sim"),
                        c("p_g1", "p_g1_sim"),

                        # Compare partis annotations to annotations
                        c("p_f1", "p_f2"),
                        c("p_f1", "p_g1"),
                        c("p_f2", "p_g1"),

                        # Compare partis simulations to simulations
                        c("p_f1_sim", "p_f2_sim"),
                        c("p_f1_sim", "p_g1_sim"),
                        c("p_f2_sim", "p_g1_sim"),

                        # Compare partis annotations using imgt database
                        # to other partis annotations using imgt database
                        c("pi_f1", "pi_f2"),
                        c("pi_f1", "pi_g1"),
                        c("pi_f2", "pi_g1"),

                        # Compare igblast annotations to annotations
                        c("i_f1", "i_f2"),
                        c("i_f1", "i_g1"),
                        c("i_f2", "i_g1")
                        )

loadNewDatasets("data/Annotations")

for(dats in data_to_compare) {
    comparison_name <- dats %>% 
        paste(collapse="_") %>%
        paste0("compare_", .)
    if(!exists(comparison_name)) {
        comparison <-
               compareRepertoires(eval(parse(text=dats[1])), 
                                  eval(parse(text=dats[2]))) %>%
               cbind(Type1=dats[1],
                     Type2=dats[2])
        saveRDS(comparison, 
                file.path("data/Comparisons", 
                          comparison_name %>% paste0('.rds')))
    }
}

