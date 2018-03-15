data_to_compare <- list(
                        # Compare partis annotations to
                        # simulations using imgt database
                        c("partis_igb_fv1", "partis_igb_fv1_sim"),
                        c("partis_igb_fv2", "partis_igb_fv2_sim"),
                        c("partis_igb_gmc1", "partis_igb_gmc1_sim"),

                        # Compare igblast annotations to partis
                        # simulations using imgt database
                        c("igb_fv1", "partis_igb_fv1_sim"),
                        c("igb_fv2", "partis_igb_fv2_sim"),
                        c("igb_gmc1", "partis_igb_gmc1_sim"),

                        # Compare partis annotations to simulations
                        c("partis_fv1", "partis_fv1_sim"),
                        c("partis_fv2", "partis_fv2_sim"),
                        c("partis_gmc1", "partis_gmc1_sim"),

                        # Compare partis annotations to annotations
                        c("partis_fv1", "partis_fv2"),
                        c("partis_fv1", "partis_gmc1"),
                        c("partis_fv2", "partis_gmc1"),

                        # Compare partis simulations to simulations
                        c("partis_fv1_sim", "partis_fv2_sim"),
                        c("partis_fv1_sim", "partis_gmc1_sim"),
                        c("partis_fv2_sim", "partis_gmc1_sim"),

                        # Compare igblast annotations to annotations
                        c("igb_fv1", "igb_fv2"),
                        c("igb_fv1", "igb_gmc1"),
                        c("igb_fv2", "igb_gmc1")
                        )

for(dats in data_to_compare) {
    comparison_name <- dats %>% 
        paste(collapse="_") %>%
        paste0("compare_", .)
    if(!exists(comparison_name)) {
        assign(comparison_name,
               compareRepertoires(eval(parse(text=dats[1])), 
                                  eval(parse(text=dats[2]))))
    }
}
