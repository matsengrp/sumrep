loadNewDatasets("data/Annotations")

data_to_compare <- list(
    # Compage IGoR annotations to simulations
    c("hsa_tra", "hsa_tra_sim"),
    c("hsa_trb", "hsa_trb_sim"),
    c("mmu_tra", "mmu_tra_sim"),
    c("mmu_trb", "mmu_trb_sim"),

    # Compare IGoR annotations to annotations
    c("hsa_tra", "mmu_tra"),
    c("hsa_trb", "mmu_trb"),

    # Compare IGoR simulations to simulations
    c("hsa_tra_sim", "mmu_tra_sim"),
    c("hsa_trb_sim", "mmu_trb_sim")
)

for(dats in data_to_compare) {
    comparison_name <- dats %>% 
        paste(collapse="_") %>%
        paste0("compare_", .)
    if(!exists(comparison_name)) {
        comparison <-
            compareRepertoires(eval(parse(text=dats[1])), 
                               eval(parse(text=dats[2])),
                               receptor_type="TCR",
                               chain_type=ifelse(grepl("tra", dats[1]) && 
                                                 grepl("tra", dats[2]),
                                                 "alpha",
                                                 "beta"
                                                ),
                               do_full_comparison=FALSE
                              ) %>%
            cbind(Type1=dats[1],
                  Type2=dats[2])
        saveRDS(comparison, 
                file.path("data/Comparisons", 
                          comparison_name %>% paste0('.rds')))
    }
}

