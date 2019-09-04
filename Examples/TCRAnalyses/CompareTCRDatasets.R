loadNewDatasets("data/Annotations",
                pattern="A")

data_to_compare <- list(
    # Compage IGoR annotations to simulations
    c("A5_S22", "A5_S22_sim"),
    c("A5_S9", "A5_S9_sim"),
    c("A5_S10", "A5_S10_sim"),
    c("A5_S15", "A5_S15_sim"),
    c("A4_i194", "A4_i194_sim"),
    c("A4_i107", "A4_i107_sim"),
                        
    # Compare IGoR annotations to annotations
    c("A5_S22", "A5_S9"),
    c("A5_S22", "A5_S10"),
    c("A5_S22", "A5_S15"),
    c("A5_S22", "A4_i194"),
    c("A5_S22", "A4_i107"),

    c("A5_S9", "A5_S10"),
    c("A5_S9", "A5_S15"),
    c("A5_S9", "A4_i194"),
    c("A5_S9", "A4_i107"),

    c("A5_S10", "A5_S15"),
    c("A5_S10", "A4_i194"),
    c("A5_S10", "A4_i107"),

    c("A5_S15", "A4_i194"),
    c("A5_S15", "A4_i107"),

    c("A4_i194", "A4_i107"),

    # Compare IGoR simulations to simulations
    c("A5_S22_sim", "A5_S9_sim"),
    c("A5_S22_sim", "A5_S10_sim"),
    c("A5_S22_sim", "A5_S15_sim"),
    c("A5_S22_sim", "A4_i194_sim"),
    c("A5_S22_sim", "A4_i107_sim"),

    c("A5_S9_sim", "A5_S10_sim"),
    c("A5_S9_sim", "A5_S15_sim"),
    c("A5_S9_sim", "A4_i194_sim"),
    c("A5_S9_sim", "A4_i107_sim"),

    c("A5_S10_sim", "A5_S15_sim"),
    c("A5_S10_sim", "A4_i194_sim"),
    c("A5_S10_sim", "A4_i107_sim"),

    c("A5_S15_sim", "A4_i194_sim"),
    c("A5_S15_sim", "A4_i107_sim"),

    c("A4_i194_sim", "A4_i107_sim")
)

# Only the IGoR simulations have a 'vj_in_frame' column since IGoR
#   doesn't include this in annotations.
# So, just check if there is such a column and subset to in-frame sequences
#   when possible.
subsetToInFrame <- function(dat) {
    if("vj_in_frame" %in% names(dat[['annotations']])) {
        dat[['annotations']] <- 
            dat[['annotations']][dat[['annotations']][["vj_in_frame"]], ]
    }

    return(dat)
}

for(dats in data_to_compare) {
    comparison_name <- dats %>% 
        paste(collapse="_") %>%
        paste0("compare_", .)
    if(!exists(comparison_name)) {
        dat_1 <- eval(parse(text=dats[1])) %>%
            subsetToInFrame
        dat_2 <- eval(parse(text=dats[2])) %>%
            subsetToInFrame
        comparison <-
            compareRepertoires(dat_1, 
                               dat_2,
                               locus="trb",
                               do_full_comparison=FALSE
                              ) %>%
            cbind(Type1=dats[1],
                  Type2=dats[2])
        saveRDS(comparison, 
                file.path("data/Comparisons", 
                          comparison_name %>% paste0('.rds')))
    }
}

