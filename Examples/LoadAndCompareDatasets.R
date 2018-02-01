data_dir <- "data/Annotations"

# Load datasets that are not already in the workspace
for(data_file in list.files(data_dir)) {
    var_name <- data_file %>%
        gsub(pattern="-", replace="_") %>%
        gsub(pattern=".rds", replace="")
    if(!exists(var_name)) {
        assign(var_name, 
            readRDS(file.path(data_dir, data_file)))
    }
}

data_to_compare <- list(
                        c("partis_fv1", "partis_fv1_sim"),
                        c("igb_fv1", "partis_fv1_sim"),
                        c("partis_fv2", "partis_fv2_sim"),
                        c("igb_fv2", "partis_fv2_sim"),
                        c("partis_gmc1", "partis_gmc1_sim"),
                        c("igb_gmc1", "partis_gmc1_sim")
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
