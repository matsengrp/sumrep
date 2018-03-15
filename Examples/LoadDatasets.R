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
