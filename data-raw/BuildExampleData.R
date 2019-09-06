# Generate example data objects

# Load data
test_dat <- readRDS(file.path("data-raw", "test_dat.rds"))
test_dat_boot <- readRDS(file.path("data-raw", "test_dat_boot.rds"))
test_simu <- readRDS(file.path("data-raw", "test_simu.rds"))

# Save data
usethis::use_data(test_dat, overwrite=TRUE)
usethis::use_data(test_dat_boot, overwrite=TRUE)
usethis::use_data(test_simu, overwrite=TRUE)