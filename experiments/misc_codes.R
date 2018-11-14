library(data.table)
# Get P1 PCs
P1pc <- fread('../LDS/Data/Ping_PCs.csv', col.names = paste0('PC', 1:12))
P1pc$year <- 1600:2005
usethis::use_data(P1pc)

P1annual <- fread('../LDS/Data/Ping_annual.csv')
usethis::use_data(P1annual)

P1model <- readRDS('experiments/lds1000.RDS')
P1model10k <- readRDS('experiments/lds10000.RDS')
P1cv <- readRDS('experiments/proper_cv_1000_restart.RDS')
usethis::use_data(P1model, P1model10k, P1cv, internal = TRUE, overwrite = TRUE)
