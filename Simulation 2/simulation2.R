# simulation2.R
source("Simulation_Setup.R")


commandArgs()
(args <- commandArgs(trailingOnly=TRUE))

assert_that( length(args) == 1
           , msg='Expecting a command line argument for the set #.' )
assert_that( file.exists("parameters.rds")
           , msg='You must run simulation2-setup.R first to create the parameters object.'
           )

set <- as.integer(args[[1]])
message("I am running set", set)

params <- readRDS("parameters.rds")[set,] %>%
    select(m, n, beta, Sigma)
results <- pmap(params, run_simulation, analyze=analyze_with_enet)
saveRDS(results, file = paste0('results-', set, '.rds'))
