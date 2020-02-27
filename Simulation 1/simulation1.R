# simulation1.R
source("Simulation_Setup.R")

p  <- 50
ps <- p/5
magnitude <- 1
rho <- 0.8

beta  <- sample(magnitude*c(rep(c(-1,1),each=ps/2),c(rep(0,p-ps))))
Sigma <- rho**abs(outer(seq_along(beta), seq_along(beta), `-`))

one_result <- simulate_one( 250, beta, Sigma
                          , analyze=analyze_with_enet
                          )

results <- run_simulation(m=1000, n=250, beta, sigma, analyze_with_enet)
