# simulation2-setup.R
library(tidyverse)

parameters <-
    list( m = 100
        , n = 500
        , magnitude = c(0.5, 1.0)
        , rho =  c(0.5, 0.8)
        ) %>% cross_df() %>%
    mutate( beta  = map(magnitude, ~sample(.x*c(rep(c(-1,1),each=5),c(rep(0,40)))))
          , Sigma = map(rho, ~.x^abs(outer( 1:50, 1:50, `-`)))
          )
saveRDS(parameters, 'parameters.rds')
