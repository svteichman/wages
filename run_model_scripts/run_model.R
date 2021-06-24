library(rstan) 

load(file = "data/data_no_s6.rda")

curr_file <- "../Model/stan files/wages_simple6-18-21.stan"

fit_test <- stan( 
   file = curr_file,  # Stan program
   data = data_no_s6,    # named list of data
   chains = 4,             # number of Markov chains
   warmup = 2,          # number of warmup iterations per chain
   iter = 10,            # total number of iterations per chain
   cores = 2              # number of cores (could use one per chain)
)

#save(fit_test, file = 'fit_test2-4.rda')
