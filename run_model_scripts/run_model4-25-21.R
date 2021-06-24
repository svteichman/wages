library(rstan)

load(file = "expanded_data1.rda")
stan_dat <- expanded_data1

curr_file <- "wages_simple4-25-21.stan"

fit_ifls <- stan(
  file = curr_file,  # Stan program
  data = stan_dat,  # named list of data
  chains = 4,       # number of Markov chains
  warmup = 500,    # number of warmup iterations per chain
  iter = 25000,    # total number of iterations per chain
  cores = 4,       # number of cores (could use one per chain)
  thin = 5        # only save every 5th draw 
)
save(fit_ifls, file = "fit_ifls_4-25-21.rda")