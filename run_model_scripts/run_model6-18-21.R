library(rstan)

load(file = "data_no_s6.rda")
stan_dat <- data_no_s6

curr_file <- "wages_simple6-18-21.stan"

fit_ifls <- stan(
  file = curr_file,  # Stan program
  data = stan_dat,  # named list of data
  chains = 4,       # number of Markov chains
  warmup = 500,    # number of warmup iterations per chain
  iter = 25000,    # total number of iterations per chain
  cores = 1,       # number of cores (could use one per chain)
  thin = 5        # only save every 5th draw 
)
save(fit_ifls, file = "fit_ifls_6-18-21.rda")