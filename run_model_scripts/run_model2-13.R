library(rstan)

load(file = "stan_dat_rm_sec.rda")
#load(file = "simulated_wage_dat2-12.rda")
stan_dat <- stan_dat_rm_sec
#stan_dat <- wage_data

curr_file <- "wages_simple1-30.stan"

fit_ifls <- stan(
  file = curr_file,  # Stan program
  data = stan_dat,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 500,          # number of warmup iterations per chain
  iter = 2500,            # total number of iterations per chain
  cores = 4              # number of cores (could use one per chain)
)
save(fit_ifls, file = "fit_ifls_2-13.rda")