library(rstan)

load(file = "deflate_y_dat.rda")
stan_dat <- deflate_y_dat

curr_file <- "wages_simple4-11.stan"

fit_ifls <- stan(
  file = curr_file,  # Stan program
  data = stan_dat,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 500,          # number of warmup iterations per chain
  iter = 2500,            # total number of iterations per chain
  cores = 4              # number of cores (could use one per chain)
)
save(fit_ifls, file = "fit_ifls_4-11.rda")