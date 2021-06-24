library(rstan)

load(file = "deflate_y_dat.rda")
deflate_y_dat <- list(deflate_y_dat$N,
                      deflate_y_dat$T,
                      deflate_y_dat$K,
                      deflate_y_dat$Xnum,
                      deflate_y_dat$y,
                      deflate_y_dat$x,
                      deflate_y_dat$a,
                      deflate_y_dat$first_class_means)
names(deflate_y_dat) <- c("N","T","K","Xnum","y","x","a","first_class_means")
stan_dat <- deflate_y_dat

curr_file <- "wages_simple7-3.stan"

fit_ifls <- stan(
  file = curr_file,  # Stan program
  data = stan_dat,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 500,          # number of warmup iterations per chain
  iter = 2500,            # total number of iterations per chain
  cores = 4              # number of cores (could use one per chain)
)
save(fit_ifls, file = "fit_ifls_7-3.rda")