library(rstan)

load(file = "deflate_y_dat.rda")
deflate_y_dat <- list(deflate_y_dat$N,
                      deflate_y_dat$T,
                      deflate_y_dat$K,
                      deflate_y_dat$y,
                      deflate_y_dat$a,
                      deflate_y_dat$r,
                      mean(log(deflate_y_dat$y)))
names(deflate_y_dat) <- c("N","T","K","y","a","r","mean_log_wage")
stan_dat <- deflate_y_dat

curr_file <- "wages_simple2-19-21.stan"

fit_ifls <- stan(
  file = curr_file,  # Stan program
  data = stan_dat,  # named list of data
  chains = 4,       # number of Markov chains
  warmup = 500,    # number of warmup iterations per chain
  iter = 25000,    # total number of iterations per chain
  cores = 4,       # number of cores (could use one per chain)
  thin = 5        # only save every 5th draw 
)
save(fit_ifls, file = "fit_ifls_2-19-21.rda")
