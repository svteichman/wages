library(geoR)
library(arm)
library(rstan)
library(raster)
library(invgamma)

load("stan_dat_rm_sec.rda")
load("init_vals_sec_rm.rda")
load("noisy_init_vals_sec_rm.rda")

dat <- stan_dat_rm_sec
N <- dat$N
K <- dat$K
Ti <- dat$T

curr_file <- "../../Model/wages_simple11-5.stan"

init_sim <- list(init_1,init_1,init_1,init_1)

fit_test <- stan(
  file = curr_file,  # Stan program
  data = stan_dat_rm_sec,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 20,          # number of warmup iterations per chain
  iter = 100,            # total number of iterations per chain
  cores = 2,              # number of cores (could use one per chain)
  init = init_sim
)

# evaluate log-likelihood by hand from initial values 

#latent space model
log(dinvchisq(init_1$sigmasq_w0, 5))
log(dinvchisq(init_1$sigmasq_w, 25))
log(dinvchisq(init_1$sigmasq_z, 5))
w_log_lik <- array(data = NA, dim = c(dat$K, dat$T, dat$Xnum))
for (i in 1:K) {
  w_log_lik[i,1,1] <- log(dnorm(init_1$w[i,1,1],0,init_1$sigmasq_w0))
  w_log_lik[i,1,2] <- log(dnorm(init_1$w[i,1,2],0,init_1$sigmasq_w0))
  for (t in 2:Ti) {
    w_log_lik[i,t,1] <- log(dnorm(init_1$w[i,t,1],init_1$w[i,t-1,1],init_1$sigmasq_w))
    w_log_lik[i,t,2] <- log(dnorm(init_1$w[i,t,2],init_1$w[i,t-1,2],init_1$sigmasq_w)) 
  }
}
z_log_lik <- array(data = NA, dim = c(N,2))
for (i in 1:N) {
  z_log_lik[i,1] <- log(dnorm(init_1$z[i,1],0,init_1$sigmasq_z))
  z_log_lik[i,2] <- log(dnorm(init_1$z[i,2],0,init_1$sigmasq_z))
}

#probabilities and multinomial job selection 
log(dinvchisq(init_1$sigmasq_bk, 5))
beta_k_log_lik <- array(data = NA, dim = K)
for (i in 1:K) {
  beta_k_log_lik[i] <- log(dnorm(init_1$beta_k[i],0,init_1$sigmasq_bk))
}
log(dinvchisq(init_1$sigmasq_a, 5))
log(dnorm(init_1$alpha,0,init_1$sigmasq_a))
log(dinvchisq(init_1$sigmasq_d,13))
log(dnorm(init_1$delta,0,init_1$sigmasq_d))
log(dinvchisq(init_1$sigmasq_b,13))
log(dnorm(init_1$beta,0,init_1$sigmasq_b))
log(dinvchisq(init_1$sigmasq_u,3))
log(dnorm(init_1$u,0,init_1$sigmasq_b))
log(dinvchisq(init_1$sigmasq_v,3))
log(dnorm(init_1$v,0,init_1$sigmasq_b))

#generate p values from data & parameters so far 
adj_alpha <- array(data = NA, dim = 2)
adj_alpha[1] <- init_1$alpha[1]/mean(dat$x[,,1])
adj_alpha[2] <- init_1$alpha[2]/mean(dat$x[,,2])
p <- array(0,dim=c(N,Ti,K))
for (n in 1:N) {
  for (t in 1:Ti) {
    for (k in 1:K) {
      p[n,t,k] <- invlogit(init_1$beta_k[k] - pointDistance(init_1$z[n,],init_1$w[k,t,],longlat = FALSE) 
                           + t(adj_alpha)%*%dat$x[n,t,]) 
    }
    p[n,t,] <- p[n,t,]/sum(p[n,t,])
  }
}

#multinomial likelihood when there is only 1 draw is just the probability of that category
sum(p < 10^(-16)) #none close to 0 

#generate eta and lambda values from data and parameters so far 
job <- array(0,dim=c(N,T))
for (n in 1:N) { 
  for (t in 1:Ti) {
    job[n,t] <- which(dat$a[n,1,]==1)
  }
}
adj_beta <- array(data = NA, dim = 2)
adj_beta[1] <- init_1$beta[1]/mean(dat$x[,,1])
adj_beta[2] <- init_1$beta[2]/mean(dat$x[,,2])
adj_delta <- array(data = NA, dim = 2)
adj_delta[1] <- init_1$delta[1]/mean(dat$x[,,1])
adj_delta[2] <- init_1$delta[2]/mean(dat$x[,,2])
eta <- array(0,dim = c(N,T-1))
lambda <- array(0,dim = c(N,T-1))
for (t in 2:T) {
  for (n in 1:N) {
    eta[n,t-1] <- invlogit(t(adj_beta)%*%x[n,t,] + v[job[n,t-1]])
    lambda[n,t-1] <- invlogit(t(adj_delta)%*%x[n,t,] + u[job[n,t-1]])
  }
}
    r[t-1,n] <- rbern(n = 1, prob = invlogit(t(beta)%*%x[n,t,] + v[job[n,t-1]]))
    s[t-1,n] <- (1-r[t-1,n])*rbern(n = 1, prob = invlogit(t(delta)%*%x[n,t,] + u[job[n,t-1]]))
    l[t-1,n] <- (1-r[t-1,n])*(1-s[t-1,n])
    if (r[t-1,n]+s[t-1,n]==1) {
      temp <- array(0,dim=K)
      temp[job[n,t-1]] <- 1
      a[n,t,] <- temp
      job[n,t] <- job[n,t-1]
    }
    else {
      a[n,t,] <- rmultinom(n = 1, size = 1, prob = p[n,,t])
      job[n,t] <- which(a[n,t,]==1)
    }
  }
}


