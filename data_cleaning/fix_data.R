library(geoR)
library(arm)
library(rstan)
library(raster)

load(file = "stan_dat.rda")

#take parameters from data (N,K,T)
K <- stan_dat$K
N <- stan_dat$N
Ti <- stan_dat$T

#hyperparameter for inverse chi square draw 
hp <- 5

# simulate data for latent space positions ----
sigmasq_w0 <- rinvchisq(1, hp)
sigmasq_w <- rinvchisq(1, 25)
sigmasq_z <- rinvchisq(1, hp)

#drawing latent positions
w <- array(data = NA, dim = c(K,Ti,2))
w[,1,] <- rnorm(K*2,0,sqrt(sigmasq_w0))
for (t in 2:Ti) {
  w[,t,1] <- rnorm(K,w[,t-1,1],sqrt(sigmasq_w))
  w[,t,2] <- rnorm(K,w[,t-1,2],sqrt(sigmasq_w))
}
z <- rnorm(N*2,0,sqrt(sigmasq_z))
dim(z) <- c(N,2)

#redo data vector to include latent positions 
stan_dat_latent <- list(
  N = N,
  T = time,
  K = K,
  Xnum = Xnum,
  y = y, 
  x = x,
  a = a,
  r = r,
  s = s,
  first_class_means = first_class_means,
  w = w,
  z = z
)

#use changed model for latent positions as data
latent_file <- "../../Model/wages_simple_latent_data.stan"

#test rstan 
fit_test <- stan(
  file = latent_file,  # Stan program
  data = stan_dat_latent,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 20,          # number of warmup iterations per chain
  iter = 100,            # total number of iterations per chain
  cores = 2              # number of cores (could use one per chain)
)



# simulate data for class probabilities ----
sigmasq_bk <- rinvchisq(1, hp)
sigmasq_a <- rinvchisq(1, hp)

#drawing probabilities
x <- stan_dat$x
w <- stan_dat_latent$w
z <- stan_dat_latent$z
beta_k <- rnorm(K,0,sqrt(sigmasq_bk))
alpha <- rnorm(2,0,sqrt(sigmasq_a))
adj_alpha <- array(data = NA, dim = 2)
adj_alpha[1] <- alpha[1]/mean(x[,,1])
adj_alpha[2] <- alpha[2]/mean(x[,,2])
p <- array(0,dim=c(N,Ti,K))
for (n in 1:N) {
  for (t in 1:Ti) {
    for (k in 1:K) {
      p[n,t,k] <- invlogit(beta_k[k] - pointDistance(z[n,],w[k,t,],longlat = FALSE) + t(adj_alpha)%*%x[n,t,]) 
    }
    p[n,t,] <- p[n,t,]/sum(p[n,t,])
  }
}

#redo data vector to include latent positions 
stan_dat_probs <- list(
  N = N,
  T = time,
  K = K,
  Xnum = Xnum,
  y = y, 
  x = x,
  a = a,
  r = r,
  s = s,
  first_class_means = first_class_means,
  w = w,
  z = z,
  p = p
)

#use changed model for latent positions as data
probs_file <- "../../Model/wages_simple_probs_data.stan"

#test rstan 
fit_test <- stan(
  file = probs_file,  # Stan program
  data = stan_dat_probs,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 20,          # number of warmup iterations per chain
  iter = 100,            # total number of iterations per chain
  cores = 2              # number of cores (could use one per chain)
)


# remove sectors with low counts (<20) - 0, 2, 4, 10 ---- 
load(file = "wave3to5_filt.rda")
rm <- c(0,2,4,10)
sec_rm <- unique(c(which(wave3to5_filt$w3_sector %in% rm), 
            which(wave3to5_filt$w4_sector %in% rm), 
            which(wave3to5_filt$w5_sector %in% rm)))
wave3to5_filt_sec <- wave3to5_filt[-sec_rm,]
save(wave3to5_filt_sec, file = "wave3to5_filt_sec.rda")

# rest copied from transform_data.R
N <- nrow(wave3to5_filt_sec)
time <- 3

#Start by using sector as class 
K <- 7
Xnum <- 2

#create X array
sex <- grep("sex", names(wave3to5_filt_sec))
age <- grep("_age", names(wave3to5_filt_sec))
x <- array(data = NA, dim = c(N,time,2))
for (i in 1:N) {
  x[i, ,1] <- as.numeric(wave3to5_filt_sec[i,age])
  x[i, ,2] <- as.numeric(wave3to5_filt_sec[i,sex])
}

#create Y array
wage <- grep("wage_yr", names(wave3to5_filt_sec))
y <- array(data = NA, dim = c(N,time))
for (i in 1:N) {
  y[i,] <- as.numeric(wave3to5_filt_sec[i,wage])
}

#convert sectors to indicator of job category 
sec <- grep("sec", names(wave3to5_filt_sec))
job <- array(data = NA, dim = c(N,time)) 
for (i in 1:N) {
  job[i,] <- as.numeric(wave3to5_filt_sec[i,sec])
}
a <- array(data = 0, dim = c(N,time,K))
inc <- c(1,3,5,6,7,8,9)
for (i in 1:N) {
  for (ti in 1:time) {
    sector <- wave3to5_filt_sec[i,sec[ti]]
    #sectors are coded {1, 3, 5, 6, 7, 8, 9}, store in values 1 to 7
    ind <- which(inc == as.numeric(sector), arr.ind = T)
    a[i,ti,ind] <- 1
  }
}

#make remaining at the same job array
r <- array(data = 0, dim = c(time-1,N))
for (i in 1:N) {
  if (wave3to5_filt_sec$w4_yr_curr_job[i] > 7) {
    r[1,i] = 1
  } 
  if (wave3to5_filt_sec$w5_yr_curr_job[i] > 7) {
    r[2,i] = 1
  }
}

#make remaining in sector but switching jobs array 
s <- 1-r
for (i in 1:N) {
  for (ti in 1: (time-1)) {
    if (s[ti,i] == 1) {
      if (!(job[i,ti] == job[i,ti+1])) {
        s[ti,i] = 0
      }
    }
  }
}

#average wages for each class in first time period
first_class_means <- array(0,dim = K)
for (k in 1:K) {
  gr <- which(job[,1]==k)
  first_class_means[k] <- mean(log(y[gr,1]))
}
fix <- which(is.na(first_class_means) | first_class_means == -Inf)
#replace classes that are NA or have mean 0 in the first time period with the mean of means
first_class_means[fix] <- mean(first_class_means[is.finite(first_class_means)],na.rm=T)

stan_dat_rm_sec <- list(
  N = N,
  T = time,
  K = K,
  Xnum = Xnum,
  y = y, 
  x = x,
  a = a,
  r = r,
  s = s,
  first_class_means = first_class_means
)

save(stan_dat_rm_sec, file = "stan_dat_rm_sec.rda")

curr_file <- "../../Model/wages_simple1-22.stan"

fit_test <- stan(
  file = curr_file,  # Stan program
  data = stan_dat_rm_sec,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 20,          # number of warmup iterations per chain
  iter = 100,            # total number of iterations per chain
  cores = 2              # number of cores (could use one per chain)
)

save(fit_test, file = "fit_test_res.rda")

# simulate data for latent space positions with filtered data (fewer sectors) ---- 

#take parameters from data (N,K,T)
K <- stan_dat_rm_sec$K
N <- stan_dat_rm_sec$N
Ti <- stan_dat_rm_sec$T

#hyperparameter for inverse chi square draw 
hp <- 5

sigmasq_w0 <- rinvchisq(1, hp)
sigmasq_w <- rinvchisq(1, 25)
sigmasq_z <- rinvchisq(1, hp)

#drawing latent positions
w <- array(data = NA, dim = c(K,Ti,2))
w[,1,] <- rnorm(K*2,0,sqrt(sigmasq_w0))
for (t in 2:Ti) {
  w[,t,1] <- rnorm(K,w[,t-1,1],sqrt(sigmasq_w))
  w[,t,2] <- rnorm(K,w[,t-1,2],sqrt(sigmasq_w))
}
z <- rnorm(N*2,0,sqrt(sigmasq_z))
dim(z) <- c(N,2)

#redo data vector to include latent positions 
stan_dat_latent_rm_sec <- stan_dat_rm_sec
stan_dat_latent_rm_sec$w <- w
stan_dat_latent_rm_sec$z <- z

#use changed model for latent positions as data
latent_file <- "../../Model/wages_simple_latent_data.stan"

#test rstan 
fit_test <- stan(
  file = latent_file,  # Stan program
  data = stan_dat_latent_rm_sec,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 20,          # number of warmup iterations per chain
  iter = 100,            # total number of iterations per chain
  cores = 2              # number of cores (could use one per chain)
)

# simulate data for class probabilities with filtered data (fewer sectors) ----

#take parameters from data (N,K,T)
K <- stan_dat_rm_sec$K
N <- stan_dat_rm_sec$N
Ti <- stan_dat_rm_sec$T

#hyperparameter for inverse chi square draw 
hp <- 5

sigmasq_bk <- rinvchisq(1, hp)
sigmasq_a <- rinvchisq(1, hp)

#drawing probabilities
x <- stan_dat_rm_sec$x
w <- stan_dat_latent_rm_sec$w
z <- stan_dat_latent_rm_sec$z
beta_k <- rnorm(K,0,sqrt(sigmasq_bk))
alpha <- rnorm(2,0,sqrt(sigmasq_a))
adj_alpha <- array(data = NA, dim = 2)
adj_alpha[1] <- alpha[1]/mean(x[,,1])
adj_alpha[2] <- alpha[2]/mean(x[,,2])
p <- array(0,dim=c(N,Ti,K))
for (n in 1:N) {
  for (t in 1:Ti) {
    for (k in 1:K) {
      p[n,t,k] <- invlogit(beta_k[k] - pointDistance(z[n,],w[k,t,],longlat = FALSE) + t(adj_alpha)%*%x[n,t,]) 
    }
    p[n,t,] <- p[n,t,]/sum(p[n,t,])
  }
}

#redo data vector to include latent positions 
stan_dat_probs_rm_sec <- stan_dat_latent_rm_sec
stan_dat_probs_rm_sec$p <- p

#use changed model for latent positions as data
probs_file <- "../../Model/wages_simple_probs_data.stan"

#test rstan 
fit_test <- stan(
  file = probs_file,  # Stan program
  data = stan_dat_probs_rm_sec,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 20,          # number of warmup iterations per chain
  iter = 100,            # total number of iterations per chain
  cores = 2              # number of cores (could use one per chain)
)
