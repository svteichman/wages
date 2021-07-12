library(Rfast)
library(extraDistr)
library(rstan)
load(file = "data/deflate_y_dat.rda")
load("fit_ifls_3-9-21.rda")
set.seed(66666)
all_parm <- rstan::extract(fit_ifls)
N <- deflate_y_dat$N 
TT <- deflate_y_dat$T
K <- deflate_y_dat$K
## Prepare hyperparameters
e <- 5
# nu <- rnorm(K, 0, sqrt(e))
tau <- rgamma(1, 1, 1)
tau_b <-  rgamma(1, 0.5, 1)
tau_gamma <- rgamma(1, 0.5, 1)
# tau_k <- rgamma(K, 1, 1)
tau_z <- rgamma(1, 0.5, 1)
tau_w0 <- rgamma(1, 0.5, 1)
# tau_w <- rgamma(1, 0.5, 1)
tau_x0 <- rgamma(1, 0.5, 1)
# tau_x <- rgamma(1, 0.5, 1)

# # Using parameters from pre-trained model
nu <- colMeans(all_parm$nu)
tau <- mean(all_parm$tau)
tau_k <- colMeans(all_parm$tau_k)
tau_w <- mean(all_parm$tau_w)
tau_w0 <- mean(all_parm$tau_w0)

L <- rep(10, K) # For simplicity, 10 firms for each sector
# L <- sample(5: 10, K, replace = TRUE)
m_bar <- mean(log(deflate_y_dat$y))

# all_beta <- rmvnorm(3, mu = rep(0, K), sigma = 1/tau_b*diag(nrow = K))
# all_beta[, 1] <- 0
# all_beta[, 1] <- rep(all_beta[1, 1], 3)

all_beta <- rbind(colMeans(all_parm$beta0), colMeans(all_parm$beta1), colMeans(all_parm$beta2))
rownames(all_beta) <- paste("beta", 0: 2, sep = "")

all_gamma <- rmvnorm(3, mu = rep(0, max(L)), sigma = 1/tau_gamma*diag(nrow = max(L)))
all_gamma[, 1] <- 0
all_gamma[, 1] <- rep(all_beta[1, 1], 3)
rownames(all_gamma) <- paste("gamma", 0: 2, sep = "")

## Sample for data
Y <- matrix(NA, N, TT)
firm_effect <- array(NA, dim = c(K, max(L), TT))
mu <- rnorm(N, m_bar, sqrt(e))
z <- rnorm(N, sqrt(1/tau_z))
w <- matrix(0, K, (TT + 1))

J <- matrix(0, N, (TT+1))
J_indicator <- array(0, dim = c(N, (TT + 1), K))
FF <- array(NA, dim = c(N, K, max(L), TT)) 
x <- array(NA, dim = c(K, max(L), (TT + 1)))

# Set some columns as NA since not all sector have the same number of firms
for (k in 1: K) {
  FF[, k, , ][, 1: L[k], ] <- 0
  x[k, , ][1: L[k], ] <- 0
}

# case <- "simple"
case <- "latent"
for (t in 1: TT) {
  w[-1, (t + 1)] <- rmvnorm(1, mu = w[-1, t], sigma = 1/(as.numeric(t == 1)*tau_w0 + as.numeric(t != 1)*tau_w)*diag(1, nrow = (K-1)))
  dist_mat1 <- sqrt(matrix(z^2, nrow = N, ncol = K) + z %*% t(w[, (t+1)]) + matrix(w[, (t+1)]^2, nrow = N, ncol = K, byrow = T))
  beta0 <- all_beta[1, ]
  beta1 <- all_beta[2, ]
  beta2 <- all_beta[3, ]
  all_q <- exp(matrix(beta0, N, K, byrow = T) + matrix(beta1, N, K, byrow = T)*as.numeric(J[, t] %in% 1: K)*J_indicator[, t, ] - 
                 matrix(beta2, N, K, byrow = T)*dist_mat1)
  # P[, t, ] <- all_q/rowSums(all_q)
  J[, (t + 1)] <- apply(all_q/rowSums(all_q), 1, rcat, n = 1)
  J_indicator[, (t + 1), ][cbind(1: N, J[, (t + 1)])] <- 1
  
  if(case == "simple"){
    current_firm_idx <- sapply(J[, (t + 1)], function(k){sample(1: L[k], 1)})
    full_idx <- cbind(1: N, J[, (t + 1)], current_firm_idx) # worker, sector, firm
    FF[, , , t][full_idx] <- 1
  }
  if(case == "latent"){
    temp_x <- x[, -1, (t + 1)] # prepare for latent space at time t, the first firm of each sector has space position zero
    temp_x_prev <- x[, -1, t] # retrieve latent space at time t-1
    # sample for all firm latent space at time t, which are all independent
    temp_x[!is.na(temp_x)] <- rmvnorm(1, mu = temp_x_prev[!is.na(temp_x_prev)], 
                                      sigma = 1/(as.numeric(t == 1)*tau_x0 + as.numeric(t != 1)*tau_x)*diag(1, nrow = sum(!is.na(temp_x))))
    x[, -1, (t + 1)] <- temp_x
    dist_mat2 <- abs(x[, , (t + 1)] - w[, (t + 1)]) # compute distance between latent space of sector and firm
    # N*L_max indicator matrix and the product between t=0 and if a worker stays in the same sector
    firm_indicator <- matrix(as.numeric(J[, t] %in% 1: K)*as.numeric(J[, (t + 1)] == J[, t]), N, max(L))
    prev_full_idx <- which(FF != 0, arr.ind = T)
    prev_full_idx <- prev_full_idx[order(prev_full_idx[, 1]), ]
    firm_idx <- prev_full_idx[prev_full_idx[, 4] == t, 3] # retrieve the firm of each worker at time t-1
    firm_indicator[cbind(1: N, firm_idx)] <- 1
    gamma0 <- all_gamma[1, ]
    gamma1 <- all_gamma[2, ]
    gamma2 <- all_gamma[3, ]
    all_q_firm <- exp(matrix(gamma0, N, max(L), byrow = T) + matrix(gamma1, N, max(L), byrow = T)*firm_indicator - 
                        matrix(gamma2, N, max(L), byrow = T)*dist_mat2[J[, (t + 1)], ])
    full_index <- cbind(1: N, J[, (t + 1)], apply(all_q_firm/rowSums(all_q_firm), 1, rcat, n = 1))
    FF[, , , t][full_index] <- 1
  }
  for (k in 1: K) {
    firm_effect[k, 1: L[k] , t] <- rnorm(L[k], nu[k], sd = sqrt(1/tau_k[k]))
  }
  worker_mean <- mu + firm_effect[cbind(full_index[, 2: 3], t)]
  Y[, t] <- sapply(1: N,function(i){rnorm(1, mean = worker_mean[i], sd = sqrt(1/tau))})
}

full_info <- which(FF != 0, arr.ind = T)
full_info <- full_info[order(full_info[, 1]), ]

job <- array(data = NA, dim = c(N,TT))
for (i in 1:N) {
  for (t in 1:TT) {
    job[i,t] <- which(deflate_y_dat$a[i,t,]==1)
  }
}
table(job)

par(mfrow = c(3, 2))
hist(log(deflate_y_dat$y[, 1]), main = "IFLS Log Wage at t=1", xlab = "wage")
hist(Y[, 1], main = "Simulated Log Wage at t=1", xlab = "wage")
hist(log(deflate_y_dat$y[, 2]), main = "IFLS log Wage at t=2", xlab = "wage")
hist(Y[, 2], main = "Simulated Log Wage at t=2", xlab = "wage")
hist(log(deflate_y_dat$y[, 3]), main = "IFLS Log Wage at t=3", xlab = "wage")
hist(Y[, 3], main = "Simulated Log Wage at t=3", xlab = "wage")


par(mfrow = c(3, 2))
sector_count_IFLS1 = table(job[, 1])
sector_count_IFLS2 = table(job[, 2])
sector_count_IFLS3 = table(job[, 3])
sector_count_sim1 = table(full_info[full_info[, 4] == 1, 2])
sector_count_sim2 = table(full_info[full_info[, 4] == 2, 2])
sector_count_sim3 = table(full_info[full_info[, 4] == 3, 2])
p1 <- barplot(sector_count_IFLS1, main = "IFLS #workers in Sectors at t=1", xlab = "sector")
text(p1, sector_count_IFLS1*0.8, as.character(sector_count_IFLS1), cex = 1)
p2 <- barplot(sector_count_sim1, main = "Simulated #workers in Sectors at t=1", xlab = "sector")
text(p2, sector_count_sim1*0.8, as.character(sector_count_sim1), cex = 1)

p3 <- barplot(sector_count_IFLS2, main = "IFLS #workers in Sectors at t=2", xlab = "sector")
text(p3, sector_count_IFLS2*0.8, as.character(sector_count_IFLS2), cex = 1)
p4 <- barplot(sector_count_sim2, main = "Simulated #workers in Sectors at t=2", xlab = "sector")
text(p4, sector_count_sim2*0.8, as.character(sector_count_sim2), cex = 1)

p5 <- barplot(sector_count_IFLS3, main = "IFLS #workers in Sectors at t=3", xlab = "sector")
text(p5, sector_count_IFLS3*0.8, as.character(sector_count_IFLS3), cex = 1)
p6 <- barplot(sector_count_sim3, main = "Simulated #workers in Sectors at t=3", xlab = "sector")
text(p6, sector_count_sim3*0.8, as.character(sector_count_sim3), cex = 1)

par(mfrow = c(2, 2))
IFLS_diff_12 <- sum(deflate_y_dat$y[, 2] >= deflate_y_dat$y[, 1])
sim_diff_12 <- sum(Y[, 2] >= Y[, 1])
compare_1 <- c(IFLS_diff_12, sim_diff_12)
names(compare_1) <- c("IFLS", "Simulated")
IFLS_diff_23 <- sum(deflate_y_dat$y[, 3] >= deflate_y_dat$y[, 2])
sim_diff_23 <- sum(Y[, 3] >= Y[, 2])
compare_2 <- c(IFLS_diff_23, sim_diff_23)
names(compare_2) <- c("IFLS", "Simulated")
p7 <- barplot(compare_1, main = "#Wage-increase from 1 to 2")
text(p7, compare_1*0.8, as.character(compare_1), cex = 1)
p8 <- barplot(compare_2, main = "#wage-increase from 2 to 3")
text(p8, compare_2*0.8, as.character(compare_2), cex = 1)

IFLS_worker_stay_12 <- sum(job[, 2] == job[, 1])
IFLS_worker_stay_23 <- sum(job[, 3] == job[, 2])
sim_worker_stay_12 <- sum(J[, 3] == J[, 2])
sim_worker_stay_23 <- sum(J[, 4] == J[, 3])

compare_3 <- c(IFLS_worker_stay_12, sim_worker_stay_12)
names(compare_3) <- c("IFLS", "Simulated")
compare_4 <- c(IFLS_worker_stay_23, sim_worker_stay_23)
names(compare_4) <- c("IFLS", "Simulated")
p9 <- barplot(compare_3, main = "#Worker Stay in Same Sector from 1 to 2")
text(p9, compare_3*0.8, as.character(compare_3), cex = 1)
p10 <- barplot(compare_4, main = "#Worker Stay in Same Sector from 2 to 3")
text(p10, compare_4*0.8, as.character(compare_4), cex = 1)


