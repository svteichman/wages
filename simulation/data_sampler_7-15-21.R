library(Rfast)
library(extraDistr)
library(rstan)
library(locfit)
library(gdata)

set.seed(8888)
#### Load the data and pre-trained model ####
load(file = "data/deflate_y_dat.rda")
# load('fit_ifls_7-13-21.rda')
load(file = 'parm_fit_ifls_7-13-21.rda')

N <- deflate_y_dat$N 
TT <- deflate_y_dat$T
K <- deflate_y_dat$K

#### Extract all parameter estiamtes from latent space model ####

# all_parm <- matrix(NA, nrow = fit_ifls@sim$chains, ncol = length(fit_ifls@sim$samples[[1]]))
# for (i in 1: nrow(all_parm)) {
#   all_parm[i, ] <- sapply(fit_ifls@sim$samples[[1]], mean)
# }
# colnames(all_parm) <- names(fit_ifls@sim$samples[[1]])
# all_parm <- colMeans(all_parm)
# save(all_parm, file = "parm_fit_ifls_7-13-21.rda")

all_parm_names <- names(all_parm)
# Extract worker latent space
z <- matrix(all_parm[grep("z", all_parm_names)], ncol = 2)

# Extract sector latent space
w <- array(all_parm[setdiff(setdiff(grep("w", all_parm_names), grep("free_w", all_parm_names)), grep("tau", all_parm_names))], 
           dim = c(K, TT, 2))

# Extract worker effect
mu <- all_parm[grep("mu", all_parm_names)]

# Extract beta 
all_beta <- matrix(all_parm[setdiff(grep("beta", all_parm_names), grep("free_beta", all_parm_names))], nrow = 3, byrow = T)
rownames(all_beta) <- paste("beta", 0: 2, sep = "")

# Extract edge probability between workers and sectors
p <- exp(array(all_parm[grep("p", all_parm_names, value = T)], dim = c(N, TT, K)))

# Extract sigma square for log wage distribution
sigmasq <- matrix(all_parm[grep("sigmasq", all_parm_names)], ncol = TT)
sigma <- matrix(all_parm[setdiff(grep("sigma", all_parm_names), grep("sigmasq", all_parm_names))], ncol = TT)
# sigma <- apply(log(deflate_y_dat$y), 2, sd)

# Extract mean of firm effect
nu <- all_parm[grep("nu", all_parm_names)]

# Extract hypeparameters
tau_w0_idx <- grep("tau_w0", all_parm_names)
tau_w0 <- all_parm[tau_w0_idx]
tau_w_idx <- setdiff(grep("tau_w", all_parm_names), grep("tau_w0", all_parm_names))
tau_w <- all_parm[tau_w_idx]
tau_k_idx <- grep("tau_k", all_parm_names)
tau_k <- all_parm[tau_k_idx]
tau_m_idx <- grep("tau_m", all_parm_names)
tau_m <- all_parm[tau_m_idx]
tau_idx <- setdiff(grep("tau", all_parm_names), c(tau_w0_idx, tau_w_idx, tau_k_idx, tau_m_idx))
tau <- all_parm[tau_idx]

# Sample for gamma (coef of edge between sector and firm)

# number of firms of each sector
L <- rep(30, K) # For simplicity, 10 firms for each sector
# L <- sample(5: 10, K, replace = TRUE)

m_bar <- mean(log(deflate_y_dat$y))
# tau_gamma <- rgamma(1, 0.5, 1)
all_gamma <- matrix(NA, nrow = 3, ncol = max(L))
# all_gamma <- Rfast::rmvnorm(3, mu = rep(0, max(L)), sigma = 1/tau_gamma*diag(nrow = max(L)))
all_gamma[, 1] <- 0
all_gamma[1, -1] <- sample(all_beta[1, -1], max(L)-1, replace = T)
all_gamma[2, -1] <- rep(all_beta[2, 2], max(L)-1) 
all_gamma[3, -1] <- rep(all_beta[3, 2], max(L)-1)
rownames(all_gamma) <- paste("gamma", 0: 2, sep = "")

# Sample for Hyperparameter of firm latent space
# tau_x0 <- rgamma(1, 0.5, 1)
# tau_x <- rgamma(1, 0.5, 1)
tau_x0 <- tau_w0
tau_x <- tau_w
#### Sample for data ####

Y <- matrix(NA, N, TT)
log_Y <- matrix(NA, N, TT)
firm_effect <- array(NA, dim = c(N, TT))
# mu <- rnorm(N, m_bar, sqrt(e))
# z <- rnorm(N, sqrt(1/tau_z))
# w <- matrix(0, K, (TT + 1))

J <- matrix(0, N, (TT+1))
FF <- matrix(0, N, TT)
J_indicator <- array(0, dim = c(N, (TT + 1), K))
F_indicator <- array(NA, dim = c(N, TT, max(L))) 
# sector_firm_idx <- c()
# for (k in 1: K) {
#   sector_firm_idx <- c(sector_firm_idx, paste(k, 1: L[k], sep = "_"))
# }
# colnames(F_indicator) <- sector_firm_idx
x <- array(NA, dim = c(K, max(L), (TT + 1), 2))

# Set some columns as NA since not all sector have the same number of firms
for (k in 1: K) {
  x[k, , , ][1: L[k], , ] <- 0
}
correction_factor <- c(2, 2.3, 2.5)
sigma_corrected <- eachrow(sigma, correction_factor, "/")
# case <- "simple"
case <- "latent"
for (t in 1: TT) {
  # t = 1
  dist_mat1 <- sqrt(matrix(rowsums(z^2), nrow = N, ncol = K) + z %*% t(w[, t, ]) + matrix(rowsums(w[, t, ]^2), nrow = N, ncol = K, byrow = T))
  beta0 <- all_beta[1, ]
  beta1 <- all_beta[2, ]
  beta2 <- all_beta[3, ]
  all_q <- exp(matrix(beta0, N, K, byrow = T) + matrix(beta1, N, K, byrow = T)*J_indicator[, t, ] - matrix(beta2, N, K, byrow = T)*dist_mat1)
  sector_of_worker <- apply(all_q/rowsums(all_q), 1, function(prob){sample(1: K, 1, prob = prob)})
  # sector_of_worker <- apply(p[, t, ]/rowsums(p[, t, ]), 1, function(prob){sample(1: K, 1, prob = prob)})
  J[, (t + 1)] <- sector_of_worker
  J_indicator[, (t + 1), ][cbind(1: N, sector_of_worker)] <- 1
  
  # if(case == "simple"){
  #   current_firm_idx <- sapply(J[, (t + 1)], function(k){sample(1: L[k], 1)})
  #   full_idx <- cbind(1: N, J[, (t + 1)], current_firm_idx) # worker, sector, firm
  #   FF[, , , t][full_idx] <- 1
  # }
  if(case == "latent"){
    temp_x <- x[, -1, (t + 1), ] # prepare for latent space at time t, the first firm of each sector has space position zero
    temp_x_prev <- x[, -1, t, ] # retrieve latent space at time t-1
    # sample for all firm latent space at time t, which are all independent
    temp_x[!is.na(temp_x)] <- Rfast::rmvnorm(1, mu = temp_x_prev[!is.na(temp_x_prev)], 
                                      sigma = 1/(as.numeric(t == 1)*tau_x0 + as.numeric(t != 1)*tau_x)*diag(sum(!is.na(temp_x))))
    x[, -1, (t + 1), ] <- temp_x
    dist_mat2 <- matrix(NA, nrow = K, ncol = max(L)) # compute distance between latent space of sector and firm
    for (i in 1: K) {
      dist_mat2[i, 1: L[i]] <- sqrt(rowsums(x[i, , (t + 1), ] - matrix(w[i, t, ], L[i], 2, byrow = T))^2)
    }
    
    # N*L_max indicator matrix and the product between t=0 and if a worker stays in the same sector
    
    firm_indicator <- matrix(as.numeric(J[, (t + 1)] == J[, t]), N, max(L))
    
    gamma0 <- all_gamma[1, ]
    gamma1 <- all_gamma[2, ]
    gamma2 <- all_gamma[3, ]
    all_q_firm <- exp(matrix(gamma0, N, max(L), byrow = T) + matrix(gamma1, N, max(L), byrow = T)*firm_indicator - 
                        matrix(gamma2, N, max(L), byrow = T)*dist_mat2[sector_of_worker, ])
    firm_of_worker <- apply(all_q_firm/rowSums(all_q_firm), 1, function(prob){sample(which(!is.na(prob)), 1, prob = prob)})
    FF[, t] <- firm_of_worker
    firm_indicator[cbind(1: N, firm_of_worker)] <- 1
    F_indicator[, t, ] <- firm_indicator
  }
  firm_effect[, t] <- sapply(1: N, function(i){rnorm(1, nu[J[i, (t + 1)]], sqrt(1/tau_k[J[i, (t + 1)]]))})
  if (t > 1){
    firm_effect[J[, (t + 1)] == J[, t], t] <- firm_effect[J[, (t + 1)] == J[, t], (t - 1)]
  }
  log_wage_IFLS <- log(deflate_y_dat$y[, t])
  fit_lm <- lm(log_wage_IFLS ~ mu + firm_effect[, t] - 1)
  log_wage_mean <- fitted(fit_lm)
  par(mfrow = c(2, 2))
  plot(fit_lm, 1: 4)
  # log_wage_mean <- mu + firm_effect[, t]
  Y[, t] <- sapply(1: N, function(i){rlnorm(1, mean = log_wage_mean[i], sd = sigma_corrected[i, t])})
  log_Y[, t] <- log(Y[, t])
}
data_sim <- list(N = N, time.num = TT, K = K, firm.num = L, y = Y, sector = J, firm = FF)
# save(data_sim, file = "data_sim7-16-2021") 

# Sector changes of IFLS
job <- array(data = NA, dim = c(N,TT))
for (i in 1:N) {
  for (t in 1:TT) {
    job[i,t] <- which(deflate_y_dat$a[i,t,]==1)
  }
}

# par(mfrow = c(3, 2))
# hist(log(deflate_y_dat$y[, 1]), breaks = 20, main = "IFLS Log Wage at t=1", xlab = "wage")
# hist(log_Y[, 1], breaks = 20, main = "Simulated Log Wage at t=1", xlab = "wage")
# hist(log(deflate_y_dat$y[, 2]), breaks = 20, main = "IFLS log Wage at t=2", xlab = "wage")
# hist(log_Y[, 2], breaks = 20, main = "Simulated Log Wage at t=2", xlab = "wage")
# hist(log(deflate_y_dat$y[, 3]), breaks = 20, main = "IFLS Log Wage at t=3", xlab = "wage")
# hist(log_Y[, 3], breaks = 20, main = "Simulated Log Wage at t=3", xlab = "wage")

# par(mfrow = c(3, 2))
# hist_max <- 1.5e8
# hist(deflate_y_dat$y[, 1], breaks = seq(0, hist_max, by = 1e7), main = "IFLS Wage at t=1", xlab = "wage")
# hist(exp(log_Y[log_Y[, 1] <= log(hist_max), 1]), breaks = seq(0, hist_max, by = 1e7), main = "Simulated exp of Log Wage at t=1", xlab = "wage")
# hist(deflate_y_dat$y[, 2], breaks = seq(0, hist_max, by = 1e7), main = "IFLS Wage at t=2", xlab = "wage")
# hist(exp(log_Y[log_Y[, 2] <= log(hist_max), 2]), breaks = seq(0, hist_max, by = 1e7), main = "Simulated exp of Log Wage at t=2", xlab = "wage")
# hist(deflate_y_dat$y[, 3], breaks = seq(0, hist_max, by = 1e7), main = "IFLS Wage at t=3", xlab = "wage")
# hist(exp(log_Y[log_Y[, 3] <= log(hist_max), 3]), breaks = seq(0, hist_max, by = 1e7), main = "Simulated exp of Log Wage at t=3", xlab = "wage")

pdf(file = "simulation/wage_comparison_2.pdf", height = 8.14, width = 7.4)
  par(mfrow = c(3, 2))
  hist_max <- 1.5e8
  hist(deflate_y_dat$y[, 1], breaks = seq(0, hist_max, by = 1e7), main = "IFLS Wage at t=1", xlab = "wage")
  hist(Y[Y[, 1] <= hist_max, 1], breaks = seq(0, hist_max, by = 1e7), main = "Simulated Wage at t=1", xlab = "wage")
  hist(deflate_y_dat$y[, 2],  breaks = seq(0, hist_max, by = 1e7), main = "IFLS Wage at t=2", xlab = "wage")
  hist(Y[Y[, 2] <= hist_max, 2], breaks = seq(0, hist_max, by = 1e7), main = "Simulated Wage at t=2", xlab = "wage")
  hist(deflate_y_dat$y[, 3], breaks = seq(0, hist_max, by = 1e7), main = "IFLS Wage at t=3", xlab = "wage")
  hist(Y[Y[, 3] <= hist_max, 3], breaks = seq(0, hist_max, by = 1e7), main = "Simulated Wage at t=3", xlab = "wage")
dev.off()

pdf(file = "simulation/sector_distr_comparison_2.pdf", height = 8.14, width = 7.4)
  par(mfrow = c(3, 2))
  sector_count_IFLS1 = table(job[, 1])
  sector_count_IFLS2 = table(job[, 2])
  sector_count_IFLS3 = table(job[, 3])
  sector_count_sim1 = table(J[, 2])
  sector_count_sim2 = table(J[, 3])
  sector_count_sim3 = table(J[, 4])
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
dev.off()


pdf(file = "simulation/wage_stay_comparison_2.pdf", height = 8.14, width = 7.4)
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
  p9 <- barplot(compare_3, main = "#Worker in Same Sector from 1 to 2")
  text(p9, compare_3*0.8, as.character(compare_3), cex = 1)
  p10 <- barplot(compare_4, main = "#Worker in Same Sector from 2 to 3")
  text(p10, compare_4*0.8, as.character(compare_4), cex = 1)
dev.off()




