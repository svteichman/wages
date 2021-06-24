load("data/expanded_data1.rda")

# first check how many observations in data include sector 6 
N <- expanded_data1$N
Ti <- expanded_data1$T
K <- expanded_data1$K
sector <- matrix(NA, nrow = N, ncol = Ti)
for (n in 1:N) {
  for (t in 1:Ti) {
    for (k in 1:K) {
      if (expanded_data1$a[n,t,k] == 1) {
        sector[n,t] <- k
      }
    }
  }
}
sector6 <- which(sector[,1] == 6 | sector[,2] == 6 | sector[,3] == 6)
length(sector6) # 150 people at some point work in sector 6
length(sector6)/N # this is ~8% of the overall dataset 

# update first parts of the dataset 
names(expanded_data1)
new_N <- N - length(sector6)
new_K <- K - 1
new_y <- expanded_data1$y[-sector6,]
new_a <- expanded_data1$a[-sector6,,-6]
new_r <- expanded_data1$r[,-sector6]
new_s <- expanded_data1$s[,-sector6]
new_mean_log_wage <- mean(log(new_y))
new_sector <- sector[-sector6, ]
new_sector[new_sector == 7] <- 6

# create an r_strict that only has a 1 if the sector before and after match 
new_r_strict <- new_r
for (n in 1:new_N) {
  for (t in 1:(Ti-1)) {
    if (new_sector[n, t] != new_sector[n, t+1]) {
      new_r_strict[t, n] <- 0 
    }
  }
}

# use r_strict to figure out how many distinct firm effects should be drawn
new_M <- sum(1-new_r_strict)
new_M1 <- sum(1-new_r_strict[1,])
new_M2 <- sum(1-new_r_strict[2,])
new_job <- vector(length = new_N+new_M)
new_job[1:new_N] <- new_sector[,1]
new_job[(new_N+1):(new_N + new_M1)] <- new_sector[which(new_r_strict[1,] == 0),2]
new_job[(new_N+1+new_M1):(new_N + new_M)] <- new_sector[which(new_r_strict[2,] == 0), 3]

# make N x T matrix which gives coordinate of unique firm effects to use
# there will be N + M unique firm effects
new_firm_effect_ind <- matrix(data = NA, nrow = new_N, ncol = Ti)
new_firm_effect_ind[,1] <- 1:new_N
count <- new_N + 1
for (t in 2:Ti) {
  for (n in 1:new_N) {
    if (new_r_strict[t-1,n] == 1) {
      new_firm_effect_ind[n,t] <- new_firm_effect_ind[n,t-1]
    } else {
      new_firm_effect_ind[n,t] <- count
      count <- count + 1
    }
  }
}

# get empirical variances for people who remain in sectors but switch firms 
diff_mat <- matrix(data = NA, nrow = new_N*Ti, ncol = 2)
count <- 1 
for (n in 1:new_N) {
  for (t in 1:(Ti - 1)) {
    if (new_sector[n,t] == new_sector[n, t+1] & new_r[t,n] == 0) {
      diff <- log(new_y[n,t]) - log(new_y[n,t+1])
      diff_mat[count,] <- c(diff, new_sector[n,t])
      count <- count + 1
    }
  }
}
diff_mat <- diff_mat[!is.na(diff_mat[,1]),]
diff_var <- var(diff_mat[,1])
new_moving_precision <- 1/diff_var

data_no_s6 <- list(new_N,
                   Ti,
                   new_K,
                   new_M, 
                   new_y,
                   new_a,
                   new_sector,
                   new_r,
                   new_s,
                   new_mean_log_wage,
                   new_job,
                   new_firm_effect_ind,
                   new_moving_precision)
names(data_no_s6) <- c("N","T","K","M","y","a","sector","r","s","mean_log_wage",
                          "new_job","firm_effect_ind","moving_precision")
save(data_no_s6, file = "data/data_no_s6.rda")

