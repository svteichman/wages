library(tidyverse)

load(file = "data/deflate_y_dat.rda")
deflate_y_dat <- list(deflate_y_dat$N,
                      deflate_y_dat$T,
                      deflate_y_dat$K,
                      deflate_y_dat$y,
                      deflate_y_dat$a,
                      deflate_y_dat$r,
                      deflate_y_dat$s,
                      mean(log(deflate_y_dat$y)))
names(deflate_y_dat) <- c("N","T","K","y","a","r","s","mean_log_wage")

sector <- matrix(NA, nrow = deflate_y_dat$N, ncol = deflate_y_dat$T)
for (n in 1:deflate_y_dat$N) {
  for (t in 1:deflate_y_dat$T) {
    for (k in 1:deflate_y_dat$K) {
      if (deflate_y_dat$a[n,t,k] == 1) {
        sector[n,t] <- k
      }
    }
  }
}

# start by checking each case where r = 1 and which sectors the job was in before and after
transition_mat_r <- matrix(data = 0, nrow = 7, ncol = 7) 
for (n in 1:deflate_y_dat$N) {
  for (t in 1:(deflate_y_dat$T-1)) {
    if (deflate_y_dat$r[t,n] == 1) {
      sec1 <- sector[n, t]
      sec2 <- sector[n, t+1]
      transition_mat_r[sec1, sec2] <- transition_mat_r[sec1, sec2] + 1
    }
  }
}
diag_prop <- sum(diag(transition_mat_r))/sum(deflate_y_dat$r)
off_diag_prop <- (sum(transition_mat_r) - diag_prop)/sum(deflate_y_dat$r)

# also look for transition mat for s 
transition_mat_s <- matrix(data = 0, nrow = 7, ncol = 7) 
for (n in 1:deflate_y_dat$N) {
  for (t in 1:(deflate_y_dat$T-1)) {
    if (deflate_y_dat$s[t,n] == 1) {
      sec1 <- sector[n, t]
      sec2 <- sector[n, t+1]
      transition_mat_s[sec1, sec2] <- transition_mat_s[sec1, sec2] + 1
    }
  }
}

# create an r_strict that only has a 1 if the sector before and after match 
r_strict <- deflate_y_dat$r
for (n in 1:deflate_y_dat$N) {
  for (t in 1:(deflate_y_dat$T-1)) {
    if (sector[n, t] != sector[n, t+1]) {
      r_strict[t, n] <- 0 
    }
  }
}

# use r_strict to figure out how many distinct firm effects should be drawn
M <- sum(1-r_strict)
M1 <- sum(1-r_strict[1,])
M2 <- sum(1-r_strict[2,])
new_job <- vector(length = deflate_y_dat$N+M)
new_job[1:deflate_y_dat$N] <- sector[,1]
new_job[(deflate_y_dat$N+1):(deflate_y_dat$N + M1)] <- sector[which(r_strict[1,] == 0),2]
new_job[(deflate_y_dat$N+1+M1):(deflate_y_dat$N + M)] <- sector[which(r_strict[2,] == 0), 3]

# # make a two dimensional list of the worker and time of each new job 
# new_job_coord <- matrix(data = NA, nrow = length(new_job), ncol = 2)
# new_job_coord[1:deflate_y_dat$N,1] <- 1:deflate_y_dat$N
# new_job_coord[1:deflate_y_dat$N,2] <- 1
# count <- 1
# for (t in 1:(deflate_y_dat$T-1)) {
#   for (n in 1:deflate_y_dat$N) {
#     if (r_strict[t,n] == 0) {
#       new_job_coord[count + deflate_y_dat$N, ] <- c(n,t+1)
#       count <- count + 1
#     }
#   }
# }

# make N x T matrix which gives coordinate of unique firm effects to use
# there will be N + M unique firm effects
firm_effect_ind <- matrix(data = NA, nrow = deflate_y_dat$N, ncol = deflate_y_dat$T)
firm_effect_ind[,1] <- 1:deflate_y_dat$N
count <- deflate_y_dat$N + 1
for (t in 2:deflate_y_dat$T) {
  for (n in 1:deflate_y_dat$N) {
    if (r_strict[t-1,n] == 1) {
      firm_effect_ind[n,t] <- firm_effect_ind[n,t-1]
    } else {
      firm_effect_ind[n,t] <- count
      count <- count + 1
    }
  }
}

# get empirical variances for people who remain in sectors but switch firms 
diff_mat <- matrix(data = NA, nrow = deflate_y_dat$N*deflate_y_dat$T, ncol = 2)
count <- 1 
for (n in 1:deflate_y_dat$N) {
  for (t in 1:(deflate_y_dat$T - 1)) {
    if (sector[n,t] == sector[n, t+1] & deflate_y_dat$r[t,n] == 0) {
      diff <- log(deflate_y_dat$y[n,t]) - log(deflate_y_dat$y[n,t+1])
      diff_mat[count,] <- c(diff, sector[n,t])
      count <- count + 1
    }
  }
}
diff_mat <- diff_mat[!is.na(diff_mat[,1]),]
sector_variances <- as.data.frame(diff_mat) %>%
  group_by(V2) %>%
  summarise(sector_var = var(V1))
sector_precisions <- 1/sector_variances$sector_var

# add new information to data
expanded_data <- list(deflate_y_dat$N,
                      deflate_y_dat$T,
                      deflate_y_dat$K,
                      M, 
                      deflate_y_dat$y,
                      deflate_y_dat$a,
                      deflate_y_dat$r,
                      deflate_y_dat$s,
                      mean(log(deflate_y_dat$y)),
                      new_job,
                      firm_effect_ind,
                      sector_precisions)
names(expanded_data) <- c("N","T","K","M","y","a","r","s","mean_log_wage",
                          "new_job","firm_effect_ind","sector_precisions")
save(expanded_data, file = "data/expanded_data.rda")

diff_var <- var(diff_mat[,1])
diff_precision <- 1/diff_var
expanded_data1 <- expanded_data
expanded_data1$sector_precisions <- NULL
expanded_data1$moving_precision <- diff_precision
save(expanded_data1, file = "data/expanded_data1.rda")
