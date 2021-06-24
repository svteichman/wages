library(tidyverse)

# transforming data using occupation as class ----
# note - we decided not to do this after discussion with Rachel 

load(file = "wave1to5.rda")
N <- nrow(wave1to5)
time <- 5
#Start by using occupation as class 
K <- 92
Xnum <- 2

sex <- grep("sex", names(wave1to5))
age <- grep("_age", names(wave1to5))
x <- array(data = NA, dim = c(N,time,2))
for (i in 1:N) {
  x[i, ,1] <- as.numeric(wave1to5[i,age])
  x[i, ,2] <- as.numeric(wave1to5[i,sex])
}

wage <- grep("wage_yr", names(wave1to5))
y <- array(data = NA, dim = c(N,time))
for (i in 1:N) {
  y[i,] <- as.numeric(wave1to5[i,wage])
}

occ <- grep("occ", names(wave1to5))
job <- array(data = NA, dim = c(N,time)) 
for (i in 1:N) {
  job[i,] <- as.numeric(wave1to5[i,occ])
}
a <- array(data = 0, dim = c(N,time,K))
for (i in 1:N) {
  for (ti in 1:time) {
    occup <- wave1to5[i,occ[ti]]
    a[i,ti,as.numeric(occup)] <- 1
  }
}

r <- array(data = 0, dim = c(time-1,N))
for (i in 1:N) {
  if (wave1to5$w2_same[i] == "1. Yes" & !is.na(wave1to5$w2_same[i])) {
    r[1,i] = 1
  }
  if (wave1to5$w3_yr_curr_job[i] > 3) {
    r[2,i] = 1
  }
  if (wave1to5$w4_yr_curr_job[i] > 7) {
    r[3,i] = 1
  } 
  if (wave1to5$w5_yr_curr_job[i] > 7) {
    r[4,i] = 1
  }
}

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

first_class_means <- array(0,dim = K)
for (k in 1:K) {
  gr <- which(job[,1]==k)
  first_class_means[k] <- mean(log(y[gr,1]))
}
fix <- which(is.na(first_class_means) | first_class_means == -Inf)
first_class_means[fix] <- mean(first_class_means[is.finite(first_class_means)],na.rm=T)
  
stan_dat <- list(
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

save(stan_dat, file = "stan_dat.rda")


# transforming data using sector as class ----

load(file = "wave3to5.rda")
#remove observations where years at current job is missing
#remove observations where sector is missing 
sec_na <- c(which(is.na(wave3to5$w5_yr_curr_job) | is.na(wave3to5$w4_yr_curr_job)),
            which(wave3to5$w3_sector == ""),
            which(wave3to5$w4_sector == ""),
            which(wave3to5$w5_sector == ""))
wave3to5_filt <- wave3to5[-sec_na,]
#remove observations where yearly wage is 0
wage_0 <- which(wave3to5_filt$w3_wage_yr == 0 | wave3to5_filt$w4_wage_yr == 0 | wave3to5_filt$w5_wage_yr == 0)
wave3to5_filt <- wave3to5_filt[-wage_0,]
#recode sectors so the values are the same across years 
wave3to5_filt <- wave3to5_filt %>%
  mutate(w4_sector = as.numeric(levels(wave3to5_filt$w4_sector))[w4_sector],
         w3_sector = as.character(w3_sector),
         w3_sector = recode(w3_sector, `X` = "10"))
#recode sex in wave 3 and 5 to match wave 4, to avoid differences in coding
wave3to5_filt$w3_sex01 <- wave3to5_filt$w4_sex01
wave3to5_filt$w5_sex01 <- wave3to5_filt$w4_sex01
#recode age in wave 5 to be age in wave 4 plus the 7 years between waves
wave3to5_filt$w5_age <- wave3to5_filt$w4_age + 7
save(wave3to5_filt, file = "wave3to5_filt.rda")

N <- nrow(wave3to5)
time <- 3
#Start by using sector as class 
K <- 11
Xnum <- 2

#create X array
sex <- grep("sex", names(wave3to5))
age <- grep("_age", names(wave3to5))
x <- array(data = NA, dim = c(N,time,2))
for (i in 1:N) {
  x[i, ,1] <- as.numeric(wave3to5[i,age])
  x[i, ,2] <- as.numeric(wave3to5[i,sex])
}

#create Y array
wage <- grep("wage_yr", names(wave3to5))
y <- array(data = NA, dim = c(N,time))
for (i in 1:N) {
  y[i,] <- as.numeric(wave3to5[i,wage])
}

#convert sectors to indicator of job category 
sec <- grep("sec", names(wave3to5))
job <- array(data = NA, dim = c(N,time)) 
for (i in 1:N) {
  job[i,] <- as.numeric(wave3to5[i,sec])
}
a <- array(data = 0, dim = c(N,time,K))
for (i in 1:N) {
  for (ti in 1:time) {
    sector <- wave3to5[i,sec[ti]]
    #sectors are coded 0 to 10, store in values 1 to 11 
    a[i,ti,as.numeric(sector)+1] <- 1
  }
}

#make remaining at the same job array
r <- array(data = 0, dim = c(time-1,N))
for (i in 1:N) {
  if (wave3to5$w4_yr_curr_job[i] > 7) {
    r[1,i] = 1
  } 
  if (wave3to5$w5_yr_curr_job[i] > 7) {
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

stan_dat <- list(
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

save(stan_dat, file = "stan_dat.rda")


