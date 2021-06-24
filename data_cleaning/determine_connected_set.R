library(tidyverse)
load("data/deflate_y_dat.rda")

dat <- deflate_y_dat
K <- 7
Ti <- 3
N <- 1973
labs <- c("Agriculture","Manufacturing","Construction","Hospitality",
          "Transportation","Business","Social Services")

job <- array(data = NA, dim = c(N,Ti))
for (i in 1:N) {
  for (t in 1:Ti) {
    job[i,t] <- which(dat$a[i,t,]==1)
  }
}
res <- matrix(data = 0, nrow = K, ncol = K)
for (i in 1:N) {
  uni <- unique(job[i,])
  if (length(uni) == 2) {
    res[uni[1],uni[2]] <- res[uni[1],uni[2]] + 1
  }
  if (length(uni) == 3) {
    res[uni[1],uni[2]] <- res[uni[1],uni[2]] + 1
    res[uni[1],uni[3]] <- res[uni[1],uni[3]] + 1
    res[uni[2],uni[3]] <- res[uni[2],uni[3]] + 1
  }
}
for (i in 1:K) {
  for (j in 1:K) {
    if (i != j) {
      res[i,j] <- res[i,j] + res[j,i]
      res[j,i] <- 0
    }
  }
}
res
colnames(res) <- labs
rownames(res) <- labs
