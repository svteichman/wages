library(rstan)
library(tidyverse)
library(splitstackshape)
library(stringi)
library(parallelDist)
#library(labdsv)
library(ape)

load("fit_ifls_3-14.rda")
load("deflate_y_dat.rda")

fit <- fit_ifls
dat <- deflate_y_dat
fit_summ <- as.data.frame(summary(fit)$summary)
post <- as.matrix(fit)

K <- 7
Ti <- 3
N <- 1973

job <- array(data = NA, dim = c(N,Ti))
for (i in 1:N) {
  for (t in 1:Ti) {
    job[i,t] <- which(dat$a[i,t,]==1)
  }
}

wcols_t1_d1 <-grep("w\\[.,1,1", names(fit))
wcols_t1_d2 <-grep("w\\[.,1,2", names(fit))
wcols_t2_d1 <-grep("w\\[.,2,1", names(fit))
wcols_t2_d2 <-grep("w\\[.,2,2", names(fit))
wcols_t3_d1 <-grep("w\\[.,3,1", names(fit))
wcols_t3_d2 <-grep("w\\[.,3,2", names(fit))
zcols_d1 <- grep("z\\[.*,1\\]", names(fit))
zcols_d2 <- grep("z\\[.*,2\\]", names(fit))
n_post <- dim(post)[1]
dist_mat1 <- parDist(x = matrix(c(post[1,wcols_t1_d1],post[1,zcols_d1],
                                  post[1,wcols_t1_d2],post[1,zcols_d2]),
                                ncol = 2), method = "euclidean") 
dist_mat2 <- parDist(x = matrix(c(post[1,wcols_t2_d1],post[1,zcols_d1],
                                  post[1,wcols_t2_d2],post[1,zcols_d2]),
                                ncol = 2), method = "euclidean") 
dist_mat3 <- parDist(x = matrix(c(post[1,wcols_t3_d1],post[1,zcols_d1],
                                  post[1,wcols_t3_d2],post[1,zcols_d2]),
                                ncol = 2), method = "euclidean") 
for (i in 2:n_post) {
  dist1 <- parDist(x = matrix(c(post[i,wcols_t1_d1],post[i,zcols_d1],
                                post[i,wcols_t1_d2],post[i,zcols_d2]),
                              ncol = 2), method = "euclidean") 
  dist_mat1 <- dist_mat1 + dist1 
  dist2 <- parDist(x = matrix(c(post[i,wcols_t2_d1],post[i,zcols_d1],
                                post[i,wcols_t2_d2],post[i,zcols_d2]),
                              ncol = 2), method = "euclidean") 
  dist_mat2 <- dist_mat2 + dist2
  dist3 <- parDist(x = matrix(c(post[i,wcols_t3_d1],post[i,zcols_d1],
                                post[i,wcols_t3_d2],post[i,zcols_d2]),
                              ncol = 2), method = "euclidean") 
  dist_mat3 <- dist_mat3 + dist3
}
dist_mat1 <- dist_mat1/n_post
dist_mat2 <- dist_mat2/n_post
dist_mat3 <- dist_mat3/n_post

# mds1 <- cmdscale(dist_mat1, k=2)
# mds2 <- cmdscale(dist_mat2, k=2)
# mds3 <- cmdscale(dist_mat3, k=2)
# #make data frame to plot
# plot_mds_z <- data.frame(mds_dim1 = c(mds1[,1],mds2[,1],mds3[,1]),
#                          mds_dim2 = c(mds1[,2],mds2[,2],mds3[,2]),
#                          time = sort(rep(1:3,N+7)),
#                          type = rep(c(rep("w",7),rep("z",N)),3)) %>%
#   filter(type == "z")
# plot_mds_z$class <- as.factor(c(job[,1],job[,2],job[,3]))
# plot_mds_w <- data.frame(mds_dim1 = c(mds1[,1],mds2[,1],mds3[,1]),
#                          mds_dim2 = c(mds1[,2],mds2[,2],mds3[,2]),
#                          time = sort(rep(1:3,N+7)),
#                          type = rep(c(rep("w",7),rep("z",N)),3)) %>%
#   filter(type == "w")
# plot_mds_w$class <- as.factor(rep(1:7,3))
# ggplot() +  
#   geom_point(data = plot_mds_z, aes(x = mds_dim1, y = mds_dim2, color = class, alpha = 0.5)) + 
#   geom_point(data = plot_mds_w, aes(x = mds_dim1, y = mds_dim2, fill = class), 
#              shape = 25, size = 4, color = "black")  +
#   facet_wrap(~time)

save(dist_mat1, file = "dist_mat1.rda")
save(dist_mat2, file = "dist_mat2.rda")
save(dist_mat3, file = "dist_mat3.rda")
