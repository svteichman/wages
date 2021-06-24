library(parallelDist)

load("data/deflate_y_dat.rda")
load("data/dist_mat1.rda")
load("data/dist_mat2.rda")
load("data/dist_mat3.rda")

dat <- deflate_y_dat

mds1 <- cmdscale(dist_mat1, k=2)
mds2 <- cmdscale(dist_mat2, k=2)
mds3 <- cmdscale(dist_mat3, k=2)

#make data frame to plot

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

plot_mds_z <- data.frame(mds_dim1 = c(mds1[,1],mds2[,1],mds3[,1]),
                         mds_dim2 = c(mds1[,2],mds2[,2],mds3[,2]),
                         time = sort(rep(1:3,N+7)),
                         type = rep(c(rep("w",7),rep("z",N)),3)) %>%
  filter(type == "z")
plot_mds_z$class <- labs[c(job[,1],job[,2],job[,3])]
plot_mds_w <- data.frame(mds_dim1 = c(mds1[,1],mds2[,1],mds3[,1]),
                         mds_dim2 = c(mds1[,2],mds2[,2],mds3[,2]),
                         time = sort(rep(1:3,N+7)),
                         type = rep(c(rep("w",7),rep("z",N)),3)) %>%
  filter(type == "w")
plot_mds_w$class <- labs[(rep(1:7,3))]
ggplot() +
  geom_point(data = plot_mds_z, aes(x = mds_dim1, y = mds_dim2, color = class, alpha = 0.5)) +
  geom_point(data = plot_mds_w, aes(x = mds_dim1, y = mds_dim2, fill = class),
             shape = 25, size = 4, color = "black")  +
  facet_wrap(~time)

ggsave("figs/mds_latent_pos.png")
