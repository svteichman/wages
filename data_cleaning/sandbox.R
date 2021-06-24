library(ggplot2)
load(file = "wave1to5.rda")

wave1to5_sector <- as.data.frame(wave1to5) %>% 
  filter(w4_yr_curr_job > 7 & w5_yr_curr_job > 14) %>%
  select(w3_sector,w4_sector,w5_sector,w4_yr_curr_job,w5_yr_curr_job) %>%
  mutate(w3_sector = as.numeric(levels(w3_sector))[w3_sector],
         w4_sector = as.numeric(levels(w4_sector))[w4_sector],
    same = ifelse(w3_sector == w4_sector & w3_sector == w5_sector, T, F))
sum(wave1to5_sector$same)/nrow(wave1to5_sector)

wave1to5df_numppl <- as.data.frame(wave1to5) %>% 
  filter(w4_yr_curr_job > 7 & w5_yr_curr_job > 14) %>%
  select(w3_num_ppl, w4_num_ppl, w5_num_ppl) 

wave1to5df_occ <- as.data.frame(wave1to5) %>% 
  filter(w2_same == "1. Yes" & w3_yr_curr_job > 3 & w4_yr_curr_job > 10 & w5_yr_curr_job > 17) %>%
  select(w1_occ, w2_occ, w3_occ, w4_occ, w5_occ) 

ggplot(wave1to5, aes(x = w1_occ)) + geom_bar()

table(wave3to5$w3_sector)
table(wave3to5$w4_sector)
table(wave3to5$w5_sector)


#looking at why mds won't work 
library(parallelDist)
library(labdsv)
library(ape)

load("data/dist1.rda")
load("data/dist2.rda")
load("data/dist3.rda")

mds1 <- cmdscale(dist1, k=2)
mds2 <- cmdscale(dist2, k=2)
mds3 <- cmdscale(dist3, k=2)
#mds1 <- pco(dist1, k = 2)$points
#mds2 <- pco(dist2, k = 2)$points
#mds3 <- pco(dist3, k = 2)$points
#mds1 <- pcoa(dist1)$vectors
#mds2 <- pcoa(dist2)$vectors
#mds3 <- pcoa(dist3)$vectors

plot_mds_z <- bind_rows(time1,time2,time3) %>%
  select(-c(`1`,`2`)) %>%
  mutate(mds_dim1 = c(mds1[,1],mds2[,1],mds3[,1]),
         mds_dim2 = c(mds1[,2],mds2[,2],mds3[,2])) %>%
  filter(type == "z")
plot_mds_w <- bind_rows(time1,time2,time3) %>%
  select(-c(`1`,`2`)) %>%
  mutate(mds_dim1 = c(mds1[,1],mds2[,1],mds3[,1]),
         mds_dim2 = c(mds1[,2],mds2[,2],mds3[,2])) %>%
  filter(type == "w")
ggplot() +  
  geom_point(data = plot_mds_z, aes(x = mds_dim1, y = mds_dim2, color = class, alpha = 0.5)) + 
  geom_point(data = plot_mds_w, aes(x = mds_dim1, y = mds_dim2, fill = class), 
             shape = 25, size = 4, color = "black")  +
  facet_wrap(~time)


# look closer at data 
wave5_TK <- wave5_TK %>% 
  mutate(wage_employment = ifelse(work_cat %in% c("1:Self-employed",
                                                  "2:Self employed with  unpaid family/temporary worker",
                                                  "3:Self-employed with employees/permanent workers",
                                                  "6:Unpaid family worker")
                                             ,0,1))
wave5 <- inner_join(wave5_K,wave5_TK, by = "id")
wave5_wage_emp <- wave5 %>% filter(wage_employment == 1)
wave5_wage <- wave5_wage_emp %>% filter(!is.na(wage_yr))
names(wave5_wage_emp) <- add_prefixes(names(wave5_wage_emp),5)
names(wave5_wage) <- add_prefixes(names(wave5_wage),5)

wave5 %>%
  group_by(sector) %>%
  summarise(count = n())

wave5_wage_emp %>% 
  group_by(w5_sector) %>% 
  summarise(count = n())

wave5_wage %>% 
  group_by(w5_sector) %>% 
  summarise(count = n())

# work through prior for variance (or persistence) for log wages
x <- seq(0, 7, by=.001)
plot(x, dgamma(x, 1/2, 1), type="l",
     ylim=c(0,2), ylab="Density",
     main="Gamma Densities: shape=.5;
 scale = 1")

plot(x, dgamma(x, 1/2, 1)^(-1), type="l",
     ylim=c(0,2), ylab="Density",
     main="Gamma Densities: shape=.5;
 scale = 1")

# check on firm changes within a sector 
load("data/deflate_y_dat.rda")
Ti <- 3
N <- 1973
job <- array(data = NA, dim = c(N,Ti))
for (i in 1:N) {
  for (t in 1:Ti) {
    job[i,t] <- which(deflate_y_dat$a[i,t,]==1)
  }
}
mat <- matrix(data = 0, nrow = 7, ncol = 2)
stay_t2 <- which(deflate_y_dat$s[1,] == 1)
stay_t3 <- which(deflate_y_dat$s[2,] == 1)
for (i in stay_t2) {
  sec <- job[i,1] 
  mat[sec,1] <- mat[sec,1] + 1
}
for (i in stay_t3) {
  sec <- job[i,2] 
  mat[sec,2] <- mat[sec,2] + 1
}
