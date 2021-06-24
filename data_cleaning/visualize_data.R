library(tidyverse)
library(xtable)

load(file = "wave3to5_filt.rda")
load(file = "stan_dat.rda")

# visualize r and s by sector ----
# make jobs array
N <- stan_dat$N
K <- stan_dat$K
time <- stan_dat$T
sec <- grep("sec", names(wave3to5_filt))
job <- array(data = NA, dim = c(time,N)) 
job <- data.frame(matrix(nrow = N, ncol = time))
for (i in 1:N) {
  job[i,] <- as.numeric(wave3to5_filt[i,sec])
}

# make array to hold r and s percentages by sector and time
# auxiliary function to add values with 0 people to jobs list 
fill_in <- function(jobs,vec_ind) {
  all <- sort(unique(unlist(job)))
  uninc <- all[which(!(all %in% vec_ind), arr.ind = T)]
  return(uninc)
}

# make dataframe 
rs_dat <- data.frame(matrix(ncol = 5, nrow = 2*K))
names(rs_dat) <- c("sector","time","count","per_r","per_s")

# set sector and time 
rs_dat$sector <- rep(0:(K-1),2)
rs_dat$time <- c(rep("2",K),rep("3",K))

# find count per sector for each time
job1 <- job %>% group_by(X1) %>% summarise(count = n())
job2 <- job %>% group_by(X2) %>% summarise(count = n())
uninc <- fill_in(job, job2$X2)
for (i in 1:length(uninc)) {
  job2 <- add_row(job2, X2= uninc[i], count = 0)
}
job2 <- job2 %>% arrange(X2)
job3 <- job %>% group_by(X3) %>% summarise(count = n())
uninc <- fill_in(job, job3$X3)
for (i in 1:length(uninc)) {
  job3 <- add_row(job3, X3 = uninc[i], count = 0)
}
job3 <- job3 %>% arrange(X3)
count_dat <- data.frame(sector = rep(0:(K-1),3),
                        time = c(rep("1",K),rep("2",K),rep("3",K)),
                        count = c(job1$count, job2$count, job3$count))
rs_dat$count <- c(job2$count, job3$count)

# find r percentage by sector 
r <- stan_dat$r
s <- stan_dat$s
r1 <- data.frame(sec = job[,2], r = r[1,]) %>%
  group_by(sec) %>%
  summarise(count = n(),
            r_count = sum(r)) %>%
  mutate(prop = r_count/count)
uninc <- fill_in(job, r1$sec)
for (i in 1:length(uninc)) {
  r1 <- add_row(r1, sec= uninc[i], count = 0, r_count = 0, prop = 0)
}
r1 <- r1 %>% arrange(sec)
r2 <- data.frame(sec = job[,3], r = r[2,]) %>%
  group_by(sec) %>%
  summarise(count = n(),
            r_count = sum(r)) %>%
  mutate(prop = r_count/count)
uninc <- fill_in(job, r2$sec)
for (i in 1:length(uninc)) {
  r2 <- add_row(r2, sec= uninc[i], count = 0, r_count = 0, prop = 0)
}
r2 <- r2 %>% arrange(sec)
rs_dat$per_r <- c(r1$prop,r2$prop)

# find s percentage by sector 
s1 <- data.frame(sec = job[,2], s = s[1,]) %>%
  group_by(sec) %>%
  summarise(count = n(),
            s_count = sum(s)) %>%
  mutate(prop = s_count/count)
uninc <- fill_in(job, s1$sec)
for (i in 1:length(uninc)) {
  s1 <- add_row(s1, sec= uninc[i], count = 0, s_count = 0, prop = 0)
}
s1 <- s1 %>% arrange(sec)
s2 <- data.frame(sec = job[,3], s = s[2,]) %>%
  group_by(sec) %>%
  summarise(count = n(),
            s_count = sum(s)) %>%
  mutate(prop = s_count/count)
uninc <- fill_in(job, s2$sec)
for (i in 1:length(uninc)) {
  s2 <- add_row(s2, sec= uninc[i], count = 0, s_count = 0, prop = 0)
}
s2 <- s2 %>% arrange(sec)
rs_dat$per_s <- c(s1$prop,s2$prop)

# make plots
ggplot(count_dat, aes(x = as.factor(sector), y = count, group = time, fill = time)) + 
  geom_bar(position = "dodge", stat = "identity") + ggtitle("Worker count by sector") +
  xlab("Sector") + ylab("Count") + scale_fill_brewer(palette="Dark2")
ggsave("count_plot.png")
ggplot(rs_dat, aes(x = as.factor(sector), y = per_r, group = time, fill = time)) + 
  geom_bar(position = "dodge", stat = "identity") + ggtitle("Proportion remaining by sector") +
  xlab("Sector") + ylab("Proportion") + scale_fill_brewer(palette="Dark2")
ggsave("r_prop_plot.png")
ggplot(rs_dat, aes(x = as.factor(sector), y = per_s, group = time, fill = time)) + 
  geom_bar(position = "dodge", stat = "identity") + 
  ggtitle("Proportion switching within sector by sector") + xlab("Sector") + ylab("Proportion") +
  scale_fill_brewer(palette="Dark2")
ggsave("s_prop_plot.png")

# visualize mistakes in r by sector ----
r_dat <- data.frame(sector = rs_dat$sector[1:K], r1_count = rep(0,K), 
                    r1_mist = rep(0,K), r2_count = rep(0,K), r2_mist = rep(0,K))
r1_ind <- which(r[1,]==1)
r2_ind <- which(r[2,]==1)

# encode mistakes for the sector in which they end i.e. job 2 if job 1 neq job 2 
# mistakes for r1
mist_r1 <- list()
for (ind in r1_ind) {
  sec <- job$X2[ind]
  r_dat$r1_count[sec+1] <- r_dat$r1_count[sec+1] + 1
  if (!(job$X1[ind] == job$X2[ind])) {
    r_dat$r1_mist[sec+1] <- r_dat$r1_mist[sec+1] + 1
    mist_r1 <- append(mist_r1, ind)
  }
}

# mistakes for r2
mist_r2 <- list()
for (ind in r2_ind) {
  sec <- job$X3[ind]
  r_dat$r2_count[sec+1] <- r_dat$r2_count[sec+1] + 1
  if (!(job$X2[ind] == job$X3[ind])) {
    r_dat$r2_mist[sec+1] <- r_dat$r2_mist[sec+1] + 1
    mist_r2 <- append(mist_r2, ind)
  }
}

# add columns for percentage wrong
r_dat <- r_dat %>% 
  mutate(r1_prop_wrong = r1_mist/r1_count,
         r2_prop_wrong = r2_mist/r2_count)

# plot proportions wrong 
r_plot_dat <- data.frame(sector = c(r_dat$sector, r_dat$sector), 
                         count_wrong = c(r_dat$r1_mist, r_dat$r2_mist),
                         prop_wrong = c(r_dat$r1_prop_wrong, r_dat$r2_prop_wrong),
                         time = c(rep("2",K),rep("3",K)))
ggplot(r_plot_dat, aes(x = as.factor(sector), y = count_wrong, group = time, fill = time)) + 
  geom_bar(position = "dodge", stat = "identity") + 
  ggtitle("Count sector incorrectly coded") + xlab("Sector") + ylab("Count") +
  scale_fill_brewer(palette="Dark2")
ggsave("r_mistake_counts.png")
ggplot(r_plot_dat, aes(x = as.factor(sector), y = prop_wrong, group = time, fill = time)) + 
  geom_bar(position = "dodge", stat = "identity") + 
  ggtitle("Proportion sector incorrectly coded") + xlab("Sector") + ylab("Proportion") +
  scale_fill_brewer(palette="Dark2")
ggsave("r_mistake_props.png")


# make transition matrices for sector ----

#make transition matrices for everyone
trans_1 <- array(data = as.integer(0), dim = c(K, K)) 
trans_2 <- array(data = as.integer(0), dim = c(K,K))
for (i in 1:N) {
  trans_1[job[i,1]+1,job[i,2]+1] <- trans_1[job[i,1]+1,job[i,2]+1] + as.integer(1)
  trans_2[job[i,2]+1,job[i,3]+1] <- trans_2[job[i,2]+1,job[i,3]+1] + as.integer(1)
}

trans <- trans_1 + trans_2

x <- xtable(trans)
print.xtable(x, type="latex")

#make transition matrices for people that seem like mistakes 
trans_mist1 <- array(data = as.integer(0), dim = c(K, K)) 
trans_mist2 <- array(data = as.integer(0), dim = c(K, K)) 
for (ind in mist_r1) {
  trans_mist1[job[ind,1]+1,job[ind,2]+1] <- trans_mist1[job[ind,1]+1,job[ind,2]+1] + as.integer(1)
}
for (ind in mist_r2) {
  trans_mist2[job[ind,2]+1,job[ind,3]+1] <- trans_mist2[job[ind,2]+1,job[ind,3]+1] + as.integer(1)
}

trans_mist <- trans_mist1 + trans_mist2

x <- xtable(trans_mist)
print.xtable(x, type="latex")

