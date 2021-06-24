library(tidyverse)
library(splitstackshape)
library(stringi)
load("data/deflate_y_dat.rda")

dat <- deflate_y_dat
K <- 7
Ti <- 3
N <- 1973

labs <- c("Agriculture","Manufacturing","Construction","Hospitality",
          "Transportation","Business","Social Services")

log_y <- log(dat$y)
job <- array(data = NA, dim = c(N,Ti))
for (i in 1:N) {
  for (t in 1:Ti) {
    job[i,t] <- which(dat$a[i,t,]==1)
  }
}


# real_trans <- data.frame(past = c(job[,1],job[,2]), current = c(job[,2],job[,3]), 
#                          time = c(rep(2,N),rep(3,N))) %>%
#   group_by(past,time) %>%
#   mutate(tot_sum = n(),
#          prop = 1/tot_sum) %>%
#   group_by(past,current,time) %>%
#   summarise(prob_sum = sum(prop)) %>% 
#   arrange(past,time) %>%
#   mutate(current1 = labs[current],
#          past1 = labs[past])

res <- matrix(0, nrow = K, ncol = 2)
for (i in 1:N) {
  col <- 1
  if (job[i,1] != job[i,2]) {col <- 2}
  res[job[i,1],col] <- res[job[i,1],col] + 1
  
  col <- 1
  if (job[i,2] != job[i,3]) {col <- 2}
  res[job[i,2],col] <- res[job[i,2],col] + 1
}
res <- as.data.frame(res) %>%
  mutate(size = V1+V2,
         sector = labs) %>% 
  arrange(size) %>%
  mutate(prop_stay = V1/size,
         prop_leave = V2/size) %>%
  select(-c(V1,V2))
ggplot(res, aes(x = prop_stay, y = prop_leave, color = size)) + geom_jitter() +
  xlim(c(0,1)) + ylim(c(0,1)) +
  xlab("Proportion of people that stay") +
  ylab("Proportion of people that leave") + 
  geom_text_repel(x = res$prop_stay, y = res$prop_leave, label = res$size, segment.color = NA)
ggsave("figs/prop_stay.png")

res$prop_size <- res$size/(N*2)
res$size
res$prop_stay
