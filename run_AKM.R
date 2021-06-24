library(tidyverse)
load("data/deflate_y_dat.rda")
theme_update(plot.title = element_text(hjust = 0.5))

dat <- deflate_y_dat
N <- dat$N
Ti <- dat$T
K <- dat$K
log_y <- log(dat$y)
job <- array(data = NA, dim = c(N,Ti))
for (i in 1:N) {
  for (t in 1:Ti) {
    job[i,t] <- which(dat$a[i,t,]==1)
  }
}

df <- data.frame(log_wage = c(log_y[,1],log_y[,2],log_y[,3]),
                  worker = as.factor(rep(1:N,Ti)),
                  sector = as.factor(c(job[,1],job[,2],job[,3])))
mod <- lm(log_wage ~ -1 + worker + sector, data = df)
tot_var <- var(df$log_wage)
worker_eff <- mod$coefficients[1:N]
sector_eff <- mod$coefficients[N:length(mod$coefficients)]
sector_eff[1] <- 0
worker_var <- var(worker_eff)/tot_var
sector_var <- var(sector_eff)/tot_var

df <- df %>%
  mutate(worker_effect = worker_eff[worker],
         sector_effect = sector_eff[sector])
worker_sector_cov <- cov(df$worker_effect, df$sector_effect)
worker_var # 56% of variance from worker
sector_var # 5% of variance from worker
2*worker_sector_cov/tot_var # 3% of variance from covariance between worker and sector 
cor(df$worker_effect, df$sector_effect)

df <- df %>%
  mutate(expected = worker_effect + sector_effect)
df_excluded <- df %>%
  filter(expected > 5)
ggplot(df_excluded, aes(x = expected)) + geom_histogram() +
  xlab("Expected Log Wage") + ylab("Count") + 
  ggtitle("AKM Results") 
ggsave("figs/AKM_expected_wage.png")
summary(df_excluded$expected)
sd(df_excluded$expected)
sd(df_excluded$expected)/sqrt(nrow(df_excluded))
