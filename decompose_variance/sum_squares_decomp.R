library(bayesplot)
library(rstan)
load("fit_ifls_7-22a.rda")
load("deflate_y_dat.rda")

dat <- deflate_y_dat
fit <- fit_ifls

ind <- grep("^mu", rownames(as.data.frame(summary(fit)$summary)))
means <- summary(fit)$summary[ind,1]

log_y <- log(dat$y)
df <- data.frame(log_wage = c(log_y[,1],log_y[,2],log_y[,3]),
                 mu = rep(means,3))
tot_var <- var(df$log_wage)
mean_var <- var(df$mu)

tot_var
mean_var
tot_var - mean_var