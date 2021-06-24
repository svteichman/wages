library(bayesplot)
library(rstan)
library(tidyverse)

K <- 7
Ti <- 3
N <- 1973

#load("fit_ifls_2-13.rda")
#fit <- fit_ifls
load("stan_res/fit_test_ifls_2-17.rda")
fit <- fit_test_ifls

#rhats
pdf("2-16_rhat.pdf") 
stan_rhat(fit)
dev.off()

#Ratio of effective sample size / Sample size 
pdf("2-16_ess.pdf") 
stan_ess(fit)
dev.off()

#ratio of Monte Carlo standard error to posterior standard deviation 
pdf("2-16_mcse.pdf") 
print(stan_mcse(fit))
dev.off()

#effective sample size 
pdf("2-16_ess.pdf") 
print(stan_ess(fit))
dev.off()

fit_summ <- as.data.frame(summary(fit)$summary)

w_medians <- fit_summ$`50%`[4:(4+K*2*Ti)]
pdf("2-16_w_medians.pdf") 
hist(w_medians)
dev.off()

z_medians <- fit_summ$`50%`[(5+K*2*Ti):(5+K*2*Ti + N*2)]
pdf("2-16_z_medians.pdf") 
hist(z_medians)
dev.off()

sigma_params <- sort(grep("sigma",row.names(fit_summ), value=T)[1:32])[1:16]
beta_k <- grep("beta_k",row.names(fit_summ), value=T)
u <- grep("u",row.names(fit_summ), value=T)[2:8]
v <- grep("v",row.names(fit_summ), value=T)[2:8]
q <- grep("q",row.names(fit_summ), value=T)[8:14]
other_1d <- c("alpha[1]","alpha[2]","delta[1]","delta[2]","beta[1]","beta[2]","zeta[1]","zeta[2]")
lambda <- grep("lambda",row.names(fit_summ), value=T)
eta <- grep("eta",row.names(fit_summ), value=T)
kappa <- grep("kappa",row.names(fit_summ), value=T)

### Individual Parameter Diagnostics 

get_single_plots <- function(fit, param) { 
  print(stan_ac(fit, pars = param))
  print(rstan::traceplot(fit, pars = param))
} 
get_aggreg_plots <- function(fit, param, trim = F, trim_amount) {
  ind <- grep(paste0("^",param), rownames(as.data.frame(summary(fit)$summary)))
  means <- data.frame(avg = as.data.frame(summary(fit)$summary)$mean[ind])
  print(summary(means))
  title <- paste0("Posterior Means of ",param)
  print(ggplot(means, aes(x = avg)) + geom_histogram(bins = 30) + ggtitle(title))
  if (trim == T) {
    lim <- quantile(abs(means$avg), probs = trim_amount)
    means_trim <- means %>% filter(abs(means$avg) < lim)
    print(ggplot(means_trim, aes(x = avg)) + geom_histogram(bins = 30) + 
            ggtitle(paste0(title, " Without Extremes")))
    means_out <- means %>% filter(abs(means$avg) > lim) 
    summary(means_out)
  }
}
plot_fit <- function(fit) {
  get_single_plots(fit, sigma_params)
  get_single_plots(fit, beta_k)
  get_single_plots(fit, other_1d)
  get_single_plots(fit, u)
  get_single_plots(fit, v)
  get_single_plots(fit, q)
  get_aggreg_plots(fit, "w")
  get_aggreg_plots(fit, "z")
  get_aggreg_plots(fit, "p", trim = T, trim_amount = .90)
  get_aggreg_plots(fit, "lambda", trim = T, trim_amount = .60)
  get_aggreg_plots(fit, "eta", trim = T, trim_amount = .60)
  get_aggreg_plots(fit, "kappa", trim = T, trim_amount = .60)
}
plot_fit(fit)

#find large Rhats
summary(fit_summ$Rhat)
big_Rhat <- fit_summ$Rhat > 6
big_Rhat_dat <- fit_summ[big_Rhat,]
