library(rstan)
library(tidyverse)
library(splitstackshape)
library(stringi)
load("fit_ifls_2-25.rda")

fit <- fit_ifls
fit_summ <- as.data.frame(summary(fit)$summary)
r_hat_z <- data.frame(par = rownames(fit_summ), r_hat = fit_summ[,10]) %>%
  filter(grepl("z\\[", par)) %>%
  cSplit("par",direction = "wide") %>%
  mutate(type = "z",
         worker = stri_extract_first_regex(par_1,"[0-9]+"),
         dim = stri_extract_first_regex(par_2,"[0-9]+")) %>%
  select(-c("par_1","par_2")) %>% 
  group_by(worker) %>%
  summarise(big_rhat = min(r_hat)) %>%
  arrange(as.numeric(worker))
r_hat_p <- data.frame(par = rownames(fit_summ), r_hat = fit_summ[,10]) %>%
  filter(grepl("p\\[", par)) %>%
  cSplit("par",direction = "wide") %>%
  mutate(type = "p",
         worker = stri_extract_first_regex(par_1,"[0-9]+"),
         time = stri_extract_first_regex(par_2,"[0-9]+"),
         class = stri_extract_first_regex(par_3,"[0-9]+")) %>%
  select(-c("par_1","par_2","par_3")) %>% 
  group_by(worker) %>%
  summarise(big_rhat = min(r_hat)) %>%
  arrange(as.numeric(worker))
plot_dat <- data.frame(worker = as.numeric(r_hat_p$worker), 
                       rhat_p = r_hat_p$big_rhat,
                       rhat_z = r_hat_z$big_rhat)
ggplot(plot_dat, aes(x = rhat_p, y = rhat_z)) + geom_point() +
  xlab("Maximum p r-hat") + ylab("Maximum z r-hat")
ggsave("figs/p_z_rhat_scatterplot.png")
