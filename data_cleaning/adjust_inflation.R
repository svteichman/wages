library(tidyverse)
library(readr)
load("data/stan_dat_rm_sec.rda") 
inflate <- read.csv(file = "inflation_percentages.csv", header = T) 

dat <- stan_dat_rm_sec
ind_inf <- inflate %>%
  select(-c(Country.Code, Indicator.Name, Indicator.Code)) %>%
  filter(Country.Name == "Indonesia",) %>%
  pivot_longer(cols = -Country.Name, names_to = "year") %>%
  mutate(year = parse_number(year),
         percent = value/100,
         inflation = 1 + percent) %>%
  filter(year %in% 2001:2014)

#deflate 2007 and 2014 wages to match 2000
y <- dat$y
new_y <- y
denom2 <- prod(ind_inf$inflation[1:7])
new_y[,2] <- y[,2]/denom2
denom3 <- prod(ind_inf$inflation[1:14])
new_y[,3] <- y[,3]/denom3
dat$y <- new_y
deflate_y_dat <- dat 
save(deflate_y_dat, file = "data/deflate_y_dat.rda")
