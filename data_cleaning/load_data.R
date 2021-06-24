library(foreign)
library(tidyverse)

add_prefixes <- function(names,num) {
  names[-1] <- paste0("w",num,"_",names[-1])
  return(names)
}

#loading in wave 5 
w5_bk_ar1 <- read.dta(file = "wave5/wave5_bookK/bk_ar1.dta")
w5_b3a_tk2 <- read.dta(file = "wave5/wave5_book3a/b3a_tk2.dta")
wave5_K <- data.frame(id = as.integer(w5_bk_ar1$pidlink), sex = as.factor(w5_bk_ar1$ar07), age = w5_bk_ar1$ar09) %>% 
  distinct(id, .keep_all = T) %>%
  mutate(sex01 = ifelse(is.na(sex),NA,ifelse(sex == 1, 0, 1))) %>%
  select(-sex)
# in original sex factor, male is 1, female is 3
# in sex01, male is 0, female is 1
wave5_TK <- data.frame(id = as.integer(w5_b3a_tk2$pidlink), sector = w5_b3a_tk2$tk19ab, num_ppl = w5_b3a_tk2$tk20aa, 
                       hrs_most_wk = w5_b3a_tk2$tk22a, wk_per_yr = w5_b3a_tk2$tk23a, 
                       yr_curr_job = w5_b3a_tk2$tk23a2y, mth_curr_job = w5_b3a_tk2$tk23a2m,
                       work_cat = w5_b3a_tk2$tk24a, wage_mth = w5_b3a_tk2$tk25a1, 
                       wage_yr = w5_b3a_tk2$tk25a2, occ = w5_b3a_tk2$occ2014_primary) %>% distinct()
wave5 <- inner_join(wave5_K,wave5_TK, by = "id")
#wave5_wage <- wave5 %>% filter(!is.na(wage_mth), !is.na(wage_yr))
wave5_wage <- wave5 %>% filter(!is.na(wage_yr))
names(wave5_wage) <- add_prefixes(names(wave5_wage),5)

#loading in wave 4 
w4_bk_ar1 <- read.dta(file = "wave4/wave4_bookK/bk_ar1.dta")
w4_b3a_tk2 <- read.dta(file = "wave4/wave4_book3a/b3a_tk2.dta")
wave4_K <- data.frame(id = as.integer(w4_bk_ar1$pidlink), sex = as.factor(w4_bk_ar1$ar07), age = w4_bk_ar1$ar09) %>% 
  distinct(id, .keep_all = T) %>%
  mutate(sex01 = ifelse(is.na(sex),NA,ifelse(sex == 1, 0, 1))) %>%
  select(-sex)
wave4_TK <- data.frame(id = as.integer(w4_b3a_tk2$pidlink), sector = w4_b3a_tk2$tk19ab, num_ppl = w4_b3a_tk2$tk20aa, 
                       hrs_most_wk = w4_b3a_tk2$tk22a, wk_per_yr = w4_b3a_tk2$tk23a, 
                       yr_curr_job = w4_b3a_tk2$tk23a2y, mth_curr_job = w4_b3a_tk2$tk23a2m,
                       work_cat = w4_b3a_tk2$tk24a, wage_mth = w4_b3a_tk2$tk25a1, 
                       wage_yr = w4_b3a_tk2$tk25a2, occ = w4_b3a_tk2$occ07tk2) %>% distinct()
wave4 <- inner_join(wave4_K,wave4_TK, by = "id")
#wave4_wage <- wave4 %>% filter(!is.na(wage_mth), !is.na(wage_yr))
wave4_wage <- wave4 %>% filter(!is.na(wage_yr))
names(wave4_wage) <- add_prefixes(names(wave4_wage),4)

#loading in wave 3 
w3_bk_ar1 <- read.dta(file = "wave3/wave3_bookK/bk_ar1.dta")
w3_b3a_tk2 <- read.dta(file = "wave3/wave3_book3a/b3a_tk2.dta")
wave3_K <- data.frame(id = as.integer(w3_bk_ar1$pidlink), sex = as.factor(w3_bk_ar1$ar07), age = w3_bk_ar1$ar09) %>% 
  distinct(id, .keep_all = T) %>%
  mutate(sex01 = ifelse(is.na(sex),NA,ifelse(sex == 1, 0, 1))) %>%
  select(-sex)
wave3_TK <- data.frame(id = as.integer(w3_b3a_tk2$pidlink), sector = w3_b3a_tk2$tk19aa, num_ppl = w3_b3a_tk2$tk20aa, 
                       hrs_most_wk = w3_b3a_tk2$tk22a, wk_per_yr = w3_b3a_tk2$tk23a, 
                       yr_curr_job = w3_b3a_tk2$tk23a2,
                       work_cat = w3_b3a_tk2$tk24a, wage_mth = w3_b3a_tk2$tk25a1, 
                       wage_yr = w3_b3a_tk2$tk25a2, occ = w3_b3a_tk2$tk20ab) %>% distinct()
wave3 <- inner_join(wave3_K,wave3_TK, by = "id")
#wave3_wage <- wave3 %>% filter(!is.na(wage_mth), !is.na(wage_yr))
wave3_wage <- wave3 %>% filter(!is.na(wage_yr))
names(wave3_wage) <- add_prefixes(names(wave3_wage),3)

#loading in wave 2 
w2_bk_ar1 <- read.dta(file = "wave2/wave2_bookK/bk_ar1.dta")
w2_b3a_tk2 <- read.dta(file = "wave2/wave2_book3/b3a_tk2.dta")
w2_b3a_tk3 <- read.dta(file = "wave2/wave2_book3/b3a_tk3.dta")
wave2_K <- data.frame(id = as.integer(w2_bk_ar1$pidlink), sex = as.factor(w2_bk_ar1$ar07), age = w2_bk_ar1$ar09) %>% 
  distinct(id, .keep_all = T) %>%
  mutate(sex01 = ifelse(is.na(sex),NA,ifelse(sex == "1. M", 0, 1))) %>% 
  select(-sex)
wave2_TK <- data.frame(id = as.integer(w2_b3a_tk2$pidlink), hrs_most_wk = w2_b3a_tk2$tk22a,
                       wks_per_yr = w2_b3a_tk2$tk23a, work_cat = w2_b3a_tk2$tk24a,
                       wage_mth = w2_b3a_tk2$tk25am, wage_yr = w2_b3a_tk2$tk25ay, 
                       occ = w2_b3a_tk2$tk20aocc)
wave2_hist <- data.frame(id = as.integer(w2_b3a_tk3$pidlink), same = w2_b3a_tk3$tk29, year = w2_b3a_tk3$tk28yr) %>%
  filter(year == "1993") %>%
  select(-year)
wave2 <- inner_join(wave2_K,wave2_TK, by = "id")
wave2 <- inner_join(wave2, wave2_hist, by = "id")
#wave2_wage <- wave2 %>% filter(!is.na(wage_mth), !is.na(wage_yr))
wave2_wage <- wave2 %>% filter(!is.na(wage_yr))
names(wave2_wage) <- add_prefixes(names(wave2_wage),2)

#loading in wave 1
w1_bk_ar2 <- read.dta(file = "wave1/wave1_bookK/bukkar2.dta")
w1_b3a_tk2 <- read.dta(file = "wave1/wave1_book3/buk3tk2.dta")
w1_b3a_tk3 <- read.dta(file = "wave1/wave1_book3/buk3tk3.dta")
wave1_K <- data.frame(id = as.integer(w1_bk_ar2$pidlink), sex = as.factor(w1_bk_ar2$ar07), age = w1_bk_ar2$ar09yr) %>%
  mutate(sex01 = ifelse(is.na(sex) | sex == 6 | sex == 9,NA,ifelse(sex == 1, 0, 1))) %>% 
  select(-sex)
wave1_TK <- data.frame(id = as.integer(w1_b3a_tk2$pidlink), hrs_most_wk = w1_b3a_tk2$tk22a,
                       wks_per_yr = w1_b3a_tk2$tk23a, work_cat = w1_b3a_tk2$tk24a,
                       wage_mth = w1_b3a_tk2$tk25r1_m, wage_yr = w1_b3a_tk2$tk25r1_y, 
                       occ = w1_b3a_tk2$occ20a) 
wave1 <- inner_join(wave1_K, wave1_TK, by = "id")
#wave1_wage <- wave1 %>% filter(!is.na(wage_mth) & !is.na(wage_yr))
wave1_wage <- wave1 %>% filter(!is.na(wage_yr))
names(wave1_wage) <- add_prefixes(names(wave1_wage),1)

wave1to2 <- inner_join(wave1_wage, wave2_wage, by = "id")
wave1to3 <- inner_join(wave1to2, wave3_wage, by = "id")
wave1to4 <- inner_join(wave1to3, wave4_wage, by = "id")
wave1to5 <- inner_join(wave1to4, wave5_wage, by = "id") %>%
  filter(w1_wage_yr != 0 & w2_wage_yr != 0 & w3_wage_yr != 0 & w4_wage_yr != 0 & w5_wage_yr != 0)
save(wave1to5, file = "wave1to5.rda")

wave3to4 <- inner_join(wave3_wage, wave4_wage, by = "id")
wave3to5 <- inner_join(wave3to4, wave5_wage, by = "id")
save(wave3to5, file = "wave3to5.rda")

w1_bk_sc1 <- read.dta(file = "wave1/wave1_bookK/bukksc1.dta")
table(w1_bk_sc1$sc01)
table(w1_bk_sc1$sc02)
table(w1_bk_sc1$sc03)

