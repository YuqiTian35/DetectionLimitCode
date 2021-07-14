library(plyr)
library(lubridate)
library(rms)
library(zoo)
library(tidyverse)
library(rlist)
library(scales)
library(lmtest)



# function for mode
mode <- function(x){
  t <- table(x)
  return(as.numeric(names(t)[which.max(t)]))
}


# data
dat <- data %>% 
  filter(year(baseline_d) >= 2000 & year(rna_d) >= 2000 
         & year(baseline_d) <= 2018 & year(rna_d) <= 2018) %>% 
  filter(center %in% c('ar-ss', 'br-ipec', 'cl-fa', 'mx-incmnsz', 'pe-imtavh')) %>% 
  filter(!is.na(cd4_baseline)) %>% 
  filter(male != 9) %>%  # exlude unknown
  mutate(baseline_d = baseline_d %>% as.yearmon, # convert date format (year-month)
         rna_d = rna_d %>% as.yearmon,
         rna_baseline_d = rna_baseline_d %>% as.yearmon) %>% 
  mutate(dl = 1) %>% # change later
  mutate(cd4_baseline = as.numeric(cd4_baseline),
         rna_baseline = as.numeric(rna_baseline),
         rna_outcome = as.numeric(rna_outcome)) %>% 
  mutate(cd4_baseline_sqrt = sqrt(cd4_baseline),
         rna_baseline_log = log(rna_baseline, 10),
         year = year(rna_d))

# mode value for DL at each site and month
rna_mode_val <- dat %>% 
  filter(!is.na(rna_l)) %>%  # only need DL information
  group_by(site, rna_d) %>% 
  dplyr::summarise(dl = mode(rna_l)) %>% 
  ungroup %>% 
  mutate(site = as.factor(site))

# time range for each site
site_date <- dat %>% 
  group_by(site) %>% 
  dplyr::summarise(min_month = min(rna_d), 
                   max_month =max(rna_d))

# DL value for each site and month
rna_dl_val <- lapply(unique(rna_mode_val$site), function(x){
  temp <- site_date %>% filter(site == x)
  return(data.frame(site = x,
                    time = seq(as_date(temp$min_month), as_date(temp$max_month), by = '1 month') %>% 
                      as.yearmon))}) %>% list.rbind %>% 
  left_join(rna_mode_val, c("site" = "site", "time" = "rna_d")) 
  

# impute DL values for outcome
colnames(rna_dl_val)[3] <- 'dl_val'

dat_obs <- dat %>% 
  mutate(site = as.factor(site)) %>% 
  left_join(rna_dl_val, c("site" = "site", "rna_d" = "time")) 

# impute DL values for baseline
colnames(rna_dl_val)[3] <- 'dl_val_baseline'

dat_obs <- dat_obs %>% 
  mutate(site = as.factor(site)) %>% 
  left_join(rna_dl_val, c("site" = "site", "rna_d" = "time")) %>% 
  mutate(rna_baseline_l = ifelse(is.na(rna_baseline_l), dl_val_baseline, rna_baseline_l),
         rna_l = ifelse(is.na(rna_l), dl_val, rna_l)) %>% 
  filter(rna_baseline >= rna_baseline_l) %>%  # exclude baseline rna < dl
  mutate(rna_outcome = ifelse(rna_outcome < dl_val, dl_val, rna_outcome),
         dl = ifelse(rna_outcome < dl_val, 0, 1))

mod <- multipleDL(rna_outcome ~ age + male + site + route + prior_aids + cd4_baseline_sqrt + 
                    rna_baseline_log + regimen + rna_time + year,
                  data = dat_obs,
                  delta_lower = dat_obs$dl,
                  link = 'logit')

dummy_X <- model.matrix(~ age + male + site + route + prior_aids + cd4_baseline_sqrt + 
                          rna_baseline_log + regimen + rna_time + year, 
                        data = dat_obs[1,])[,-1]

##### age ####
new.data <- matrix(NA, ncol=ncol(dummy_X), nrow=length(18:82))
colnames(new.data) <- colnames(dummy_X)

new <- 18:82
new.data[,1] <- new
new.data[,2] <- 1
for(i in 3:6){
  new.data[,i] <- dummy_X[which(dat_obs$site == "peru")[1], i]
}
for(i in 7:8){
  new.data[,i] <- dummy_X[which(dat_obs$route == "homo/bisexual")[1], i]
}
new.data[,9] <- FALSE
new.data[,10] <- median(dat_obs %>%  pull(cd4_baseline_sqrt))
new.data[,11] <- median(dat_obs %>%  pull(rna_baseline_log))
for(i in 12:14){
  new.data[,i] <- dummy_X[which(dat_obs$regimen == "NNRTI-based")[1], i]
}
new.data[,15] <- median(dat_obs %>%  pull(rna_time))
new.data[,16] <- median(dat_obs %>%  pull(year))

# quantiles
q50.est <- quantile_dl(mod, new.data, 0.5)
q90.est <- quantile_dl(mod, new.data, 0.9)

# cdf
c1000.est <- cdf_dl(mod, new.data, 1000)
c20.est <- cdf_dl(mod, new.data, 20)

###### Prior AIDS  #####
new.data <- matrix(NA, ncol=ncol(dummy_X), nrow=2)
colnames(new.data) <- colnames(dummy_X)

new <- c(0, 1)
new.data[,1] <- 35
new.data[,2] <- 1
for(i in 3:6){
  new.data[,i] <- dummy_X[which(dat_obs$site == "peru")[1], i]
}
for(i in 7:8){
  new.data[,i] <- dummy_X[which(dat_obs$route == "homo/bisexual")[1], i]
}
new.data[,9] <- new
new.data[,10] <- median(dat_obs %>%  pull(cd4_baseline_sqrt)) 
new.data[,11] <- median(dat_obs %>%  pull(rna_baseline_log))
for(i in 12:14){
  new.data[,i] <- dummy_X[which(dat_obs$regimen == "NNRTI-based")[1], i]
}
new.data[,15] <- median(dat_obs %>%  pull(rna_time)) 
new.data[,16] <- median(dat_obs %>%  pull(year)) 

# quantiles
q50.est <- quantile_dl(mod, new.data, 0.5)
q90.est <- quantile_dl(mod, new.data, 0.9)

# cdf
c1000.est <- cdf_dl(mod, new.data, 1000)
c20.est <- cdf_dl(mod, new.data, 20)

##### Baseline viral load #####
new <- seq(2, 6, 0.1)
new.data <- matrix(NA, ncol=ncol(dummy_X), nrow=length(new))
colnames(new.data) <- colnames(dummy_X)

new.data[,1] <- 35
new.data[,2] <- 1
for(i in 3:6){
  new.data[,i] <- dummy_X[which(dat_obs$site == "peru")[1], i]
}
for(i in 7:8){
  new.data[,i] <- dummy_X[which(dat_obs$route == "homo/bisexual")[1], i]
}
new.data[,9] <- FALSE
new.data[,10] <- median(dat_obs %>%  pull(cd4_baseline_sqrt))
new.data[,11] <- new
for(i in 12:14){
  new.data[,i] <- dummy_X[which(dat_obs$regimen == "NNRTI-based")[1], i]
}
new.data[,15] <- median(dat_obs %>%  pull(rna_time))
new.data[,16] <- median(dat_obs %>%  pull(year))

# quantiles
q50.est <- quantile_dl(mod, new.data, 0.5)
q90.est <- quantile_dl(mod, new.data, 0.9)

# cdf
c1000.est <- cdf_dl(mod, new.data, 1000)
c20.est <- cdf_dl(mod, new.data, 20)

####### Logistic regression #####
dat_logistic <- dat_obs %>% 
  mutate(binary_outcome = ifelse(rna_outcome < 400, 0, 1))

mod_logistic <- glm(binary_outcome ~ age + male + site + route + prior_aids + cd4_baseline_sqrt + 
                      rna_baseline_log + regimen + rna_time + year,
                    data = dat_logistic,
                    family = 'binomial')
summary(mod_logistic)

####### Likelihood approach #######
dat_mle <- dat_obs %>% 
  mutate(outcome = ifelse(rna_outcome < dl_val, dl_val, rna_outcome),
         delta = ifelse(rna_outcome < dl_val, 0, 1),
         outcome_log10 = log(outcome, 10))


mod_mle <- survreg(Surv(dat_mle$outcome_log10, dat_mle$delta, type='left') ~ age + male + site + 
                     route + prior_aids + cd4_baseline_sqrt +  rna_baseline_log + regimen + rna_time + year, 
                   dist = 'gaussian',
                   data = dat_mle)

summary(mod_mle)