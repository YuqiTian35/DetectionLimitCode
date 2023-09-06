library(plyr)
library(lubridate)
library(rms)
library(zoo)
library(tidyverse)
library(rlist)
library(scales)
library(lmtest)
library(multipleDL)



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

# appropriate baseline for factor
dat_recode$male <- as.factor(dat_recode$male)
dat_recode$male <- relevel(dat_recode$male, ref = '1')

dat_recode$site <- relevel(dat_recode$site, ref = 'peru')

dat_recode$center <- as.factor(dat_recode$center)
dat_recode$center <- relevel(dat_recode$center, ref = 'pe-imtavh')

dat_recode$route <- as.factor(dat_recode$route)
dat_recode$route <- relevel(dat_recode$route, ref = 'homo/bisexual')

dat_recode$regimen <- as.factor(dat_recode$regimen)
dat_recode$regimen <- relevel(dat_recode$regimen, ref = 'NNRTI-based')


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

######## Residual Plots #########
mod_logit <- mod
cof <- mod_logit$coef # original model
k <- length(cof) - mod_logit$p
N <- dim(dat_recode)[1]
X <- mod$x 
Y <- dat_recode$rna_outcome
if(!is.factor(Y))Y <- factor(Y)
Y <- unclass(Y) - 1
cumprob <- plogis
# cumprob = pnorm (probit link), function(x) exp(-exp(x)) (loglog link), function(x) 1 - exp(-exp(x)) (cloglog link)

px <- cdf_dl(mod_logit, X, at.y = mod$yunique, se=F)$est %>% t
low.x <- rbind(0, px)[cbind(Y + 1L, 1:N)]
hi.x  <- 1 - rbind(px, 1)[cbind(Y + 1L, 1:N)]
cn_res_logit <- low.x - hi.x
# add - sign to flip right-censoring to left-censoring
fitpsr_logit <- survfit(Surv(-cn_res_logit, dat_recode$dl) ~ 1)
# use 1-() to flip
plot(1-summary(fitpsr_logit)$time, 1-qunif(1-summary(fitpsr_logit)$surv, -1, 1),
     xlab="", ylab="")
abline(0,1,col=gray(.5))
mtext("Quantiles of Unif(-1,1)",side=1,line=2,cex=.7)
mtext("Quantiles of Cox-Snell-like PSR",side=2,line=2.5,cex=.7,padj=1)
title('Logit link')

# same code for model with probit, loglog, and cloglog links

######### Residual-over-covariate plots #####
cof <- mod_logit$coef # original model
k <- length(cof) - mod$p
N <- dim(dat_recode)[1]
X <- mod_logit$x 
Y <- dat_recode$rna_outcome
if(!is.factor(Y))Y <- factor(Y)
Y <- unclass(Y) - 1
cumprob <- plogis

px <- cdf_dl(mod_logit, X, at.y = mod$yunique, se=F)$est %>% t
low.x <- rbind(0, px)[cbind(Y + 1L, 1:N)]
hi.x  <- 1 - rbind(px, 1)[cbind(Y + 1L, 1:N)]
observed <- low.x - hi.x
censored <- -hi.x

psr <- observed
for(i in 1:N){
  if(dat_recode$dl[i] == 0){
    psr[i] <- censored[i]
  }
}

plot(dat_recode$age, psr, type="n",xlim=c(18,85),xlab="",ylab="",ylim=c(-1,1))
points(dat_recode$age[dat_recode$dl==0],psr[dat_recode$dl==0],pch=4,cex=.5,col=gray(.6))
points(dat_recode$age[dat_recode$dl==1],psr[dat_recode$dl==1],pch=1,cex=.5,col=gray(.6))
abline(h=0,lty=2)
lines(supsmu(dat_recode$age, psr, bass=8),col=1,lwd=4)
mtext("Age",side=1,line=2,cex=.7 )
mtext("PSR",side=2,line=2.5,cex=.7,padj=1)

plot(dat_recode$cd4_baseline_sqrt^2, psr, type="n",xlim=c(0, 500),xlab="",ylab="",ylim=c(-1,1))
points(dat_recode$cd4_baseline_sqrt[dat_recode$dl==0]^2, psr[dat_recode$dl==0],pch=4,cex=.5,col=gray(.6))
points(dat_recode$cd4_baseline_sqrt[dat_recode$dl==1]^2, psr[dat_recode$dl==1],pch=1,cex=.5,col=gray(.6))
abline(h=0,lty=2)
lines(supsmu(dat_recode$cd4_baseline_sqrt^2, psr, bass=1),col=1,lwd=4)
mtext("Baseline CD4",side=1,line=2,cex=.7 )
mtext("PSR",side=2,line=2.5,cex=.7,padj=1)

plot(exp(dat_recode$rna_baseline_log), psr, type="n",xlim=c(0, 500),xlab="",ylab="",ylim=c(-1,1))
points(exp(dat_recode$rna_baseline_log[dat_recode$dl==0]),psr[dat_recode$dl==0],pch=1,cex=.5,col=gray(.6))
points(exp(dat_recode$rna_baseline_log[dat_recode$dl==1]),psr[dat_recode$dl==1],pch=1,cex=.5,col=gray(.6))
abline(h=0,lty=2)
lines(supsmu(exp(dat_recode$rna_baseline_log), psr, bass=1),col=1,lwd=4)
mtext("Baseline Viral Load",side=1,line=2,cex=.7 )
mtext("PSR",side=2,line=2.5,cex=.7,padj=1)

plot(dat_recode$rna_time/30, psr, type="n",xlim=c(2, 15),xlab="",ylab="",ylim=c(-1,1))
points(dat_recode$rna_time[dat_recode$dl==0]/30,psr[dat_recode$dl==0],pch=4,cex=.5,col=gray(.6))
points(dat_recode$rna_time[dat_recode$dl==1]/30,psr[dat_recode$dl==1],pch=1,cex=.5,col=gray(.6))
abline(h=0,lty=2)
lines(supsmu(dat_recode$rna_time/30, psr, bass=1),col=1,lwd=4)
mtext("Months to VL measure",side=1,line=2,cex=.7 )
mtext("PSR",side=2,line=2.5,cex=.7,padj=1)

######### Mexico vs. Non-Mexico
dat_recode <- dat_recode %>% 
  mutate(is_mexico = ifelse(site == 'mexico', TRUE, FALSE))

dat_recode_mexico <- dat_recode %>% dplyr::filter(is_mexico == TRUE)
dat_recode_rest <- dat_recode %>% dplyr::filter(is_mexico == FALSE)

mod_mexico <- multipleDL(rna_outcome ~ age + male + route + prior_aids + cd4_baseline_sqrt + 
                           rna_baseline_log + regimen + rna_time + year_since_baseline,
                         data = dat_recode_mexico,
                         delta_lower = dat_recode_mexico$dl,
                         link = 'logit')


beta <- mod_mexico$coef[(length(mod_mexico$coef)-11):length(mod_mexico$coef)]
beta.se <- (mod_mexico$var %>% diag %>% sqrt)[(length(mod_mexico$coef)-11):length(mod_mexico$coef)]
names(beta.se) <- names(beta) <- variable_name

# p-value
pval <- (1 - pnorm(abs(beta / beta.se))) * 2

tab <- data.frame(coef=beta, se = beta.se, waldz = beta / beta.se, pval)
tab %>% 
  kable(digits=4) %>% 
  kable_styling(full_width = F)

# odds ratio
or <- exp(beta)
# odds ratio lower bound
or_lb <- exp(beta - qnorm(0.975) * beta.se)
# odds ratio upper bound
or_ub <- exp(beta + qnorm(0.975) * beta.se)
# p-value
pval <- (1 - pnorm(abs(beta / beta.se))) * 2

# probabilistic index
p_ind_func <- function(z){
  return(exp(z) * (exp(z) - z - 1) / (exp(z) - 1)^2)
}
p_ind <- p_ind_func(beta)

# table for odds rato
tab_or <- data.frame(or = or, or_lb = or_lb, or_ub = or_ub, pval, p_ind)

# age (per 10 years)
tab_or[1, 'or'] <- exp(beta[1] * 10)
tab_or[1, 'or_lb'] <- exp((beta[1] - qnorm(0.975) * beta.se[1]) * 10)
tab_or[1, 'or_ub'] <- exp((beta[1] + qnorm(0.975) * beta.se[1]) * 10)
tab_or[1, 'p_ind'] <- p_ind_func(beta[1] * 10)

# rna_time (per month)
tab_or[11, 'or'] <- exp(beta[11] * 30)
tab_or[11, 'or_lb'] <- exp((beta[11] - qnorm(0.975) * beta.se[11]) * 30)
tab_or[11, 'or_ub'] <- exp((beta[11] + qnorm(0.975) * beta.se[11]) * 30)
tab_or[11, 'p_ind'] <- p_ind_func(beta[11] * 30)

mod_rest <- multipleDL(rna_outcome ~ age + male + route + prior_aids + cd4_baseline_sqrt + 
                         rna_baseline_log + regimen + rna_time + year_since_baseline,
                       data = dat_recode_rest,
                       delta_lower = dat_recode_rest$dl,
                       link = 'logit')


beta <- mod_rest$coef[(length(mod_rest$coef)-11):length(mod_rest$coef)]
beta.se <- (mod_rest$var %>% diag %>% sqrt)[(length(mod_rest$coef)-11):length(mod_rest$coef)]
names(beta.se) <- names(beta) <- variable_name

# p-value
pval <- (1 - pnorm(abs(beta / beta.se))) * 2

tab <- data.frame(coef=beta, se = beta.se, waldz = beta / beta.se, pval)

# odds ratio
or <- exp(beta)
# odds ratio lower bound
or_lb <- exp(beta - qnorm(0.975) * beta.se)
# odds ratio upper bound
or_ub <- exp(beta + qnorm(0.975) * beta.se)
# p-value
pval <- (1 - pnorm(abs(beta / beta.se))) * 2

# probabilistic index
p_ind_func <- function(z){
  return(exp(z) * (exp(z) - z - 1) / (exp(z) - 1)^2)
}
p_ind <- p_ind_func(beta)

# table for odds rato
tab_or <- data.frame(or = or, or_lb = or_lb, or_ub = or_ub, pval, p_ind)

# age (per 10 years)
tab_or[1, 'or'] <- exp(beta[1] * 10)
tab_or[1, 'or_lb'] <- exp((beta[1] - qnorm(0.975) * beta.se[1]) * 10)
tab_or[1, 'or_ub'] <- exp((beta[1] + qnorm(0.975) * beta.se[1]) * 10)
tab_or[1, 'p_ind'] <- p_ind_func(beta[1] * 10)

# rna_time (per month)
tab_or[11, 'or'] <- exp(beta[11] * 30)
tab_or[11, 'or_lb'] <- exp((beta[11] - qnorm(0.975) * beta.se[11]) * 30)
tab_or[11, 'or_ub'] <- exp((beta[11] + qnorm(0.975) * beta.se[11]) * 30)
tab_or[11, 'p_ind'] <- p_ind_func(beta[11] * 30)


####### 2000-2009 vs. 2010-2018
dat_recode <- dat_recode %>% 
  dplyr::mutate(year_period = ifelse(year <= 2009, '2000-2009', '2010-2018'))
dat_recode_2000 <- dat_recode %>% dplyr::filter(year_period == '2000-2009')
dat_recode_2010 <- dat_recode %>% dplyr::filter(year_period == '2010-2018')

mmod_2000 <- multipleDL(rna_outcome ~ age + male + site + route + prior_aids + cd4_baseline_sqrt + 
                          rna_baseline_log + regimen + rna_time + year_since_baseline,
                        data = dat_recode_2000,
                        delta_lower = dat_recode_2000$dl,
                        link = 'logit')

beta <- mod_2000$coef[(length(mod_2000$coef)-(length(variable_name)-1)):length(mod_2000$coef)]
beta.se <- (mod_2000$var %>% diag %>% sqrt)[(length(mod_2000$coef)-(length(variable_name)-1)):length(mod_2000$coef)]
names(beta.se) <- names(beta) <- variable_name

# p-value
pval <- (1 - pnorm(abs(beta / beta.se))) * 2

tab <- data.frame(coef=beta, se = beta.se, waldz = beta / beta.se, pval)

# odds ratio
or <- exp(beta)
# odds ratio lower bound
or_lb <- exp(beta - qnorm(0.975) * beta.se)
# odds ratio upper bound
or_ub <- exp(beta + qnorm(0.975) * beta.se)
# p-value
pval <- (1 - pnorm(abs(beta / beta.se))) * 2

# probabilistic index
p_ind_func <- function(z){
  return(exp(z) * (exp(z) - z - 1) / (exp(z) - 1)^2)
}
p_ind <- p_ind_func(beta)

# table for odds rato
tab_or <- data.frame(or = or, or_lb = or_lb, or_ub = or_ub, pval, p_ind)

# age (per 10 years)
tab_or[1, 'or'] <- exp(beta[1] * 10)
tab_or[1, 'or_lb'] <- exp((beta[1] - qnorm(0.975) * beta.se[1]) * 10)
tab_or[1, 'or_ub'] <- exp((beta[1] + qnorm(0.975) * beta.se[1]) * 10)
tab_or[1, 'p_ind'] <- p_ind_func(beta[1] * 10)

# rna_time (per month)
tab_or[15, 'or'] <- exp(beta[15] * 30)
tab_or[15, 'or_lb'] <- exp((beta[15] - qnorm(0.975) * beta.se[15]) * 30)
tab_or[15, 'or_ub'] <- exp((beta[15] + qnorm(0.975) * beta.se[15]) * 30)
tab_or[15, 'p_ind'] <- p_ind_func(beta[15] * 30)

tab_or_2000 <- tab_or 

mod_2010 <- multipleDL(rna_outcome ~ age + male + site + route + prior_aids + cd4_baseline_sqrt + 
                         rna_baseline_log + regimen + rna_time + year_since_baseline,
                       data = dat_recode_2010,
                       delta_lower = dat_recode_2010$dl,
                       link = 'logit')

beta <- mod_2010$coef[(length(mod_2010$coef)-(length(variable_name)-1)):length(mod_2010$coef)]
beta.se <- (mod_2010$var %>% diag %>% sqrt)[(length(mod_2010$coef)-(length(variable_name)-1)):length(mod_2010$coef)]
names(beta.se) <- names(beta) <- variable_name

# p-value
pval <- (1 - pnorm(abs(beta / beta.se))) * 2

tab <- data.frame(coef=beta, se = beta.se, waldz = beta / beta.se, pval)

# odds ratio
or <- exp(beta)
# odds ratio lower bound
or_lb <- exp(beta - qnorm(0.975) * beta.se)
# odds ratio upper bound
or_ub <- exp(beta + qnorm(0.975) * beta.se)
# p-value
pval <- (1 - pnorm(abs(beta / beta.se))) * 2

# probabilistic index
p_ind_func <- function(z){
  return(exp(z) * (exp(z) - z - 1) / (exp(z) - 1)^2)
}
p_ind <- p_ind_func(beta)

# table for odds rato
tab_or <- data.frame(or = or, or_lb = or_lb, or_ub = or_ub, pval, p_ind)

# age (per 10 years)
tab_or[1, 'or'] <- exp(beta[1] * 10)
tab_or[1, 'or_lb'] <- exp((beta[1] - qnorm(0.975) * beta.se[1]) * 10)
tab_or[1, 'or_ub'] <- exp((beta[1] + qnorm(0.975) * beta.se[1]) * 10)
tab_or[1, 'p_ind'] <- p_ind_func(beta[1] * 10)

# rna_time (per month)
tab_or[15, 'or'] <- exp(beta[15] * 30)
tab_or[15, 'or_lb'] <- exp((beta[15] - qnorm(0.975) * beta.se[15]) * 30)
tab_or[15, 'or_ub'] <- exp((beta[15] + qnorm(0.975) * beta.se[15]) * 30)
tab_or[15, 'p_ind'] <- p_ind_func(beta[15] * 30)

tab_or_2010 <- tab_or 

######## With splines ########
# data with splines
dat_rcs <- dat_recode 

# age
age_knots <- rcspline.eval(dat_recode$age, nk=4, knots.only=TRUE)
age_splines <- rcspline.eval(dat_recode$age, inclx=TRUE, knots=age_knots)
colnames(age_splines) <- paste0('age', c('','1','2'))
dat_rcs <- cbind(dat_rcs, age_splines[,-1])

# cd4_baseline_sqrt
cd4_knots <- rcspline.eval(dat_recode$cd4_baseline_sqrt, nk=4, knots.only=TRUE)
cd4_splines <- rcspline.eval(dat_recode$cd4_baseline_sqrt, inclx=TRUE, knots=cd4_knots)
colnames(cd4_splines) <- paste0('cd4_baseline_sqrt', c('','1','2'))
dat_rcs <- cbind(dat_rcs, cd4_splines[,-1])

# rna_baseline_log
rna_knots <- rcspline.eval(dat_recode$rna_baseline_log, nk=4, knots.only=TRUE)
rna_splines <- rcspline.eval(dat_recode$rna_baseline_log, inclx=TRUE, knots=rna_knots)
colnames(rna_splines) <- paste0('rna_baseline_log', c('','1', '2'))
dat_rcs <- cbind(dat_rcs, rna_splines[,-1])

# rna_time
time_knots <- rcspline.eval(dat_recode$rna_time, nk=4, knots.only=TRUE)
time_splines <- rcspline.eval(dat_recode$rna_time, inclx=TRUE, knots=time_knots)
colnames(time_splines) <- paste0('rna_time', c('','1','2'))
dat_rcs <- cbind(dat_rcs, time_splines[,-1])

mod_rcs <- multipleDL(rna_outcome ~ age + age1 + age2 + male + site + route + prior_aids + 
                        cd4_baseline_sqrt + cd4_baseline_sqrt1 + cd4_baseline_sqrt2 + 
                        rna_baseline_log + rna_baseline_log1 + rna_baseline_log2 +
                        regimen + rna_time + rna_time1 + rna_time2 + year,
                      data = dat_rcs,
                      delta_lower = dat_rcs$dl,
                      link = 'logit')

