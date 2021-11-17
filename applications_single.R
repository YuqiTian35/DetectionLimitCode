library(plyr)
library(lubridate)
library(rms)
library(zoo)
library(tidyverse)
library(rlist)
library(scales)
library(lmtest)

dat <- all_dat %>% 
  dplyr::rename(id = HATIM.ID,
                age = Age,
                sex = Sex,
                bmi = Body.mass.index..BMI.,
                status = Final.Study.Group,
                il_4 = IL.4..pg.ml.,
                il_12p70 = IL.12p70..pg.ml.) %>% 
  dplyr::select(id, age, sex, bmi, status, il_4, il_12p70) %>% 
  mutate(il_4 = replace_na(il_4, 0.018),
         il_12p70 = replace_na(il_12p70, 0.080))

# transformation
dat <- dat %>% 
  mutate(log_il_4 = log(il_4),
         log_il_12p70 = log(il_12p70))

# CPM
mod1 <- orm(il_4 ~ age + sex + bmi + status,
            data = dat)
p <- length(mod1$coefficients) - (length(mod1$yunique) - 1)


# pdds ratio table
beta <- mod1$coefficients[length(mod1$yunique):length(mod1$coefficients)]
beta.se <- sqrt(diag(mod1$var)[-1])

### odds ratio
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

# BMI (5 unit)
tab_or[3, 'or'] <- exp(beta[3] * 5)
tab_or[3, 'or_lb'] <- exp((beta[3] - qnorm(0.975) * beta.se[3]) * 5)
tab_or[3, 'or_ub'] <- exp((beta[3] + qnorm(0.975) * beta.se[3]) * 5)
tab_or[3, 'p_ind'] <- p_ind_func(beta[3] * 30)

### Varying BMI ####
new.data <- matrix(NA, ncol = p, nrow=length(22:58))
colnames(new.data) <- names(mod1$coefficients[length(mod1$yunique):length(mod1$coefficients)])

new <- 22:58
new.data[,1] <- median(dat %>%  pull(age))
new.data[,2] <- 1 # male
new.data[,3] <- new
new.data[,4] <- 0 # group 1
new.data[,5] <- 0 # group 1
new.data[,6] <- 0 # group 1

q50.est <- quantile.orm(mod1, new.data = new.data, probs = 0.5, se=TRUE)
q90.est <- quantile.orm(mod1, new.data = new.data, probs = 0.9, se=TRUE)

c003.est <- cdf.orm(mod1, new.data, at.y=0.019, se = TRUE)
c005.est <- cdf.orm(mod1, new.data, at.y=0.05, se = TRUE)


### Median Regresion ### 
dat_med_dl <- dat %>% 
  mutate(il_4 = ifelse(il_4 < 0.019, 0.019, il_4),
         il_12p70 = ifelse(il_12p70 < 0.090, 0.090, il_12p70))

#model
mod_rq <- rq(il_4  ~ age + sex + bmi + status,
             data = dat_med_dl)
mod_rq %>% tidy() %>% kable(digits=3) %>% kable_styling(full_width = F)

new.data <- matrix(NA, ncol = p, nrow=length(22:58))
colnames(new.data) <- names(mod1$coefficients[length(mod1$yunique):length(mod1$coefficients)])

new <- 22:58
new.data[,1] <- median(dat %>%  pull(age))
new.data[,2] <- 1 # male
new.data[,3] <- new
new.data[,4] <- 0 # group 1
new.data[,5] <- 0 # group 1
new.data[,6] <- 0 # group 1

X <-  as.matrix(cbind(1, new.data))
pred <-drop(X %*% mod_rq$coef)
V <- summary(mod_rq, cov = TRUE)
sdpred <- sqrt(diag(X %*% V$cov %*% t(X)))
med_est <- pred
med_lower <- pred - qnorm(0.975) * sdpred
med_upper <- pred + qnorm(0.975) * sdpred


### Likelihood Approach (assume log-normal) ###
#### data for MLE (survival) ####
delta <- dat$il_4 >= 0.019
dat_mle <- dat
dat_mle$il_4 <- ifelse(dat_mle$il_4 < 0.019, 0.019, dat_mle$il_4)

mod_mle <- survreg(Surv(dat_mle$il_4, delta, type='left') ~ age + sex + bmi + status, 
                   data = dat_mle, dist = 'lognormal')
new <- 22:58
new.data <- data.frame(age = median(dat %>% pull(age)),
                       sex = 'Male',
                       bmi = new,
                       status = levels(dat$status)[1])

# median (censored quantile regression)
med_mle <- predict(mod_mle, newdata=new.data, 
                   type = "quantile", se=T, p=0.5)

med_est2 <- med_mle$fit
med_lower2 <- med_mle$fit - qnorm(0.975) * med_mle$se.fit
med_upper2 <- med_mle$fit + qnorm(0.975) * med_mle$se.fit




