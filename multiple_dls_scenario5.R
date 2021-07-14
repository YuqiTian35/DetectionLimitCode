# multiple dls - scenario 1

library(tidyverse)
library(rms)
library(multipleDL)

reps <- 1e4
new.data <- data.frame(x = c(0,1))

# true value
beta.true <- 1

set.seed(35)
beta.true <- 1
med.true <- c(exp(0), exp(1))
cdf.true <- c(mean(exp(rnorm(1e5, 0, 1) + 0) <= 1.5),
              mean(exp(rnorm(1e5, 0, 1) + 1) <= 1.5))

# generate the true data and the observed data
data_gen <- function(n, num_site, detect_lim_lower, detect_lim_upper){
  
  site_name <- 1:num_site # name of site
  x <- rnorm(sum(n), 0, 1) 
  e <- rnorm(sum(n), 0, 1)
  y <- exp(x + e)
  
  # true data
  dat <- data.frame(site = rep(site_name, n),
                    y = y,
                    x = x)
  
  
  # the observed data
  dat_obs <- dat
  
  # dl = indicator for observed value
  dat_obs$dl_lower <- unlist(lapply(1:length(site_name), function(i) {
    temp <- dat_obs %>% filter(site == site_name[i])
    return(ifelse(temp$y < detect_lim_lower[i], 0, 1))
  }))
  
  dat_obs$dl_upper <- unlist(lapply(1:length(site_name), function(i) {
    temp <- dat_obs %>% filter(site == site_name[i])
    return(ifelse(temp$y > detect_lim_upper[i], 0, 1))
  }))
  
  # observed value
  dat_obs$y <- unlist(lapply(1:length(site_name), function(i) {
    temp <- dat_obs %>% filter(site == site_name[i])
    temp$y <- ifelse(temp$dl_lower == 0, detect_lim_lower[i], 
                     ifelse(temp$dl_upper == 0, detect_lim_upper[i], temp$y))
    return(temp$y)
  }))
  
  
  return(dat_obs)
}


n_each <- 50
# n_each <- 300
num_site <- 3
detect_lim_lower <- c(0.2, 0.3, -1e6)
detect_lim_upper <- c(4, 3.5, 1e6)


# store results
beta.est <- se.est <- 
  q50.0.est <- q50.0.lb <- q50.0.ub <- 
  q50.1.est <- q50.1.lb <- q50.1.ub <- 
  cdf.0.est <- cdf.0.se <- cdf.0.lb <- cdf.0.ub <- 
  cdf.1.est <- cdf.1.se <- cdf.1.lb <- cdf.1.ub <-
  rep(NA, reps)



for(i in 1:reps){
  set.seed(i)
  data <- data_gen(n = rep(n_each, num_site), num_site = num_site, 
                   detect_lim_lower = detect_lim_lower, detect_lim_upper = detect_lim_upper)
  
  mod <- multipleDL(formula = y ~ x, 
                    data = data, 
                    delta_lower = data$dl_lower,
                    delta_upper = data$dl_upper,
                    link='probit')
  
  coef <- mod$coef
  se <- mod$var %>% diag %>% sqrt
  q.est <- quantile_dl(mod, new.data, probs = 0.5)
  cdf.est <- cdf_dl(mod, new.data, at.y=1.5)
  
  q50.0.est[i] <- q.est$est[1,1]
  q50.0.lb[i] <- q.est$lb[1,1]
  q50.0.ub[i] <- q.est$ub[1,1]
  
  q50.1.est[i] <- q.est$est[2,1]
  q50.1.lb[i] <- q.est$lb[2,1]
  q50.1.ub[i] <- q.est$ub[2,1]
  
  cdf.0.est[i] <- cdf.est$est[1,1]
  cdf.0.se[i] <- cdf.est$se[1,1]
  cdf.0.lb[i]  <- cdf.est$lb[1,1]
  cdf.0.ub[i]  <- cdf.est$ub[1,1]
  
  cdf.1.est[i] <- cdf.est$est[2,1]
  cdf.1.se[i] <- cdf.est$se[2,1]
  cdf.1.lb[i]  <- cdf.est$lb[2,1]
  cdf.1.ub[i]  <- cdf.est$ub[2,1]
} 
```