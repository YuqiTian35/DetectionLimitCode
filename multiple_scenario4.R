library(tidyverse)
library(rms)
library(kableExtra)
library(multipleDL)
library(rlist)
library(xtable)
library(knitr)

# truth
set.seed(35)
beta.true <- 1
med.true <- c(exp(0), exp(1))
cdf.true <- c(mean(exp(rnorm(1e5, 0, 1) + 0) <= 1.5),
              mean(exp(rnorm(1e5, 0, 1) + 1) <= 1.5))

# generate the true data and the observed data
# generate the true data and the observed data
data_gen <- function(n_each, num_site, lower_dl = c(0.2, 0.3, -Inf), upper_dl = c(Inf, 4, 3.5)){
  
  site_name <- 1:num_site # name of site
  x <- rnorm(n_each * num_site, 0, 1) 
  e <- rnorm(n_each * num_site, 0, 1)
  y <- exp(x + e)
  
  # true data
  dat <- data.frame(site = rep(site_name, each = n_each),
                    y = y,
                    x = x)
  
  if(is.null(lower_dl)){
    lower_dl <- rep(-Inf, num_site)
  }
  if(is.null(upper_dl)){
    upper_dl <- rep(Inf, num_site)
  }
  
  dat$lower_dl <- rep(lower_dl, each = n_each)
  dat$upper_dl <- rep(upper_dl, each = n_each)
  
  # the observed data
  dat_obs <- dat
  dat_obs$y_obs <- ifelse(dat$y < dat$lower_dl, dat$lower_dl,
                          ifelse(dat$y > dat$upper_dl, dat$upper_dl, dat$y))
  # dl = indicator for observed value
  dat_obs$dl_lower <- unlist(lapply(1:length(site_name), function(i) {
    temp <- dat_obs %>% filter(site == site_name[i])
    return(ifelse(temp$y < lower_dl[i], 0, 1))
  }))
  
  dat_obs$dl_upper <- unlist(lapply(1:length(site_name), function(i) {
    temp <- dat_obs %>% filter(site == site_name[i])
    return(ifelse(temp$y > upper_dl[i], 0, 1))
  }))
  
  return(list(dat = dat, dat_obs = dat_obs))
}


reps <- 1e3
new.data <- data.frame(X = c(0,1))

n_each <- 50 # 300
num_site <- 3

# store results
beta.est  <- se.est <- 
  q50.0.est <- q50.0.lb <- q50.0.ub <- 
  q50.1.est <- q50.1.lb <- q50.1.ub <- 
  cdf.0.est <- cdf.0.se <- cdf.0.lb <- cdf.0.ub <- 
  cdf.1.est <- cdf.1.se <- cdf.1.lb <- cdf.1.ub <-
  rep(NA, reps)


for(i in 1:reps){
  set.seed(i)
  data <- data_gen(n_each = n_each, num_site = num_site)
  data.obs <- data$dat_obs
  
  mod <-  multipleDL(y_obs ~ x, data = data.obs, link = 'probit')
  # quantiles
  q.est <- quantile_dl(mod, new.data, probs=0.5)
  # cdf
  cdf.est <- cdf_dl(mod, new.data, at.y = 1.5)
  
  beta.est[i] <- mod$coef['x']
  se.est[i] <- sqrt(mod$var['x', 'x'])
  
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