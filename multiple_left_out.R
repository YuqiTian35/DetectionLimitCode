# Table S8: left-out variable 
library(tidyverse)
library(rms)
library(kableExtra)
library(multipleDL)
library(rlist)
library(xtable)
library(knitr)

# truth
beta.true <- 1

# generate the true data and the observed data
data_gen <- function(n_each, num_site, lower_dl = c(0.16, 0.30, 0.5), upper_dl = NULL){
  
  site_name <- 1:num_site # name of site
  # x dependent of site, z independent of site
  x1 <- rnorm(n_each, -0.5, 1)
  z1 <- rnorm(n_each, 0, 1)
  x2 <- rnorm(n_each, 0, 1)
  z2 <- rnorm(n_each, 0, 1)
  x3 <- rnorm(n_each, 0.5, 1) 
  z3 <- rnorm(n_each, 0, 1)
  
  # same transformation for each site
  y1 <- exp(1*x1 + 0.5 * z1 + rnorm(n_each, 0, 1))
  y2 <- exp(1*x2 + 0.5 * z2 + rnorm(n_each, 0, 1))
  y3 <- exp(1*x3 + 0.5 * z3 + rnorm(n_each, 0, 1))
  
  x <- c(x1, x2, x3)
  y <- c(y1, y2, y3)
  
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
n_each <- 50 # 300
num_site <- 3

# store results
est  <- se <- mse <-
  rep(NA, reps)


for(i in 1:reps){
  set.seed(i)
  data <- data_gen(n_each = n_each, num_site = num_site)
  data.obs <- data$dat_obs
  
  mod <-  multipleDL(y_obs ~ x, data = data.obs, link = 'probit', delta_lower = data.obs$dl_lower)
  
  est[i] <- mod$coef['x']
  se[i] <- sqrt(mod$var['x', 'x'])
  mse[i] <- (est[i] - beta.true)^2
} 
