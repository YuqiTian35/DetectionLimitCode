# Figure S5 - p-values
library(tidyverse)
library(multipleDL)


# truth
set.seed(35)
beta.true <- 1
med.true <- c(exp(0), exp(1))
cdf.true <- c(mean(exp(rnorm(1e5, 0, 1) + 0) <= 1.5),
              mean(exp(rnorm(1e5, 0, 1) + 1) <= 1.5))

# detection limits
lower_dl <- c(0.16, 0.30, 0.5)

# generate the true data and the observed data
data_gen <- function(n_each, num_site, lower_dl = NULL, upper_dl = NULL){
  
  site_name <- 1:num_site # name of site
  # x dependent of site
  x1 <- rnorm(n_each, -0.5, 1)
  x2 <- rnorm(n_each, 0, 1)
  x3 <- rnorm(n_each, 0.5, 1) 
  
  # same transformation for each site
  y1 <- exp(0*x1 + rnorm(n_each, 0, 1))
  y2 <- exp(0*x2 + rnorm(n_each, 0, 1))
  y3 <- exp(0*x3 + rnorm(n_each, 0, 1))
  # combine data from all 3 sites
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
  dat$y_obs <- ifelse(dat$y < dat$lower_dl, dat$lower_dl,
                          ifelse(dat$y > dat$upper_dl, dat$upper_dl, dat$y))
  # dl = indicator for observed value
  dat$dl_lower <- unlist(lapply(1:length(site_name), function(i) {
    temp <- dat %>% filter(site == site_name[i])
    return(ifelse(temp$y < lower_dl[i], 0, 1))
  }))
  
  dat$dl_upper <- unlist(lapply(1:length(site_name), function(i) {
    temp <- dat %>% filter(site == site_name[i])
    return(ifelse(temp$y > upper_dl[i], 0, 1))
  }))
  
  return(dat)
}

# number of replications
reps <- 1000

# store results
pval_150 <- pval_900 <- rep(NA, reps)

for(i in 1:reps){
  set.seed(i)
  # n=50*3
  dat <- data_gen(n_each = 50, num_site = 3, lower_dl = lower_dl)
  mod <-  multipleDL(y_obs ~ x, data = dat, link = 'probit', delta_lower = dat$dl_lower)
  est <- mod$coef['x']
  se  <- sqrt(mod$var['x', 'x'])
  pval_150[i] <-  2 * (1 - pnorm(abs(est / se)))
  
  # n=300*3
  dat <- data_gen(n_each = 300, num_site = 3)
  mod <-  multipleDL(y_obs ~ x, data = dat, link = 'probit', delta_lower = dat$dl_lower)
  est <- mod$coef['x']
  se  <- sqrt(mod$var['x', 'x'])
  pval_900[i] <-  2 * (1 - pnorm(abs(est / se)))
} 

# uncomment the code to plot histograms
hist(pval_150, breaks = 50, xlab = 'p-value', main = 'n = 50*3')
hist(pval_900, breaks = 50, xlab = 'p-value', main = 'n = 300*3')

print(hist(pval_150, breaks = 50, xlab = 'p-value', main = 'n = 50*3'))
print(hist(pval_900, breaks = 50, xlab = 'p-value', main = 'n = 300*3'))
