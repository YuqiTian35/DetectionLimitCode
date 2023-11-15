# simulation - multiple DLs scenario 2 (left-out variable) [Table S8]
# lower DLs at 0.16, 0.30, 0.50

library(tidyverse)
library(multipleDL)
library(rlist)

# truth
set.seed(35)
beta.true <- 1
q50.true <- c(exp(0), exp(1))
cdf.true <- c(mean(exp(rnorm(1e5, 0, 1) + 0) <= 1.5),
              mean(exp(rnorm(1e5, 0, 1) + 1) <= 1.5))

# use 10 as an example. Please set to 1000
reps <- 10
# new data for prediction
new.data <- data.frame(X = c(0,1))

# generate the true data and the observed data
## n_each: sample size at each site
## num_site: the number of sites
## lower_dl: lower DLs at each site
## upper_dl: upper DLs at each site
data_gen <- function(n_each, num_site, lower_dl = NULL, upper_dl = NULL){
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
  
  # if no DL specified, use -Inf/Inf
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
  # dl: indicator for observed value
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


# function to run simulation each iteration
func_res <- function(n_each, num_site=3){
  res <- sapply(1:reps, function(i){
    
    # data 
    set.seed(i*35)
    
    # data generation
    data <- data_gen(n_each = n_each, num_site = num_site, lower_dl = c(0.16, 0.30, 0.5))
    
    # model
    mod <-  multipleDL(y_obs ~ x, data = data, link = 'probit')
    
    ## beta
    beta.est <- mod$coef["x"]
    beta.se <- sqrt(mod$var["x","x"])
    beta.mse <- (beta.est - beta.true)^2
    
    return(c(beta.est,
             beta.se,
             beta.mse)) 
  }) %>% t %>% as.data.frame
  
  # names for results
  colnames(res) <- Cs(beta.est,
                      beta.se,
                      beta.mse)
  return(res)
}

# function: summary of result
sum_tab <- function(result){
  return(c(
    # bias of beta
    mean(result$beta.est) - beta.true,
    # percent bias
    (mean(result$beta.est) - beta.true)/beta.true*100,
    # sd of beta est
    sd(result$beta.est),
    # mean of sd(beta)
    mean(result$beta.se),
    # mse of beta
    sqrt(mean(result$beta.mse)),
    # ci coverage of beta
    mean((result$beta.est - qnorm(0.975) * result$beta.se <= beta.true) &
           (result$beta.est + qnorm(0.975) * result$beta.se >= beta.true))
  ))
}

# function: create a table combining results
func_tab <- function(result1, result2){
  tab <- as.data.frame(cbind(sum_tab(result1), 
                             sum_tab(result2)))
  
  colnames(tab) <- c("n=50*3","n=300*3")
  rownames(tab) <- c("bias of beta",  "percent bias of beta", 
                     "empirical estimates of sd(beta)", "mean of estimated sd(beta)",
                     "rmse of beta", "ci of beta")
  
  return(tab)
}

# results for n=50*3 and n=300*3
result_150 <- func_res(n_each=50, num_site=3)
result_900 <- func_res(n_each=300, num_site=3)
print(func_tab(result_150, result_900))