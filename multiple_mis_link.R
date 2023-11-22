# simulation - multiple DLs scenario 2 (misspecification of link functions) [Table S7]
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

reps <- 1000
# new data for prediction
new.data <- data.frame(X = c(0,1))

# generate the true data and the observed data
## n_each: sample size at each site
## num_site: the number of sites
## lower_dl: lower DLs at each site
## upper_dl: upper DLs at each site
data_gen <- function(n_each, num_site, lower_dl = NULL, upper_dl = NULL){
  site_name <- 1:num_site # name of site
  # x dependent of site
  x1 <- rnorm(n_each, -0.5, 1)
  x2 <- rnorm(n_each, 0, 1)
  x3 <- rnorm(n_each, 0.5, 1) 
  
  # same transformation for each site
  y1 <- exp(x1 + rnorm(n_each, 0, 1))
  y2 <- exp(x2 + rnorm(n_each, 0, 1))
  y3 <- exp(x3 + rnorm(n_each, 0, 1))
  
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


# function to run simulation each iteration - logit link
func_res_logit <- function(n_each, num_site=3){
  res <- sapply(1:reps, function(i){
    
    # data 
    set.seed(i*35)
    
    # data generation
    data <- data_gen(n_each = n_each, num_site = num_site, lower_dl = c(0.16, 0.30, 0.5))
    
    # model
    mod <-  multipleDL(y_obs ~ x, data = data, link = 'logit')
    
    ## beta
    beta.est <- mod$coef["x"] * sqrt(3 / pi^2) # scale back 
    beta.se <- sqrt(mod$var["x","x"]) * sqrt(3 / pi^2)  # scale back
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

# function to run simulation each iteration - loglog link
func_res_loglog <- function(n_each, num_site=3){
  res <- sapply(1:reps, function(i){
    
    # data 
    set.seed(i*35)
    
    # data generation
    data <- data_gen(n_each = n_each, num_site = num_site, lower_dl = c(0.16, 0.30, 0.5))
    
    # model
    mod <-  multipleDL(y_obs ~ x, data = data, link = 'loglog')
    
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

# results for logit and loglog link functions
result_logit_150 <- func_res_logit(n_each=50, num_site=3)
result_logit_900 <- func_res_logit(n_each=300, num_site=3)
print("Logit link")
print(func_tab(result_logit_150, result_logit_900))

result_loglog_150 <- func_res_loglog(n_each=50, num_site=3)
result_loglog_900 <- func_res_loglog(n_each=300, num_site=3)
print("Loglog link")
print(func_tab(result_loglog_150, result_loglog_900))