# simulation - multiple DLs scenario 5 [Table 2]
# lower DLs at 0.4, 0.1, 2.5

library(tidyverse)
library(multipleDL)
library(rlist)

# truth
set.seed(35)
beta.true <- 1
q50.true <- c(quantile(exp(rnorm(1e5, 0, 1) + 0), 0.9), quantile(exp(rnorm(1e5, 0, 1) + 1), 0.9))
cdf.true <- c(mean(exp(rnorm(1e5, 0, 1) + 0) <= 3),
              mean(exp(rnorm(1e5, 0, 1) + 1) <= 3))


# number of replications
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


# function to run simulation each iteration
func_res <- function(n_each, num_site=3){
  res <- sapply(1:reps, function(i){
    
    # data 
    set.seed(i*35)
    
    # data generation
    data <- data_gen(n_each = n_each, num_site = num_site, lower_dl = c(0.4, 1.0, 2.5))
    
    # model
    mod <-  multipleDL(y_obs ~ x, data = data, link = 'probit')
    
    ## beta
    beta.est <- mod$coef["x"]
    beta.se <- sqrt(mod$var["x","x"])
    beta.mse <- (beta.est - beta.true)^2
    
    ## conditional quantile
    q50.res <- quantile_dl(mod, new.data, probs=0.9)
    q50.est <- q50.res[["est"]]
    q50.lb <- q50.res[["lb"]]
    q50.ub <- q50.res[["ub"]]
    q50.mse <- (q50.est - q50.true)^2
    
    ## conditional cdf
    cdf.res <- cdf_dl(mod, new.data, at.y = 3)
    cdf.est <- cdf.res[["est"]] 
    cdf.se <- cdf.res[["se"]]
    cdf.lb <- cdf.res[["lb"]]
    cdf.ub <- cdf.res[["ub"]]
    cdf.mse <- (cdf.est - cdf.true)^2
    
    return(c(beta.est,
             beta.se,
             beta.mse,
             q50.est,
             q50.mse,
             q50.lb,
             q50.ub,
             cdf.est,
             cdf.se,
             cdf.mse,
             cdf.lb,
             cdf.ub)) 
  }) %>% t %>% as.data.frame
  
  # names for results
  colnames(res) <- Cs(beta.est,
                      beta.se,
                      beta.mse,
                      q50.est.0,
                      q50.est.1,
                      q50.mse.0,
                      q50.mse.1,
                      q50.lb.0,
                      q50.lb.1,
                      q50.ub.0,
                      q50.ub.1,
                      cdf.est.0,
                      cdf.est.1,
                      cdf.se.0,
                      cdf.se.1,
                      cdf.mse.0,
                      cdf.mse.1,
                      cdf.lb.0,
                      cdf.lb.1,
                      cdf.ub.0,
                      cdf.ub.1)
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
           (result$beta.est + qnorm(0.975) * result$beta.se >= beta.true)),
    
    # bias of q50.0 est
    mean(result$q50.est.0) - q50.true[1],
    # bias of q50.1 est
    mean(result$q50.est.1) - q50.true[2],
    # percent bias of q50.0 est
    (mean(result$q50.est.0) - q50.true[1]) / q50.true[1] * 100,
    # percent bias of q50.1 est
    (mean(result$q50.est.1) - q50.true[2]) / q50.true[2] * 100, 
    # mse of q50.0
    sqrt(mean(result$q50.mse.0)),
    # mse of q50.1
    sqrt(mean(result$q50.mse.1)),
    # ci coverage of q50.0
    mean((result$q50.lb.0 <= q50.true[1]) & (result$q50.ub.0 >= q50.true[1])),
    # ci coverage of q50.1
    mean((result$q50.lb.1 <= q50.true[2]) & (result$q50.ub.1 >= q50.true[2])),
    
    # bias of cdf.0 est
    mean(result$cdf.est.0) - cdf.true[1],
    # bias of cdf.1 est
    mean(result$cdf.est.1) - cdf.true[2],
    # percent bias of cdf.0 est
    (mean(result$cdf.est.0) - cdf.true[1]) / cdf.true[1] * 100,
    # percent bias of cdf.1 est
    (mean(result$cdf.est.1) - cdf.true[2]) / cdf.true[2] * 100,
    # sd of cdf.0 est
    sd(result$cdf.est.0),
    # sd of cdf.1 est
    sd(result$cdf.est.1),
    # mean of sd(cdf.est)
    mean(result$cdf.se.0),
    # mean of sd(cdf.est)
    mean(result$cdf.se.1),
    # mse of cdf.0
    sqrt(mean(result$cdf.mse.0)),
    # mse of cdf.1
    sqrt(mean(result$cdf.mse.1)),
    # ci coverage of cdf.0
    mean((result$cdf.lb.0 <= cdf.true[1]) & (result$cdf.ub.0 >= cdf.true[1])),
    # ci coverage of cdf.1
    mean((result$cdf.lb.1 <= cdf.true[2]) & (result$cdf.ub.1 >= cdf.true[2]))
  ))
}

# function: create a table combining results
func_tab <- function(result1, result2){
  tab <- as.data.frame(cbind(sum_tab(result1), 
                             sum_tab(result2)))
  
  colnames(tab) <- c("n=50*3","n=300*3")
  rownames(tab) <- c("bias of beta",  "percent bias of beta", 
                     "empirical estimates of sd(beta)", "mean of estimated sd(beta)",
                     "rmse of beta", "ci of beta",
                     
                     "bias of Q(0.9|X=0)", "bias of Q(0.9|X=1)",
                     "percent bias of Q(0.9|X=0)", "percentbias of Q(0.9|X=1)",
                     "rmse of Q(0.9|X=0)", "rmse of Q(0.9|X=1)",
                     "ci of Q(0.9|X=0)", "ci of Q(0.9|X=1)",
                     
                     "bias of F(3|X=0)", "bias of F(3|X=1)",
                     "percent bias of F(3|X=0)", "percent bias of F(3|X=1)",
                     "empirical estimates of sd(F(3|X=0))", "empirical estimates of sd(F(3|X=1))",
                     "mean of estimated sd(F(3|X=0))", "mean of estimated sd(F(3|X=1))",
                     "rmse of F(3|X=0)", "rmse of F(3|X=1)",
                     "ci of F(3|X=0)", "ci of F(3|X=1)")
  
  return(tab)
}

# results for n=50*3 and n=300*3
result_150 <- func_res(n_each=50, num_site=3)
result_900 <- func_res(n_each=300, num_site=3)
print(func_tab(result_150, result_900))