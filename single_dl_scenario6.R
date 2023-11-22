# single DL - scenario 6
# lower DL at 0.0625 with a complex transformation

library(tidyverse)
library(multipleDL)

reps <- 1000

# new data for prediction
new.data <- data.frame(x=c(0,1))

# true value
set.seed(35)
e <- rnorm(1e5, 0, 1)
y0 <- exp(0 + e)
y0.weird <- ifelse((y0 < 2 & y0 >= 0.25), sqrt(y0),
                   ifelse(y0 < 0.25, y0^2, y0))
y1 <- exp(1 + e)
y1.weird <- ifelse((y1 < 2 & y1 >= 0.25), sqrt(y1),
                   ifelse(y1 < 0.25, y1^2, y1))
med.true <- c(median(y0.weird), median(y1.weird))
cdf.true <- c(mean(y0.weird <= 1.5), mean(y1.weird <= 1.5))


# add DL indicator based on detection limit specified
# allow for 3 types of censoring
data_one_dl <- function(dat, type=c("lower", "upper", "both")){
  if(type == "lower"){
    dat$dl <- ifelse(dat$y < lower_dl, 0, 1)
    dat$y <- ifelse(dat$dl == 0, lower_dl, dat$y)
  }else if(type == "upper"){
    dat$dl <- ifelse(dat$y > upper_dl, 0, 1)
    dat$y <- ifelse(dat$dl == 0, upper_dl, dat$y)
  }else{
    dat$lower_dl <- ifelse(dat$y < lower_dl, 0, 1)
    dat$upper_dl <- ifelse(dat$y > upper_dl, 0, 1)
    dat$y <- ifelse(dat$lower_dl == 0, lower_dl,
                    ifelse(dat$upper_dl == 0, upper_dl, dat$y))
  }
  return(dat)
}


# function to run simulation each iteration
func_res <- function(n){
  res <- sapply(1:reps, function(i){
    
    # data 
    set.seed(i*35)
    
    # data generation
    x <- rnorm(n, 0, 1) 
    e <- rnorm(n, 0, 1)
    y1 <- exp(x + e)
    y <- ifelse((y1 < 2 & y1 >= 0.25), sqrt(y1),
                ifelse(y1 < 0.25, y1^2, y1))
    dat <- data.frame(y = y,
                      x = x) 
    
    # detection limit 
    lower_dl <- 0.0625
    dat <- data_one_dl(dat, 'lower')
    
    # model
    mod <- multipleDL(y ~ x, data = dat, delta_lower = dat$dl, link = 'probit')
    
    ## beta
    beta.est <- mod$coef["x"]
    beta.se <- sqrt(mod$var["x","x"])
    beta.mse <- (beta.est - beta.true)^2
    
    ## median
    med.res <- quantile_dl(mod, new.data, 0.5)
    med.est <- med.res[["est"]]
    med.lb <- med.res[["lb"]]
    med.ub <- med.res[["ub"]]
    med.mse <- (med.est - med.true)^2
    
    ## cdf
    cdf.res <- cdf_dl(mod, new.data, at.y=1.5)
    cdf.est <- cdf.res[["est"]] 
    cdf.se <- cdf.res[["se"]]
    cdf.lb <- cdf.res[["lb"]]
    cdf.ub <- cdf.res[["ub"]]
    cdf.mse <- (cdf.est - cdf.true)^2
    
    return(c(beta.est,
             beta.se,
             beta.mse,
             med.est,
             med.mse,
             med.lb,
             med.ub,
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
                      med.est.0,
                      med.est.1,
                      med.mse.0,
                      med.mse.1,
                      med.lb.0,
                      med.lb.1,
                      med.ub.0,
                      med.ub.1,
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
    
    # bias of med.0 est
    mean(result$med.est.0) - med.true[1],
    # bias of med.1 est
    mean(result$med.est.1) - med.true[2],
    # percent bias of med.0 est
    (mean(result$med.est.0) - med.true[1]) / med.true[1] * 100,
    # percent bias of med.1 est
    (mean(result$med.est.1) - med.true[2]) / med.true[2] * 100, 
    # mse of med.0
    sqrt(mean(result$med.mse.0)),
    # mse of med.1
    sqrt(mean(result$med.mse.1)),
    # ci coverage of med.0
    mean((result$med.lb.0 <= med.true[1]) & (result$med.ub.0 >= med.true[1])),
    # ci coverage of med.1
    mean((result$med.lb.1 <= med.true[2]) & (result$med.ub.1 >= med.true[2])),
    
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
  
  colnames(tab) <- c("n=100","n=500")
  rownames(tab) <- c("bias of beta",  "percent bias of beta", 
                     "empirical estimates of sd(beta)", "mean of estimated sd(beta)",
                     "rmse of beta", "ci of beta",
                     
                     "bias of med(Y|X=0)", "bias of med(Y|X=1)",
                     "percent bias of med(Y|X=0)", " percentbias of med(Y|X=1)",
                     "rmse of med(Y|X=0)", "rmse of med(Y|X=1)",
                     "ci of med(Y|X=0)", "ci of med(Y|X=1)",
                     
                     "bias of F(1.5|X=0)", "bias of F(1.5|X=1)",
                     "percent bias of F(1.5|X=0)", "percent bias of F(1.5|X=1)",
                     "empirical estimates of sd(F(1.5|X=0))", "empirical estimates of sd(F(1.5|X=1))",
                     "mean of estimated sd(F(1.5|X=0))", "mean of estimated sd(F(1.5|X=1))",
                     "rmse of F(1.5|X=0)", "rmse of F(1.5|X=1)",
                     "ci of F(1.5|X=0)", "ci of F(1.5|X=1)")
  
  return(tab)
}

# results for n=100 and n=500
result_100 <- func_res(n=100)
result_500 <- func_res(n=500)
print(func_tab(result_100, result_500))
