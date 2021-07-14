# single DL - scenario 1

library(tidyverse)
library(rms)
library(parallel)

source('conditional_func')

reps <- 10000

new.data <- data.frame(x=c(0,1))

# true value
set.seed(35)
beta.true <- 1
med.true <- c(exp(0), exp(1))
cdf.true <- c(mean(exp(rnorm(1e6, 0, 1) + 0) <= 1.5),
              mean(exp(rnorm(1e6, 0, 1) + 1) <= 1.5))


func_res <- function(n){
  res <- parSapply(cl, 1:reps, function(i){
    
    # data 
    set.seed(i*35)
   
    # data generation
    x <- rnorm(n, 0, 1) 
    e <- rnorm(n, 0, 1)
    y <- exp(x + e)
    dat <- data.frame(y = y,
                      x = x)
    
    mod <- orm(y ~ x,
               data = dat,
               family = probit)
    
    ## beta
    beta.est <- mod$coefficients["x"]
    beta.se <- sqrt(mod$var["x","x"])
    beta.mse <- (beta.est - beta.true)^2
    
    ## median
    med.res <- quantile.orm(mod, new.data, se=T)
    med.est <- med.res[["quantile"]]
    med.lb <- med.res[["lb"]]
    med.ub <- med.res[["ub"]]
    # if(type == "upper" | type=="both"){
    #   med.est <- sapply(med.est, function(x) ifelse(is.finite(x), x, upper_dl)) # for upper dl
    #   med.ub <- sapply(med.ub, function(x) ifelse(is.finite(x), x, upper_dl)) # for upper dl
    # }
    med.mse <- (med.est - med.true)^2
    
    ## cdf
    cdf.res <- cdf.orm(mod, new.data, at.y=1.5, se=T)
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



cl <- makeCluster(detectCores())
clusterEvalQ(cl, library(rms, quietly = TRUE))
clusterExport(cl,varlist=ls())

result_100 <- func_res(n=100)
result_500 <- func_res(n=500)