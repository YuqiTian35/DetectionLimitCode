# misspecification - logit link

library(tidyverse)
library(rms)
library(parallel)

source('conditional_func')

reps <- 10000

new.data <- data.frame(x=c(0,1))

# true value
beta.true <- 1
med.true <- c(exp(0), exp(1))
cdf.true <- c(mean(exp(rnorm(1e5, 0, 1) + 0) <= 1.5),
              mean(exp(rnorm(1e5, 0, 1) + 1) <= 1.5))

lower_dl <- 0.25

data_one_dl <- function(n, type=c("lower", "upper", "both")){
  # data generation
  x <- rnorm(n, 0, 1)
  e <- rnorm(n, 0, 1)
  y <- exp(x + e)
  dat <- data.frame(y = y,
                    x = x)
  # recoded data
  dat_dl <- dat
  
  if(type == "lower"){
    dat_dl$y <- ifelse(dat_dl$y < lower_dl, lower_dl - 1e-5, dat_dl$y)
  }else if(type == "upper"){
    dat_dl$y <- ifelse(dat_dl$y > upper_dl, upper_dl + 1e-5, dat_dl$y)
  }else{
    dat_dl$y <- ifelse(dat_dl$y < lower_dl, lower_dl - 1e-5,
                       ifelse(dat_dl$y > upper_dl, upper_dl + 1e-5, dat_dl$y))
  }
  return(list(dat = dat, dat_dl = dat_dl))
}


func_res <- function(n, type){
  for(i in 1:reps){
    
    # data 
    set.seed(i+8000)
    
    data <- data_one_dl(n, type)
    dat <- data$dat # true data
    dat_dl <- data$dat_dl # censored data
    order.y <- unique(dat_dl$y[order(dat_dl$y)])
    
    # model
    tryCatch({
      mod <- orm(y ~ x,
                 data = dat_dl,
                 family = cloglog)
      
      
      ## beta
      
      beta.est[i] <- mod$coefficients["x"]
      beta.se[i] <- sqrt(mod$var["x","x"])
      beta.mse[i] <- (beta.est - beta.true)^2
      
      ## median
      med.res <- quantile.orm(mod, new.data, se=T)
      med.est.0[i] <- med.res[["quantile"]][1]
      med.est.1[i] <- med.res[["quantile"]][2]
      med.lb.0[i] <- med.res[["lb"]][1]
      med.lb.1[i] <- med.res[["lb"]][2]
      med.ub.0[i] <- med.res[["ub"]][1]
      med.ub.1[i] <- med.res[["ub"]][2]
      med.mse.0[i] <- (med.est.0 - med.true[1])^2
      med.mse.1[i] <- (med.est.1 - med.true[2])^2
      
      ## cdf
      cdf.res <- cdf.orm(mod, new.data, at.y=1.5, se=T)
      cdf.est.0[i] <- cdf.res[["est"]][1]
      cdf.est.1[i] <- cdf.res[["est"]][2]
      cdf.lb.0[i] <- cdf.res[["lb"]][1]
      cdf.lb.1[i] <- cdf.res[["lb"]][2]
      cdf.ub.0[i] <- cdf.res[["ub"]][1]
      cdf.ub.1[i] <- cdf.res[["ub"]][2]
      cdf.mse.0[i] <- (cdf.est.0 - cdf.true[1])^2
      cdf.mse.1[i] <- (cdf.est.1 - cdf.true[2])^2
    },
    error = function(e){
      print(i)
      beta.est <- beta.se <- beta.mse <-
        med.est.0 <- med.est.1 <- med.mse.0 <- med.mse.1 <- med.lb.0 <- med.lb.1 <- med.ub.0 <- med.ub.1 <-
        cdf.est.0 <- cdf.est.1 <- cdf.mse.0 <- cdf.mse.1 <- cdf.lb.0 <- cdf.lb.1 <- cdf.ub.0 <- cdf.ub.1 <- NA
    })
  }
  
  res <- cbind(beta.est,
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
               cdf.mse.0,
               cdf.mse.1,
               cdf.lb.0,
               cdf.lb.1,
               cdf.ub.0,
               cdf.ub.1) %>% 
    as.data.frame() %>% 
    filter(!is.na(beta.est))
  
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

result_100 <- func_res(n=100, type="lower")
result_500 <- func_res(n=500, type="lower")
