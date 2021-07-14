# single DL - comparison

library(tidyverse)
library(rms)
library(parallel)

source('conditional_func')

reps <- 1e4
lower_dl <- 0.25
n <- 10000


# multiple imputation
dl_mi <- function(b=10, data=dat){
  dat_fit <- data %>% 
    mutate(left_val = ifelse(y < lower_dl, NA, y),
           right_val = y)
  
  ## get lognormal distribution parameters
  lnorm_par <- fitdistcens(dat_fit %>% dplyr::select(left_val, right_val), "lnorm")
  
  mi_res <- sapply(1:b, function(i){
    # data for this iteration
    data_b <- data
    
    ## get lognormal distribution parameters by sample
    lnorm_par_b <- rmvnorm(1, mean = lnorm_par$estimate, sigma = lnorm_par$vcov)
    
    # number of observations need to be imputed
    n_imp <- sum(data_b$y < lower_dl)
    # imputed values
    val_imp <- rep(NA, n_imp) # imputed values
    temp <- 1 
    while(temp <= n_imp){
      val <- rlnorm(1, lnorm_par_b[1], lnorm_par_b[2]) # random 
      if(val < lower_dl){ # below DL
        val_imp[temp] <- val
        temp <- temp + 1
      }else{
        next
      }
    }
    
    # insert imputed values
    data_b$y[data_b$y < lower_dl] <- val_imp
    
    
    mod_mi <- lm(log(y) ~ x, # correct transformation
                 data = data_b)
    
    mod_rq <- rq(log(y) ~ x,
                 data = data_b)
    med_mi <- predict(mod_rq, new.data, interval = 'confidence')
    
    return(c(mod_mi$coefficients["x"],
             vcov(mod_mi)['x','x'],
             med_mi[,1], med_mi[,2], med_mi[,3]))
    
  }) 
  
  
  b_est_mi <- mean(mi_res[1,])
  b_se_mi <- sqrt(mean(mi_res[2,]) + (1 + 1 / b) * var(mi_res[1,]))
  b_ci_mi <- (b_est_mi - qnorm(0.975) * b_se_mi <= beta_true) &
    (b_est_mi + qnorm(0.975) * b_se_mi >= beta_true)
  
  med_est_mi_0 <- exp(mean(mi_res[3,]))
  med_ci_mi_0 <- (exp(mean(mi_res[5,])) <= med_0_true) &
    (exp(mean(mi_res[7,])) >= med_0_true)
  
  med_est_mi_1 <- exp(mean(mi_res[4,]))
  med_ci_mi_1 <- (exp(mean(mi_res[6,])) <= med_1_true) &
    (exp(mean(mi_res[8,])) >= med_1_true)
  
  return(c(b_est_mi, b_ci_mi, 
           med_est_mi_0, med_ci_mi_0,
           med_est_mi_1, med_ci_mi_1))
  
}


#### correctly specified ####

func_res <- function(n=1000){
  res <- parSapply(cl, 1:reps, function(i){
    
    # data 
    set.seed(i)
    
    # data generation
    x <- rnorm(n, 0, 1) 
    e <- rnorm(n, 0, 1)
    y <- exp(x + e)
    dat <- data.frame(y = y,
                      x = x)
    
    
    ################## data ##################
    ##### data for CPM ####
    dat_cpm <- dat
    dat_cpm$y <- ifelse(dat_cpm$y < lower_dl, lower_dl - 1e-5, dat_cpm$y)
    
    ##### data for single imputation - sqrt(2) ####
    dat_imp <- dat
    dat_imp$y <- ifelse(dat_imp$y < lower_dl, lower_dl/sqrt(2), dat_imp$y)
    
    ##### data for single imputation - 2 ####
    dat_imp2 <- dat
    dat_imp2$y <- ifelse(dat_imp2$y < lower_dl, lower_dl/2, dat_imp2$y)
    
    
    #### data for MLE (survival) ####
    delta <- dat$y > lower_dl
    dat_mle <- dat
    dat_mle$y <- ifelse(dat_mle$y < lower_dl, lower_dl, dat_mle$y)
    
    ################## models ##################
    ## CPM ##
    mod_cpm <- orm(y ~ x,
                   data = dat_cpm,
                   family = probit)
    # beta
    b_est_cpm <- mod_cpm$coefficients["x"]
    b_se_cpm <- sqrt(mod_cpm$var["x","x"])
    b_ci_cpm <- (b_est_cpm - qnorm(0.975) * b_se_cpm <= beta_true) & 
      (b_est_cpm + qnorm(0.975) * b_se_cpm >= beta_true)
    b_mse_cpm <- (b_est_cpm - beta_true)^2
    # median 0
    med_cpm <- quantile.orm(mod_cpm, new.data, 0.5)
    med_est_cpm_0 <- med_cpm$quantile[1]
    med_ci_cpm_0 <- (med_cpm$lb[1] <= med_0_true) & 
      (med_cpm$ub[1] >= med_0_true)
    med_mse_cpm_0 <- (med_est_cpm_0 - med_0_true)^2
    # median 1
    med_est_cpm_1 <- med_cpm$quantile[2]
    med_ci_cpm_1 <- (med_cpm$lb[2] <= med_1_true) & 
      (med_cpm$ub[2] >= med_1_true)
    med_mse_cpm_1 <- (med_est_cpm_1 - med_1_true)^2
    
    ## single imputation - sqrt(2) ## 
    mod_imp <- lm(log(y) ~ x, # correct transformation
                  data = dat_imp)
    b_est_imp <- mod_imp$coefficients["x"]
    b_se_imp <- vcov(mod_imp)['x','x'] %>% sqrt
    b_ci_imp <- (b_est_imp - qnorm(0.975) * b_se_imp <= beta_true) & 
      (b_est_imp + qnorm(0.975) * b_se_imp >= beta_true)
    b_mse_imp <- (b_est_imp - beta_true)^2
    # medan
    mod_rq <- rq(log(y) ~ x,
                 data = dat_imp)
    med_imp <- predict(mod_rq, new.data, interval = 'confidence')
    med_est_imp_0 <- exp(med_imp[1,1])
    med_ci_imp_0 <- (exp(med_imp[1,2]) <= med_0_true) & 
      (exp(med_imp[1,3]) >= med_0_true)
    med_mse_imp_0 <- (med_est_imp_0 - med_0_true)^2
    med_est_imp_1 <- exp(med_imp[2,1])
    med_ci_imp_1 <- (exp(med_imp[2,2]) <= med_1_true) & 
      (exp(med_imp[2,3]) >= med_1_true)
    med_mse_imp_1 <- (med_est_imp_1 - med_1_true)^2
    
    ## single imputation - 2 ## 
    mod_imp2 <- lm(log(y) ~ x, # correct transformation
                   data = dat_imp2)
    b_est_imp2 <- mod_imp2$coefficients["x"]
    b_se_imp2 <- vcov(mod_imp2)['x','x'] %>% sqrt
    b_ci_imp2 <- (b_est_imp2 - qnorm(0.975) * b_se_imp2 <= beta_true) & 
      (b_est_imp2 + qnorm(0.975) * b_se_imp2 >= beta_true)
    b_mse_imp2 <- (b_est_imp2 - beta_true)^2
    # medan
    mod_rq <- rq(log(y) ~ x,
                 data = dat_imp2)
    med_imp2 <- predict(mod_rq, new.data, interval = 'confidence')
    med_est_imp2_0 <- exp(med_imp2[1,1])
    med_ci_imp2_0 <- (exp(med_imp2[1,2]) <= med_0_true) & 
      (exp(med_imp2[1,3]) >= med_0_true)
    med_mse_imp2_0 <- (med_est_imp2_0 - med_0_true)^2
    med_est_imp2_1 <- exp(med_imp2[2,1])
    med_ci_imp2_1 <- (exp(med_imp2[2,2]) <= med_1_true) & 
      (exp(med_imp2[2,3]) >= med_1_true)
    med_mse_imp2_1 <- (med_est_imp2_1 - med_1_true)^2
    
    ## MLE (lognormal) ##
    mod_mle <- survreg(Surv(dat_mle$y, delta, type='left') ~ x, 
                       dist = 'lognormal')
    b_est_mle <- mod_mle$coefficients["x"]
    b_se_mle <- vcov(mod_mle)['x','x'] %>% sqrt
    b_ci_mle <- (b_est_mle - qnorm(0.975) * b_se_mle <= beta_true) & 
      (b_est_mle + qnorm(0.975) * b_se_mle >= beta_true)
    b_mse_mle <- (b_est_mle - beta_true)^2
    # median (censored quantile regression)
    med_mle <- predict(mod_mle, newdata=new.data, 
                       type = "quantile", se=T, p=0.5)
    med_est_mle_0 <- med_mle$fit[1]
    med_ci_mle_0 <- (med_mle$fit[1] - qnorm(0.975) * med_mle$se.fit[1] <= med_0_true) & (med_mle$fit[1] + qnorm(0.975) * med_mle$se.fit[1] >= med_0_true)
    med_mse_mle_0 <- (med_est_mle_0 - med_0_true)^2
    med_est_mle_1 <- med_mle$fit[2]
    med_ci_mle_1 <- (med_mle$fit[2] - qnorm(0.975) * med_mle$se.fit[2] <= med_1_true) & (med_mle$fit[2] + qnorm(0.975) * med_mle$se.fit[2] >= med_1_true)
    med_mse_mle_1 <- (med_est_mle_1 - med_1_true)^2
    
    
    ## multiple imputation ##
    mi_result <- dl_mi(b=10, data = dat)
    b_est_mi <- mi_result[1]
    b_ci_mi <- mi_result[2]
    b_mse_mi <- (b_est_mi - beta_true)^2
    med_est_mi_0 <- mi_result[3]
    med_ci_mi_0 <- mi_result[4]
    med_mse_mi_0 <- (med_est_mi_0 - med_0_true)^2
    med_est_mi_1 <- mi_result[5]
    med_ci_mi_1 <- mi_result[6]
    med_mse_mi_1 <- (med_est_mi_1 - med_1_true)^2
    
    
    return(c(b_est_cpm,
             b_ci_cpm,
             b_mse_cpm,
             med_est_cpm_0,
             med_ci_cpm_0,
             med_mse_cpm_0,
             med_est_cpm_1,
             med_ci_cpm_1,
             med_mse_cpm_1,
             
             b_est_imp,
             b_ci_imp,
             b_mse_imp,
             med_est_imp_0,
             med_ci_imp_0,
             med_mse_imp_0,
             med_est_imp_1,
             med_ci_imp_1,
             med_mse_imp_1,
             
             b_est_imp2,
             b_ci_imp2,
             b_mse_imp2,
             med_est_imp2_0,
             med_ci_imp2_0,
             med_mse_imp2_0,
             med_est_imp2_1,
             med_ci_imp2_1,
             med_mse_imp2_1,
             
             b_est_mle,
             b_ci_mle,
             b_mse_mle,
             med_est_mle_0,
             med_ci_mle_0,
             med_mse_mle_0,
             med_est_mle_1,
             med_ci_mle_1,
             med_mse_mle_1,
             
             b_est_mi,
             b_ci_mi,
             b_mse_mi,
             med_est_mi_0,
             med_ci_mi_0,
             med_mse_mi_0,
             med_est_mi_1,
             med_ci_mi_1,
             med_mse_mi_1
    )) 
  }) %>% t %>% as.data.frame
  
  colnames(res) <- Cs(
    b_est_cpm,
    b_ci_cpm,
    b_mse_cpm,
    med_est_cpm_0,
    med_ci_cpm_0,
    med_mse_cpm_0,
    med_est_cpm_1,
    med_ci_cpm_1,
    med_mse_cpm_1,
    
    b_est_imp,
    b_ci_imp,
    b_mse_imp,
    med_est_imp_0,
    med_ci_imp_0,
    med_mse_imp_0,
    med_est_imp_1,
    med_ci_imp_1,
    med_mse_imp_1,
    
    b_est_imp2,
    b_ci_imp2,
    b_mse_imp2,
    med_est_imp2_0,
    med_ci_imp2_0,
    med_mse_imp2_0,
    med_est_imp2_1,
    med_ci_imp2_1,
    med_mse_imp2_1,
    
    b_est_mle,
    b_ci_mle,
    b_mse_mle,
    med_est_mle_0,
    med_ci_mle_0,
    med_mse_mle_0,
    med_est_mle_1,
    med_ci_mle_1,
    med_mse_mle_1,
    
    b_est_mi,
    b_ci_mi,
    b_mse_mi,
    med_est_mi_0,
    med_ci_mi_0,
    med_mse_mi_0,
    med_est_mi_1,
    med_ci_mi_1,
    med_mse_mi_1
  )
  return(res)
}

cl <- makeCluster(detectCores())
clusterEvalQ(cl, library(rms, quietly = TRUE))
clusterEvalQ(cl, library(tidyverse, quietly = TRUE))
clusterEvalQ(cl, library(fitdistrplus, quietly = TRUE))
clusterEvalQ(cl, library(survival, quietly = TRUE))
clusterEvalQ(cl, library(mvtnorm, quietly = TRUE))
clusterEvalQ(cl, library(quantreg, quietly = TRUE))
clusterExport(cl,varlist=ls())

result_1000 <- func_res(n=1000)

#### misspecification ####
lower_dl <- 13.12
beta_true <- 1

set.seed(35)
x <- rnorm(1e6, 5, 1) 
e <- rnorm(1e6, 0, 1)
y <- (x + e)^2
# quantile(y, 0.165)
med_0_true <- median((e+5)^2)
med_1_true <- median((6+e)^2)

new.data <- data.frame(x = c(5, 6))
              
# multiple imputation
dl_mi <- function(b=10, data=dat){
  dat_fit <- data %>% 
    mutate(left_val = ifelse(y < lower_dl, NA, y),
           right_val = y)
  
  ## get lognormal distribution parameters
  lnorm_par <- fitdistcens(dat_fit %>% dplyr::select(left_val, right_val), "lnorm")
  
  mi_res <- sapply(1:b, function(i){
    # data for this iteration
    data_b <- data
    
    ## get lognormal distribution parameters by sample
    lnorm_par_b <- rmvnorm(1, mean = lnorm_par$estimate, sigma = lnorm_par$vcov)
    
    # number of observations need to be imputed
    n_imp <- sum(data_b$y < lower_dl)
    # imputed values
    val_imp <- rep(NA, n_imp) # imputed values
    temp <- 1 
    while(temp <= n_imp){
      val <- rlnorm(1, lnorm_par_b[1], lnorm_par_b[2]) # random 
      if(val < lower_dl){ # below DL
        val_imp[temp] <- val
        temp <- temp + 1
      }else{
        next
      }
    }
    
    # insert imputed values
    data_b$y[data_b$y < lower_dl] <- val_imp
    
    
    mod_mi <- lm(log(y) ~ x, # correct transformation
                 data = data_b)
    
    mod_rq <- rq(log(y) ~ x,
                 data = data_b)
    X <-  as.matrix(cbind(1, new.data))
    pred <- drop(X %*% mod_rq$coef)
    V <- summary(mod_rq, cov = TRUE)=
    sdpred <- sqrt(diag(X %*% V$cov %*% t(X)))
    
    return(c(mod_mi$coefficients["x"],
             vcov(mod_mi)['x','x'],
             pred, sdpred))
    
  }) 
  
  b_est_mi <- sum(mi_res[1,] / mi_res[2,]) / 
    sum(1 / mi_res[2,])
  b_se_mi <- sqrt(mean(mi_res[2,]) + (1 + 1 / b) * var(mi_res[1,])) 
  b_ci_mi <- (b_est_mi - qnorm(0.975) * b_se_mi <= beta_true) & 
    (b_est_mi + qnorm(0.975) * b_se_mi >= beta_true)
  
  med_est_mi_0 <- sum(mi_res[3,] / mi_res[5,]) / 
    sum(1 / mi_res[5,])
  med_se_mi_0 <- sqrt(mean(mi_res[5,]) + (1 + 1 / b) * var(mi_res[3,])) 
  med_ci_mi_0 <- (exp(med_est_mi_0 - qnorm(0.975) * med_se_mi_0) <= med_0_true) & 
    (exp(med_est_mi_0 + qnorm(0.975) * med_se_mi_0) >= med_0_true)
  
  med_est_mi_1 <- sum(mi_res[4,] / mi_res[6,]) / 
    sum(1 / mi_res[6,])
  med_se_mi_1 <- sqrt(mean(mi_res[6,]) + (1 + 1 / b) * var(mi_res[4,])) 
  med_ci_mi_1 <- (exp(med_est_mi_1 - qnorm(0.975) * med_se_mi_1) <= med_1_true) & 
    (exp(med_est_mi_1 + qnorm(0.975) * med_se_mi_1) >= med_1_true)
  
  return(c(b_est_mi, b_ci_mi, 
           exp(med_est_mi_0), med_ci_mi_0,
           exp(med_est_mi_1), med_ci_mi_1))
  
} 

func_res <- function(n=1000){
  res <- parSapply(cl, 1:reps, function(i){
    # res <- sapply(1:reps, function(i){
    # data 
    set.seed(i)
    # data generation
    x <- rnorm(n, 5, 1) 
    e <- rnorm(n, 0, 1)
    y <- (x + e)^2
    dat <- data.frame(y = y,
                      x = x)
    
    
    ################## data ##################
    ##### data for CPM ####
    dat_cpm <- dat
    dat_cpm$y <- ifelse(dat_cpm$y < lower_dl, lower_dl - 1e-5, dat_cpm$y)
    
    ##### data for single imputation - sqrt(2) ####
    dat_imp <- dat
    dat_imp$y <- ifelse(dat_imp$y < lower_dl, lower_dl/sqrt(2), dat_imp$y)
    
    ##### data for single imputation - 2 ####
    dat_imp2 <- dat
    dat_imp2$y <- ifelse(dat_imp2$y < lower_dl, lower_dl/2, dat_imp2$y)
    
    
    #### data for MLE (survival) ####
    delta <- dat$y > lower_dl
    dat_mle <- dat
    dat_mle$y <- ifelse(dat_mle$y < lower_dl, lower_dl, dat_mle$y)
    
    ################## models ##################
    ## CPM ##
    mod_cpm <- orm(y ~ x,
                   data = dat_cpm,
                   family = probit)
    # beta
    b_est_cpm <- mod_cpm$coefficients["x"]
    b_se_cpm <- sqrt(mod_cpm$var["x","x"])
    b_ci_cpm <- (b_est_cpm - qnorm(0.975) * b_se_cpm <= beta_true) & 
      (b_est_cpm + qnorm(0.975) * b_se_cpm >= beta_true)
    b_mse_cpm <- (b_est_cpm - beta_true)^2
    # median 0
    med_cpm <- quantile.orm(mod_cpm, new.data, 0.5)
    med_est_cpm_0 <- med_cpm$quantile[1]
    med_ci_cpm_0 <- (med_cpm$lb[1] <= med_0_true) & 
      (med_cpm$ub[1] >= med_0_true)
    med_mse_cpm_0 <- (med_est_cpm_0 - med_0_true)^2
    # median 1
    med_est_cpm_1 <- med_cpm$quantile[2]
    med_ci_cpm_1 <- (med_cpm$lb[2] <= med_1_true) & 
      (med_cpm$ub[2] >= med_1_true)
    med_mse_cpm_1 <- (med_est_cpm_1 - med_1_true)^2
    
    ## dichotomization ##
    # mod_dic <- glm(y ~ x, 
    #                data = dat_dic,
    #                family = 'binomial')
    # # beta
    # est_dic <- mod_dic$coefficients["x"]
    # se_dic <- vcov(mod_dic)['x','x'] %>% sqrt
    # mse_dic <- (est_dic - beta_true)^2
    
    ## single imputation - sqrt(2) ## 
    mod_imp <- lm(log(y) ~ x, # correct transformation
                  data = dat_imp)
    b_est_imp <- mod_imp$coefficients["x"]
    b_se_imp <- vcov(mod_imp)['x','x'] %>% sqrt
    b_ci_imp <- (b_est_imp - qnorm(0.975) * b_se_imp <= beta_true) & 
      (b_est_imp + qnorm(0.975) * b_se_imp >= beta_true)
    b_mse_imp <- (b_est_imp - beta_true)^2
    # medan
    mod_rq <- rq(log(y) ~ x,
                 data = dat_imp)
    X <-  as.matrix(cbind(1, new.data))
    pred <-drop(X %*% mod_rq$coef)
    V <- summary(mod_rq, cov = TRUE)
    # df <- V$rdf
    # tfrac <- qt((1 - 0.95)/2, df)
    sdpred <- sqrt(diag(X %*% V$cov %*% t(X)))
    # med_imp <- cbind(pred, pred + tfrac * sdpred %o% 
    # c(1, -1)) %>% exp
    # colnames(med_imp) <- c("fit", "lower", "upper")
    med_est_imp_0 <- exp(pred[1])
    med_ci_imp_0 <- (exp(pred[1] - qnorm(0.975) * sdpred[1]) <= med_0_true) & 
      (exp(pred[1] + qnorm(0.975) * sdpred[1]) >= med_0_true)
    med_mse_imp_0 <- (med_est_imp_0 - med_0_true)^2
    med_est_imp_1 <- exp(pred[2])
    med_ci_imp_1 <- (exp(pred[2] - qnorm(0.975) * sdpred[2]) <= med_1_true) & 
      (exp(pred[2] + qnorm(0.975) * sdpred[2]) >= med_1_true)
    med_mse_imp_1 <- (med_est_imp_1 - med_1_true)^2
    
    ## single imputation - 2 ## 
    mod_imp2 <- lm(log(y) ~ x, # correct transformation
                   data = dat_imp2)
    b_est_imp2 <- mod_imp2$coefficients["x"]
    b_se_imp2 <- vcov(mod_imp2)['x','x'] %>% sqrt
    b_ci_imp2 <- (b_est_imp2 - qnorm(0.975) * b_se_imp2 <= beta_true) & 
      (b_est_imp2 + qnorm(0.975) * b_se_imp2 >= beta_true)
    b_mse_imp2 <- (b_est_imp2 - beta_true)^2
    # medan
    mod_rq <- rq(log(y) ~ x,
                 data = dat_imp2)
    X <-  as.matrix(cbind(1, new.data))
    pred <- drop(X %*% mod_rq$coef)
    V <- summary(mod_rq, cov = TRUE)
    # df <- V$rdf
    # tfrac <- qt((1 - 0.95)/2, df)
    sdpred <- sqrt(diag(X %*% V$cov %*% t(X)))
    # med_imp2 <- cbind(pred, pred + tfrac * sdpred %o% 
    #                        c(1, -1)) %>% exp
    # colnames(med_imp2) <- c("fit", "lower", "upper")
    med_est_imp2_0 <- exp(pred[1])
    med_ci_imp2_0 <- (exp(pred[1] - qnorm(0.975) * sdpred[1]) <= med_0_true) & 
      (exp(pred[1] + qnorm(0.975) * sdpred[1]) >= med_0_true)
    med_mse_imp2_0 <- (med_est_imp2_0 - med_0_true)^2
    med_est_imp2_1 <- exp(pred[2])
    med_ci_imp2_1 <- (exp(pred[2] - qnorm(0.975) * sdpred[2]) <= med_1_true) & 
      (exp(pred[2] + qnorm(0.975) * sdpred[2]) >= med_1_true)
    med_mse_imp2_1 <- (med_est_imp2_1 - med_1_true)^2
    
    ## MLE (lognormal) ##
    mod_mle <- survreg(Surv(dat_mle$y, delta, type='left') ~ x, 
                       dist = 'lognormal')
    b_est_mle <- mod_mle$coefficients["x"]
    b_se_mle <- vcov(mod_mle)['x','x'] %>% sqrt
    b_ci_mle <- (b_est_mle - qnorm(0.975) * b_se_mle <= beta_true) & 
      (b_est_mle + qnorm(0.975) * b_se_mle >= beta_true)
    b_mse_mle <- (b_est_mle - beta_true)^2
    # median (censored quantile regression)
    med_mle <- predict(mod_mle, newdata=new.data, 
                       type = "quantile", se=T, p=0.5)
    med_est_mle_0 <- med_mle$fit[1]
    med_ci_mle_0 <- (med_mle$fit[1] - qnorm(0.975) * med_mle$se.fit[1] <= med_0_true) & (med_mle$fit[1] + qnorm(0.975) * med_mle$se.fit[1] >= med_0_true)
    med_mse_mle_0 <- (med_est_mle_0 - med_0_true)^2
    med_est_mle_1 <- med_mle$fit[2]
    med_ci_mle_1 <- (med_mle$fit[2] - qnorm(0.975) * med_mle$se.fit[2] <= med_1_true) & (med_mle$fit[2] + qnorm(0.975) * med_mle$se.fit[2] >= med_1_true)
    med_mse_mle_1 <- (med_est_mle_1 - med_1_true)^2
    
    
    ## multiple imputation ##
    mi_result <- dl_mi(b=10, data = dat)
    b_est_mi <- mi_result[1]
    b_ci_mi <- mi_result[2]
    b_mse_mi <- (b_est_mi - beta_true)^2
    med_est_mi_0 <- mi_result[3]
    med_ci_mi_0 <- mi_result[4]
    med_mse_mi_0 <- (med_est_mi_0 - med_0_true)^2
    med_est_mi_1 <- mi_result[5]
    med_ci_mi_1 <- mi_result[6]
    med_mse_mi_1 <- (med_est_mi_1 - med_1_true)^2
    
    
    return(c(b_est_cpm,
             b_ci_cpm,
             b_mse_cpm,
             med_est_cpm_0,
             med_ci_cpm_0,
             med_mse_cpm_0,
             med_est_cpm_1,
             med_ci_cpm_1,
             med_mse_cpm_1,
             
             b_est_imp,
             b_ci_imp,
             b_mse_imp,
             med_est_imp_0,
             med_ci_imp_0,
             med_mse_imp_0,
             med_est_imp_1,
             med_ci_imp_1,
             med_mse_imp_1,
             
             b_est_imp2,
             b_ci_imp2,
             b_mse_imp2,
             med_est_imp2_0,
             med_ci_imp2_0,
             med_mse_imp2_0,
             med_est_imp2_1,
             med_ci_imp2_1,
             med_mse_imp2_1,
             
             b_est_mle,
             b_ci_mle,
             b_mse_mle,
             med_est_mle_0,
             med_ci_mle_0,
             med_mse_mle_0,
             med_est_mle_1,
             med_ci_mle_1,
             med_mse_mle_1,
             
             b_est_mi,
             b_ci_mi,
             b_mse_mi,
             med_est_mi_0,
             med_ci_mi_0,
             med_mse_mi_0,
             med_est_mi_1,
             med_ci_mi_1,
             med_mse_mi_1
    )) 
  }) %>% t %>% as.data.frame
  
  colnames(res) <- Cs(
    b_est_cpm,
    b_ci_cpm,
    b_mse_cpm,
    med_est_cpm_0,
    med_ci_cpm_0,
    med_mse_cpm_0,
    med_est_cpm_1,
    med_ci_cpm_1,
    med_mse_cpm_1,
    
    b_est_imp,
    b_ci_imp,
    b_mse_imp,
    med_est_imp_0,
    med_ci_imp_0,
    med_mse_imp_0,
    med_est_imp_1,
    med_ci_imp_1,
    med_mse_imp_1,
    
    b_est_imp2,
    b_ci_imp2,
    b_mse_imp2,
    med_est_imp2_0,
    med_ci_imp2_0,
    med_mse_imp2_0,
    med_est_imp2_1,
    med_ci_imp2_1,
    med_mse_imp2_1,
    
    b_est_mle,
    b_ci_mle,
    b_mse_mle,
    med_est_mle_0,
    med_ci_mle_0,
    med_mse_mle_0,
    med_est_mle_1,
    med_ci_mle_1,
    med_mse_mle_1,
    
    b_est_mi,
    b_ci_mi,
    b_mse_mi,
    med_est_mi_0,
    med_ci_mi_0,
    med_mse_mi_0,
    med_est_mi_1,
    med_ci_mi_1,
    med_mse_mi_1
  )
  return(res)
}

result_1000 <- func_res(n=1000)
              