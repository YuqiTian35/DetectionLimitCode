# Table 3 
library(tidyverse)
library(rms)
library(survival)
library(fitdistrplus)
library(mvtnorm)
library(quantreg)
library(multipleDL)

# multiple imputation function
## b: number of iterations combined for multiple imputation 
dl_mi <- function(b=10, data){
  lower_dl <- data$lower_dl
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
    ind_imp <- which(data_b$y < lower_dl)
    temp <- 1 
    while(temp <= n_imp){
      val <- rlnorm(1, lnorm_par_b[1], lnorm_par_b[2]) # random 
      if(val < lower_dl[ind_imp[temp]]){ # below DL
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

####### Correct Transformation #######
# use 10 as an example. Please set to 1000
reps <-10

# true values 
beta_true <- 1
med_0_true <- 1
med_1_true <- exp(1)

# new data
new.data <- data.frame(x = c(0, 1))

# detection limits
lower_dl <- c(0.16, 0.30, 0.5)

# function for data generation
data_gen <- function(n_each, num_site, lower_dl = c(0.16, 0.30, 0.5), upper_dl = NULL){
  
  site_name <- 1:num_site # name of site
  # x dependent of site
  x1 <- rnorm(n_each, -0.5, 1)
  x2 <- rnorm(n_each, 0, 1)
  x3 <- rnorm(n_each, 0.5, 1) 
  
  # same transformation for each site
  y1 <- exp(x1 + rnorm(n_each, 0, 1))
  y2 <- exp(x2 + rnorm(n_each, 0, 1))
  y3 <- exp(x3 + rnorm(n_each, 0, 1))
  # combine all 3 sites
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
  
  # DL values 
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

# function to generate results
func_res <- function(n_each, num_site = 3){
  res <- sapply(1:reps, function(i){
    # data 
    set.seed(i)
    dat <- data_gen(n_each = n_each, num_site = num_site)
    
    ################## data ##################
    ##### data for CPM ####
    dat_cpm <- dat
  
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
    mod_cpm <- multipleDL(y_obs ~ x, data = dat_cpm, link = 'probit', delta_lower = dat_cpm$dl_lower)
    
    # beta
    b_est_cpm <- mod_cpm$coef["x"]
    b_se_cpm <- sqrt(mod_cpm$var["x","x"])
    b_ci_cpm <- (b_est_cpm - qnorm(0.975) * b_se_cpm <= beta_true) & 
      (b_est_cpm + qnorm(0.975) * b_se_cpm >= beta_true)
    b_mse_cpm <- (b_est_cpm - beta_true)^2
    # median 0
    med_cpm <- quantile_dl(mod_cpm, new.data, 0.5)
    med_est_cpm_0 <- med_cpm$est[1]
    med_ci_cpm_0 <- (med_cpm$lb[1] <= med_0_true) & 
      (med_cpm$ub[1] >= med_0_true)
    med_mse_cpm_0 <- (med_est_cpm_0 - med_0_true)^2
    # median 1
    med_est_cpm_1 <- med_cpm$est[2]
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
    # median
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
    # median
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
    mod_mle <- survreg(Surv(dat_mle$y, delta, type='left') ~ x,  data=dat_mle, 
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

# summary of each iteration
sum_tab <- function(result){
  return(matrix(c(
    ### CPM
    # bias of beta
    mean(result$b_est_cpm) - beta_true,
    # percent bias
    (mean(result$b_est_cpm) - beta_true) / beta_true * 100,
    # sd of beta est
    sd(result$b_est_cpm),
    # rmse of beta
    sqrt(mean(result$b_mse_cpm)),
    # ci coverage of beta
    mean(result$b_ci_cpm),
    # bias of med 0
    mean(result$med_est_cpm_0) - med_0_true,
    # percent bias
    (mean(result$med_est_cpm_0) - med_0_true) / med_0_true * 100,
    # sd of med0 est
    sd(result$med_est_cpm_0),
    # rmse of beta
    sqrt(mean(result$med_mse_cpm_0)),
    # ci coverage of beta
    mean(result$med_ci_cpm_0),
    # bias of med 1
    mean(result$med_est_cpm_1) - med_1_true,
    # percent bias
    (mean(result$med_est_cpm_1) - med_1_true) / med_1_true * 100,
    # sd of med0 est
    sd(result$med_est_cpm_1),
    # rmse of beta
    sqrt(mean(result$med_mse_cpm_1)),
    # ci coverage of beta
    mean(result$med_ci_cpm_1),
    
    ### single imputation - sqrt(2)
    # bias of beta
    mean(result$b_est_imp) - beta_true,
    # percent bias
    (mean(result$b_est_imp) - beta_true) / beta_true * 100,
    # sd of beta est
    sd(result$b_est_imp),
    # rmse of beta
    sqrt(mean(result$b_mse_imp)),
    # ci coverage of beta
    mean(result$b_ci_imp),
    # bias of med 0
    mean(result$med_est_imp_0) - med_0_true,
    # percent bias
    (mean(result$med_est_imp_0) - med_0_true) / med_0_true * 100,
    # sd of med0 est
    sd(result$med_est_imp_0),
    # rmse of beta
    sqrt(mean(result$med_mse_imp_0)),
    # ci coverage of beta
    mean(result$med_ci_imp_0),
    # bias of med 1
    mean(result$med_est_imp_1) - med_1_true,
    # percent bias
    (mean(result$med_est_imp_1) - med_1_true) / med_1_true * 100,
    # sd of med0 est
    sd(result$med_est_imp_1),
    # rmse of beta
    sqrt(mean(result$med_mse_imp_1)),
    # ci coverage of beta
    mean(result$med_ci_imp_1),
    
    ### single imputation - 2
    # bias of beta
    mean(result$b_est_imp2) - beta_true,
    # percent bias
    (mean(result$b_est_imp2) - beta_true) / beta_true * 100,
    # sd of beta est
    sd(result$b_est_imp2),
    # rmse of beta
    sqrt(mean(result$b_mse_imp2)),
    # ci coverage of beta
    mean(result$b_ci_imp2),
    # bias of med 0
    mean(result$med_est_imp2_0) - med_0_true,
    # percent bias
    (mean(result$med_est_imp2_0) - med_0_true) / med_0_true * 100,
    # sd of med0 est
    sd(result$med_est_imp2_0),
    # rmse of beta
    sqrt(mean(result$med_mse_imp2_0)),
    # ci coverage of beta
    mean(result$med_ci_imp2_0),
    # bias of med 1
    mean(result$med_est_imp2_1) - med_1_true,
    # percent bias
    (mean(result$med_est_imp2_1) - med_1_true) / med_1_true * 100,
    # sd of med0 est
    sd(result$med_est_imp2_1),
    # rmse of beta
    sqrt(mean(result$med_mse_imp2_1)),
    # ci coverage of beta
    mean(result$med_ci_imp2_1),
    
    ### MLE
    # bias of beta
    mean(result$b_est_mle) - beta_true,
    # percent bias
    (mean(result$b_est_mle) - beta_true) / beta_true * 100,
    # sd of beta est
    sd(result$b_est_mle),
    # rmse of beta
    sqrt(mean(result$b_mse_mle)),
    # ci coverage of beta
    mean(result$b_ci_mle),
    # bias of med 0
    mean(result$med_est_mle_0) - med_0_true,
    # percent bias
    (mean(result$med_est_mle_0) - med_0_true) / med_0_true * 100,
    # sd of med0 est
    sd(result$med_est_mle_0),
    # rmse of beta
    sqrt(mean(result$med_mse_mle_0)),
    # ci coverage of beta
    mean(result$med_ci_mle_0),
    # bias of med 1
    mean(result$med_est_mle_1) - med_1_true,
    # percent bias
    (mean(result$med_est_mle_1) - med_1_true) / med_1_true * 100,
    # sd of med0 est
    sd(result$med_est_mle_1),
    # rmse of beta
    sqrt(mean(result$med_mse_mle_1)),
    # ci coverage of beta
    mean(result$med_ci_mle_1),
    
    ### multiple imputation
    # bias of beta
    mean(result$b_est_mi) - beta_true,
    # percent bias
    (mean(result$b_est_mi) - beta_true) / beta_true * 100,
    # sd of beta est
    sd(result$b_est_mi),
    # rmse of beta
    sqrt(mean(result$b_mse_mi)),
    # ci coverage of beta
    mean(result$b_ci_mi),
    # bias of med 0
    mean(result$med_est_mi_0) - med_0_true,
    # percent bias
    (mean(result$med_est_mi_0) - med_0_true) / med_0_true * 100,
    # sd of med0 est
    sd(result$med_est_mi_0),
    # rmse of beta
    sqrt(mean(result$med_mse_mi_0)),
    # ci coverage of beta
    mean(result$med_ci_mi_0),
    # bias of med 1
    mean(result$med_est_mi_1) - med_1_true,
    # percent bias
    (mean(result$med_est_mi_1) - med_1_true) / med_1_true * 100,
    # sd of med0 est
    sd(result$med_est_mi_1),
    # rmse of beta
    sqrt(mean(result$med_mse_mi_1)),
    # ci coverage of beta
    mean(result$med_ci_mi_1)),
    ncol=5))
}

# create a table to demonstrate results
func_tab <- function(result){
  tab <- as.data.frame(result)
  
  colnames(tab) <- c("CPM",
                     "Single imputation - DL/sqrt(2)",
                     "Single imputation - DL/2",
                     "MLE - lognormal",
                     "Multiple imputation - lognormal")
  
  rownames(tab) <- c("bias of beta",  
                     "percent bias of beta", 
                     "empirical estimates of sd(beta)",
                     "rmse of beta",
                     "ci of beta",
                     
                     "bias of med_0",  
                     "percent bias of med_0", 
                     "empirical estimates of sd(med_0)",
                     "rmse of med_0",
                     "ci of med_0",
                     
                     "bias of med_1",  
                     "percent bias of med_1", 
                     "empirical estimates of sd(med_1)",
                     "rmse of med_1",
                     "ci of med_1")
  
  return(tab)
}

result_correct <- func_res(n_each=300, num_site=3)
print("Correct Transformation")
print(func_tab(sum_tab(result_correct)))


######### Incorrect Transformation ##########
# true values
beta_true <- 1
set.seed(35)
y_sim_0 <-  qchisq(pnorm((0 + rnorm(1e6))/2), 5)
y_sim_1 <-  qchisq(pnorm((1 + rnorm(1e6))/2), 5)
# median
med_0_true <- median(y_sim_0)
med_1_true <- median(y_sim_1)

# detection limilts
lower_dl <- c(2, 3, 4)

# function for data generation
data_gen_incorrect <- function(n_each, num_site, lower_dl = c(2, 3, 4), upper_dl = NULL){
  
  site_name <- 1:num_site # name of site
  # x dependent of site
  x1 <- rnorm(n_each, 4.5, 1)
  x2 <- rnorm(n_each, 5, 1)
  x3 <- rnorm(n_each, 5.5, 1) 
  
  # same transformation for each site
  qchisq(pnorm((0 + rnorm(1e6))/2), 5)
  y1 <- qchisq(pnorm((x1 + rnorm(n_each, 0, 1))/2), 5)
  y2 <- qchisq(pnorm((x2 + rnorm(n_each, 0, 1))/2), 5)
  y3 <- qchisq(pnorm((x3 + rnorm(n_each, 0, 1))/2), 5)
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

func_res_incorrect <- function(n_each=300, num_site =3){
  res <- sapply(1:reps, function(i){
    # data 
    set.seed(i)
    # data generation
    dat <- data_gen_incorrect(n_each = n_each, num_site = num_site)
    
    ################## data ##################
    ##### data for CPM ####
    dat_cpm <- dat
    
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
    mod_cpm <- multipleDL(y_obs ~ x, data = dat_cpm, link = 'probit', delta_lower = dat_cpm$dl_lower)
    
    # beta
    b_est_cpm <- mod_cpm$coef["x"]
    b_se_cpm <- sqrt(mod_cpm$var["x","x"])
    b_ci_cpm <- (b_est_cpm - qnorm(0.975) * b_se_cpm <= beta_true) & 
      (b_est_cpm + qnorm(0.975) * b_se_cpm >= beta_true)
    b_mse_cpm <- (b_est_cpm - beta_true)^2
    # median 0
    med_cpm <- quantile_dl(mod_cpm, new.data, 0.5)
    med_est_cpm_0 <- med_cpm$est[1]
    med_ci_cpm_0 <- (med_cpm$lb[1] <= med_0_true) & 
      (med_cpm$ub[1] >= med_0_true)
    med_mse_cpm_0 <- (med_est_cpm_0 - med_0_true)^2
    # median 1
    med_est_cpm_1 <- med_cpm$est[2]
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
    # median
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
    # median
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
    mod_mle <- survreg(Surv(dat_mle$y, delta, type='left') ~ x,  data = dat_mle,
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

result_incorrect <- func_res_incorrect(n_each=300, num_site=3)
print("Incorrect Transformation")
print(func_tab(sum_tab(result_incorrect)))
