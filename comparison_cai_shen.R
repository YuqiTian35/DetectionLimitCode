# Comparison: Table 4
library(MASS)
library(survival)
library(multipleDL)


cai_function <-
  function(formula = formula,data = data,r=0,subset,dx=0.001,iter.max=100,num.sim=200, 
           left_censoring, right_censoring){

    obs.t <- data$y_obs
    delta_left <- ifelse(obs.t > data$lower_dl, 1, 0)
    delta_right <- ifelse(obs.t < data$upper_dl, 1, 0)
    delta <- delta_left * delta_right 
    z <- model.matrix(formula, data)[, -1]
    if(class(z) == 'numeric'){
      p <- 1
      n <- length(z)
    }else{
      p = dim(z)[2]
      n = dim(z)[1]
    }
    
    left_censoring <- left_censoring
    right_censoring <- right_censoring
    
    
    if(p>1) {zc = t(t(z) - apply(z,2,mean))}
    if(p==1) {zc = as.matrix(z-mean(z))}
    varnms<-colnames(z)
    if(p>0){
      ix = order(obs.t)
      ord.t <- obs.t[ix]
      ord.delta = delta[ix]
      ord.delta_right = delta_right[ix]
      ord.delta_left = delta_left[ix]
      ord.z = as.matrix(zc[ix,])
      colnames(ord.z)<- colnames(z)<-varnms
      
      beta.ini = rep(0,p)
      Rt.ini = cumsum(1/(n:1))
      
      ###solving estimating equations for beta and Rt
      wt<-rep(1,n)
      temp.b<-solve.beta.cai(beta.ini, Rt.ini, obs.t, delta, z, wt, r, dx, iter.max, left_censoring, right_censoring)
      
      beta<-temp.b$Beta
      Rt<-temp.b$data$Rt
      converged<-temp.b$converged
      iter<-temp.b$iter

      names(beta) <- colnames(z)
      
      output<-list(coefficients=-beta, formula=formula,
                   r=r,converged=c(converged,iter),df=n-1,call=match.call())
      
      output$p = p
      output$ord.time<-sort(obs.t)
      output$ord.delta=ord.delta
      output$z=z
      # output$ord.z=as.matrix(z[order(obs.t),])
      output$ord.z=as.matrix(z[order(obs.t)])
      output$Rt=temp.b$data$Rt
    }
    
    return(output)
  }


solve.beta.cai <-
  function(beta.ini,Rt.ini,obs.t,delta,z,wt,r,dx,iter.max, left_censoring, right_censoring){
    p = length(beta.ini)
    n = length(Rt.ini)
    
    if(p>1) zc = t(t(z) - apply(z,2,mean))
    if(p==1) zc = as.matrix(z-mean(z))
    ix = order(obs.t)
    ord.delta = delta[ix]
    ord.z = as.matrix(zc[ix,])
    colnames(ord.z)= colnames(z)
    ord.t <- obs.t[ix]
    ord.wt = wt[ix]
    ord.left_censoring =left_censoring[ix]
    ord.right_censoring = right_censoring[ix]
    
    # loglog
    g <- function(x){
      # return(1 - exp(-exp(x)))
      return(exp(-exp(x)))
    }
    
    ###solving estimating equations for beta and Rt
    dif = 1
    iter = 0
    while(dif > dx & iter <= iter.max){
      if(r == 0){
        # ee3 
        ai <- rep(0, n)
        for(i in 1:n){
          ai[i] <- sum(sapply(1:n, function(j){
            (ord.left_censoring[i] < ord.t[j]) & (ord.t[j] <= ord.t[i])
          }))
        }
        
        j_indicator <- vector("list", length = n)
        for(i in 1:n){
          j_indicator[[i]] <- sapply(1:n, function(j){
            (ord.left_censoring[i] < ord.t[j]) & (ord.t[j] <= ord.right_censoring[i])
          })
        }
        
        
        ee3b <- function(b){
          g_part_i <- sapply(1:n, function(i){
            return(sum(g(as.numeric(ord.z[i,] %*% b) + Rt.ini)[j_indicator[[i]]]))
          })
          result <- sum(ord.z * (ai - g_part_i))
          # result <- sum(result^2)  # Square the column sums and sum them
          return(result)
        }

        beta <- uniroot(ee3b, interval = c(-20, 20))$root 
        
        # ee4
        Rt <- sapply(1:n, function(k){
          ee4_k <- function(ht){
            left_side <- sum(sapply(1:n, function(i){
              (ord.left_censoring[i] < ord.t[k]) & (ord.t[k] <= ord.t[i])
            }))
            
            right_side <- sum(sapply(1:n, function(i){
              ((ord.left_censoring[i] < ord.t[k]) & (ord.t[k] <= ord.right_censoring[i])) * 
                g(as.numeric(ord.z[i,] %*% beta) + ht)
            }))
            
            return(left_side - right_side)
          }
          return(uniroot(ee4_k, interval = c(-20, 20))$root)
        })
        
        
        dif = max(abs(beta-beta.ini))
        iter = iter + 1
        beta.ini = beta
        Rt.ini = Rt
      }
    }
    converged = as.numeric(dif > dx)  
    Hbeta = 'not yet'
    return(list(Beta=beta,converged=converged,iter=iter,data=data.frame(ord.time=sort(obs.t),ord.delta=ord.delta,Rt=Rt),Hbeta=Hbeta))
  }


shen_function <-
  function(formula = formula,data = data,r=0,subset,dx=0.001,iter.max=100,num.sim=200, left_censoring){
    obs.t <- data$y_obs
    delta_left <- ifelse(obs.t > data$lower_dl, 1, 0)
    delta_right <- ifelse(obs.t < data$upper_dl, 1, 0)
    delta <- delta_left * delta_right 
    z <- model.matrix(formula, data)[, -1]
    if(class(z) == 'numeric'){
      p <- 1
      n <- length(z)
    }else{
      p = dim(z)[2]
      n = dim(z)[1]
    }
    
    left_censoring <- left_censoring
    
    
    if(p>1) {zc = t(t(z) - apply(z,2,mean))}
    if(p==1) {zc = as.matrix(z-mean(z))}
    varnms<-colnames(z)
    if(p>0){
      ix = order(obs.t)
      ord.t <- obs.t[ix]
      ord.delta = delta[ix]
      ord.delta_right = delta_right[ix]
      ord.delta_left = delta_left[ix]
      ord.z = as.matrix(zc[ix,])
      colnames(ord.z)<- colnames(z)<-varnms
      
      beta.ini = rep(0,p)
      Rt.ini = cumsum(1/(n:1))
      
      ###solving estimating equations for beta and Rt
      wt<-rep(1,n)
      temp.b<-solve.beta.shen(beta.ini, Rt.ini, obs.t, delta, z, wt, r, dx, iter.max, left_censoring)
      
      beta<-temp.b$Beta
      Rt<-temp.b$data$Rt
      converged<-temp.b$converged
      iter<-temp.b$iter
      ###compute the variance-covariance matrix of beta
      
      
      output<-list(coefficients=-beta, formula=formula,
                   r=r,converged=c(converged,iter),df=n-1,call=match.call())
      
      output$p = p
      output$ord.time<-sort(obs.t)
      output$ord.delta=ord.delta
      output$z=z
      # output$ord.z=as.matrix(z[order(obs.t),])
      output$ord.z=as.matrix(z[order(obs.t)])
      output$Rt=temp.b$data$Rt
    }
    
    return(output)
  }


solve.beta.shen <-
  function(beta.ini,Rt.ini,obs.t,delta,z,wt,r,dx,iter.max, left_censoring){
    p = length(beta.ini)
    n = length(Rt.ini)
    
    if(p>1) zc = t(t(z) - apply(z,2,mean))
    if(p==1) zc = as.matrix(z-mean(z))
    ix = order(obs.t)
    ord.delta = delta[ix]
    ord.z = as.matrix(zc[ix,])
    colnames(ord.z)= colnames(z)
    ord.t <- obs.t[ix]
    ord.wt = wt[ix]
    ord.left_censoring = left_censoring[ix]
    
    ###solving estimating equations for beta and Rt
    dif = 1
    iter = 0
    while(dif > dx & iter <= iter.max){
      if(r == 0){
        # H(Li)
        H_Li <- sapply(ord.left_censoring, function(x){
          if(x < min(ord.t)){
            return(-9999)
          }else{
            return(Rt.ini[max(which(ord.t <= x))])
          }})
        
        # EE1 (estimating equation for solving beta)
        ee1 <- function(b){
          result <- sum(
            sapply(1:n, function(i){
              if(ord.t[i] == ord.left_censoring[i]){
                return(0)
              }else if(ord.left_censoring[i] == -Inf){
                return(ord.z[i,] * (ord.delta[i] - exp(as.numeric(ord.z[i,] %*% b) + Rt.ini[i])))
              }else{
                return(ord.z[i,] * (ord.delta[i] - exp(as.numeric(ord.z[i,] %*% b) + Rt.ini[i]) + 
                                      exp(as.numeric(ord.z[i,] %*% b) + H_Li[i])))
              }
            })
          )
          return(result)
          
        }
  
        beta <- uniroot(ee1, interval = c(-20, 20))$root 
        
        # EE2 (estimating equation for solving H)
        Rt <- rep(0, n)
        ord.bz <- as.numeric(ord.z %*% beta)
        if(any(ord.left_censoring < ord.t[1])){
          ee2_1 <- function(h){
            return(ord.delta[1] - sum(exp(ord.bz + h)[ord.left_censoring < ord.t[1]]))
          }
          Rt[1] <- uniroot(ee2_1, interval = c(-20, 20))$root
          
          for(k in 2:n){
            if(ord.t[k-1] == ord.t[k]){
              Rt[k] <- Rt[k-1]
            }else{
              denominator_indicator <- sapply(1:n, function(i){
                return((ord.left_censoring[i] < ord.t[k]) &  (ord.t[k] <= ord.t[i]))
              })
              Rt[k] <- Rt[k-1] + ord.delta[k] / sum(exp(ord.bz[denominator_indicator] + Rt[k-1]))
            }
          }
        }else{
          Rt[1] <- -Inf
          
          for(k in 2:n){
            if(ord.t[k-1] == ord.t[k]){
              Rt[k] <- Rt[k-1]
            }else{
              ee2_i <- function(h){
                ord.delta[k] - sum((exp(ord.bz + h) - exp(ord.bz + Rt[k-1]))[(ord.left_censoring < ord.t[k]) &
                                                                               (ord.t[k] <= ord.t)])
              }
              Rt[k] <- uniroot(ee2_i, interval = c(-20, 20))$root
            }
          }
        }
        
        dif = max(abs(beta-beta.ini))
        iter = iter + 1
        beta.ini = beta
        Rt.ini = Rt
      }
    }
    converged = as.numeric(dif > dx)  
    Hbeta = 'not yet'
    return(list(Beta=beta,converged=converged,iter=iter,data=data.frame(ord.time=sort(obs.t),ord.delta=ord.delta,Rt=Rt),Hbeta=Hbeta))
  }


########## Scenario 1 #######
data_gen <- function(n_each, num_site, lower_dl, upper_dl){
  
  site_name <- 1:num_site # name of site
  x <- rnorm(n_each * num_site, 0, 1) 
  e <-  rgumbel(sum(n_each), max=F)
  y <- exp(x + e)
  
  # true data
  dat <- data.frame(site = rep(site_name, each = n_each),
                    y = y,
                    x = x)
  dat$lower_dl <- rep(lower_dl, each = n_each)
  dat$upper_dl <- rep(upper_dl, each = n_each)
  
  
  # the observed data
  dat_obs <- dat
  dat_obs$y_obs <- ifelse(dat$y < dat$lower_dl, dat$lower_dl,
                          ifelse(dat$y > dat$upper_dl, dat$upper_dl, dat$y))
  dat_obs$left <- ifelse(dat$y < dat$lower_dl, -Inf, 
                         ifelse(dat$y > dat$upper_dl, dat$upper_dl, dat$y))
  dat_obs$right <- ifelse(dat$y > dat$upper_dl, Inf, 
                          ifelse(dat$y < dat$lower_dl, dat$lower_dl, dat$y))
  # dat_obs <- dat_obs[,-2]
  
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


simulation <- function(iter){
  
  reps <- 1
  res_all <- lapply(1:reps, function(i){
    set.seed(i + iter * reps)
    
    lower_dl <- c(-Inf, -Inf, -Inf)
    upper_dl <- c(Inf, Inf, Inf)
    beta.positive <- 1
    num_site <- 3
    
    # n=50*3
    data <- data_gen(n_each = 50, num_site = num_site, lower_dl = lower_dl, upper_dl = upper_dl)
    data.obs <- data$dat_obs
    
    tian_50 <- tryCatch({
      multipleDL(y_obs ~ x, data = data.obs, delta_lower = data.obs$dl_lower,
                 delta_upper = data.obs$dl_upper, link = 'cloglog')$coef[['x']]
    }
    , error = function(e) return(NA))
    
    cai_50 <- tryCatch({
      cai_function(formula=Surv(left, right, type = 'interval2') ~ x, data = data.obs, r=0,
                   left_censoring = data.obs$lower_dl, right_censoring = data.obs$upper_dl)$coefficients
    }, error = function(e) return(NA))
    
    shen_50 <- tryCatch({
      shen_function(formula=Surv(left, right, type = 'interval2') ~ x, data = data.obs, r=0, left_censoring = data.obs$lower_dl)$coefficients
    }, error = function(e) return(NA))
    
    tian_50.mse <- (tian_50 - beta.positive)^2
    cai_50.mse <- (cai_50 - beta.positive)^2
    shen_50.mse <- (shen_50 - beta.positive)^2
    
    # n=300*3
    data <- data_gen(n_each = 300, num_site = num_site, lower_dl = lower_dl, upper_dl = upper_dl)
    data.obs <- data$dat_obs
    
    tian_300 <- tryCatch({
      multipleDL(y_obs ~ x, data = data.obs, delta_lower = data.obs$dl_lower,
                 delta_upper = data.obs$dl_upper, link = 'cloglog')$coef[['x']]
    }
    , error = function(e) return(NA))
    
    cai_300 <- tryCatch({
      cai_function(formula=Surv(left, right, type = 'interval2') ~ x, data = data.obs, r=0,
                   left_censoring = data.obs$lower_dl, right_censoring = data.obs$upper_dl)$coefficients
    }, error = function(e) return(NA))
    
    shen_300 <- tryCatch({
      shen_function(formula=Surv(left, right, type = 'interval2') ~ x, data = data.obs, r=0, left_censoring = data.obs$lower_dl)$coefficients
    }, error = function(e) return(NA))
    
    tian_300.mse <- (tian_300 - beta.positive)^2
    cai_300.mse <- (cai_300 - beta.positive)^2
    shen_300.mse <- (shen_300 - beta.positive)^2
    
    return(list(tian_50 = tian_50, tian_50.mse = tian_50.mse, 
                tian_300 = tian_300, tian_300.mse = tian_300.mse,
                shen_50 = shen_50, shen_50.mse = shen_50.mse, 
                cai_50 = cai_50, cai_50.mse = cai_50.mse,
                shen_300 = shen_300, shen_300.mse = shen_300.mse,
                cai_300 = cai_300, cai_300.mse = cai_300.mse))
  })
  
  save(res_all, file=paste('output/sim_',iter,'.Rdata', sep=""))
  
}


######### Simulation 2 #########
# generate the true data and the observed data
data_gen <- function(n_each, num_site, lower_dl, upper_dl){
  
  site_name <- 1:num_site # name of site
  x1 <- rnorm(n_each, -0.5, 1) 
  y1 <- exp(x1 + rgumbel(n_each, max=F))
  x2 <- rnorm(n_each, 0, 1) 
  y2 <- exp(x2 + rgumbel(n_each, max=F))
  x3 <- rnorm(n_each, 0.5, 1) 
  y3 <- exp(x3 + rgumbel(n_each, max=F))
  
  # true data
  dat <- data.frame(site = rep(site_name, each = n_each),
                    y = c(y1, y2, y3),
                    x = c(x1, x2, x3))
  dat$lower_dl <- rep(lower_dl, each = n_each)
  dat$upper_dl <- rep(upper_dl, each = n_each)
  
  
  # the observed data
  dat_obs <- dat
  dat_obs$y_obs <- ifelse(dat$y < dat$lower_dl, dat$lower_dl,
                          ifelse(dat$y > dat$upper_dl, dat$upper_dl, dat$y))
  dat_obs$left <- ifelse(dat$y < dat$lower_dl, -Inf, 
                         ifelse(dat$y > dat$upper_dl, dat$upper_dl, dat$y))
  dat_obs$right <- ifelse(dat$y > dat$upper_dl, Inf, 
                          ifelse(dat$y < dat$lower_dl, dat$lower_dl, dat$y))
  # dat_obs <- dat_obs[,-2]
  
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


simulation <- function(iter){
  
  reps <- 1
  res_all <- lapply(1:reps, function(i){
    set.seed(i + iter * reps)
    
    lower_dl <- c(0.08, 0.16, 0.26)
    upper_dl <- c(Inf, Inf, Inf)
    beta.positive <- 1
    num_site <- 3
    
    # n=50*3
    data <- data_gen(n_each = 50, num_site = num_site, lower_dl = lower_dl, upper_dl = upper_dl)
    data.obs <- data$dat_obs
    
    tian_50 <- tryCatch({
      multipleDL(y_obs ~ x, data = data.obs, delta_lower = data.obs$dl_lower,
                 delta_upper = data.obs$dl_upper, link = 'cloglog')$coef[['x']]
    }
    , error = function(e) return(NA))
    
    cai_50 <- tryCatch({
      cai_function(formula=Surv(left, right, type = 'interval2') ~ x, data = data.obs, r=0,
                   left_censoring = data.obs$lower_dl, right_censoring = data.obs$upper_dl)$coefficients
    }, error = function(e) return(NA))
    
    shen_50 <- tryCatch({
      shen_function(formula=Surv(left, right, type = 'interval2') ~ x, data = data.obs, r=0, 
                    left_censoring = data.obs$lower_dl)$coefficients
    }, error = function(e) return(NA))
    
    tian_50.mse <- (tian_50 - beta.positive)^2
    cai_50.mse <- (cai_50 - beta.positive)^2
    shen_50.mse <- (shen_50 - beta.positive)^2
    
    # n=300*3
    data <- data_gen(n_each = 300, num_site = num_site, lower_dl = lower_dl, upper_dl = upper_dl)
    data.obs <- data$dat_obs
    
    tian_300 <- tryCatch({
      multipleDL(y_obs ~ x, data = data.obs, delta_lower = data.obs$dl_lower,
                 delta_upper = data.obs$dl_upper, link = 'cloglog')$coef[['x']]
    }
    , error = function(e) return(NA))
    
    cai_300 <- tryCatch({
      cai_function(formula=Surv(left, right, type = 'interval2') ~ x, data = data.obs, r=0,
                   left_censoring = data.obs$lower_dl, right_censoring = data.obs$upper_dl)$coefficients
    }, error = function(e) return(NA))
    
    shen_300 <- tryCatch({
      shen_function(formula=Surv(left, right, type = 'interval2') ~ x, data = data.obs, r=0, left_censoring = data.obs$lower_dl)$coefficients
    }, error = function(e) return(NA))
    
    tian_300.mse <- (tian_300 - beta.positive)^2
    cai_300.mse <- (cai_300 - beta.positive)^2
    shen_300.mse <- (shen_300 - beta.positive)^2
    
    return(list(tian_50 = tian_50, tian_50.mse = tian_50.mse, 
                tian_300 = tian_300, tian_300.mse = tian_300.mse,
                shen_50 = shen_50, shen_50.mse = shen_50.mse, 
                cai_50 = cai_50, cai_50.mse = cai_50.mse,
                shen_300 = shen_300, shen_300.mse = shen_300.mse,
                cai_300 = cai_300, cai_300.mse = cai_300.mse))
  })
  
  save(res_all, file=paste('output/sim_',iter,'.Rdata', sep=""))
  
}


########## Simulation 3 ###########
# generate the true data and the observed data
data_gen <- function(n_each, num_site, lower_dl, upper_dl){
  
  site_name <- 1:num_site # name of site
  
  x <- rnorm(n_each * num_site, 0, 1) 
  e <-  rgumbel(n_each * num_site, max=F)
  y <- exp(x + e)
  # true data
  dat <- data.frame(site = rep(site_name, each = n_each),
                    y = y,
                    x = x)
  
  dat$lower_dl <- rep(lower_dl, each = n_each)
  dat$upper_dl <- rep(upper_dl, each = n_each)
  
  
  # the observed data
  dat_obs <- dat
  dat_obs$y_obs <- ifelse(dat$y < dat$lower_dl, dat$lower_dl,
                          ifelse(dat$y > dat$upper_dl, dat$upper_dl, dat$y))
  dat_obs$left <- ifelse(dat$y < dat$lower_dl, -Inf, 
                         ifelse(dat$y > dat$upper_dl, dat$upper_dl, dat$y))
  dat_obs$right <- ifelse(dat$y > dat$upper_dl, Inf, 
                          ifelse(dat$y < dat$lower_dl, dat$lower_dl, dat$y))
  
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


simulation <- function(iter){
  
  reps <- 1
  res_all <- lapply(1:reps, function(i){
    set.seed(i + iter * reps)
    
    lower_dl <- -c(Inf, Inf, Inf)
    upper_dl <- c(0.07, 0.16, 0.28)
    beta.positive <- 1
    num_site <- 3
    
    # n=50*3
    data <- data_gen(n_each = 50, num_site = num_site, lower_dl = lower_dl, upper_dl = upper_dl)
    data.obs <- data$dat_obs
    
    tian_50 <- tryCatch({
      multipleDL(y_obs ~ x, data = data.obs, delta_lower = data.obs$dl_lower,
                 delta_upper = data.obs$dl_upper, link = 'cloglog')$coef[['x']]
    }
    , error = function(e) return(NA))
    
    cai_50 <- tryCatch({
      cai_function(formula=Surv(left, right, type = 'interval2') ~ x, data = data.obs, r=0,
                   left_censoring = data.obs$lower_dl, right_censoring = data.obs$upper_dl)$coefficients
    }, error = function(e) return(NA))
    
    shen_50 <- tryCatch({
      shen_function(formula=Surv(left, right, type = 'interval2') ~ x, data = data.obs, r=0, 
                    left_censoring = data.obs$lower_dl)$coefficients
    }, error = function(e) return(NA))
    
    tian_50.mse <- (tian_50 - beta.positive)^2
    cai_50.mse <- (cai_50 - beta.positive)^2
    shen_50.mse <- (shen_50 - beta.positive)^2
    
    # n=300*3
    data <- data_gen(n_each = 300, num_site = num_site, lower_dl = lower_dl, upper_dl = upper_dl)
    data.obs <- data$dat_obs
    
    tian_300 <- tryCatch({
      multipleDL(y_obs ~ x, data = data.obs, delta_lower = data.obs$dl_lower,
                 delta_upper = data.obs$dl_upper, link = 'cloglog')$coef[['x']]
    }
    , error = function(e) return(NA))
    
    cai_300 <- tryCatch({
      cai_function(formula=Surv(left, right, type = 'interval2') ~ x, data = data.obs, r=0,
                   left_censoring = data.obs$lower_dl, right_censoring = data.obs$upper_dl)$coefficients
    }, error = function(e) return(NA))
    
    shen_300 <- tryCatch({
      shen_function(formula=Surv(left, right, type = 'interval2') ~ x, data = data.obs, r=0, left_censoring = data.obs$lower_dl)$coefficients
    }, error = function(e) return(NA))
    
    tian_300.mse <- (tian_300 - beta.positive)^2
    cai_300.mse <- (cai_300 - beta.positive)^2
    shen_300.mse <- (shen_300 - beta.positive)^2
    
    return(list(tian_50 = tian_50, tian_50.mse = tian_50.mse, 
                tian_300 = tian_300, tian_300.mse = tian_300.mse,
                shen_50 = shen_50, shen_50.mse = shen_50.mse, 
                cai_50 = cai_50, cai_50.mse = cai_50.mse,
                shen_300 = shen_300, shen_300.mse = shen_300.mse,
                cai_300 = cai_300, cai_300.mse = cai_300.mse))
  })
  
  save(res_all, file=paste('output/sim_',iter,'.Rdata', sep=""))
  
}


####### Simulation 4 #########
# generate the true data and the observed data
data_gen <- function(n_each, num_site, lower_dl, upper_dl){
  
  site_name <- 1:num_site # name of site
  
  x <- rnorm(n_each * num_site, 0, 1) 
  e <- rgumbel(n_each * num_site, max=F)
  y <- exp(x + e)
  # true data
  dat <- data.frame(site = rep(site_name, each = n_each),
                    y = y,
                    x = x)
  
  dat$lower_dl <- rep(lower_dl, each = n_each)
  dat$upper_dl <- rep(upper_dl, each = n_each)
  
  
  # the observed data
  dat_obs <- dat
  dat_obs$y_obs <- ifelse(dat$y < dat$lower_dl, dat$lower_dl,
                          ifelse(dat$y > dat$upper_dl, dat$upper_dl, dat$y))
  dat_obs$left <- ifelse(dat$y < dat$lower_dl, -Inf, 
                         ifelse(dat$y > dat$upper_dl, dat$upper_dl, dat$y))
  dat_obs$right <- ifelse(dat$y > dat$upper_dl, Inf, 
                          ifelse(dat$y < dat$lower_dl, dat$lower_dl, dat$y))
  
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


simulation <- function(iter){
  
  reps <- 1
  res_all <- lapply(1:reps, function(i){
    set.seed(i + iter * reps)
    
    lower_dl <- c(0.10, 0.16, -Inf)
    upper_dl <- c(Inf, 2.3, 2.6)
    beta.positive <- 1
    num_site <- 3
    
    # n=50*3
    data <- data_gen(n_each = 50, num_site = num_site, lower_dl = lower_dl, upper_dl = upper_dl)
    data.obs <- data$dat_obs
    
    tian_50 <- tryCatch({
      multipleDL(y_obs ~ x, data = data.obs, delta_lower = data.obs$dl_lower,
                 delta_upper = data.obs$dl_upper, link = 'cloglog')$coef[['x']]
    }
    , error = function(e) return(NA))
    
    cai_50 <- tryCatch({
      cai_function(formula=Surv(left, right, type = 'interval2') ~ x, data = data.obs, r=0,
                   left_censoring = data.obs$lower_dl, right_censoring = data.obs$upper_dl)$coefficients
    }, error = function(e) return(NA))
    
    shen_50 <- tryCatch({
      shen_function(formula=Surv(left, right, type = 'interval2') ~ x, data = data.obs, r=0, 
                    left_censoring = data.obs$lower_dl)$coefficients
    }, error = function(e) return(NA))
    
    tian_50.mse <- (tian_50 - beta.positive)^2
    cai_50.mse <- (cai_50 - beta.positive)^2
    shen_50.mse <- (shen_50 - beta.positive)^2
    
    # n=300*3
    data <- data_gen(n_each = 300, num_site = num_site, lower_dl = lower_dl, upper_dl = upper_dl)
    data.obs <- data$dat_obs
    
    tian_300 <- tryCatch({
      multipleDL(y_obs ~ x, data = data.obs, delta_lower = data.obs$dl_lower,
                 delta_upper = data.obs$dl_upper, link = 'cloglog')$coef[['x']]
    }
    , error = function(e) return(NA))
    
    cai_300 <- tryCatch({
      cai_function(formula=Surv(left, right, type = 'interval2') ~ x, data = data.obs, r=0,
                   left_censoring = data.obs$lower_dl, right_censoring = data.obs$upper_dl)$coefficients
    }, error = function(e) return(NA))
    
    shen_300 <- tryCatch({
      shen_function(formula=Surv(left, right, type = 'interval2') ~ x, data = data.obs, r=0, left_censoring = data.obs$lower_dl)$coefficients
    }, error = function(e) return(NA))
    
    tian_300.mse <- (tian_300 - beta.positive)^2
    cai_300.mse <- (cai_300 - beta.positive)^2
    shen_300.mse <- (shen_300 - beta.positive)^2
    
    return(list(tian_50 = tian_50, tian_50.mse = tian_50.mse, 
                tian_300 = tian_300, tian_300.mse = tian_300.mse,
                shen_50 = shen_50, shen_50.mse = shen_50.mse, 
                cai_50 = cai_50, cai_50.mse = cai_50.mse,
                shen_300 = shen_300, shen_300.mse = shen_300.mse,
                cai_300 = cai_300, cai_300.mse = cai_300.mse))
  })
  
  save(res_all, file=paste('output/sim_',iter,'.Rdata', sep=""))
  
}

########## Simulation 5 ############
data_gen <- function(n_each, num_site, lower_dl, upper_dl){
  
  site_name <- 1:num_site # name of site
  x1 <- rnorm(n_each, -0.5, 1) 
  y1 <- exp(x1 + rgumbel(n_each, max=F))
  x2 <- rnorm(n_each, 0, 1) 
  y2 <- exp(x2 + rgumbel(n_each, max=F))
  x3 <- rnorm(n_each, 0.5, 1) 
  y3 <- exp(x3 + rgumbel(n_each, max=F))
  
  # true data
  dat <- data.frame(site = rep(site_name, each = n_each),
                    y = c(y1, y2, y3),
                    x = c(x1, x2, x3))
  dat$lower_dl <- rep(lower_dl, each = n_each)
  dat$upper_dl <- rep(upper_dl, each = n_each)
  
  
  # the observed data
  dat_obs <- dat
  dat_obs$y_obs <- ifelse(dat$y < dat$lower_dl, dat$lower_dl,
                          ifelse(dat$y > dat$upper_dl, dat$upper_dl, dat$y))
  dat_obs$left <- ifelse(dat$y < dat$lower_dl, -Inf, 
                         ifelse(dat$y > dat$upper_dl, dat$upper_dl, dat$y))
  dat_obs$right <- ifelse(dat$y > dat$upper_dl, Inf, 
                          ifelse(dat$y < dat$lower_dl, dat$lower_dl, dat$y))

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


simulation <- function(iter){
  
  reps <- 1
  res_all <- lapply(1:reps, function(i){
    set.seed(i + iter * reps)
    
    lower_dl <- c(0.24, 0.65, 1.7)
    upper_dl <- c(Inf, Inf, Inf)
    beta.positive <- 1
    num_site <- 3
    
    # n=50*3
    data <- data_gen(n_each = 50, num_site = num_site, lower_dl = lower_dl, upper_dl = upper_dl)
    data.obs <- data$dat_obs
    
    tian_50 <- tryCatch({
      multipleDL(y_obs ~ x, data = data.obs, delta_lower = data.obs$dl_lower,
                 delta_upper = data.obs$dl_upper, link = 'cloglog')$coef[['x']]
    }
    , error = function(e) return(NA))
    
    cai_50 <- tryCatch({
      cai_function(formula=Surv(left, right, type = 'interval2') ~ x, data = data.obs, r=0,
                   left_censoring = data.obs$lower_dl, right_censoring = data.obs$upper_dl)$coefficients
    }, error = function(e) return(NA))
    
    shen_50 <- tryCatch({
      shen_function(formula=Surv(left, right, type = 'interval2') ~ x, data = data.obs, r=0, 
                    left_censoring = data.obs$lower_dl)$coefficients
    }, error = function(e) return(NA))
    
    tian_50.mse <- (tian_50 - beta.positive)^2
    cai_50.mse <- (cai_50 - beta.positive)^2
    shen_50.mse <- (shen_50 - beta.positive)^2
    
    # n=300*3
    data <- data_gen(n_each = 300, num_site = num_site, lower_dl = lower_dl, upper_dl = upper_dl)
    data.obs <- data$dat_obs
    
    tian_300 <- tryCatch({
      multipleDL(y_obs ~ x, data = data.obs, delta_lower = data.obs$dl_lower,
                 delta_upper = data.obs$dl_upper, link = 'cloglog')$coef[['x']]
    }
    , error = function(e) return(NA))
    
    cai_300 <- tryCatch({
      cai_function(formula=Surv(left, right, type = 'interval2') ~ x, data = data.obs, r=0,
                   left_censoring = data.obs$lower_dl, right_censoring = data.obs$upper_dl)$coefficients
    }, error = function(e) return(NA))
    
    shen_300 <- tryCatch({
      shen_function(formula=Surv(left, right, type = 'interval2') ~ x, data = data.obs, r=0, left_censoring = data.obs$lower_dl)$coefficients
    }, error = function(e) return(NA))
    
    tian_300.mse <- (tian_300 - beta.positive)^2
    cai_300.mse <- (cai_300 - beta.positive)^2
    shen_300.mse <- (shen_300 - beta.positive)^2
    
    return(list(tian_50 = tian_50, tian_50.mse = tian_50.mse, 
                tian_300 = tian_300, tian_300.mse = tian_300.mse,
                shen_50 = shen_50, shen_50.mse = shen_50.mse, 
                cai_50 = cai_50, cai_50.mse = cai_50.mse,
                shen_300 = shen_300, shen_300.mse = shen_300.mse,
                cai_300 = cai_300, cai_300.mse = cai_300.mse))
  })
  
  save(res_all, file=paste('output/sim_',iter,'.Rdata', sep=""))
  
}



