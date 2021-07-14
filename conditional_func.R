# functions for conditional CDF and conditional quantile 


### CDF
cdf.orm <- function(mod, new.data, at.y=0,se=TRUE){
  if(is.null(mod$yunique)) {
    stop("Need to set x=TRUE and y=TRUE for orm") 
  } else{
    order.y <- mod$yunique
    xb <- as.matrix(new.data)%*%matrix(coef(mod)[colnames(new.data)])
    
    index <- sapply(at.y, FUN=function(x) {if(x<min(order.y)[1]) result <- Inf 
    else if (x==min(order.y)[1]) result <- 1
    else if(x >= max(order.y)[1]) result <- -Inf
    else which(order.y>=x)[1]-1})
    
    m.alpha <- mod$coef[index]
    m.alpha <- ifelse(is.infinite(index), index, m.alpha)
    if(length(at.y)==1){
      lb <- as.matrix(outer(m.alpha, xb, "+")[,,1])
    } else lb <- t(outer(m.alpha, xb, "+")[,,1])
    m.cdf <- 1- mod$trans$cumprob(lb)
    
    
    if(se){
      if(mod$family=="logistic") mod$trans$deriv <- function(x) exp(-x)/(1+exp(-x))^2
      if(mod$family=='probit') mod$trans$deriv <- function(x) dnorm(x)
      if(mod$family=='loglog') mod$trans$deriv <- function(x) exp(-x - exp(-x))
      if(mod$family=='cloglog') mod$trans$deriv <- function(x) exp(x - exp(x))
      
      cdf.se <- matrix(NA, ncol=length(at.y), nrow=dim(new.data)[1])
      lb.se <- matrix(NA, ncol=length(at.y), nrow=dim(new.data)[1])
      
      var <- as.matrix(solve(mod$info.matrix))
      
      for(i in 1:length(at.y)) {
        
        var.i <- var[c(index[i], which(names(coef(mod)) %in% colnames(new.data))), 
                     c(index[i], which(names(coef(mod)) %in% colnames(new.data)))]
        dcdf.dtheta <- cbind(-mod$trans$deriv(lb[,i]),  
                             -mod$trans$deriv(lb[,i])*as.matrix(new.data) )
        dlb.dtheta <- as.matrix(cbind(1, new.data))
        cdf.se[,i] <- sqrt(diag(dcdf.dtheta %*% var.i%*% t(dcdf.dtheta)))
        lb.se[, i] <- sqrt(diag(dlb.dtheta%*%var.i%*% t(dlb.dtheta)))
        
      }
      ci.lb <- sapply(1:length(at.y), FUN=function(i) 
      { 1- mod$trans$cumprob(lb[, i] +qnorm(0.975)*lb.se[, i])})
      ci.ub <- sapply(1:length(at.y), FUN=function(i) 
      { 1- mod$trans$cumprob(lb[, i] -qnorm(0.975)*lb.se[, i])})
      
      
      result <- list(est=m.cdf,
                     se=cdf.se,
                     lb=ci.lb,
                     ub=ci.ub)
    } else{
      result <- list(est=m.cdf)
    }
    
    
    return(result)
    
  } 
  
}


### Quantile
quantile.orm <- function(mod, new.data, probs=0.5, se=TRUE){
  
  quantile <- matrix(NA, nrow=dim(new.data)[1], ncol=length(probs))
  order.y <- mod$yunique
  # Qi's (left)
  order.y_left <- c(order.y, order.y[length(order.y)])
  # Frank's (right)
  order.y_right <- c(order.y[1], order.y)
  
  #n.alpha <- length(order.y)-1
  xb <- as.matrix(new.data)%*%matrix(coef(mod)[colnames(new.data)])
  alpha <- mod$coef[1:(length(unique(order.y))-1)]
  lb <- t(outer(alpha, xb, "+")[,,1])
  m.cdf <- 1 - mod$trans$cumprob(lb)
  m.cdf <- cbind(0, m.cdf, 1)
  for(i in 1: length(probs)){
    try({
      index.1 <- apply(m.cdf, 1, FUN=function(x){ max(which(x<=probs[i]))[1]} )
      index.2 <- apply(m.cdf, 1, FUN=function(x){ min(which(x>=probs[i]))[1]} )
      
      # index - Qi's
      index.y1_left <- ifelse(index.1>length(order.y_left), Inf, order.y_left[index.1])
      index.y2_left <- ifelse(index.2>length(order.y_left), Inf, order.y_left[index.2])
      # index - Frank's 
      index.y1_right <- ifelse(index.1>length(order.y_right), Inf, order.y_right[index.1])
      index.y2_right <- ifelse(index.2>length(order.y_right), Inf, order.y_right[index.2])
      
      # CDF - Qi's
      index.y1.cdf_left <- ifelse(index.1 == 0, 0, m.cdf[cbind(1:dim(new.data)[1], index.1)])
      index.y2.cdf_left <- ifelse(index.2 > length(order.y_left), 1, m.cdf[cbind(1:dim(new.data)[1], index.2)])
      # CDF - Frank's
      index.y1.cdf_right <- ifelse(index.1 == 0, 0, m.cdf[cbind(1:dim(new.data)[1], index.1)])
      index.y2.cdf_right <- ifelse(index.2 > length(order.y_right), 1, m.cdf[cbind(1:dim(new.data)[1], index.2)])
      
      # quantile - Qi's
      quantile_left <- ifelse(index.1==index.2, index.y1_left, 
                              (index.y2_left - index.y1_left)/(index.y2.cdf_left - index.y1.cdf_left)*
                                (probs[i] - index.y1.cdf_left) + index.y1_left) 
      quantile_left <- ifelse(is.na(quantile_left), max(order.y), quantile_left)
      # quantile - Frank's
      quantile_right <- ifelse(index.1==index.2, index.y1_right, 
                               (index.y2_right - index.y1_right)/(index.y2.cdf_right - index.y1.cdf_right)*
                                 (probs[i] - index.y1.cdf_right) + index.y1_right) 
      quantile_right <- ifelse(is.na(quantile_right), max(order.y), quantile_right)
      
      # weight
      all_weights <- c(0, (m.cdf[-length(m.cdf)] - m.cdf[1]) / (m.cdf[length(m.cdf) - 1] - m.cdf[1]), 1)
      weight <- all_weights[index.1]
      # weighted average 
      quantile[,i] <- weight * quantile_left + (1 - weight) * quantile_right
      
    })
    
  }
  result <- quantile
  
  if(se){
    if(mod$family=="logistic") mod$trans$deriv <- function(x) exp(-x)/(1+exp(-x))^2
    if(mod$family=='probit') mod$trans$deriv <- function(x) dnorm(x)
    if(mod$family=='loglog') mod$trans$deriv <- function(x) exp(-x - exp(-x))
    if(mod$family=='cloglog') mod$trans$deriv <- function(x) exp(x - exp(x))
    
    quantile.lb <- quantile.ub <- matrix(NA, nrow=dim(new.data)[1], ncol=length(probs))
    lb.se <- matrix(NA, ncol=dim(lb)[2], nrow=dim(new.data)[1])
    var <- as.matrix(solve(mod$info.matrix))
    
    for(i in 1:dim(lb)[2]){
      var.i <- var[c(i, which(names(coef(mod)) %in% colnames(new.data))), 
                   c(i, which(names(coef(mod)) %in% colnames(new.data)))]
      
      dcdf.dtheta <- cbind(-mod$trans$deriv(lb[,i]),  
                           -mod$trans$deriv(lb[,i])*as.matrix(new.data) )
      dlb.dtheta <- as.matrix(cbind(1, new.data))
      lb.se[,i] <- sqrt(diag(dlb.dtheta%*%var.i%*% t(dlb.dtheta)))
    }
    
    ci.lb <- sapply(1:dim(lb)[2], FUN=function(i) 
    { 1- mod$trans$cumprob(lb[, i] +qnorm(0.975)*lb.se[, i])})
    ci.ub <- sapply(1:dim(lb)[2], FUN=function(i) 
    { 1- mod$trans$cumprob(lb[, i] -qnorm(0.975)*lb.se[, i])})
    ci.lb <- matrix(ci.lb, nrow=dim(new.data)[1])
    ci.ub <- matrix(ci.ub, nrow=dim(new.data)[1])
    
    ci.lb <- cbind(0, ci.lb, 1)
    ci.ub <- cbind(0, ci.ub, 1)
    
    for(i in 1: length(probs)){
      try({
        #### lower bound ####
        index.1 <- apply(ci.lb, 1, FUN=function(x){ max(which(x<=probs[i]))[1]} )
        index.2 <- apply(ci.lb, 1, FUN=function(x){ min(which(x>=probs[i]))[1]} )
        
        # index - Qi's
        index.y1_left <- ifelse(index.1>length(order.y_left), Inf, order.y_left[index.1])
        index.y2_left <- ifelse(index.2>length(order.y_left), Inf, order.y_left[index.2])
        # index - Frank's 
        index.y1_right <- ifelse(index.1>length(order.y_right), Inf, order.y_right[index.1])
        index.y2_right <- ifelse(index.2>length(order.y_right), Inf, order.y_right[index.2])
        
        # CDF - Qi's
        index.y1.cdf_left <- ifelse(index.1 == 0, 0, ci.lb[cbind(1:dim(new.data)[1], index.1)])
        index.y2.cdf_left <- ifelse(index.2 > length(order.y_left), 1, ci.lb[cbind(1:dim(new.data)[1], index.2)])
        # CDF - Frank's
        index.y1.cdf_right <- ifelse(index.1 == 0, 0, ci.lb[cbind(1:dim(new.data)[1], index.1)])
        index.y2.cdf_right <- ifelse(index.2 > length(order.y_right), 1, ci.lb[cbind(1:dim(new.data)[1], index.2)])
        
        # quantile - Qi's
        quantile_left <- ifelse(index.1==index.2, index.y1_left, 
                                (index.y2_left - index.y1_left) / (index.y2.cdf_left - index.y1.cdf_left) *
                                  (probs[i] - index.y1.cdf_left) + index.y1_left) 
        quantile_left <- ifelse(is.infinite(quantile_left), max(order.y), quantile_left)
        # quantile - Frank's
        quantile_right <- ifelse(index.1==index.2, index.y1_right, 
                                 (index.y2_right - index.y1_right)/(index.y2.cdf_right - index.y1.cdf_right) *
                                   (probs[i] - index.y1.cdf_right) + index.y1_right) 
        quantile_right <- ifelse(is.infinite(quantile_right), max(order.y), quantile_right)
        
        # weight
        all_weights <- c(0, (ci.lb[-length(ci.lb)] - ci.lb[1]) / (ci.lb[length(ci.lb) - 1] - ci.lb[1]), 1)
        weight <- all_weights[index.1]
        # weighted average 
        quantile.lb[,i] <- weight * quantile_left + (1 - weight) * quantile_right
        
        
        #### upper bound ####
        index.1 <- apply(ci.ub, 1, FUN=function(x){ max(which(x<=probs[i]))[1]} )
        index.2 <- apply(ci.ub, 1, FUN=function(x){ min(which(x>=probs[i]))[1]} )
        
        # index - Qi's
        index.y1_left <- ifelse(index.1>length(order.y_left), Inf, order.y_left[index.1])
        index.y2_left <- ifelse(index.2>length(order.y_left), Inf, order.y_left[index.2])
        # index - Frank's 
        index.y1_right <- ifelse(index.1>length(order.y_right), Inf, order.y_right[index.1])
        index.y2_right <- ifelse(index.2>length(order.y_right), Inf, order.y_right[index.2])
        
        # CDF - Qi's
        index.y1.cdf_left <- ifelse(index.1 == 0, 0, ci.ub[cbind(1:dim(new.data)[1], index.1)])
        index.y2.cdf_left <- ifelse(index.2 > length(order.y_left), 1, ci.ub[cbind(1:dim(new.data)[1], index.2)])
        # CDF - Frank's
        index.y1.cdf_right <- ifelse(index.1 == 0, 0, ci.ub[cbind(1:dim(new.data)[1], index.1)])
        index.y2.cdf_right <- ifelse(index.2 > length(order.y_right), 1, ci.ub[cbind(1:dim(new.data)[1], index.2)])
        
        # quantile - Qi's
        quantile_left <- ifelse(index.1==index.2, index.y1_left, 
                                (index.y2_left - index.y1_left) / (index.y2.cdf_left - index.y1.cdf_left) *
                                  (probs[i] - index.y1.cdf_left) + index.y1_left) 
        quantile_left <- ifelse(is.infinite(quantile_left), max(order.y), quantile_left)
        # quantile - Frank's
        quantile_right <- ifelse(index.1==index.2, index.y1_right, 
                                 (index.y2_right - index.y1_right)/(index.y2.cdf_right - index.y1.cdf_right) *
                                   (probs[i] - index.y1.cdf_right) + index.y1_right) 
        quantile_right <- ifelse(is.infinite(quantile_right), max(order.y), quantile_right)
        
        # weight
        all_weights <- c(0, (ci.ub[-length(ci.ub)] - ci.ub[1]) / (ci.ub[length(ci.ub) - 1] - ci.ub[1]), 1)
        weight <- all_weights[index.1]
        # weighted average 
        quantile.ub[,i] <- weight * quantile_left + (1 - weight) * quantile_right
        
      })
      
    }
    
    result <- list(quantile=quantile,
                   lb=quantile.ub,
                   ub=quantile.lb)
  }
  return(result)
}
