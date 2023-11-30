# Load packages
?library(glmnet)

#' @param n: length of historical data
#' @param h: length of forecasting horizon (size of validation set)
#' @param pm: >1, maximum number of column used to build the design matrix. Upper bound of the model order p. Size of the sliding window for a training set.
#' @param k in [2pm+2, ..., n], size of the development set
#' @param gamma suggested to be 1 as tuned parameter in Adaptive Lasso.
#' @param alpha.cv = TRUE if CV is applied to select the alpha, otherwise, alpha = alpha.default.
#' @param alpha.default = 1 is L-1 LASSO penalty; = 0 if for L-2 Ridge penalty. 
AR_Elas <- function(data_input, n=100, h = 3, pm = 10, k = 20, gamma = 1, alpha.cv = TRUE, alpha.default = 1){ 
  # # Manually attemp.
  # data_input = seq(1, 1200);
  # n = 1000; p = 5 ; ar = c(p:1) / (2*p) * rep(c(1, -1), p)[1:p]; ar
  # N_simu = 1000; pm0=pm = 10;gamma = 1;l=3
  # data_input = arima.sim(list(order=c(5,0,0), ar=ar), n=1200, sd=1)
  # 
  # n=1000;h=3;pm=10;k=10;gamma=1;alpha.cv = TRUE;alpha.default = 1
  
  
  
  # Initialization.
  # Let the maiximum model order (pm + 1) be less than the sample size k. Otherwise overfitting.
  
  km = k - pm
  # Check the region of k .
  if(k < 2 * pm + 1){
    stop(paste("Please input the development size k no less than", 2*pm+1))
  }
  
  # Check pm.
  if(pm < 2){
    stop("Please input the pm no less than 2.")
  }
 
  if(k > (n)){
  # y_train_all <- data_input[c((pm+1):(n-l))]
  # (n-l-k+1) > (pm+1) => k <= n-l-pm
    stop(paste("Please try a smaller development size k that is less than", (n-1)))
  }
  
  # Convert time series data into matrix form with most recent k observations.
  Z = matrix(0, nrow=km+1, ncol=pm)
  for (i in 1:(km+1)){
    for (j in 1:pm){
      Z[i,j] = data_input[n-k+pm+i-j]
    }
  }
  colnames(Z) = c(paste("lag", 1:pm, sep="")); colnames(Z)
  Z_train <- Z[1:km,]
  Z_test = matrix( Z[km+1,], nrow=1)
  # all.equal(km, nrow(Z))
  y_train <- data_input[c((n-k+pm+1):(n))] # length = km = k - pm
  y_test <- data_input[c((n+1):(n+h))] # length = h
  
  # Estimate the adaptive weight using ridge
  
  # ols
  # fit.ols  <- lm(y_train ~ Z_train-1) # Z_train has no extra intercept.
  # summary(fit.ols)
  # coef_ols <- fit.ols$coefficients
  
  fit.ridge.cv <- cv.glmnet(Z_train, y_train, alpha=0, intercept = FALSE, nfolds = 3, grouped=FALSE)
  lambda_ridge = fit.ridge.cv$lambda.min
  fit.ridge <- glmnet(Z_train, y_train, alpha=0, lambda=lambda_ridge, intercept = FALSE, nfolds = 3, grouped=FALSE)
  coef_ridge = coef(fit.ridge)[-1]
  
  # # Estimate the adaptive weight using elasticnet
  # fit.elasini.cv <- cv.glmnet(Z_train, y_train, alpha=0, intercept = FALSE, nfolds = 3, grouped=FALSE)
  # lambda_ridge = fit.ridge.cv$lambda.min
  # fit.ridge <- glmnet(Z_train, y_train, alpha=0, lambda=lambda_ridge, intercept = FALSE, nfolds = 3, grouped=FALSE)
  # coef_ridge = coef(fit.ridge)[-1]
  
  if (alpha.cv){
    alphalist <- seq(0.5, 1, by=0.05) # This range is defined to deal with sparsity.
    elasticnet <- lapply(alphalist, function(a){
      cv.glmnet(Z_train, y_train,  alpha = a, intercept = FALSE, grouped=FALSE) # , lambda = lambdalist
    })
    cv.err <- (sapply(elasticnet, "[",'cvm'))
    min.index = which.min(sapply(cv.err, min))
    alpha_cv = alphalist[[min.index]]
  } else {
    alpha_cv = alpha.default
  }

  # Define adaptive weights
  weight_ridge = (1/abs(coef_ridge))^gamma
  
  # Sorted Adaptive weight -> cut ridge into the decreasing order
  # Define sorting strategy in sortb()
  b = abs(coef_ridge)
  coef_ridge_sort = sortb(b)
  weight_ridge_sort = (coef_ridge_sort+1/k)^(-gamma)

  # Select lambda cv of elasticnet
  fit_elase_cv = cv.glmnet(Z_train, y_train, alpha = alpha_cv, intercept = FALSE, grouped=FALSE) # , lambda = lambdalist
  lambda_elas_cv = fit_elase_cv$lambda.min
  # Fit adaptive elasticm net with selected lambda
  fit.adaelas = glmnet(Z_train, y_train, alpha = alpha_cv, intercept = FALSE, grouped=FALSE,
                    penalty.factor = weight_ridge_sort, lambda = lambda_elas_cv)
  # Define the order as the location of the last non-zero element.
  order = 1
  for (i in 1:pm){
    if ((fit.adaelas$beta)[pm-i+1] != 0){
      order = pm-i+1
      break
    }
  }
  multiplier = (1+lambda_elas_cv*(1-alpha_cv)/(2*k))
  coef_adaelas = multiplier * fit.adaelas$beta[1:order]
  # This is not for fitting the model, but to use the sliding window forecast strategy built-in within ar()
  mod = ar(data_input[(n-k+1):n], aic = FALSE, order.max = order, demean=FALSE) 
  mod$ar = coef_adaelas # We change the coefficients into adaptive elastic net estimators or whatever other \hat{\phi}
  y_hat = predict(mod, n.ahead = h, ci = 0.95)$pred
  # mod$x.mean == 0 # demean=FALSE
  res = ( y_hat - y_test )
  MSE = (mean( res^2 ) )
  result <- list(MSE = MSE, order = order, coef = coef_adaelas, paracv = c(alpha=alpha_cv, lambda=lambda_elas_cv), yhat=y_hat)
  return(result)
}

#' @param alpha.cv = TRUE if CV is applied to select the alpha, otherwise, alpha = alpha.default.
#' @param alpha.default = 1 is L-1 LASSO penalty; = 0 if for L-2 Ridge penalty. 
AR_Lasso <- function(data_input, n=100, k=80, h = 7, pm = 10, gamma = 1){
  AR_Elas(data_input = data_input, n=n, h = h, pm = pm, k = k, gamma = gamma, alpha.cv = FALSE, alpha.default = 1)
}

#' @param alpha.cv = TRUE if CV is applied to select the alpha, otherwise, alpha = alpha.default.
#' @param alpha.default = 1 is L-1 LASSO penalty; = 0 if for L-2 Ridge penalty. 
AR_Elas_5 <- function(data_input, n=100, k=80, h = 7, pm = 10, gamma = 1){
  AR_Elas(data_input = data_input, n=n, h = h, pm = pm, k = k, gamma = gamma, alpha.cv = FALSE, alpha.default = 0.5)
}

AR_True <- function(data_input,  n=100, k=80, h = 3, pm = 10, p = 5, gamma = 1, Fuchini=FALSE , no=TRUE){
  # set.seed(1)
  # n = length(data_input)
  # k=n-h;k
  if(Fuchini==TRUE){
    MSE_list = numeric(n-h-k+1)
    if(no==TRUE){
      # Fuchini subsmapling -> non-overlapping
      l_list = seq(from=1, by=k+h, length.out=floor(n/(k+h)))
    } else{
      # Fuchini subsampling -> overlapping
      l_list = seq(n-h-k+1)
    }
  } else{
    MSE_list = numeric()
    l_list = n-k+1
  }
  
  order.true = min(k-1, p, pm)
  for(i in 1:(length(l_list))){
    l=l_list[i]
    data_train = ts( data_input[l:(l+k-1)] ) # length = k
    data_test = data_input[(l+k):(l+k+h-1)] # length = h
    mod = ar(data_train, aic = FALSE, order.max = order.true) # k=4, p=5, order should < # observations
    coef = mod$ar
    est = as.numeric(predict(mod, n.ahead = h)$pred)
    res = ( est - data_test )
    MSE_list[i] = (mean( res^2 ) )
  }
  names(MSE_list) <- paste("l=",l_list,sep="")
  MSE_tail = tail(MSE_list,1)
  MSE_mean = mean(MSE_list)
  result <- list(MSE = MSE_mean, MSE_list=MSE_list, coef=coef,  order = order.true)
  return(result)
}


ARpre_True <- function(data_input, n=100, k=80, h = 7, pm = 10, p = 5, gamma = 1){
  # set.seed(1)
  n = length(data_input)
  # k=n-l;k
  data_train = ts( data_input[(n-h-k+1):(n-h)] ) # length = k
  data_test = data_input[(n-h+1):(n)] # length = l
  order.true = min(k-1, p, pm)
  mod = ar(ts(data_train), aic = FALSE, order.max = order.true, demean=FALSE) # k=4, p=5, order should < # observations
  est = predict(mod, n.ahead = h, ci = 0.95)$pred
  

  # (mod$ar%*%tail(data_train,order.true)[order.true:1]) == est[1]
  
  if (h>1){
    a = mod$ar
    tail = tail(data_train,order.true)[order.true:1]
    for (step in 2:h){
      xinv = c(est[step-1], tail[-length(tail)])
      est[step] = a%*% xinv
      tail = xinv
    }
  }
  res = ( as.numeric( est) - data_test )
  MSE = (mean( res^2 ) )
  result <- list(RMSE = RMSE, coef = mod$ar, order = mod$order)
  return(result)
}



sortb <- function(b){
  bcopy = b
  nb = length(b)
  index1 = 1
  pointer = index1
  index2 = 0
  count = 0 
  while(pointer < nb){
    pointer = index1
    upper = b[index1]
    # Check first point that is less than max
    for( i in 2:nb){
      if(b[i]< upper){
        index2 = i
        lower = b[i]
        # print(index2)
        break
      }
    }
    if(index2 == 0 ){
      bcopy = rep( max(b[1], 0), nb)
      return(bcopy)
    }
    if( (index2==index1) && (index2 < nb) ){
      for (i in (index2+1):(nb)){
        slope = (b[index2] - b[1])/(index2-1)
        bcopy[i] = bcopy[index2] + slope * i
      }
      break
    }
    # if b is only increasing then set all afterwards beta into zero
    # change interpolating values
    if((index2-index1) > 1){
      for (i in 1:(index2-index1-1)){
        slope = (b[index2] - b[index1])/(index2-index1)
        bcopy[index1+i] = bcopy[index1] + slope * i
      }
    } 
    index1 = index2
  } 
  
  for(i in 1:nb){
    if(bcopy[i]<0){
      bcopy[i]=0
    }
  }
  
  return(bcopy)
  
}

model_fullname <- function(model_index){ 
  result = switch(model_index,
                  "1" = "1. Oracle Model",
                  "2" = "2. Lasso Model",
                  "3" = "3. Fixed Elastic",
                  "4" = "4. Tuned Elastic")
  return(result)
}

method_namef <- function(model_index){ 
  result = switch(model_index,
                  "1" = AR_True,
                  "2" = AR_Lasso,
                  "3" = AR_Elas_5,
                  "4" = AR_Elas)
  return(result)
}

data_fullname <- function(model_index){ 
  result = switch(model_index,
                  "1" = "Amazon",
                  "2" = "Google",
                  "3" = "Facebook",
                  "4" = "Apple",
                  "5" = "Ozone")
  return(result)
}

data_namef <- function(model_index){ 
  result = switch(model_index,
                  "1" = "AMZN",
                  "2" = "GOOG",
                  "3" = "META",
                  "4" = "AAPL",
                  "5" = "ozone")
  return(result)
}









