rm(list=ls())

# Set working directory----
main_path = getwd();main_path # project directory
setwd(main_path)

# Load functions
# In local PC
source("./codes/AR_Elas.R")

# Global Variable ----
#' @param n: sample size of the time series.
#' @param l: size of testing data (prediction step).
#' @param k: size of training data.
#' @param p: AR(p) process.
#' @param ar: ar coefficients with decreasing order.
#' @param pmax0: maximum p to be selected
#' @param ar: a vector of length p containing the AR(p) coefficients
#' @param N_simu: Monte Carlo simulation sample size.
#' @param gamma: Coefficient in elastic net.

n = 1000; p = 5 ; ar = c(p:1) / (2*p) * rep(c(1, -1), p)[1:p]; ar
N_simu = 1000; pmax0=pmax = 10;gamma = 1;l=3
k=20

set.seed(1000)
ar_sim = arima.sim(list(order=c(p,0,0), ar=ar, sd = sqrt(0.0025)), n=n)
data_input = ar_sim - mean(ar_sim) # Here we consider the ts with zero mean. Reason literature review.
data_input = ar_sim # Here we consider the ts with zero mean. Reason literature review.

# data_input = c(1:n)



# Try to call the function first
k = 20
AR_Elas(data_input, h = 7, k = k, gamma = gamma, pmax = pmax0)
AR_Elas_5(data_input, h = 7, k = k, gamma = gamma, pmax = pmax0)
AR_Lasso(data_input, h = 7, k = k, gamma = gamma, pmax = pmax0)
AR_True(data_input, h = 7, k = k, gamma = gamma, pmax = pmax0)


# Define functions ----
# getRMSEOnePred -----
#' @description Compute the RMSE(k) for k in {(h+1), (h+2), ..., (n-h)}
#' @param data_input: Simulated data or real data input into the model.
#' @param h: Steps ahead to predict.
#' @param model: The type of prediction model to be used. Input data_input
#' @param model(data_input, h = h, pmax = pmax0, k = k, gamma = gamma, )
getRMSEOnePred <- function(data_input, h = 7, model=AR_Elas){
  # Experiment
  # data_input= arima.sim(list(order=c(p,0,0), ar=ar), n=n);  model="ar"; k = h+1
  # Assume n/2 < n-h-pmax0 <=> n/2 > h+pmax0
  k_list = c((h+1):(n/2))
  RMSE_k <- numeric(length(k_list))
  
  getRMSEoneK<- function(k, h){
    model_fit = model(data_input, h = h, k = k, gamma = gamma, pmax = pmax0)
    RMSE = model_fit$RMSE
    cat("==============================================================\n")
    cat("When k=", k, " and h=",h ,", the fitted model is:\n", sep="")
    print(model_fit$coef)
    cat("  while the order selected is", model_fit$order,".\n")
    return(RMSE)
  }
  
  getRMSEoneKV = Vectorize(getRMSEoneK, vectorize.args = "k")
  RMSE_k = getRMSEoneKV(k_list, h)
  index = which.min(RMSE_k)
  k_opt = k_list[index] # --> k=50
  # plot(k_list, RMSE_k, type="h", xlab="k")
  result = list(RMSE = RMSE_k, k_opt=k_opt, k_list = k_list)
  return(result)
}

# oneTraj = getRMSEOnePred(data_input = data_input, h = 7)

# oneSimuAR -----
#' @description A simulation with given h and N_simu
#' @param N_simu: # of replicates for a simulation.
#' @param h: Steps ahead to predict.
#' @export version_log.txt: Log of the fitted models for each replication .
#' @export version.RDS: A list with simulated result for further analysis.
oneSimuAR <- function(simun=1, h=7, method=AR_Elas){
  # Experiment
  # N_simu=10; h=7
  cat("==============================================================\n")
  cat(simun, "th simulation:\n")
  set.seed(simun)
  ar_sim = arima.sim(list(order=c(p,0,0), ar=ar), n=n)
  pred = getRMSEOnePred(data_input=ar_sim, h, model=method)
  # result[[simun]] = pred
  # cat("==============================================================\n")
  # cat("The summary Statistics:\n")
  return(pred)
}

# a = oneSimuAR(1, 7)
# Run the Monte Carlo Simulation ----
library(doParallel)
library(foreach)
# library(doSNOW)
# library(Rmpi)
# cl<-makeCluster(detectCores(),type="SOCK")
# registerDoSNOW(cl)
N_simu = 10
n_simu_list = c(1:N_simu)
n_cores = detectCores();n_cores

# Customize in script. 
method_name = AR_True

# Run the simulation in parallel
n_cores = 4; h=3
version_name = paste("ar",p,"_h",h,"_N",N_simu, sep="")
simu.log.path = paste(main_path, "/simu/", version_name, "_log.txt", sep = "")

cl <- parallel::makeCluster(n_cores, outfile = simu.log.path)
registerDoParallel(cl)
start_time <- Sys.time()
final_result <- foreach(s = n_simu_list, .packages=c("glmnet")) %dopar% {
  oneSimuAR(s, h=h, method = method_name)
}
names(final_result) = paste("Simu",n_simu_list, sep="")
end_time <- Sys.time()
runtime = end_time - start_time; runtime
stopCluster(cl)

saveRDS(final_result, paste("./simu/",version_name,".RDS", sep=""))












