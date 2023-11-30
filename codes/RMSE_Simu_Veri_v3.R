rm(list=ls())

# Set working directory----
main_path = getwd();main_path # project directory
setwd(main_path)

# Global Variable ----
#' @param n: size of the time series samples=
#' @param p: AR(p) process
#' @param ar: a vector of length p containing the AR(p) coefficients
#' @param N_simu: a variable identifying the number of replicates
#' @param k: k=80 past k observations. k=80 for the finite sample property
n = 1000; p = 5; ar = c(p:1) / (2*p) * rep(c(1, -1), p)[1:p]; ar
n.total = as.integer(1.2*n)
N_simu = 1000; l=7; k = 100; gamma = 1; pmax=10
sigma=1
# In local PC
# source("./codes/AR_Elas.R")
# In HPC Server
source("./codes/AR_Elas.R")
set.seed(1)


# # Exporation of the variance
# ar_sim = arima.sim(list(order=c(p,0,0), ar=ar), n=n.total, sd=sigma)
# mod = ar(ar_sim, aic = FALSE, order.max = p, method = c("yw"))
# mod$var.pred # predicted sigma.hat.2 supposed to be 1
# 
ar(ar_sim, aic = FALSE, order.max = p, method = c("yule-walker")) # "yule"="yw"
ar(ar_sim, aic = FALSE, order.max = p, method = c("yw"))
ar(ar_sim, aic = FALSE, order.max = p, method = c("burg"))
ar(ar_sim, aic = FALSE, order.max = p, method = c("ols"))

#' Compare SSE for past k observations with prediction horizontal l
compareSSE <- function(l=7){
  # Get a longer data n+h
  ar_sim = arima.sim(list(order=c(p,0,0), ar=ar), n=n.total, sd=sigma)
  data_input = ar_sim
  data_train = ts( data_input[(n-k+1):n] ) # length = k
  data_test = data_input[(n+1):(n+l)] # length = l
  # demean=F: a mean be estimated during fitting? Nope. Zero mean
  mod = ar(data_train, aic = FALSE, order.max = p, demean=F, method = c("yule-walker")) # k=4, p=5, order should < # observations
  # true coef
  coef_true = ar; coef_true
  coef_yw = mod$ar; coef_yw
  est = predict(mod, n.ahead = l, method="yule", se.fit = F); est
  est = as.numeric(est)
  res = ( as.numeric( est) - data_test )
  simu_SSE = sum(res^2)
  R = ar - mod$ar
  
  # #' prediction function
  # xkl <- function(ll=l, arc = coef_yw){
  #   pred = numeric()
  #   dat = data_train[k:(k-p+1)] # notice the reverse sequence
  #   for (i in 1:ll){
  #     pred[i] = arc %*% dat
  #     dat = c(pred[i], dat[-p])
  #   }
  #   return(pred)
  # }
  # est_theo = xkl(l=l, arc = coef_yw);est_theo
  # est_simu = est;est_simu
  # all.equal(est_theo, est_simu, tol=.2)

  
  #' Iterative function
  #' @param phi estimated phik that depends on the previous k observations
  #' @param truep the true order p of AR model
  phik <- function(h, truep = p, phi=mod$ar){
    if (h <= 0){
      return(0)
    } else if (h == 1) {
      return(1)
    } else{
      temp = 0
      for(i in 1:truep){
        temp = temp + phi[i] * phik(h-i) 
      }
    }
    return(temp)
  }
  phikvalue = sapply(c(1:l), phik)
  

  # Define matrix F
  FlMat <- matrix(0, nrow = l, ncol = l)
  fl = sapply(c(1:l), phik)
  for(i in 1:l){
    FlMat[i,i:l] =  fl[1:(l-i+1)]
  }
  
  #' Define function e_phi^(k)(h)
  #' @param h interger
  ephi <- function(h, pp=p){
    if (h <= 0){
      return(0)
    } else{
      x = data_input[(n-l+h-1):(n-l+h-pp)]
      temp = R %*% x
      return(temp)
    }
  }
  el = sapply(c(1:l), ephi)
  
  ekt <- function(t){
    temp = 0
    for (i in 1:t){
      temp = temp + phik(i) * (ephi(t+1-i) + rnorm(1))
    }
    return(temp)
  }
  ek = sapply(c(1:l), ekt)
  
  # This step can be revised to improve matrix multiplication time O(l^2)
  # all.equal(3^2+4^2, norm(c(3,4), type="2")^2)
  left = norm(x = t(FlMat) %*% phikvalue, type="2");left
  lowertri = lower.tri(matrix(1, nrow = l, ncol = l), diag = 1)
  right = norm(x = lowertri %*% phikvalue, type="2")
  theo_SSE = left^2 + right^2 * mod$var.pred^2
  
  diff = (theo_SSE - simu_SSE)
  # theo_RMSE = sqrt(theo_SSE/l);  simu_RMSE = sqrt(simu_SSE/l)
  return(list(theo_SSE=theo_SSE, simu_SSE=simu_SSE, diff=diff))
}

sapply(c(3, 5, 7, 10), compareSSE)

# Run the Monte Carlo Simulation ----
library(doParallel)
library(foreach)
set.seed(1)
N_simu = 1000
n_simu_list = c(1:N_simu)
n_cores = detectCores();n_cores
n_cores = 5 
h=7
cl <- makeCluster(n_cores)
registerDoParallel(cl)
start_time <- Sys.time()
final_result <- foreach(N_simu = n_simu_list, .packages=c("glmnet")) %dopar% {
  compareSSE(h)
}
names(final_result) = paste("Simu",n_simu_list, sep="")
end_time <- Sys.time()
end_time - start_time
stopCluster(cl)


library(ggplot2)
library(ggpubr)
library(plyr)
simu <- rep(1:N_simu, 2)
theo_SSE = as.numeric(unlist( sapply(final_result, '[', 1)  ))
simu_SSE = as.numeric(unlist( sapply(final_result, '[', 2)  ))
SSE = c(theo_SSE, simu_SSE)
label=c(rep(c("approx","empirical"), each=N_simu))
data = data.frame(s=simu, SSE, label)
mu <- ddply(data, "label", summarise, grp.mean=mean(SSE))
g0 <- ggplot(data, aes(x = SSE, color=label, linetype=label)) +
  geom_density( size=1, alpha=0.9) + # color="#69b3a2"
  geom_vline(data=mu, aes(xintercept=grp.mean, color=label, linetype=label)) + 
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  ggtitle("Distribution of SSE") +
  scale_color_manual(values=c("#E69F00", "#999999",  "#56B4E9"))

g1 <- ggplot(data, aes(y = SSE, x=label, color=label, linetype=label)) +
  geom_boxplot( size=1, alpha=0.9) + # color="#69b3a2"
  geom_hline(data=mu, aes(yintercept=grp.mean, color=label, linetype=label)) + 
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  )  + coord_flip() +
  ggtitle(paste("Distribution of SSE when h=", l, sep="")) +
  labs(x="") +
  scale_color_manual(values=c("#E69F00", "#999999",  "#56B4E9"))

g2 <- ggplot(data, aes(x = simu, y=SSE, color=label, linetype=label), size=.1, alpha=0.9) +
  geom_line( ) + # color="#69b3a2"
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  ggtitle("SSE for each simulation") +
  scale_color_manual(values=c("#E69F00", "#999999",  "#56B4E9"))

g = ggarrange(g1, g2, nrow = 2, 
              legend = "right",
              label.x = 0,
              label.y = 0,
              font.label = list(size = 15, face = "bold"),
              common.legend = T, 
              align = "v")
g1
g

mean(theo_SSE);sd(theo_SSE)
mean(simu_SSE);sd(simu_SSE)
width = 10
height = 5.5
figurepath = paste("./figure/sse_h7_2.png", sep="")
ggsave( filename = figurepath, plot = g, width = width,  height = height,  units = c("in"), dpi = 500)

width = 10
height = 3
figurepath = paste("./figure/sse_h7_1.png", sep="")
ggsave( filename = figurepath, plot = g1, width = width,  height = height,  units = c("in"), dpi = 500)










