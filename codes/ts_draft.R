

n = 100; p = 5 ; ar = c(p:1) / (2*p) * rep(c(1, -1), p)[1:p] ;ar
N=N_simu = 100; pmax0 = 10;gamma = 1;l=3
sub_path = paste(main_path,"/simu/n",n,"_N",N,"/",sep="");sub_path
sub_path = "./"
set.seed(1)
ar_sim = arima.sim(list(order=c(p,0,0), ar=ar), n=n)
ar_sim

ar_sim2 = tail(ar_sim)
ar_sim
ar_sim2

order.true = 2
ar2 = ar(ar_sim2, aic = FALSE, order.max = order.true)
ar1 = ar(ar_sim, aic = FALSE, order.max = order.true)
ar1
ar2





