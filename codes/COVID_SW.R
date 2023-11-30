

library(tseries)
library(forecast)
library(Metrics)
library(ggplot2)
library(readr)
library(WaveletArima)
library(caret)
library(nnfor)
library(tsDyn)
library(fracdiff)
library(bsts)
library(forecastHybrid)
library(e1071)
library(tseriesChaos)
library(pracma)
library(Kendall)
library(nonlinearTseries)
library(GeneCycle)
library(fpp2)


# data <- read_csv("location of the USA data")
data <- read_csv("USA_data.csv")
datats=ts(data$USA)
train = window(datats, start = 1, end=210) #define training data
test = window(datats, start= 211,end=240) #define test data

#plot time series with training and test period
autoplot(train)+ggtitle('COVID-19 cases in USA')+
  ylab('Cases')+xlab('Days')+autolayer(test)

#ACF & PACF plots
#diffcase=ndiffs(train)
diffset = diff(train, differences = ndiffs(train))

ggAcf(diffset) +
  geom_point(color = 'navy blue') +
  ggtitle("ACF plot")

ggPacf(train) +
  geom_point(color = 'navy blue') +
  ggtitle("PACF plot")

#performing tests
fracdiff(train)
hurstexp(train)
Box.test(train, lag = 1, type = c("Box-Pierce"))
skewness(train)
kurtosis(train)
max_lyapunov_expo <- lyap_k(train, m=1, d=2, s=1, t=4, ref=length(train), k=2, eps=4)

#Stationarity tests
kpss.test(train) 

#nonlinearity tests
nonlinearityTest(train, verbose = TRUE)
test=data$USA[211:240] #load test data as a column vector

#### Sliding window parameters
n = 210
h = 30
kmin = h*2
nk = n-kmin+1
sw_results = list()




#fitting arima model
fitARIMA = auto.arima(train) 
summary(fitARIMA)
predARIMA = forecast::forecast(fitARIMA, h=30)
plot(predARIMA)
#accuracy metrics
a1<-forecast::accuracy(predARIMA, test);a1
smape(test,predARIMA$mean)


criteria = data.frame(k=c(kmin:n), RMSE=numeric(nk), SMAPE = numeric(nk))
for (k in c(kmin:n)){
  train_k = window(datats, start = n-k+1, end=n) 
  predARIMA = forecast::forecast(auto.arima(train_k), h=30)
  a1<-forecast::accuracy(predARIMA, test)
  criteria[k-kmin+1,'RMSE'] = a1[2,2]
  criteria[k-kmin+1,'SMAPE'] = smape(test,predARIMA$mean)
}

plot(criteria$k, criteria$RMSE, type="l")
sw_results$arima = criteria

#fitting ETS model
fitETS=ets(train)
summary(fitETS)
predETS=forecast::forecast(fitETS, h=30)
plot(predETS)
a2<-forecast::accuracy(predETS, test);a2
smape(test,predETS$mean)


criteria = data.frame(k=c(kmin:n), RMSE=numeric(nk), SMAPE = numeric(nk))
for (k in c(kmin:n)){
  train_k = window(datats, start = n-k+1, end=n) 
  predETS = forecast::forecast(ets(train_k), h=30)
  a1<-forecast::accuracy(predETS, test)
  criteria[k-kmin+1,'RMSE'] = a1[2,2]
  criteria[k-kmin+1,'SMAPE'] = smape(test,predETS$mean)
}

plot(criteria$k, criteria$RMSE, type="l")
sw_results$ets = criteria

#SETAR model
fit_SETAR = setar(train, m = 4)
fc_SETAR = predict(fit_SETAR, n.ahead = 30)
a3<-forecast::accuracy(fc_SETAR, test);a3

smape(test, fc_SETAR)


criteria = data.frame(k=c(kmin:n), RMSE=numeric(nk), SMAPE = numeric(nk))
for (k in c(kmin:n)){
  train_k = window(datats, start = n-k+1, end=n) 
  predSETAR = predict(setar(train_k, m=4), n.ahead = 30)
  a1<-forecast::accuracy(predSETAR, test)
  criteria[k-kmin+1,'RMSE'] = a1[,2]
  criteria[k-kmin+1,'SMAPE'] = smape(test,predSETAR)
}

plot(criteria$k, criteria$RMSE, type="l")
sw_results$setar = criteria



##fitting tbats
fit_tbats = tbats(train)
summary(fit_tbats)
predTBATS=forecast::forecast(fit_tbats, h=30)
autoplot(predTBATS)
a4<-forecast::accuracy(predTBATS, test);a4
smape(test,predTBATS$mean)


criteria = data.frame(k=c(kmin:n), RMSE=numeric(nk), SMAPE = numeric(nk))
for (k in c(kmin:n)){
  train_k = window(datats, start = n-k+1, end=n) 
  predTBATS= forecast::forecast(tbats(train_k), h=30)
  a1<-forecast::accuracy(predTBATS, test)
  criteria[k-kmin+1,'RMSE'] = a1[2,2]
  criteria[k-kmin+1,'SMAPE'] = smape(predTBATS$mean, test)
  print(k-kmin+1)
}

plot(criteria$k, criteria$RMSE, type="l")
sw_results$tbats = criteria


#Theta model
fit_theta=thetaf(train, h=30)
autoplot(fit_theta)
a5<-forecast::accuracy(fit_theta$mean, test);a5
smape(test,fit_theta$mean)

criteria = data.frame(k=c(kmin:n), RMSE=numeric(nk), SMAPE = numeric(nk))
for (k in c(kmin:n)){
  train_k = window(datats, start = n-k+1, end=n) 
  predTHETA= thetaf(train_k, h=30)
  a<-forecast::accuracy(predTHETA, test)
  criteria[k-kmin+1,'RMSE'] = a[2,2]
  criteria[k-kmin+1,'SMAPE'] = smape(predTHETA$mean, test)
}

plot(criteria$k, criteria$RMSE, type="l")
sw_results$theta = criteria

str(sw_results)
#fitting MLP/ANN
fit_ANN = mlp(train)
predANN = forecast::forecast(fit_ANN, h=30)
autoplot(predANN)
a6<-forecast::accuracy(predANN$mean, test);a6
smape(test,predANN$mean)

criteria = data.frame(k=c(kmin:n), RMSE=numeric(nk), SMAPE = numeric(nk))
for (k in c(kmin:n)){
  train_k = window(datats, start = n-k+1, end=n) 
  predANN = forecast::forecast(mlp(train_k), h=30)
  a<-forecast::accuracy(predANN, test)
  criteria[k-kmin+1,'RMSE'] = a[2,2]
  criteria[k-kmin+1,'SMAPE'] = smape(predANN$mean, test)
  print(k-kmin+1)
}
plot(criteria$k, criteria$RMSE, type="l")
sw_results$ann = criteria


#fitting ARNN model
fit_ARNN = nnetar(train)
predARNN=forecast::forecast(fit_ARNN, h= 30)
plot(predARNN)
a7<-forecast::accuracy(predARNN$mean, test);a7
smape(test, predARNN$mean)

criteria = data.frame(k=c(kmin:n), RMSE=numeric(nk), SMAPE = numeric(nk))
for (k in c(kmin:n)){
  train_k = window(datats, start = n-k+1, end=n) 
  predARNN = forecast::forecast(nnetar(train_k), h=30)
  a<-forecast::accuracy(predARNN, test)
  criteria[k-kmin+1,'RMSE'] = a[2,2]
  criteria[k-kmin+1,'SMAPE'] = smape(predARNN$mean, test)
  print(k-kmin+1)
}
plot(criteria$k, criteria$RMSE, type="l")
sw_results$arnn = criteria

#fitting of Wavelet ARIMA
fit_wa <- WaveletFittingarma(train, Waveletlevels = floor(log(length(train))), boundary = 'periodic', FastFlag = TRUE, MaxARParam = 5,
                             MaxMAParam = 5, NForecast = 30)

a8<-forecast::accuracy(fit_wa$Finalforecast, test);a8
smape(test,fit_wa$Finalforecast)

criteria = data.frame(k=c(kmin:n), RMSE=numeric(nk), SMAPE = numeric(nk))
for (k in c(kmin:n)){
  train_k = window(datats, start = n-k+1, end=n) 
  predWVARIMA = WaveletFittingarma(train_k, Waveletlevels = floor(log(length(train_k))), boundary = 'periodic', FastFlag = TRUE, MaxARParam = 5,
                              MaxMAParam = 5, NForecast = 30)
  a<-forecast::accuracy(predWVARIMA$Finalforecast, test)
  criteria[k-kmin+1,'RMSE'] = a[2,2]
  criteria[k-kmin+1,'SMAPE'] = smape(redWVARIMA$Finalforecast, test)
}
plot(criteria$k, criteria$RMSE, type="l")
sw_results$wvarima = criteria



#### summary results
l = length(sw_results)
method_names = names(sw_results)

for(i in 1:l){
  plot(criteria$k, criteria$RMSE, type="l", xlab="sample size k", ylab="Test RMSE", 
       main=method_names[l], col="dark blue")
}





#### models run by now
#fitting bsts model
ss <- AddLocalLinearTrend(list(), train)
fit_bsts=bsts(train,state.specification = ss, niter = 1000)
predBSTS <- predict(fit_bsts, horizon = 30)
plot(predBSTS, plot.original = 211)

burn <- SuggestBurn(0.1, fit_bsts)
fitted_bsts=as.numeric(-colMeans(fit_bsts$one.step.prediction.errors[-(1:burn),])+train)

a9<-forecast::accuracy(predBSTS$mean, test);a9
smape(test,predBSTS$mean)

#fitting ARFIMA model
fit_ARFIMA=arfima(train)
predARFIMA = forecast::forecast(fit_ARFIMA, h=30)
autoplot(predARFIMA)
a10<-forecast::accuracy(predARFIMA$mean, test);a10

smape(test,predARFIMA$mean)


#fitting arima+ann model
fit_res_ANN=mlp(fitARIMA$residuals)
pred_res_ANN = forecast::forecast(fit_res_ANN, h=30)
pred_arima_ann=predARIMA$mean+pred_res_ANN$mean
a11<-forecast::accuracy(pred_arima_ann, test);a11

smape(test,pred_arima_ann)


#fitting arima+arnn model
fit_res_ARNN=nnetar(fitARIMA$residuals)
pred_res_ARNN = forecast::forecast(fit_res_ARNN, h=30)
pred_arima_arnn=predARIMA$mean+pred_res_ARNN$mean
a12<-forecast::accuracy(pred_arima_arnn, test);a12

smape(test,pred_arima_arnn)

#fitting arima+wbf model
fit_res_wbf=WaveletFittingarma(fitARIMA$residuals, Waveletlevels = floor(log(length(train))), boundary = 'periodic', FastFlag = TRUE, MaxARParam = 5,
                               MaxMAParam = 5, NForecast = 30)

pred_arima_wbf=predARIMA$mean+fit_res_wbf$Finalforecast
a13<-forecast::accuracy(pred_arima_wbf, test);a13

smape(test,pred_arima_wbf)

#fitting warima+ann model
res_wa = train - fit_wa$FinalPrediction
fit_wa_ANN=mlp(res_wa)
pred_wa_ANN = forecast::forecast(fit_wa_ANN, h=30)
pred_wa_ann=fit_wa$Finalforecast+pred_wa_ANN$mean
a14<-forecast::accuracy(pred_wa_ann, test);a14

smape(test,pred_wa_ann)


#fitting warima+arnn model
fit_wa_ARNN=nnetar(res_wa)
pred_wa_ARNN = forecast::forecast(fit_wa_ARNN, h=30)
pred_wa_arnn=fit_wa$Finalforecast+pred_wa_ARNN$mean
a15<-forecast::accuracy(pred_wa_arnn, test);a15

smape(test,pred_wa_arnn)

#fitting Avg(ARIMA, ETS, Theta) model using forecast hybrid
avg_aef=hybridModel(train, weights="equal",errorMethod = "MASE", models = "aef")
pred_aef = forecast::forecast(avg_aef, h=30)
autoplot(pred_aef)
a16<-forecast::accuracy(pred_aef$mean, test);a16
smape(test, pred_aef$mean)

#fitting Avg(ARIMA, ETS, nnar) model using forecast hybrid
avg_aen=hybridModel(train, weights="equal",errorMethod = "MASE", models = "aen")
pred_aen = forecast::forecast(avg_aen, h=30)
autoplot(pred_aen)
a17<-forecast::accuracy(pred_aen$mean, test);a17
smape(test, pred_aen$mean)

#fitting Avg(ARIMA, theta, nnar) model using forecast hybrid
avg_afn=hybridModel(train, weights="equal",errorMethod = "MASE", models = "afn")
pred_afn = forecast::forecast(avg_afn, h=30)
autoplot(pred_afn)
a18<-forecast::accuracy(pred_afn$mean, test);a18
smape(test, pred_afn$mean)

#fitting Avg(ETS, theta, nnar) model using forecast hybrid
avg_efn=hybridModel(train, weights="equal",errorMethod = "MASE", models = "efn")
pred_efn = forecast::forecast(avg_efn, h=30)
autoplot(pred_efn)
a19<-forecast::accuracy(pred_efn$mean, test);a19
smape(test, pred_efn$mean)

#fitting Avg(ANN,NNAR,WA) manually
pred_hybrid=0.333*predANN$mean+0.333*predARNN$mean+0.333*fit_wa$Finalforecast
a20<-forecast::accuracy(pred_hybrid, test);a20
smape(test, pred_hybrid)


#fitting AR Elas
source("~/Desktop/Project1_TS/codes/AR_Elas.R")
data = window(datats, start = 1, end=240) #define training data
ark = AR_Elas(data, n=210, h = 30, k = 210, gamma = 1)
a21 <- forecast::accuracy(ark$yhat, test);a21
smape(test, ark$yhat)

armod = ar(train, aic=FALSE, order.max=ark$order, demean=FALSE)
ar_pred = predict(armod, n.ahead = 30)$pred
a22 <- forecast::accuracy(ar_pred, test)
smape(test, ar_pred)


ark$order
a20
a21
a22

