rm(list=ls())





#####

#' 
#' 
#' 
#' 
#' 

getwd()
library(quantmod)
getSymbols("AAPL",  # "AMZN" Amazon
           from = "2019/12/31",
           to = "2020/12/31",
           periodicity = "daily")
head(AAPL)
myStocks <-lapply(c("AAPL", "GOOG"), function(x) {getSymbols(x, 
                                                             from = "2016/12/31", 
                                                             to = "2018/12/31",
                                                             periodicity = "daily",
                                                             auto.assign=FALSE)} )
names(myStocks) <- c("AAPL", "GOOG")
head(myStocks$AAPL)
adjustedPrices <- lapply(myStocks, Ad)
adjustedPrices <- do.call(merge, adjustedPrices)
adjprice = adjustedPrices$GOOG.Adjusted[1:100]


date <- index(adjprice)
data <- data.frame(date,GOOG = adjprice)
data$year <- format(data$date, "%Y") 
data$month <- format(data$date, "%m")
data$ym <- format(data$date, "%Y-%m")
ts.plot(data$GOOG.Adjusted, main = "TS plot of GOOG's adjusted price")




   
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 

library("quantmod") #https://www.quantmod.com/
#https://bookdown.org/kochiuyu/Technical-Analysis-with-R/downloading-data.html
library(TSA) # eacf
if(!require(tseries)) install.packages("tseries") # adf.test
library(forecast) #auto.arima



getSymbols("AAPL", from = "2020/1/1",  to = "2021/3/6",  periodicity = "daily", "getSymbols.warning4.0"=FALSE)
company = AAPL
AdjClose=Ad(company) #adjusted closing price
n = dim(AdjClose)



#' chartSeries() provides OHLC chart, where OHLC stands
#' for Open-high-low-close chart. It is a typical finantial plot.
chartSeries(company,type="bar",theme=chartTheme("white"))

install.packages("PerformanceAnalytics")
library(PerformanceAnalytics)

return = dailyReturn(AAPL)
str(return)
str(logReturn)
# Performance Summary
charts.PerformanceSummary(return, main="Naive Buy Rule")

# color bars
library(scales) # show_col
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
show_col(cbp1)
clb2 <- unlist(chartTheme('white'))
show_col(clb2)



AR1 <- arima(logReturn,c(1,0,0)) 
AR1


autoarma <-auto.arima (logReturn) 
autoarma

l=5

theForecast <-predict(AR1,
                      n.ahead=2,
                      se.fit=TRUE)

plot.ts(theForecast$pred)
plot.ts(theForecast$se)







# ozone data
#' chartSeries() provides OHLC chart, where OHLC stands
#' for Open-high-low-close chart. It is a typical finantial plot.
chartSeries(AAPL,type="bar",theme=chartTheme("white"))
addSMA(n=20,on=1,col = "lightblue")
addBBands(n=20,sd=1)
