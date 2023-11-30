###############################################################################
#' @source https://kourentzes.com/forecasting/wp-content/uploads/2017/06/Forecasting-with-R-notes.pdf
#' @description Xinyu's exploration of time series in R
###############################################################################
#load the necessary libraries
library(forecast)
library(tstools)
library(tsutils)
library(stats)

# Load data
library(tseries)
ts = get.hist.quote(instrument = "AAPL", start = "2018-01-02", end="2019-12-27", quote = "Close")
summary(ts)
y <- ts(ts)
y <- AirPassengers
plot(y)
summary(y)

tail(y,1)
frequency(y)
start(y)
# Let us look for trend in the data, by calculating the Centred Moving
# Average
cma <- cmav(y, outplot=1)
print(cma)

# Cox-Stuart on the estimated CMA to test for significant trend
coxstuart(cma)
par(mfrow=c(2,2))
# We can test for seasonality visually by producing a seasonal plot
seasplot(y, m=12)
# This functions removes the trend automatically but we can control this.
# It also provides other interesting visualisations of the seasonal component
seasplot(y,m=12,outplot=2)
seasplot(y,m=12,outplot=3)
seasplot(y,m=12,outplot=4)

# The equivalent function in the forecast package is:
seasonplot(AirPassengers, year.labels=TRUE)

# Decomposition
par(mfrow=c(1,1))
# We can perform classical decomposition using the decomp function
decomp(y,outplot=1)
# or using the pure.seasonal for estimating the seasonal indices
y.dc <- decomp(y,outplot=1, type="pure.seasonal", h=12)
# Control chart of residuals
residout(y.dc$irregular)

# STL time series decomposition
y.stl <- stl(y,s.window=7)
plot(y.stl)



# Load the necessary libraries
library(forecast)
library(MAPA)
library(TStools)

fit.ets <- ets(y)
print(fit.ets)
# And then we use this to forecast
f.ets <- forecast(fit.ets,h=1)
# Notice that since now we have a model we can produce analytical prediction
intervals
print(f.ets)
# Plot the resulting forecast and prediction intervals
plot(f.ets)
