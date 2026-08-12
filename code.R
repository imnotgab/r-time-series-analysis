library(forecast)
data <- read.table("data/dane_as.txt", header = TRUE, fill = TRUE)
data.ts <- ts(data=data, start=c(2015,7), frequency = 365.25)

head(data.ts)

#aggregation -----------

library(zoo)
#aggregate from daily to monthly
data$date <- as.Date(data$date, format="%d.%m.%Y")
data.zoo <- zoo(data$pm25, order.by = data$date)
data.monthly <- aggregate(data.zoo, as.yearmon, mean, na.rm = TRUE)

plot(data.monthly)


#initial plots ------

data.ts <- as.ts(data.monthly)

seasonplot(data.ts, year.labels = T, year.labels.left = T, col=rainbow(6))

monthplot(data.ts)

ts.plot(data.ts)

lag.plot(data.ts, lags=12, do.lines = FALSE, pch=20)

Acf(data.ts, main='ACF', lag.max=120)
#seasonality present
Pacf(data.ts, main='PACF', lag.max=120)
#trend present
#last statistically significant lag value is 12

plot(decompose(data.ts))

tsoutliers(data.ts)
#no outliers


#split dataset into train and test -------
data.ts.train <- window(data.ts, end=c(2023,6))
data.ts.test <- window(data.ts, start=c(2023,7))
ts.plot(data.ts.train)
ts.plot(data.ts.test)



#Box-Cox transformation ------
(lambda <- BoxCox.lambda(data.ts.train))
data.ts.train.bc <- BoxCox(data.ts.train, lambda)
tsoutliers(data.ts.train.bc)
data.ts.train.bc = tsclean(data.ts.train.bc) # cleaning
ts.plot(data.ts.train.bc) # ready train set with stabilized variance


#AR12 MA12 differencing ------
data.ts.train.bc.diff12 <- diff(data.ts.train.bc, lag = 12)
#removing seasonality
ts.plot(data.ts.train.bc.diff12)
Acf(data.ts.train.bc.diff12, lag.max = 120)
#last significant lag is 12 - q AR

Pacf(data.ts.train.bc.diff12, lag.max = 120)
Pacf(data.ts.train.bc.diff12)
#last significant lag is 19 or 12 - p MA

lag.plot(data.ts.train.bc.diff12, lags=12, do.lines = FALSE, pch=20)

#data.diff.residuals ------
data.diff = Arima(data.ts.train.bc,
                    order=c(0,1,0),
                    seasonal = c(0,1,0))
data.diff.residuals = data.diff$residuals

summary(data.diff)
Acf(data.diff.residuals)
Pacf(data.diff.residuals)
shapiro.test(data.diff.residuals)
wilcox.test(data.diff.residuals)
Box.test(data.diff.residuals)
Box.test(data.diff.residuals,
         type = 'Ljung-Box')


#tslm ----
data.tslm = tslm(data.ts.train.bc ~ trend + season)
data.tslm.residuals = data.tslm$residuals
Pacf(data.tslm.residuals, lag.max = 120)


#tslm MA12 -----
data.tslm.MA12 = Arima(data.tslm.residuals,
                         order = c(0,0,12),
                         seasonal = c(0,0,0),
                         include.mean = F)
summary(data.tslm.MA12)

#coefficient significance
coefs = data.tslm.MA12$coef
coefs.se = sqrt(diag(data.tslm.MA12$var.coef))
ind = abs(coefs/(1.96*coefs.se))
signif = which(ind >= 1) # ma1 and 7 are significant
temp.fixed = numeric(length(coefs))
temp.fixed[signif] = NA # preparing vector for arima function
temp.fixed

#remodeling
data.tslm.MA12.fixed = Arima(data.tslm.residuals,
                       order = c(0,0,12),
                       seasonal = c(0,0,0),
                       include.mean = F,
                       fixed = temp.fixed)
summary(data.tslm.MA12.fixed)

#checking residuals
data.tslm.MA12.fixed.residuals = data.tslm.MA12.fixed$residuals
summary(data.tslm.MA12.fixed.residuals)
Acf(data.tslm.MA12.fixed.residuals)
Pacf(data.tslm.MA12.fixed.residuals)
shapiro.test(data.tslm.MA12.fixed.residuals)
wilcox.test(data.tslm.MA12.fixed.residuals)
Box.test(data.tslm.MA12.fixed.residuals)
Box.test(data.tslm.MA12.fixed.residuals,
         type = 'Ljung-Box')


#tslm AR12 -----
data.tslm.AR12 = Arima(data.tslm.residuals,
                              order = c(12,0,0),
                              seasonal = c(0,0,0),
                              include.mean = F)
summary(data.tslm.AR12)
# coefficient significance
coefs = data.tslm.AR12$coef
coefs.se = sqrt(diag(data.tslm.AR12$var.coef))
ind = abs(coefs/(1.96*coefs.se))
(signif = which(ind >= 1))
temp.fixed = numeric(length(coefs))
temp.fixed[signif] = NA
temp.fixed

#remodeling
data.tslm.AR12.fixed = Arima(data.tslm.residuals,
                             order = c(12,0,0),
                             seasonal = c(0,0,0),
                             include.mean = F,
                             fixed = temp.fixed)
summary(data.tslm.AR12.fixed)

#checking residuals
data.tslm.AR12.fixed.residuals = data.tslm.AR12.fixed$residuals
summary(data.tslm.AR12.fixed.residuals)
Acf(data.tslm.AR12.fixed.residuals)
Pacf(data.tslm.AR12.fixed.residuals)
shapiro.test(data.tslm.AR12.fixed.residuals)
wilcox.test(data.tslm.AR12.fixed.residuals)
Box.test(data.tslm.AR12.fixed.residuals)
Box.test(data.tslm.AR12.fixed.residuals,
         type = 'Ljung-Box')


#autoarima ------
data.tslm.auto = auto.arima(data.tslm.residuals)
summary(data.tslm.auto) # suggests white noise

coefs = data.tslm.auto$coef
coefs.se = sqrt(diag(data.tslm.auto$var.coef))
ind = abs(coefs/(1.96*coefs.se))
signif = which(ind >= 1)
temp.fixed = numeric(length(coefs))
temp.fixed[signif] = NA
temp.fixed

#check residuals
Acf(data.tslm.auto$residuals, lag.max=120)
Pacf(data.tslm.auto$residuals, lag.max=120)
shapiro.test(data.tslm.auto$residuals)
t.test(data.tslm.auto$residuals)
Box.test(data.tslm.auto$residuals)
Box.test(data.tslm.auto$residuals,
         type='Ljung-Box')


#forecasts -----
length(data.ts.test)
#24
data.tslm.forecast <- forecast(data.tslm, h=24)

#tslm
data.tslm.forecast
plot(data.tslm.forecast)
plot(data.tslm.forecast$mean)

data.tslm.forecast.inv <- InvBoxCox(data.tslm.forecast$mean, lambda)


ts.plot(data.ts.test,
        data.tslm.forecast.inv, col=c('red','blue'))

#MA12
data.tslm.MA12.forecast <- forecast(data.tslm.MA12.fixed, h=24)
data.tslm.MA12.forecast.total <- data.tslm.forecast$mean + data.tslm.MA12.forecast$mean

data.tslm.MA12.forecast.total.inv <- InvBoxCox(data.tslm.MA12.forecast.total, lambda)

ts.plot(data.ts.test, 
        data.tslm.MA12.forecast.total.inv, col=c('red','blue'))

#AR12
data.tslm.AR12.forecast <- forecast(data.tslm.AR12.fixed, h=24)
data.tslm.AR12.forecast.total <- data.tslm.forecast$mean + data.tslm.AR12.forecast$mean

data.tslm.AR12.forecast.total.inv <- InvBoxCox(data.tslm.AR12.forecast.total, lambda)

ts.plot(data.ts.test, 
        data.tslm.AR12.forecast.total.inv, col=c('red','blue'))

ts.plot(data.ts.test, data.tslm.forecast.inv, data.tslm.MA12.forecast.total.inv,
        data.tslm.AR12.forecast.total.inv, col=c('red','blue','green', 'purple'))

#autoarima
data.tslm.auto.forecast <- forecast(data.tslm.auto, h=24)
data.tslm.auto.forecast.total <- data.tslm.forecast$mean + data.tslm.auto.forecast$mean

data.tslm.auto.forecast.total.inv <- InvBoxCox(data.tslm.auto.forecast.total, lambda)


ts.plot(data.ts.test, data.tslm.forecast.inv, data.tslm.MA12.forecast.total.inv,
        data.tslm.AR12.forecast.total.inv, data.tslm.auto.forecast.total.inv, col=c('red','blue','green', 'purple','yellow'))

#differencing
data.diff.forecast <- forecast(data.diff, h=24)
plot(data.diff.forecast)
data.diff.forecast.inv = InvBoxCox(data.diff.forecast$mean, lambda)
accuracy(data.diff)

ts.plot(data.ts.test, data.tslm.forecast.inv, data.tslm.MA12.forecast.total.inv,
        data.tslm.AR12.forecast.total.inv, data.diff.forecast.inv,
        col=c('red','blue','green', 'purple','yellow', 'pink'))

#autoarima on baseline
data.auto <- auto.arima(data.ts.train.bc)
summary(data.auto)
data.auto.forecast <- forecast(data.auto, h=24)

data.auto.forecast.inv <- InvBoxCox(data.auto.forecast$mean, lambda)

ts.plot(data.ts.test, data.tslm.forecast.inv, data.tslm.MA12.forecast.total.inv,
        data.tslm.AR12.forecast.total.inv, data.diff.forecast.inv,
        data.auto.forecast.inv,
        col=c('red','blue','green', 'purple','yellow', 'pink','magenta'))

#exponential smoothing ----
hw(data.ts.train.bc)
data.hw.forecast.bc = hw(data.ts.train.bc, h=24)
summary(data.hw.forecast.bc)
plot(data.hw.forecast.bc)

data.hw.forecast.inv = InvBoxCox(data.hw.forecast.bc$mean, lambda)

#normality
shapiro.test(data.hw.forecast.bc$residuals) # assumption is met
Acf(data.hw.forecast.bc$residuals, lag.max = 120)
Pacf(data.hw.forecast.bc$residuals, lag.max=120)
#mean
t.test(data.hw.forecast.bc$residuals) # do not reject
#randomness
Box.test(data.hw.forecast.bc$residuals)
Box.test(data.hw.forecast.bc$residuals, type=c('Ljung-Box'))

#final plot ----
ts.plot(data.ts.test, data.tslm.forecast.inv, data.tslm.MA12.forecast.total.inv,
        data.tslm.AR12.forecast.total.inv, data.diff.forecast.inv,
        data.auto.forecast.inv, data.hw.forecast.inv,
        col=c('red','blue','green', 'purple','yellow', 'pink','magenta','cyan'))
legend('topleft',
       col = c('red','blue','green', 'purple','yellow', 'pink','magenta','cyan'),
       c('actual values',
         'ARIMA(insert)[insert]',
         'tslm + MA(12)', 'tslm + AR(12)',
         'tslm + insert',
        'ARIMA(insert)',
        'exponential smoothing'),
       lty = c(1,1,1,1,1,1,1))

#forecast accuracy ----
#calculates accuracy with respect to the test model
#best to make a table and compare everything
accuracy(data.diff.forecast.inv, data.ts.test)
accuracy(data.tslm.forecast.inv, data.ts.test)
accuracy(data.tslm.AR12.forecast.total.inv, data.ts.test)
accuracy(data.tslm.MA12.forecast.total.inv, data.ts.test)
accuracy(data.auto.forecast.inv, data.ts.test)
accuracy(data.hw.forecast.inv, data.ts.test)

plot(data.tslm.MA12.forecast)
plot(data.tslm.auto.forecast)

#assuming we have:
#data.tslm - tslm model with trend and seasonality
#data.tslm.MA12.fixed - MA(12) model fitted to tslm residuals
#lambda - Box-Cox parameter

#generate forecast
forecast.tslm <- forecast(data.tslm, h = 24)
forecast.MA12 <- forecast(data.tslm.MA12.fixed, h = 24)

#sum: trend + seasonality + residuals
forecast.total <- forecast.tslm$mean + forecast.MA12$mean

#invert Box-Cox transformation
forecast.total.inv <- InvBoxCox(forecast.total, lambda)

#plot
ts.plot(data.ts.test, forecast.total.inv,
        col = c("red", "blue"),
        main = "PM2.5 - TSLM + MA(12) Forecast",
        ylab = "PM2.5 [µg/m³]")
legend("topleft", legend = c("Actual values", "Forecast"), col = c("red", "blue"), lty = 1)

pchisq(2.25, df=2, lower.tail = FALSE)
