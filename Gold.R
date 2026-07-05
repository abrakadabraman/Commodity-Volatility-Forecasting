
library(dplyr)
library(forecast)
library(tsutils)
library(tseries)
library(xts)
library(smooth)
library(zoo)
library(ggplot2)
library(urca)
library(rugarch)
library(FinTS)
library(moments)
library(MASS)
library(foreach)
library(doParallel)
library(MCS)


oil <- read.csv("Resources_Dataset.csv")

oil_wti <- oil[, c("Date", "GC.F_Close")]
oil_wti$Date <- as.Date(oil_wti$Date, format = "%d/%m/%Y")

oil_xts <- xts(oil_wti$GC.F_Close, order.by = oil_wti$Date)

### Elininate top NAs
oil_xts <- oil_xts[which(!is.na(oil_xts))[1]:NROW(oil_xts)]


sum(is.na(oil_xts))

# drop na
oil_xts <- oil_xts[!is.na(oil_xts)]

sum(is.na(oil_xts))

mean(oil_xts, na.rm = TRUE)
median(oil_xts, na.rm = TRUE)
min(oil_xts, na.rm = TRUE)
max(oil_xts, na.rm = TRUE)
sd(oil_xts, na.rm = TRUE)

oil_xts[which.min(oil_xts)]
oil_xts[which.max(oil_xts)]


# Stationarity
adf.test(oil_xts)
kpss.test(oil_xts)

#  Log returns 
r_xts <- na.omit(diff(log(oil_xts)))
plot(r_xts, main = "WTI Prices - Log Returns", ylab = "Log return")

# ADF/KPSS 
adf.test(r_xts)
kpss.test(r_xts)


mean(r_xts, na.rm = TRUE)
median(r_xts, na.rm = TRUE)
min(r_xts, na.rm = TRUE)
max(r_xts, na.rm = TRUE)
sd(r_xts, na.rm = TRUE)
skewness(r_xts, na.rm = TRUE)
kurtosis(r_xts, na.rm = TRUE)


# Centered residuals
# Residuals 
r <- as.numeric(r_xts)
eps <- r                       
abs_returns <- abs(eps)
squared_returns <- eps^2      


# arch test
ArchTest(eps, lags = 12)

# Normality
jarque.bera.test(eps)

# Histogram of log returns
returns_df <- data.frame(log_returns = r)
ggplot(returns_df, aes(x = log_returns)) +
  geom_histogram(bins = 100, fill = "steelblue", color = "white") +
  labs(title = "Histogram of Gold Log Returns",
       x = "Log Returns", y = "Frequency") +
  theme_minimal()

# Histogram + Normal overlay
mu <- mean(returns_df$log_returns)
sigma <- sd(returns_df$log_returns)
ggplot(returns_df, aes(x = log_returns)) +
  geom_histogram(aes(y = ..density..), bins = 100, fill = "steelblue", color = "white", alpha = 0.6) +
  stat_function(fun = dnorm, args = list(mean = mu, sd = sigma), size = 1.2) +
  labs(title = "Histogram of Gold Log Returns with Normal Curve",
       x = "Log Returns", y = "Density") +
  theme_minimal()


#Skewness & Kurtosis
skewness(eps)     
kurtosis(eps)    

#QQ plots
qqnorm(eps, main = "QQ Plot of Daily Gold Log Returns"); qqline(eps, col = "red")



# Autocorrelation
par(mfrow =c(2,2))
acf(r_xts,    main = "ACF of Daily Gold Returns")
pacf(r_xts,   main = "PACF of Daily Gold Returns")


acf(r_xts^2,  main = "ACF of Squared Daily Gold Returns")
pacf(r_xts^2, main = "PACF of Squared Daily Gold Returns")


par(mfrow = c(1,1))
acf(abs(r_xts),  main = "ACF of Absolute Gold Returns")   
pacf(abs(r_xts), main = "PACF of Absolute Gold Returns")  

# Ljung–Box
Box.test(r_xts, lag = 40, type = "Ljung-Box")

# Train/Test 
n <- length(eps)
train_size <- floor(0.8 * n)

# GARCH uses returns
garch_train <- eps[1:train_size]
garch_test  <- eps[(train_size + 1):n]
n_test <- length(garch_test)

# STES/ES uses variance proxy 
stes_train <- squared_returns[1:train_size]
stes_test  <- squared_returns[(train_size + 1):n]

# Rolling window 
window <- 1000   

# plot
plot(ts(squared_returns), main = "Daily Squared Residuals (Variance Proxy)", ylab = expression(epsilon[t]^2))
plot(oil_xts, main = "Gold Daily Closing Price", ylab = "Price")
