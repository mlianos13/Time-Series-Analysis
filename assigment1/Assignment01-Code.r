library(ggplot2)
library(diagram)
library(emmeans)
library(dplyr)
library(MASS)
library(car)
library(multcomp)
library(readxl)
library(tidyverse)
library(fpp2)
library(tidyr)
library(lme4)
library(lmtest)
library(languageserver)


data <- read_excel("DST_BIL54_train.xlsx")
colnames(data)[1] <- "Category"

dt1 <- filter(data, Category == "Drivmidler i alt")
dt1 <- dt1[, -1]
# Create the time variable
months <- seq(as.Date("2018-01-01"), as.Date("2022-11-01"), by="month")
# months <- seq(as.Date("2018-01-01"), as.Date("2022-11-01"), by="month") creates a sequence of dates from January 2018 to November 2022, one for each month.

# iterate through months and assign to x 
x <- as.numeric(format(months, "%Y")) + (as.numeric(format(months, "%m")) - 1) / 12 
print(x)
dt <- pivot_longer(dt1, 
                          cols = everything(),
                          names_to = "Time", 
                          values_to = "Value")

dt$Time <- x

# Now, let's plot
ggplot(dt, aes(x = x, y = Value)) + 
  geom_point() + 
  labs(title = "Training Data over Time", x = "Time", y = "Value") +
  theme_minimal()



##--------------------------------------------------------------------------------------------------------------------------##


##--------------------------------------------------------------------------------------------------------------------------##
## 2. OLS Model


nparams <- 2 # number of parameters
n <- length(dt$Time) # number of observations
X <- cbind(1, dt$Time) # X is the "design matrix"
print(X)
y <- cbind(dt$Value) # y is vector with observations:
print(y)

## 2.1 Estimate the parameters β0 and β1 using the training data (OLS model). Describe how you estimate the parameters.

# Estimate the parameters β0 and β1 using the training data (OLS model). Describe how you estimate the parameters.
# The OLS estimate of the parameters is given by βˆ = (X'X)^(-1)X'y
# where X is the "design matrix" and y is the vector with observations.

# to estimate parameters we solve the "normal equations":
(OLS <- solve(t(X)%*%X)%*%t(X)%*%y)
# these are the parameter estimates!
theta_0 <- OLS[1]
theta_1 <- OLS[2]



## 2.2 Present the values of the parameter estimates βˆ0 and βˆ1 and their estimated standard errors ˆσβˆ0 and ˆσβˆ1

# now we compute y_hat values for all observations:
yhat_ols <- X%*%OLS
print(yhat_ols)

# plot:
ggplot(dt, aes(x=x, y=Value)) +
  geom_point() + 
  geom_line(aes(y=yhat_ols), col="#ff0000", size=.5) 

# we will now calculate the standard errors on the parameters beta_0 and beta_1:

# first compute residuals:
e_ols <- y - yhat_ols
print(e_ols)

# calculate sum of squared residuals:
RSS_ols <- t(e_ols)%*%e_ols
print(RSS_ols)

# calculate sigma^2:
sigma2_ols <- as.numeric(RSS_ols/(n - nparams)) # = "Residual variance"
sqrt(sigma2_ols) # = "Residual standard error"

# calculate variance-covariance matrix of _parameters_:
V_ols <- sigma2_ols * solve(t(X) %*% X)
print(V_ols)

# the variances of the parameters are the values in the diagonal:
diag(V_ols)

# and the standard errors are given by:
(sqrt(diag(V_ols))) 

se_theta_0 <- (sqrt(diag(V_ols)))[1]
se_theta_1 <- (sqrt(diag(V_ols)))[2]

# now we have both point estimates and standard errors:
# intercept:
theta_0
se_theta_0
# slope:
theta_1
se_theta_1



## 2.3 Make a forecast for the next 12 months - i.e., compute predicted values with corresponding prediction intervals. Present these values in a table.

# Step 1: Identify the last time point
last_x <- max(dt$Time)  # Assuming dt$Time contains your 'x' values

# Step 2: Generate future time points for the next 12 months
future_x <- seq(from = last_x + 1/12, by = 1/12, length.out = 12)


# If you need these as a design matrix for prediction (with an intercept):
Xtest <- cbind(1, future_x)  

# compute predictions 
y_pred <- Xtest%*%OLS
print(y_pred)

# compute prediction variance-covariance matrix:
Vmatrix_pred <- sigma2_ols*(1+(Xtest%*%solve(t(X)%*%X))%*%t(Xtest))

# the variances of individual predictions are in the diagonal of the matrix above

# compute "prediction intervals"
y_pred_lwr <- y_pred - 1.96*sqrt(diag(Vmatrix_pred))
y_pred_upr <- y_pred + 1.96*sqrt(diag(Vmatrix_pred))



## 2.4 Plot the fitted model together with the training data and the forecasted values (also plot the prediction intervals of the forecasted values).

# Create a data frame for the forecasted values
forecast_df <- data.frame(
  Time = future_x, 
  y_pred = as.vector(y_pred), 
  y_pred_lwr = as.vector(y_pred_lwr), 
  y_pred_upr = as.vector(y_pred_upr)
)

# Convert Time to numeric since ggplot2's geom_smooth works with numeric x-axis for prediction bands
dt$TimeNumeric <- as.numeric(dt$Time)
forecast_df$TimeNumeric <- as.numeric(forecast_df$Time)

# Add the forecasted points and prediction intervals
ggplot() +
  geom_point(data = dt, aes(x = TimeNumeric, y = Value)) +  # Plot the historical data points
  geom_line(data = dt, aes(x = TimeNumeric, y = yhat_ols), color = "red") +  # Plot the fitted line
  geom_point(data = forecast_df, aes(x = TimeNumeric, y = y_pred), color = "blue", shape = 1) +  # Plot the forecasted points
  geom_ribbon(data = forecast_df, aes(x = TimeNumeric, ymin = y_pred_lwr, ymax = y_pred_upr), fill = "blue", alpha = 0.2) +  # Prediction intervals
  labs(title = "Historical Data and Forecast", x = "Time", y = "Value") +
  theme_minimal()

# 2.6 Investigate the residuals from the OLS model. Are there any patterns in the residuals? If so, describe them.

# Summary of residuals
summary(e_ols)

# Plot of residuals vs. fitted values
plot(yhat_ols, e_ols, xlab = "Fitted Values", ylab = "Residuals", main = "Residuals vs. Fitted Values")
abline(h = 0, col = "red", lty = 2)

# QQ plot
qqnorm(e_ols)
qqline(e_ols, col = "red")
## A Q-Q plot compares the quantiles of the residuals to the quantiles of a theoretical normal distribution.

# Histogram of the residuals
hist(e_ols, breaks = 20, main = "Histogram of Residuals", xlab = "Residuals")

##--------------------------------------------------------------------------------------------------------------------------##


##--------------------------------------------------------------------------------------------------------------------------##
## 3 WLS-Local Linear Model

# 3.2 plot the "λ-weights" vs time. 
n <- 59 
lambda = 0.9
weights <- lambda^((n-1):0)
# plot the weights:
barplot(weights, names=1:59)


# 3.3 Calculate the sum of all the λ-weights. 
sum <- sum(weights[1:59])
print(sum)

# 3.4 Estimate and present the hat_beta0 and hat_beta1 using the training data (WLS model).

SIGMA <- diag(n)
for (i in 1:n) {
  SIGMA[i,i] <- 1/lambda^(n-i)
}
print(SIGMA[52:59,52:59])

# estimate parameters with WLS (using only first 26 observations)
WLS <- solve(t(X[1:59,])%*%solve(SIGMA)%*%X[1:59,])%*%(t(X[1:59,])%*%solve(SIGMA)%*%y[1:59])
yhat_wls <- X[1:59,]%*%WLS

theta_0 <- WLS[1]
theta_1 <- WLS[2]
print(theta_0)
print(theta_1)


# 3.5 Make a forecast for the next 12 months - i.e., compute predicted values with corresponding prediction intervals

## λ = 0.9
e_wls <- y - yhat_wls
RSS_wls <- t(e_wls)%*%e_wls
sigma2_wls <- as.numeric(RSS_wls/(n - nparams))

y_pred <- Xtest%*%WLS
Vmatrix_pred <- sigma2_wls*(1+(Xtest%*%solve(t(X)%*%X))%*%t(Xtest))
y_pred_lwr <- y_pred - 1.96*sqrt(diag(Vmatrix_pred))
y_pred_upr <- y_pred + 1.96*sqrt(diag(Vmatrix_pred))
forecast_df <- data.frame(
  Time = future_x, 
  y_pred = as.vector(y_pred), 
  y_pred_lwr = as.vector(y_pred_lwr), 
  y_pred_upr = as.vector(y_pred_upr)
)
dt$TimeNumeric <- as.numeric(dt$Time)
forecast_df$TimeNumeric <- as.numeric(forecast_df$Time)

ggplot() +
  geom_point(data = dt, aes(x = TimeNumeric, y = Value)) +  # Plot the historical data points
  geom_line(data = dt, aes(x = TimeNumeric, y = yhat_wls), color = "red") +  # Plot the fitted line
  geom_point(data = forecast_df, aes(x = TimeNumeric, y = y_pred), color = "blue", shape = 1) +  # Plot the forecasted points
  geom_ribbon(data = forecast_df, aes(x = TimeNumeric, ymin = y_pred_lwr, ymax = y_pred_upr), fill = "blue", alpha = 0.2) +  # Prediction intervals
  labs(title = "Historical Data and Forecast λ=0.9", x = "Time", y = "Value") +
  theme_minimal()

## λ = 0.8
 
n <- 59 
lambda = 0.8
weights <- lambda^((n-1):0)
SIGMA <- diag(n)
for (i in 1:n) {
  SIGMA[i,i] <- 1/lambda^(n-i)
}
WLS <- solve(t(X[1:59,])%*%solve(SIGMA)%*%X[1:59,])%*%(t(X[1:59,])%*%solve(SIGMA)%*%y[1:59])
yhat_wls <- X[1:59,]%*%WLS
theta_0 <- WLS[1]
theta_1 <- WLS[2]
print(theta_0)
print(theta_1)

e_wls <- y - yhat_wls
RSS_wls <- t(e_wls)%*%e_wls
sigma2_wls <- as.numeric(RSS_wls/(n - nparams))

y_pred <- Xtest%*%WLS
print(y_pred)
Vmatrix_pred <- sigma2_wls*(1+(Xtest%*%solve(t(X)%*%X))%*%t(Xtest))
y_pred_lwr <- y_pred - 1.96*sqrt(diag(Vmatrix_pred))
y_pred_upr <- y_pred + 1.96*sqrt(diag(Vmatrix_pred))
forecast_df <- data.frame(
  Time = future_x, 
  y_pred = as.vector(y_pred), 
  y_pred_lwr = as.vector(y_pred_lwr), 
  y_pred_upr = as.vector(y_pred_upr)
)
dt$TimeNumeric <- as.numeric(dt$Time)
forecast_df$TimeNumeric <- as.numeric(forecast_df$Time)

ggplot() +
  geom_point(data = dt, aes(x = TimeNumeric, y = Value)) +  # Plot the historical data points
  geom_line(data = dt, aes(x = TimeNumeric, y = yhat_wls), color = "red") +  # Plot the fitted line
  geom_point(data = forecast_df, aes(x = TimeNumeric, y = y_pred), color = "blue", shape = 1) +  # Plot the forecasted points
  geom_ribbon(data = forecast_df, aes(x = TimeNumeric, ymin = y_pred_lwr, ymax = y_pred_upr), fill = "blue", alpha = 0.2) +  # Prediction intervals
  labs(title = "Historical Data and Forecast λ=0.8", x = "Time", y = "Value") +
  theme_minimal()

## λ = 0.7

n <- 59 
lambda = 0.7
weights <- lambda^((n-1):0)
SIGMA <- diag(n)
for (i in 1:n) {
  SIGMA[i,i] <- 1/lambda^(n-i)
}
WLS <- solve(t(X[1:59,])%*%solve(SIGMA)%*%X[1:59,])%*%(t(X[1:59,])%*%solve(SIGMA)%*%y[1:59])
yhat_wls <- X[1:59,]%*%WLS
theta_0 <- WLS[1]
theta_1 <- WLS[2]
print(theta_0)
print(theta_1)

e_wls <- y - yhat_wls
RSS_wls <- t(e_wls)%*%e_wls
sigma2_wls <- as.numeric(RSS_wls/(n - nparams))

y_pred <- Xtest%*%WLS
print(y_pred)
Vmatrix_pred <- sigma2_wls*(1+(Xtest%*%solve(t(X)%*%X))%*%t(Xtest))
y_pred_lwr <- y_pred - 1.96*sqrt(diag(Vmatrix_pred))
y_pred_upr <- y_pred + 1.96*sqrt(diag(Vmatrix_pred))
forecast_df <- data.frame(
  Time = future_x, 
  y_pred = as.vector(y_pred), 
  y_pred_lwr = as.vector(y_pred_lwr), 
  y_pred_upr = as.vector(y_pred_upr)
)
dt$TimeNumeric <- as.numeric(dt$Time)
forecast_df$TimeNumeric <- as.numeric(forecast_df$Time)

ggplot() +
  geom_point(data = dt, aes(x = TimeNumeric, y = Value)) +  # Plot the historical data points
  geom_line(data = dt, aes(x = TimeNumeric, y = yhat_wls), color = "red") +  # Plot the fitted line
  geom_point(data = forecast_df, aes(x = TimeNumeric, y = y_pred), color = "blue", shape = 1) +  # Plot the forecasted points
  geom_ribbon(data = forecast_df, aes(x = TimeNumeric, ymin = y_pred_lwr, ymax = y_pred_upr), fill = "blue", alpha = 0.2) +  # Prediction intervals
  labs(title = "Historical Data and Forecast λ=0.7", x = "Time", y = "Value") +
  theme_minimal()

## λ = 0.6

n <- 59 
lambda = 0.6
weights <- lambda^((n-1):0)
SIGMA <- diag(n)
for (i in 1:n) {
  SIGMA[i,i] <- 1/lambda^(n-i)
}
WLS <- solve(t(X[1:59,])%*%solve(SIGMA)%*%X[1:59,])%*%(t(X[1:59,])%*%solve(SIGMA)%*%y[1:59])
yhat_wls <- X[1:59,]%*%WLS
theta_0 <- WLS[1]
theta_1 <- WLS[2]
print(theta_0)
print(theta_1)

e_wls <- y - yhat_wls
RSS_wls <- t(e_wls)%*%e_wls
sigma2_wls <- as.numeric(RSS_wls/(n - nparams))
 
y_pred <- Xtest%*%WLS
print(y_pred)
Vmatrix_pred <- sigma2_wls*(1+(Xtest%*%solve(t(X)%*%X))%*%t(Xtest))
y_pred_lwr <- y_pred - 1.96*sqrt(diag(Vmatrix_pred))
y_pred_upr <- y_pred + 1.96*sqrt(diag(Vmatrix_pred))
forecast_df <- data.frame(
  Time = future_x, 
  y_pred = as.vector(y_pred), 
  y_pred_lwr = as.vector(y_pred_lwr), 
  y_pred_upr = as.vector(y_pred_upr)
)
dt$TimeNumeric <- as.numeric(dt$Time)
forecast_df$TimeNumeric <- as.numeric(forecast_df$Time)

ggplot() +
  geom_point(data = dt, aes(x = TimeNumeric, y = Value)) +  # Plot the historical data points
  geom_line(data = dt, aes(x = TimeNumeric, y = yhat_wls), color = "red") +  # Plot the fitted line
  geom_point(data = forecast_df, aes(x = TimeNumeric, y = y_pred), color = "blue", shape = 1) +  # Plot the forecasted points
  geom_ribbon(data = forecast_df, aes(x = TimeNumeric, ymin = y_pred_lwr, ymax = y_pred_upr), fill = "blue", alpha = 0.2) +  # Prediction intervals
  labs(title = "Historical Data and Forecast λ=0.6", x = "Time", y = "Value") +
  theme_minimal()

##--------------------------------------------------------------------------------------------------------------------------##




##--------------------------------------------------------------------------------------------------------------------------##
## 4 Iterative update and optimal λ

lambda <- 0.9

# 4.1
f <- function(j) rbind(1, j)
L <- matrix(c(1.,0., 1.,1.),
            byrow=TRUE, nrow=2)
Linv <- solve(L) 

print(f(0))


# 4.2
i <- 1
(F_N <-  (lambda^0) * f(0)%*%t(f(0)))
(h_N <-  (lambda^0) * f(0) * y[i])


print(F_N)
print(h_N)

# 4.3 
for (i in 2:10){
  F_N <- F_N + lambda^(i-1) * f(-(i-1)) %*% t(f(-(i-1)))  
  h_N <- lambda * Linv %*% h_N + f(0)*y[i]
  theta_N <- solve(F_N)%*%h_N
}

# 4.4

dt$onestep_lamb_09  <- NA
dt$sixstep_lamb_09  <- NA
dt$twelvestep_lamb_09 <- NA

# Assuming lambda and y are defined above

for (i in 11:58){
  # 1 month ahead prediction

    F_N <- F_N + lambda^(i-1) * f(-(i-1)) %*% t(f(-(i-1)))  
    h_N <- lambda * Linv %*% h_N + f(0)*y[i]
    theta_N <- solve(F_N)%*%h_N
    yhat_1 <- t(f(1))%*%theta_N
    dt$onestep_lamb_09[i+1] <- yhat_1
    if (i + 6 <= 59){
      # 6 months ahead prediction
      yhat_6 <- t(f(6))%*%theta_N
      dt$sixstep_lamb_09[i+6] <- yhat_6
    }
    if (i + 12 <= 59){
      # 12 months ahead prediction
      yhat_12 <- t(f(12))%*%theta_N
      dt$twelvestep_lamb_09[i+12] <- yhat_12
    }
}

# 4.5

ggplot(dt, aes(x = TimeNumeric, y = Value)) +
  geom_point() +  # Plot the historical data points with default color
  geom_line(aes(y = yhat_ols, color = "OLS Fit"), linetype = "solid") +  # Plot the fitted line with color aesthetic
  geom_point(aes(y = onestep_lamb_09, color = "1-Month Forecast")) +  # Plot the forecasted points with color aesthetic
  geom_line(aes(y = onestep_lamb_09, color = "1-Month Forecast")) +  # Plot the forecasted points with color aesthetic
  geom_point(aes(y = sixstep_lamb_09, color = "6-Months Forecast")) +
  geom_line(aes(y = sixstep_lamb_09, color = "6-Months Forecast")) +
  geom_point(aes(y = twelvestep_lamb_09, color = "12-Months Forecast")) +
  geom_line(aes(y = twelvestep_lamb_09, color = "12-Months Forecast")) +
  scale_color_manual(values = c("OLS Fit" = "red", "1-Month Forecast" = "blue", "6-Months Forecast" = "green", "12-Months Forecast" = "#eea1ae"),
                     name = "Forecast Horizon", 
                     labels = c("OLS Fit" = "OLS Fit", "1-Month Forecast" = "1-Month", "6-Months Forecast" = "6-Months", "12-Months Forecast" = "12-Months")) +
  labs(title = "Historical Data and Forecast λ=0.9", x = "Time", y = "Value") +
  theme_minimal()



 #plot_N <- ggplot(dt[1:59,], aes(x=TimeNumeric, y=Value)) +
 #   geom_point() + 
 #   geom_point(data=dt[1:i,], col="blue") + 
 #   geom_line(data=dt[1:i,], aes(y=yhat_1[1:i]), col="blue")+
 #   geom_line(aes(y=yhat_ols), col="red", linetype=2) +  
 #   ggtitle(paste0("N = ", i))
 # 
 # print(plot_N)



# 4.6 


rmse_results <- list()


for (lambda in seq(0.55, 0.95, by = 0.01)) {
  
  rse_1 <- as.numeric()
  rse_6 <- as.numeric()
  rse_12 <- as.numeric()
  
  
  j <- 1
  (F_N <-  (lambda^0) * f(0)%*%t(f(0)))
  (h_N <-  (lambda^0) * f(0) * y[j])

  for (j in 2:10){
    F_N <- F_N + lambda^(j-1) * f(-(j-1)) %*% t(f(-(j-1)))  
    h_N <- lambda * Linv %*% h_N + f(0)*y[j]
    theta_N <- solve(F_N)%*%h_N
  }

  for (k in 11:58){
    # 1 month ahead prediction
    F_N <- F_N + lambda^(k-1) * f(-(k-1)) %*% t(f(-(k-1)))  
    h_N <- lambda * Linv %*% h_N + f(0)*y[k]
    theta_N <- solve(F_N)%*%h_N
    yhat_1 <- t(f(1))%*%theta_N
    dt$onestep_lamb_09[k+1] <- yhat_1
    if (k + 6 <= 59){
      # 6 months ahead prediction
      yhat_6 <- t(f(6))%*%theta_N
      dt$sixstep_lamb_09[k+6] <- yhat_6
    }
    if (k + 12 <= 59){
      # 12 months ahead prediction
      yhat_12 <- t(f(12))%*%theta_N
      dt$twelvestep_lamb_09[k+12] <- yhat_12
    }
  }

  for (i in 12:59){  
    rse_1 <- c(rse_1,(dt$Value[i] - dt$onestep_lamb_09[i])^2)
  }

  for (i in 17:59){
    rse_6 <- c(rse_6,(dt$Value[i] - dt$sixstep_lamb_09[i])^2)
  }

  for (i in 23:59){
    rse_12 <- c(rse_12,(dt$Value[i] - dt$twelvestep_lamb_09[i])^2)
  }

  rmse_1 <- sqrt(mean(rse_1))
  rmse_6 <- sqrt(mean(rse_6))
  rmse_12 <- sqrt(mean(rse_12))
 


  # Store RMSE values for the current lambda
    rmse_results[[as.character(lambda)]] <- list(
        "1 Month" = sqrt(mean(rse_1)),
        "6 Months" = sqrt(mean(rse_6)),
        "12 Months" = sqrt(mean(rse_12))
    )
}

# Convert RMSE results to a data frame for easier plotting
rmse_df <- do.call(rbind, lapply(rmse_results, function(x) as.data.frame(t(x))))
rmse_df$Lambda <- as.numeric(row.names(rmse_df))

#print('Idx, Optimal Lambda, RMSE for 1 Month Forecast Horizon: ', min_lambda_1_idx, rmse_df$'Lambda'[min_lambda_1_idx] min_lambda_1_rmse)

plot(seq(0.55, 0.95, by = 0.01), rmse_df$'1 Month', type = "l", main = "lambdas for 1month", xlab = "lambda", ylab = "S(lambda)")

plot(seq(0.55, 0.95, by = 0.01), rmse_df$'6 Month', type = "l", main = "lambda for 6months", xlab = "lambda", ylab = "S(lambda)")

plot(seq(0.55, 0.95, by = 0.01), rmse_df$'12 Month', type = "l", main = "lambda for 12months", xlab = "lambda", ylab = "S(lambda)")

# 4.7 Optimal Lambda for 1 Month Forecast Horizon

min_lambda_1_idx <- which.min(rmse_df$`1 Month`)
min_lambda_1_rmse <- rmse_df$'1 Month'[min_lambda_1_idx]

# 4.8 Optimal Lambda for 6 Month Forecast Horizon

min_lambda_6_idx <- which.min(rmse_df$`6 Month`)
min_lambda_6_rmse <- rmse_df$'6 Month'[min_lambda_6_idx]

# 4.9 Optimal Lambda for 12 Month Forecast Horizon

min_lambda_12_idx <- which.min(rmse_df$`12 Month`)
min_lambda_12_rmse <- rmse_df$'12 Month'[min_lambda_6_idx]

# 4.11 It would be problematic to make λ as small as 0.5. Why is that? (hint: compare the sum of weights to the number of parameters).


for (i in 12:58){  
  rmse_1 <- c(rse_1,(dt$Value[i] - dt$Value[i+1])^2)
}
rmse_n <- sqrt(mean(rmse_1))

print(rmse_n)


# Initialize a vector to store the naive_pred
naive_pred <- numeric()

for (i in 12:58) {
  naive_pred <- c(naive_pred, y[i])
}
rmse_naive <- sqrt(mean((y[13:59] - naive_pred)^2))

print(rmse_naive)




# 4.12 

data <- read_excel("DST_BIL54_test.xlsx")
colnames(data)[1] <- "Category"

dt1 <- filter(data, Category == "Drivmidler i alt")
dt1 <- dt1[, -1]

dt_test <- pivot_longer(dt1, 
                          cols = everything(),
                          names_to = "Time", 
                          values_to = "Value")

dt$Time <- x


