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


data <- read_excel("assigment1/DST_BIL54_train.xlsx")
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
  geom_line() + 
  labs(title = "Training Data over Time", x = "Time", y = "Value") +
  theme_minimal()



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

# estimate parameters with unweighted OLS (for comparison)
OLS <- solve(t(X[1:59,])%*%X[1:59,])%*%(t(X[1:59,])%*%y[1:59])
yhat_ols <- X[1:59,]%*%OLS

theta_0 <- WLS[1]
theta_1 <- WLS[2]
print(theta_0)
print(theta_1)


# 3.5 Make a forecast for the next 12 months - i.e., compute predicted values with corresponding prediction intervals

## λ = 0.9

y_pred <- Xtest%*%WLS
Vmatrix_pred <- sigma2_ols*(1+(Xtest%*%solve(t(X)%*%X))%*%t(Xtest))
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

y_pred <- Xtest%*%WLS
print(y_pred)
Vmatrix_pred <- sigma2_ols*(1+(Xtest%*%solve(t(X)%*%X))%*%t(Xtest))
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

y_pred <- Xtest%*%WLS
print(y_pred)
Vmatrix_pred <- sigma2_ols*(1+(Xtest%*%solve(t(X)%*%X))%*%t(Xtest))
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
 
y_pred <- Xtest%*%WLS
print(y_pred)
Vmatrix_pred <- sigma2_ols*(1+(Xtest%*%solve(t(X)%*%X))%*%t(Xtest))
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




# 4 Iterative update and optimal λ

# 4.1. Provide L and f(0) for the model

# Define the transition matrix L
#L <- matrix(c(1, 0, 1, 1), nrow = 2, byrow = TRUE)
L <- matrix(c(1.,0., 1.,1.),
            byrow=TRUE, nrow=2)
Linv <- solve(L) 

f <- function(j) rbind(1, j)
f_0 <- f(0)

# Print the transition matrix L and initial value of f(0)
cat("Transition matrix L:\n")
print(L)

cat("\nInitial value of f(0):\n")
print(f_0)

# 4.2  Provide F1 and h1.

# '''
#F: The matrix 
#F represents the parameter estimates for the local trend model. It is updated recursively based on the current parameter estimates and the values of 
#f(0) for each time step. It essentially captures the evolving trend structure of the time series data.
#
#
#h: The vector 
#h represents the intercept and slope of the local trend model. It is updated recursively based on the current values of 
#h and f(0) for each time step. It essentially captures the instantaneous trend direction and magnitude at each time point.
#'''

# Initialize lambda
lambda <- 0.9

# Calculate F1 using the updating formula
i <- 1
(F_1 <-  (lambda^0) * f(0)%*%t(f(0)))

# Calculate h1 using the updating formula
(h_1 <-  (lambda^0) * f(0) * y[i])
#h_1 <- lambda * Linv %*% f_0 * y_1


# Print the values of F1 and h1
cat("F1:\n")
print(F_1)

cat("\nh1:\n")
print(h_1)

#theta_1 <- solve(F_1)%*%h_1 # does not work for i=1, is not invertible!
#print(theta_1)

## 4.3  Using λ = 0.9 update FN and hN recursively and provide F10 and h10. We will not calculate predictions for these first 10 steps.


# Initialize F and h vectors/matrices
F <- list()
h <- list()
theta <- list()

F[[1]] <- F_1
h[[1]] <- h_1
# intitialize theta as double[2x1] value with 0
theta[[1]] <- matrix(c(0,0), nrow=2, ncol=1)
#theta[[1]] <- 0

# Iterate to update FN and hN up to 10 steps
for (i in 2:59) {
  # Update FN recursively
  F[[i]] <- F[[i-1]] + lambda^(i-1) * f(-(i-1)) %*% t(f(-(i-1)))
  
  # Update hN recursively
  h[[i]] <- lambda * Linv %*% h[[i-1]] + f(0) * y[i]
}

# Print the values of F10 and h10
cat("F10:\n")
print(F[[10]])

cat("\nh10:\n")
print(h[[10]])

#str(F)
#str(h)

# 4.4 Now update the model recursively up to F59 and h59, while also calculating predictions at each step. You should calculate predictions for 1 month ahead, 6 months ahead and 12 months ahead.


# want to save onestep predictions in dataframe (for lambda = 0.6):
dt$onestep_lamb_09  <- NA
dt$sixstep_lamb_09  <- NA
dt$twelvestep_lamb_09  <- NA

# initialize yhat:
yhat <- list()
# update iteratively:
for (i in 2:59){
  
  theta[[i]] <- solve(F[[i]])%*%h[[i]]
  yhat[[i]] <- t(f(-(i-1):(59-i)))%*%theta[[i]]

  # Assign predictions to the dataframe
  dt$onestep_lamb_09[i] <- yhat[[i]][1]
  dt$sixstep_lamb_09[i] <- yhat[[i]][6]
  dt$twelvestep_lamb_09[i] <- yhat[[i]][12]

}



## 4.5 Plot the predictions for 1 month ahead, 6 months ahead, and 12 months ahead over time.
#
## Prepare the predictions for plotting
#predictions <- data.frame(Time = future_x)
#for (i in 1:num_steps) {
#  for (j in c(1, 6, 12)) {
#    pred <- F[[i]] %*% h[[i]] + j * f(0)
#    predictions[paste0("Month_", j, "_ahead")] <- pred
#  }
#}
#
## Convert predictions to long format for easier plotting
#predictions_long <- pivot_longer(predictions, 
#                                  cols = starts_with("Month"), 
#                                  names_to = "Prediction_Type", 
#                                  values_to = "Prediction_Value")
#
## Combine the training data and predictions
#combined_data <- rbind(dt, predictions_long)
#
## Plot the training data and predictions
#ggplot(combined_data, aes(x = Time, y = Value, color = Prediction_Type, linetype = Prediction_Type)) +
#  geom_line() +
#  labs(title = "Training Data and Predictions", x = "Time", y = "Value") +
#  theme_minimal()
