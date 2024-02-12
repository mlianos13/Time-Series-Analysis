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

data <- read_excel("assigment1/DST_BIL54_train.xlsx")
colnames(data)[1] <- "Category"

dt1 <- filter(data, Category == "Drivmidler i alt")
dt1 <- dt1[, -1]
# Create the time variable
months <- seq(as.Date("2018-01-01"), as.Date("2022-11-01"), by="month")
# months <- seq(as.Date("2018-01-01"), as.Date("2022-11-01"), by="month") creates a sequence of dates from January 2018 to November 2022, one for each month.

x <- 2018 + (as.numeric(format(months, "%m")) - 1) / 12 

dt <- pivot_longer(dt1, 
                          cols = everything(),
                          names_to = "Time", 
                          values_to = "Value")

dt$Time <- as.Date(months)

# Now, let's plot
ggplot(dt, aes(x = Time, y = Value)) + 
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
ggplot(dt, aes(x=Time, y=Value)) +
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
last_time_point <- as.Date("2022-11-01")

# Step 2: Generate future time points for the next 12 months
future_time_points <- seq(from = last_time_point, by = "month", length.out = 13)[-1]  # Generate one extra month and exclude the first

# If you need these as a design matrix for prediction (with an intercept):
Xtest <- cbind(1, future_time_points)


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
  Time = future_time_points, 
  y_pred = as.vector(y_pred), 
  y_pred_lwr = as.vector(y_pred_lwr), 
  y_pred_upr = as.vector(y_pred_upr)
)

# Plot the historical data and fitted line
p <- ggplot(dt, aes(x = Time, y = Value)) +
  geom_point() +
  geom_line(aes(y = yhat_ols), color = "blue") +
  labs(title = "Original Data and Forecast", x = "Time", y = "Value")

# Convert Time to numeric since ggplot2's geom_smooth works with numeric x-axis for prediction bands
dt$TimeNumeric <- as.numeric(dt$Time)
forecast_df$TimeNumeric <- as.numeric(forecast_df$Time)

# Add the forecasted points and prediction intervals
ggplot() +
  geom_point(data = dt, aes(x = Time, y = Value)) +  # Plot the historical data points
  geom_line(data = dt, aes(x = Time, y = yhat_ols), color = "red") +  # Plot the fitted line
  geom_point(data = forecast_df, aes(x = Time, y = y_pred), color = "red", shape = 1) +  # Plot the forecasted points
  geom_ribbon(data = forecast_df, aes(x = Time, ymin = y_pred_lwr, ymax = y_pred_upr), fill = "pink", alpha = 0.2) +  # Prediction intervals
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +  # Set x-axis labels to years only
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
